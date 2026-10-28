#!/usr/bin/env python3
"""plot_det_C_vgeneral_phi_v2.py -- the eigenvalue search of the negative
mode (sector n = 0) and of the zero mode (sector n = 1): |det C| against s^2,
C the constant Wronskian built at r_max, with its located zero crossing.

THE OBJECT
    det C(n, s^2) = 0   <=>   omega^2_{n,0} = -s^2 is an eigenvalue of sector n
    C[a,b] = r^3 W = delta_ab + h^-[a,b](r_max)
                     + r_max I_nu K_nu(kappa_a r_max) h^-'[a,b](r_max)
    (at (n, s^2); the stored array, untransposed)
C is the constant (r-independent) Wronskian matrix of the two Bessel-dressed
solution branches, so its r = 0 value equals the value
computed at the outer boundary; the pipeline builds it there from the value
and slope of h^- = h_1 + h_2 + h_geq3 (exact, all Born orders).  The radial
Green's function carries C^{-1}, so where det C vanishes
it has the pole of an eigenvalue.  The production census (eigenvalue_neg /
eigenvalue_zero_vgeneral_phi_v2) locates the crossings by walking the
eigenvalue curves of C and roots them with brentq
(mode_census._crossing_in_bracket).

THE PICTURE (eigenvalues omega^2_{n,r}: the negative mode omega^2_{0,0}, the
zero mode omega^2_{1,0}; the pipeline keys are omega2_00 / omega2_10), two
panels, one right under the other:
  TOP     "Negative mode eigenvalue search": |det[C^(n=0)_ab(s^2)]| at its
          scan points (log-spaced in s^2 > 0), joined by a thin line, and the
          located zero crossing s^2 = -omega^2_{0,0} as a dashed vertical line;
  BOTTOM  "Zero mode eigenvalue search": the same for n = 1, whose scan is
          evenly spaced across s^2 = 0 (from just above the continuum
          threshold -min m_i^2); the crossing s^2 = -omega^2_{1,0} sits at ~0.
Both panels share ONE linear s^2 window (no log axis can show the n = 1
scan below s^2 = 0), s2_window: (-hi/9, hi), hi = the top of the stored
n = 0 scan (the variational ceiling) rounded up to a multiple of
0.5 x 10^k, 10^k <= top < 10^(k+1); the n = 1 points below its left edge
are cut, and the n = 0 points piled up near s^2 = 0 by the log spacing
overlap there.
Each legend sits below its panel; nothing else is drawn.

PURE READ: the curves are the s2_scan / det_scan arrays the eigenvalue stages
store in eig_n0/eig_n1_<tag>.npz -- no h-ODE is re-solved.  The marked roots
are the stored production eigenvalues (omega2_00, s2_c).

GATES (fail-loud)  det_scan == prod(lambda_scan) on every point (skipped,
with a printed note, when C has a complex eigenvalue pair: lambda_scan
stores real parts); every stored root inside its scan window; finite
arrays; bounce-fingerprint provenance via load_run.

OUTPUT   plot_det_C_vgeneral_phi_v2_<bounce>.png   in --out-dir
(default the run's pictures/pictures_<potential>_<bounce>/computation/);
<bounce> is the bounce name from the bounce npz (pair_name, e.g. M2_TV0)
-- figure names never carry the run tag/date.
"""
import argparse
import os
import sys

import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.dont_write_bytecode = True
_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

from plot_common_vgeneral_phi_v2 import (POINT_SIZE, add_standard_plot_cli,
                                         bounce_label, legend_outside, load_run,
                                         resolve_run, save, tex_number)

# one colour per role
COLOR_DET_C, COLOR_ROOT = '#2a78d6', '#eb6834'


def s2_window(scan_n0):
    """The ONE linear s^2 window both panels show, (-hi/9, hi): hi is the
    top of the stored n = 0 scan (the variational ceiling, above both
    crossings) rounded up to a multiple of step = 0.5 x 10^k, 10^k <= top
    < 10^(k+1); the n = 1 scan below -hi/9 is cut off (a top of 4.46 gives
    (-0.5, 4.5))."""
    s2_top = float(scan_n0['s2'][-1])
    step = 0.5 * 10.0 ** np.floor(np.log10(s2_top))
    hi = step * np.ceil(s2_top / step)
    return -hi / 9.0, hi


def sector_scan(sector_npz, n):
    """Everything one eig npz says about its det scan, unpacked once."""
    s2 = np.asarray(sector_npz['s2_scan'], float)
    det = np.asarray(sector_npz['det_scan'], float)
    lambda_scan = np.asarray(sector_npz['lambda_scan'], float)
    if not (np.all(np.isfinite(s2)) and np.all(np.isfinite(det))
            and np.all(np.isfinite(lambda_scan))):
        raise RuntimeError(f'[ABORT] non-finite scan arrays in eig_n{n} '
                           f'(fail-loud, no nansum).')
    # det C = prod of its eigenvalue curves -- the two stored arrays must
    # agree, else the npz mixes scans from different engine builds.  The
    # identity holds only while every eigenvalue is REAL: lambda_scan stores
    # real parts, so a complex conjugate pair (legitimate -- C is not
    # symmetric) breaks it without anything being wrong.  Skip the gate there
    # and say so; such a pair carries no mode anyway (it cannot make det C
    # vanish), which is why the census never tests it.
    dev = np.max(np.abs(det - np.prod(lambda_scan, axis=1))
                 / np.maximum(np.abs(det), 1e-30))
    if dev > 1e-9:
        imag = float(sector_npz['max_abs_imag']) if 'max_abs_imag' in sector_npz.files else None
        if imag is None:
            raise RuntimeError(
                f'[ABORT] eig_n{n}: det_scan != prod(lambda_scan) (max rel dev '
                f'{dev:.3e}).  This npz stores no max_abs_imag key, so a '
                f'legitimate complex eigenvalue pair cannot be told from a '
                f'mismatched npz -- re-run eigenvalue_neg / eigenvalue_zero_'
                f'vgeneral_phi_v2 for this tag.')
        if imag > 0.0:
            print(f'[plot_det_C] n={n}: complex eigenvalue pair '
                  f'present '
                  f'(max |Im| = {imag:.3e}), so det != prod(real parts) -- '
                  f'consistency gate skipped (such a pair carries no mode).')
        else:
            raise RuntimeError(f'[ABORT] eig_n{n}: det_scan != prod(lambda_scan) '
                               f'(max rel dev {dev:.3e}) -- inconsistent npz.')
    omega2 = float(sector_npz['omega2_00']) if n == 0 else float(sector_npz['omega2_10'])
    s2_c = -omega2                     # omega^2_{n,0} = -s_c^2 in BOTH sectors
    if not (s2[0] <= s2_c <= s2[-1]):
        raise RuntimeError(f'[ABORT] eig_n{n}: stored root s_c^2 = {s2_c:g} '
                           f'lies outside the stored scan window '
                           f'[{s2[0]:g}, {s2[-1]:g}].')
    sgn = np.sign(det)
    det_C_sign_changes = int(np.sum(sgn[:-1] * sgn[1:] < 0.0))
    return dict(n=n, s2=s2, det=det, omega2=omega2, s2_c=s2_c,
                n_crossings=len(np.atleast_1d(
                    np.asarray(sector_npz[f'omega2_{n}0_census'], float))),
                resid=float(sector_npz['crossing_resid']),
                det_C_sign_changes=det_C_sign_changes,
                title=('Negative mode eigenvalue search' if n == 0
                       else 'Zero mode eigenvalue search'),
                # the numbering omega^2_{n,r} of the labels: wave n first
                sym=(r'\omega^2_{0,0}' if n == 0 else r'\omega^2_{1,0}'))


def panel(ax, f, window):
    """One sector: |det C| at its scan points and the located root, on the
    s^2 window of both panels."""
    det_C = rf'$|\det[C^{{(n={f["n"]})}}_{{\alpha\beta}}(s^2)]|$'
    # the data points and lines in the one set of sizes (plot_common)
    ax.plot(f['s2'], np.abs(f['det']), color=COLOR_DET_C, marker='o',
            ms=POINT_SIZE, mew=0, label=det_C)
    ax.axvline(f['s2_c'], color=COLOR_ROOT, ls='--',
               label=(rf'zero-crossing location at $s^2_z=-{f["sym"]}'
                      rf'={tex_number(f["s2_c"])}$'))
    ax.set_xlim(*window)
    ax.set_yscale('log')
    ax.set_xlabel(r'$s^2$')
    ax.set_ylabel(det_C)
    ax.set_title(f['title'])
    ax.grid(True, which='both')     # the one grid style of plot_common
    legend_outside(ax, 'below', below=-0.25)     # clear of the s^2 label


def main():
    ap = argparse.ArgumentParser(
        description='det C (the constant, r-independent Wronskian, built '
                    'at r_max) vs s^2 for sectors n=0 and n=1, roots marked')
    add_standard_plot_cli(ap)
    args = ap.parse_args()
    resolve_run(args)      # bare invocation: auto-discover tag + bounce npz

    run = load_run(args.data_dir, args.tag, args.bounce_npz,
                   need=('eig_n0', 'eig_n1'))
    scans_n0_n1 = {n: sector_scan(run[f'eig_n{n}'], n) for n in (0, 1)}
    for n, scan in scans_n0_n1.items():
        print(f'[det n={n}] scan {len(scan["s2"])} points on '
              f'[{scan["s2"][0]:g}, {scan["s2"][-1]:g}]: {scan["det_C_sign_changes"]} sign '
              f'change(s); root s_c^2 = {scan["s2_c"]:.12g} -> '
              f'omega2_{n}0 = {scan["omega2"]:+.12e} '
              f'({scan["n_crossings"]} located crossing(s), '
              f'|det C| = {scan["resid"]:.1e})')

    # the negative mode on top, the zero mode right under it, each legend
    # below its panel (a picture <= 8 in wide)
    fig, axes = plt.subplots(2, 1, figsize=(6.4, 8.4),
                             gridspec_kw={'hspace': 0.75})
    window = s2_window(scans_n0_n1[0])
    for ax, n in zip(axes, (0, 1)):
        panel(ax, scans_n0_n1[n], window)
    save(fig, args.out_dir,
         f'plot_det_C_vgeneral_phi_v2_{bounce_label(args)}.png')
    plt.close(fig)


if __name__ == '__main__':
    main()
