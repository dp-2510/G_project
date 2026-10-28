#!/usr/bin/env python3
"""plot_lnD_ren_high_n_fit_vgeneral_phi_v2.py -- the fit behind the high-n
tail Zcal of the final result ln D': the complete partial-wave integrals of
the computed waves, scaled by (n+1)^3, against n, with the fitted form
through them and the two n -> infinity limits -- the fitted leading
coefficient a and the analytic one, A3.

THE OBJECTS (lnD_ren npz; the tail assembly of tail_high_n_vgeneral_phi_v2)
    I_n(s^2 = inf) = I_n(s2_max) + T_n(s2_max)      every computed wave n = 2..n_max
    fit            a/nu^3 + b/nu^5 + c/nu^7, nu = n + 1, fitted on the waves
                   n = fit_lo..n_max (fit_lo: the stored start of the window)
    Zcal           sum_{n > n_max} of the fit = a zeta(3, n_max+2)
                   + b zeta(5, n_max+2) + c zeta(7, n_max+2)
    A3             the leading coefficient predicted by the third Born order,
                   A3/nu^3 for nu -> inf: the analytic check of a

THE PICTURE  one panel, "High n-fitfunction": (n+1)^3 I_n(s^2 = inf) at the
computed waves as points; the fit times (n+1)^3 = a + b/nu^2 + c/nu^4 as a
line from fit_lo to N_CONTINUED * n_max (beyond n_max it is the continuation
Zcal sums), the start and the end of the fit window (dotted, dashed); the two
n -> infinity limits as horizontal lines, the fitted a (dotted) and the
predicted A3 (dashed).  The legend sits below the panel and carries every
number: the three fitted coefficients in the fit's label, a and A3 in the
labels of their lines.

PURE READ: band_n, I_n_inf, fit_a, fit_b, fit_c, fit_lo, A3, n_max, Zcal of
the lnD_ren npz; nothing is refitted.
GATES (fail-loud)  finite arrays; the stored Zcal equals the zeta sum of the
stored coefficients; the fit window lies inside the computed waves.
OUTPUT   plot_lnD_ren_high_n_fit_vgeneral_phi_v2_<bounce>.png   in --out-dir
(default the run's pictures/pictures_<potential>_<bounce>/computation/);
<bounce> is the bounce name from the bounce npz (pair_name, e.g. M2_TV0)
-- figure names never carry the run tag/date.
"""
import argparse
import os
import sys

import numpy as np
from scipy.special import zeta

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.dont_write_bytecode = True
_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

from plot_common_vgeneral_phi_v2 import (DOTTED_LINE_WIDTH, POINT_SIZE,
                                         add_standard_plot_cli, bounce_label,
                                         legend_outside, load_run, resolve_run,
                                         save, tex_number)

#: the fit is drawn on to N_CONTINUED * n_max: the continuation Zcal sums
N_CONTINUED = 3
#: the stored Zcal must be the zeta sum of the stored coefficients
CHECK_TOL = 1e-9
# one colour per role
COLOR_DATA, COLOR_FIT, COLOR_A3 = 'tab:blue', 'tab:orange', 'tab:green'


def tail_fit(lnD_ren_npz):
    """Everything the panel needs, unpacked and gated once."""
    n = np.asarray(lnD_ren_npz['band_n'], int)
    I_inf = np.asarray(lnD_ren_npz['I_n_inf'], float)
    a, b, c = (float(lnD_ren_npz[k]) for k in ('fit_a', 'fit_b', 'fit_c'))
    A3, Zcal = float(lnD_ren_npz['A3']), float(lnD_ren_npz['Zcal'])
    n_max, fit_lo = int(lnD_ren_npz['n_max']), int(lnD_ren_npz['fit_lo'])
    if not (np.all(np.isfinite(I_inf)) and np.isfinite([a, b, c, A3, Zcal]).all()):
        raise RuntimeError('[ABORT] non-finite tail data in the lnD_ren npz '
                           '(fail-loud, no nansum).')
    if not (n[0] <= fit_lo <= n_max == n[-1]):
        raise RuntimeError(f'[ABORT] fit window n = {fit_lo}..{n_max} is not '
                           f'inside the computed waves n = {n[0]}..{n[-1]}.')
    q = n_max + 2.0
    Zcal_zeta = a * zeta(3.0, q) + b * zeta(5.0, q) + c * zeta(7.0, q)
    if abs(Zcal_zeta - Zcal) > CHECK_TOL * max(1.0, abs(Zcal)):
        raise RuntimeError(f'[ABORT] stored Zcal = {Zcal:.10g} is not the zeta '
                           f'sum of the stored coefficients, {Zcal_zeta:.10g}.')
    nu = n + 1.0
    on_window = n >= fit_lo
    fit_on_data = a / nu**3 + b / nu**5 + c / nu**7
    dev = np.max(np.abs(fit_on_data[on_window] - I_inf[on_window])
                 / np.abs(I_inf[on_window]))
    print(f'[plot_high_n_fit] fit a/nu^3 + b/nu^5 + c/nu^7 on n = {fit_lo}..{n_max}: '
          f'a = {a:.6g}, b = {b:.6g}, c = {c:.6g}; predicted A3 = {A3:.6g} '
          f'(a/A3 = {a / A3:.4f}); largest relative deviation of the fit from '
          f'I_n(inf) on the window {dev:.2e}; Zcal = {Zcal:+.6f}')
    n_fit = np.arange(fit_lo, N_CONTINUED * n_max + 1)
    nu_fit = n_fit + 1.0
    return dict(n=n, nu=nu, I_inf=I_inf, a=a, b=b, c=c, A3=A3, Zcal=Zcal,
                n_max=n_max, fit_lo=fit_lo, n_fit=n_fit,
                scaled_fit=a + b / nu_fit**2 + c / nu_fit**4)


def signed_term(x, power):
    """' + |x|/(n+1)^power' or ' - |x|/(n+1)^power' for the fit's label."""
    sign = '+' if x >= 0.0 else '-'
    return rf' {sign} {tex_number(abs(x))}/(n+1)^{power}'


def panel(ax, f):
    """(n+1)^3 I_n(s^2 = inf) against n: the computed waves (points), the
    scaled fit on its window and beyond (line), the window (dotted, dashed)
    and the two n -> infinity limits (horizontal lines)."""
    ax.plot(f['n'], f['nu']**3 * f['I_inf'], 'o', color=COLOR_DATA,
            ms=POINT_SIZE, mew=0, zorder=4, label=r'$(n+1)^3\, I_n(s^2 = \infty)$')
    ax.plot(f['n_fit'], f['scaled_fit'], '-', color=COLOR_FIT, zorder=3,
            label=(rf'fit: ${tex_number(f["a"])}' + signed_term(f['b'], 2)
                   + signed_term(f['c'], 4) + '$'))
    ax.axvline(f['fit_lo'], color='0.35', ls=':', lw=DOTTED_LINE_WIDTH,
               label=rf'fit window starts at: $n = {f["fit_lo"]}$')
    ax.axvline(f['n_max'], color='0.35', ls='--', lw=DOTTED_LINE_WIDTH,
               label=rf'fit window ends at: $n_{{max}} = {f["n_max"]}$')
    ax.axhline(f['a'], color=COLOR_FIT, ls=':', lw=DOTTED_LINE_WIDTH, zorder=2,
               label=(r'$n \to \infty$ limit of the fitfunction: '
                      rf'$a = {tex_number(f["a"])}$'))
    ax.axhline(f['A3'], color=COLOR_A3, ls='--', zorder=2,
               label=(r'$n \to \infty$ limit of the analytical coefficient: '
                      rf'$A_3 = {tex_number(f["A3"])}$'))
    ax.set_xlim(0, f['n_fit'][-1])
    # the plateau and the approach to it: from 1.5 % below the lower of the
    # two coefficients to 17 % above it; the low waves, far above, are cut
    lowest = min(f['a'], f['A3'])
    ax.set_ylim(1.015 * lowest, 0.83 * lowest)
    ax.set_xlabel(r'$n$')
    ax.set_ylabel(r'$(n+1)^3\, I_n(s^2 = \infty)$')
    ax.set_title(r'High $n$-fitfunction')
    ax.grid(True, which='both')     # the one grid style of plot_common
    legend_outside(ax, 'below', below=-0.25)     # clear of the n label


def main():
    ap = argparse.ArgumentParser(
        description='vgeneral_phi_v2: the fit behind the high-n tail Zcal -- '
                    '(n+1)^3 I_n(s^2 = inf) against n, the fit and its two limits')
    add_standard_plot_cli(ap)
    args = ap.parse_args()
    resolve_run(args)      # bare invocation: auto-discover tag + bounce npz
    run = load_run(args.data_dir, args.tag, args.bounce_npz, need=('lnD_ren',))
    f = tail_fit(run['lnD_ren'])
    # a picture <= 8 in wide (sizes: plot_common, the one set)
    fig, ax = plt.subplots(figsize=(7.4, 3.8))
    panel(ax, f)
    save(fig, args.out_dir,
         f'plot_lnD_ren_high_n_fit_vgeneral_phi_v2_{bounce_label(args)}.png')
    plt.close(fig)


if __name__ == '__main__':
    main()
