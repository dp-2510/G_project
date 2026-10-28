#!/usr/bin/env python3
"""traces_vs_s2_n_common_vgeneral_phi_v2.py -- the machinery behind the
stored-data trace picture plot_core_eq_fictitious_vs_s2_vgeneral_phi_v2.py
(the core equation with the dangerous modes subtracted and the fictitious
modes added, against the upper limit Lambda^2 of its s^2-integral).

Not a plot code itself: it reads + CERTIFIES the band tables and sums them
into the trace curves, and attaches the raw n = 0, 1 traces with their two
poles and the integrand of that core equation.

THE OBJECTS (all deg-weighted, deg_n = (n+1)^2).  The residual BAND covers
the safe waves n = 2..n_max; the n = 0 negative-mode and n = 1 zero-mode
sectors are the DANGEROUS ones (each has a resolvent pole in s^2), handled
in the assembly by their own pole-subtracted sector integrals.  The curves
attached below still start at n = 0: the RAW traces [Gbar - Gbar_FV]_{n=0,1}
-- finite away from the poles, stored in every BASE band slice (keys
Gbar_minus_fv_n0/_n1, add_trace_G_minus_G_fv_n0_n1 below; a trace that is
NaN -- a slice exactly ON a pole -- aborts) -- are added to the band waves:

    trace[G - G_FV](s^2)   = sum_{n=0}^{n_max} [Gbar_bounce - Gbar_FV]_n
                             the UNSUBTRACTED trace difference (raw n = 0, 1
                             traces + band waves): box-independent but
                             UV-divergent -- its full n-sum grows without
                             bound (the Born-1,2 content IS the UV divergence);
    core_eq_fict_all(s^2)    the same sum with the two pole terms subtracted
                             and the fictitious modes added
                             (add_core_eq_fictitious below).

The Born traces Gbar_1, Gbar_2 and the Born-subtracted remainder delta_geq3
of the band waves are summed as well: they certify the slice set (below).

EVERY CURVE STOPS AT n_max.  The waves beyond n_max enter lnD_ren only at
the determinant level, through the one fit of the s2-integrated,
Born-subtracted I_n(inf) and its Hurwitz-zeta sum
(computation/tail_high_n_zeta_vgeneral_phi_v2); no completion of a
fixed-s^2 object is made anywhere in the pipeline.

DATA -- everything from the pipeline npz, nothing refitted here:
    residual band slices       the per-wave traces (plot_common.band_tables_
                               multikey, certified slice set) summed here
                               into the n = 2..n_max trace curves, and the
                               raw n = 0, 1 traces (plot_common.trace_G_
                               minus_G_fv_n0_n1_on_band_grid, slice keys
                               Gbar_minus_fv_n0/_n1);
    band_integrals_<tag>.npz   n_max, s2_max, the certified band rectangle
                               (slices above s2_max are not read); ABSENT -> the
                               plot code aborts (load_run);
    eig_n0/_n1_<tag>.npz       omega2_00 and omega2_10, the two pole positions;
    lnD_ren_<tag>.npz          deg1 (the n = 1 residue) and m2fict
                               (the fictitious mass^2 of the assembly);
                               a missing eig or lnD_ren file aborts, naming
                               the stage that writes it.
CERTIFICATION before anything is drawn (load_stored_traces): the Born-split
identity trace_delta_geq3 = trace_Gbar_minus_fv - trace_Gbar_1 -
trace_Gbar_2 must hold on the summed curves (fail-loud).

PICTURE NAME: <code>_<bounce>.png -- the plot code's name and the bounce
name, never the run tag or a date.
"""
import os
import sys

import numpy as np

import matplotlib
matplotlib.use('Agg')

sys.dont_write_bytecode = True
_HERE = os.path.dirname(os.path.abspath(__file__))
_COMP = os.path.normpath(os.path.join(_HERE, '..', 'computation'))
for _p in (_HERE, _COMP):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from plot_common_vgeneral_phi_v2 import (add_standard_plot_cli, bounce_label,
                                      trace_G_minus_G_fv_n0_n1_on_band_grid, band_tables_multikey,
                                      load_run, resolve_run)

# one colour per role, uniform across the figures of this family
COLOR_TRUNCATED = 'tab:blue'     # the n_max-truncated object
# tolerance of the Born-split identity on the summed curves (the same slice
# files, the same ascending-n summation -> float rounding only)
REPRO_RTOL = 1e-12


# --------------------------------------------------------------------------- #
#  The CLI of a thin plot code                                                 #
# --------------------------------------------------------------------------- #
def make_cli(description):
    """The shared CLI: the standard plot selectors.  Returns the
    parsed-and-resolved args."""
    import argparse
    ap = argparse.ArgumentParser(description=description)
    add_standard_plot_cli(ap)
    args = ap.parse_args()
    resolve_run(args)      # bare invocation: auto-discover tag + bounce npz
    return args


# --------------------------------------------------------------------------- #
#  The stored trace curves + band tables, certified                            #
# --------------------------------------------------------------------------- #
def load_stored_traces(args):
    """Read the certified band tables, sum them into the trace curves, verify
    the Born-split identity, and return one dict traces:

        s2                   the band s^2 grid
        trace_Gbar_minus_fv  sum_{n=2..n_max} [Gbar - Gbar_FV]_n   per s^2
        trace_Gbar_1,
        trace_Gbar_2         the summed Born-1/2 traces
        trace_delta_geq3     the summed Born-subtracted remainder
        n_max                the certified band top (band_integrals npz)
        run, background      the loaded run (band_integrals guaranteed) + background
        tab, ns              the certified band tables + their wave axis

    ABORTS (fail-loud) when band_integrals_<tag>.npz is absent (load_run),
    belongs to a different bounce (load_run's provenance check), or when the
    summed curves violate the Born-split identity."""
    # the whole run first: every stage npz is opened ONCE (by load_run); the
    # certified band rectangle (n_max, s2_max) is the one band_integrals stored
    run = load_run(args.data_dir, args.tag, args.bounce_npz,
                   need=('band_integrals',))
    n_max = int(run['band_integrals']['n_max'])
    s2_max = float(run['band_integrals']['s2_max'])
    # ---- the per-wave traces (certified slice set), summed over n = 2..n_max --
    tab = band_tables_multikey(args.data_dir, args.tag, args.bounce_npz,
                               n_max=n_max, s2_cap=s2_max)
    band_n = tab['n']
    s2 = np.asarray(tab['s2'], float)
    curves = {}
    for name, key in (('trace_Gbar_minus_fv', 'Gbar_minus_fv'),
                      ('trace_Gbar_1', 'Gbar_1'),
                      ('trace_Gbar_2', 'Gbar_2'),
                      ('trace_delta_geq3', 'delta_geq3')):
        arr = np.asarray(tab[key], float).sum(axis=0)      # sum n = 2..n_max
        if not np.all(np.isfinite(arr)):
            raise RuntimeError(f'[ABORT] non-finite {name} after summing the '
                               f'band tables (fail-loud).')
        curves[name] = arr
    # ---- certification: the Born-split identity must hold -----------------
    dev = float(np.max(np.abs(curves['trace_delta_geq3']
                              - (curves['trace_Gbar_minus_fv']
                                 - curves['trace_Gbar_1']
                                 - curves['trace_Gbar_2']))))
    scale = float(np.max(np.abs(curves['trace_delta_geq3'])))
    if dev > REPRO_RTOL * max(scale, 1e-30):
        raise RuntimeError(f'[ABORT] Born-split identity violated on the '
                           f'summed curves (dev {dev:.2e}, scale {scale:.2e}) '
                           f'-- stale or foreign slice set.')
    print(f'[traces_vs_s2_n] band tables summed over n = {band_n[0]}..{band_n[-1]} '
          f'on {len(s2)} stored s^2 values; Born-split identity '
          f'trace_delta_geq3 = trace_Gbar_minus_fv - trace_Gbar_1 - '
          f'trace_Gbar_2 holds (dev {dev:.1e})')
    return dict(s2=s2, n_max=n_max, run=run, background=run['bg'], tab=tab, band_n=band_n,
                **curves)


def add_trace_G_minus_G_fv_n0_n1(traces, args):
    """Attach the RAW n = 0 and n = 1 deg-weighted traces
    [Gbar - Gbar_FV]_n on the band s^2 grid to traces (they come from the BASE
    band slices, keys Gbar_minus_fv_n0/_n1, read through
    plot_common.trace_G_minus_G_fv_n0_n1_on_band_grid on the certified slice set):

        Gbar_minus_fv_n0,
        Gbar_minus_fv_n1         the two traces, aligned with traces['s2']
        trace_Gbar_minus_fv_all  Gbar_minus_fv_n0 + Gbar_minus_fv_n1
                                 + trace_Gbar_minus_fv  -- the n = 0..n_max sum

    WHY THEY BELONG IN THE G - G_FV SUM: n = 0, 1 are the DANGEROUS
    sectors -- each has a resolvent pole in s^2, so their s^2-INTEGRALS into
    lnD_ren go through the pole-subtracted sector stages, never through the
    band.  But the raw TRACES are perfectly finite away from the poles,
    and the trace of [G - G_FV] is a sum over ALL waves.

    GUARD (fail-loud): a NaN trace entry means that stored slice sits exactly
    ON a dangerous-mode pole (det C = 0) -- an n >= 0 trace figure cannot be
    drawn through a pole trace, so this aborts NAMING the offending s^2."""
    Gbar_minus_fv_n0, Gbar_minus_fv_n1 = trace_G_minus_G_fv_n0_n1_on_band_grid(traces['tab'])
    for name, trace_n in (('n = 0', Gbar_minus_fv_n0), ('n = 1', Gbar_minus_fv_n1)):
        bad = np.asarray(traces['s2'], float)[~np.isfinite(trace_n)]
        if bad.size:
            raise RuntimeError(
                f'[ABORT] the raw {name} trace is NaN at s^2 = '
                f'{", ".join(f"{v:g}" for v in bad)} -- that slice sits '
                f'exactly ON the dangerous-mode pole (det C = 0), so the '
                f'core-equation figure cannot include it.  Re-run the band '
                f'off the pole.')
    traces['Gbar_minus_fv_n0'], traces['Gbar_minus_fv_n1'] = Gbar_minus_fv_n0, Gbar_minus_fv_n1
    traces['trace_Gbar_minus_fv_all'] = ((Gbar_minus_fv_n0 + Gbar_minus_fv_n1)
                                    + traces['trace_Gbar_minus_fv'])
    print(f'[traces_vs_s2_n] raw n = 0, 1 traces attached (finite on all '
          f'{len(traces["s2"])} stored s^2 values): the G - G_FV sums run '
          f'over n = 0..{traces["n_max"]}')
    add_sector_poles(traces)
    return traces


def add_sector_poles(traces):
    """Attach the two RESOLVENT POLES the unsubtracted n = 0, 1 traces carry,
    their two most prominent features:

        s2_c_n0       the n = 0 crossing s^2 = -omega^2_{0,0} > 0 (det C = 0),
                      the negative-mode pole: the unsubtracted n = 0 trace
                      diverges there and CHANGES SIGN across it (the
                      eigenvalue is negative, so the pole sits at a
                      positive s^2, inside the plotted range);
        omega2_10,
        deg1          the n = 1 zero-mode pole term deg1/(omega^2_{1,0}+s^2),
                      deg1 read from the lnD_ren npz
                      (= (n+1)^2 = 4 at n = 1, the sector stage's deg):
                      since omega^2_{1,0} = 0 to numerical accuracy, this
                      is the ~ 4/s^2 rise of the unsubtracted n = 1 trace
                      towards small s^2, not a feature of the band.

    Both eigenvalues come from the eig_n0/eig_n1 stages, which run BEFORE the
    band, so this adds no new dependency: the stage npz are already loaded
    (and provenance-checked) in traces['run'].  add_core_eq_fictitious
    subtracts both pole terms."""
    omega2_n0_n1 = {}
    for stage, key in (('eig_n0', 'omega2_00'), ('eig_n1', 'omega2_10')):
        eig_npz = traces['run'].get(stage)
        if eig_npz is None:
            raise SystemExit(
                f'[ABORT] {os.path.basename(traces["run"]["paths"][stage])} not '
                f'found -- the pole positions come from it; re-run '
                f'{"eigenvalue_neg" if stage == "eig_n0" else "eigenvalue_zero"}'
                f'_vgeneral_phi_v2 for this tag.')
        omega2_n0_n1[key] = float(eig_npz[key])
    if traces['run'].get('lnD_ren') is None:
        raise SystemExit(
            f'[ABORT] {os.path.basename(traces["run"]["paths"]["lnD_ren"])} not '
            f'found -- the n = 1 residue (key deg1) is read from '
            f'it; re-run assemble_lnD_vgeneral_phi_v2 for this tag.')
    traces['s2_c_n0'] = -omega2_n0_n1['omega2_00']
    traces['omega2_10'] = omega2_n0_n1['omega2_10']
    traces['deg1'] = float(traces['run']['lnD_ren']['deg1'])
    # warn when a stored slice sits so close to the n = 0 pole that its
    # unsubtracted trace is dominated by it
    s2 = np.asarray(traces['s2'], float)
    j = int(np.argmin(np.abs(s2 - traces['s2_c_n0'])))
    gap = abs(s2[j] - traces['s2_c_n0']) / traces['s2_c_n0']
    print(f'[traces_vs_s2_n] poles of the unsubtracted traces: n = 0 at '
          f's^2 = {traces["s2_c_n0"]:.6f} (= -omega^2_00), n = 1 term '
          f'{traces["deg1"]:g}/({traces["omega2_10"]:+.3e} + s^2) ~ '
          f'{traces["deg1"]:g}/s^2 at small s^2')
    if gap < 0.05:
        print(f'[WARN DANGEROUS-MODE] the stored slice s^2 = {s2[j]:g} sits '
              f'{100 * gap:.1f}% from the n = 0 pole s^2 = '
              f'{traces["s2_c_n0"]:.6f}: its raw n = 0 trace is pole-dominated, '
              f'so the integrand there is the difference of two large '
              f'numbers (the trace and the pole term).')
    return traces


def add_core_eq_fictitious(traces):
    """Attach the integrand of the core equation with the dangerous modes
    subtracted and the fictitious modes added -- wave by wave and summed, on
    the band s^2 grid:

        Gbar(s^2) - Gbar^v(s^2) - 1/(omega^2_{0,0} + s^2)
                  - deg1/(omega^2_{1,0} + s^2) + (1 + deg1)/(m_f^2 + s^2)

    (omega^2_{n,r} with the wave n first: omega^2_{1,0} is the pipeline's
    omega2_10), the n-sums of Gbar and Gbar^v stopped at n_max.  One
    fictitious mode replaces one removed mode (assemble_lnD), so the pole and
    fictitious terms belong to the two dangerous waves:

        core_eq_fict_n0   [Gbar - Gbar_FV]_0 - 1/(omega2_00 + s^2) + 1/(m_f^2 + s^2)
        core_eq_fict_n1   [Gbar - Gbar_FV]_1 - deg1/(omega2_10 + s^2)
                          + deg1/(m_f^2 + s^2)
        (n >= 2: the band traces [Gbar - Gbar_FV]_n themselves)
        core_eq_fict_all  core_eq_fict_n0 + core_eq_fict_n1 + trace_Gbar_minus_fv,
                          the sum n = 0..n_max

    m_f^2 = m2fict, the fictitious mass^2 the assembly used (lnD_ren npz;
    every m_f^2 > 0 gives the same lnD_ren, the + 5 ln(m_f^2) outside the
    integral compensates).  Needs add_trace_G_minus_G_fv_n0_n1 first: the
    raw n = 0, 1 traces at the stored band nodes and the two poles.  The pole
    terms are subtracted AT the stored nodes -- nothing is interpolated."""
    s2 = np.asarray(traces['s2'], float)
    m2_f = float(traces['run']['lnD_ren']['m2fict'])
    omega2_00, deg1 = -traces['s2_c_n0'], traces['deg1']
    core_eq_fict_n0 = (traces['Gbar_minus_fv_n0'] - 1.0 / (omega2_00 + s2)
                       + 1.0 / (m2_f + s2))
    core_eq_fict_n1 = (traces['Gbar_minus_fv_n1'] - deg1 / (traces['omega2_10'] + s2)
                       + deg1 / (m2_f + s2))
    for name, trace_n in (('n = 0', core_eq_fict_n0), ('n = 1', core_eq_fict_n1)):
        if not np.all(np.isfinite(trace_n)):
            raise RuntimeError(f'[ABORT] non-finite core-equation integrand '
                               f'at {name} (fail-loud).')
    traces['m2_f'] = m2_f
    traces['core_eq_fict_n0'], traces['core_eq_fict_n1'] = core_eq_fict_n0, core_eq_fict_n1
    traces['core_eq_fict_all'] = (core_eq_fict_n0 + core_eq_fict_n1
                                  + traces['trace_Gbar_minus_fv'])
    print(f'[traces_vs_s2_n] core-equation integrand: poles '
          f'1/(omega2_00 + s^2) and {deg1:g}/(omega2_10 + s^2) subtracted, '
          f'{1 + deg1:g} fictitious modes at m_f^2 = {m2_f:.6g} added; summed '
          f'over n = 0..{traces["n_max"]}')
    return traces


def picture_name(code_stem, args):
    """The picture name of a thin code: <code>_<bounce>.png."""
    return f'{code_stem}_{bounce_label(args)}.png'
