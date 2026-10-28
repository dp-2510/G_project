#!/usr/bin/env python3
"""plot_common_vgeneral_phi_v2.py -- the SHARED read layer of the vgeneral_phi_v2
diagnostic plots.  Not a physics stage: it reads the production npz set and
re-exposes the production objects.  The one formula it re-spells is the
assembly sum of lnD_from_parts, the seven terms of assemble_lnD_vgeneral_phi_v2
(Sigma_dang = Delta_0 + Delta_1 + fict_block +
ln|omega2_00|; lnD_ren = (T1) Sigma_dang + (T2) Born-order-1 add-back on
[0, s2_max] + (T3) sum I_n(s2_max) + (T4) Born-order-1 add-back beyond s2_max + (T5) sum
T_n(s2_max) + (T6) Zcal + (T7) (-A2_fin/2)), written term for term and in the
same order as assemble_lnD_vgeneral_phi_v2.main writes it, so the stored
scalar can be rebuilt and checked; fict_block_terms is imported.

WHY THIS FILE EXISTS
--------------------
The production reader pipeline_helpers_vgeneral_phi_v2.read_band_tables merges the
residual-band slice files for ONE key (delta_geq3, the object that enters
lnD_ren).  The diagnostic plots also need the UNSUBTRACTED per-wave objects
that live in the same slice files:

    Gbar_minus_fv = Gbar_bounce - Gbar_FV   (deg-weighted, per wave, per s2)
    Gbar_1   = Gbar^(1)  first Born trace   (tadpole insertion)
    Gbar_2   = Gbar^(2)  second Born trace  (fish insertion)
    delta_geq3 = Gbar_minus_fv - Gbar_1 - Gbar_2   (all Born orders >= 3)

so this module merges the SAME certified slice set for all four keys.  The
certification (bounce-fingerprint stamp, MANIFEST set equality,
duplicate-wave consistency) is NOT re-implemented: read_band_tables is called
first and its file list is the only one this module is allowed to read.

CONTENTS
    band_tables_multikey(...)   the (n, s2) tables of all four band keys
    trace_G_minus_G_fv_n0_n1_on_band_grid(...)       the raw n = 0, 1 waves on the band s2 grid
                                (the slice keys Gbar_minus_fv_n0/_n1 -- no
                                side cache, no re-solve)
    load_run(...)               every production npz of a tag, in one object
    lnD_from_parts(...)         reproduce lnD_ren from the stored parts
    certify_lnD(...)            lnD_from_parts checked against the stored scalar
    bounce_label(...)           the bounce name of a picture name
    add_standard_plot_cli / resolve_run(...)  the uniform CLI + auto-discovery
    figure_folder(...)          THE folder a figure is saved in (one rule for
                                every plot code of the pipeline)
    legend_outside(...)         THE legend rule: never on top of a curve
    save(...)                   the one writer of a figure (bbox tight)

Every helper here is fail-loud: a missing file, a non-finite value or a
reproduction mismatch raises, it never silently degrades a figure.
"""
import hashlib
import os
import sys

import numpy as np
import matplotlib
import matplotlib.patheffects as patheffects

# ONE grid for every figure (2-D and 3-D axes, major and minor lines alike):
# opaque light grey, 1.3 pt -- thinner or transparent lines wash out unevenly
# when a PDF viewer shrinks the picture, so some lines would vanish
matplotlib.rcParams.update({'grid.color': '0.72', 'grid.linewidth': 1.3,
                            'grid.alpha': 1.0, 'grid.linestyle': '-'})
# ONE set of sizes for every figure.  Every picture is meant to be printed
# at the same scale, \includegraphics[scale=0.69] (so a picture must stay
# <= 8 in wide), and these sizes then print alike in every figure: tick labels
# 7.6 pt, legends 8.3 pt, axis labels and every title 9 pt (titles thicker:
# thick_titles, called by save), lines 1.0 pt, data points 2.8 pt
matplotlib.rcParams.update({'font.size': 11, 'axes.labelsize': 13,
                            'xtick.labelsize': 11, 'ytick.labelsize': 11,
                            'legend.fontsize': 12, 'legend.title_fontsize': 12,
                            'axes.titlesize': 13, 'figure.titlesize': 13,
                            'lines.linewidth': 1.5})
#: every line of a figure (solid, dashed, dash-dotted) and every data point
#: (marker size, drawn without an edge: mew=0); a DOTTED line is drawn
#: thicker, its dots would otherwise read lighter than the other lines
LINE_WIDTH, POINT_SIZE, DOTTED_LINE_WIDTH = 1.5, 4.0, 2.2

sys.dont_write_bytecode = True
_HERE = os.path.dirname(os.path.abspath(__file__))
_COMP = os.path.normpath(os.path.join(_HERE, '..', 'computation'))
if _COMP not in sys.path:
    sys.path.insert(0, _COMP)

from assemble_lnD_vgeneral_phi_v2 import fict_block_terms
from pipeline_helpers_vgeneral_phi_v2 import (bounce_fingerprint,
                                           bounce_name, load_background,
                                           load_stage_npz, read_band_tables,
                                           stage_paths)
from potential_vgeneral_phi_v2 import monomial_order

#: the band slice keys this module merges (delta_geq3 LAST: it is the check)
BAND_KEYS = ('Gbar_minus_fv', 'Gbar_1', 'Gbar_2', 'delta_geq3')

#: the pipeline's own pictures/ folder: every figure is saved below it
PICTURES = os.path.join(os.path.dirname(_HERE), 'pictures')


def bounce_label(args):
    """The BOUNCE NAME for figure filenames (e.g. 'M2_TV0'): the run's
    bounce npz (args.bounce_npz) read through pipeline_helpers.bounce_name
    (key pair_name = <falseLabel>_<trueLabel> from the vacuum catalogue).
    Fail-loud: no bounce npz, or one without pair_name, aborts.

    THE NAMING RULE this serves: a figure name is the PLOT CODE's name plus
    the bounce name and the (n, s^2) parameter tags (s2_filetag) -- enough to
    tell two parameter choices apart -- NEVER the full run tag or a date (the
    input npz lookup keeps the full --tag; only the picture drops it)."""
    p = getattr(args, 'bounce_npz', None)
    if not (p and os.path.isfile(p)):
        raise SystemExit(f'[ABORT] bounce npz {p} not found -- the figure '
                         f'name needs its pair_name.')
    with np.load(p, allow_pickle=True) as bounce_npz_data:
        return bounce_name(bounce_npz_data)


# --------------------------------------------------------------------------- #
#  Band slice tables, all keys                                                 #
# --------------------------------------------------------------------------- #
def band_tables_multikey(data_dir, tag, bounce_npz, n_max, s2_cap,
                         check_tol=1e-10):
    """Merge the residual-band slices of `tag` into dense (n_waves, n_s2)
    tables for every key of BAND_KEYS, on the certified band rectangle:
    waves up to n_max and slices up to s2_cap (the n_max and s2_max of
    band_integrals).

    The slice SET is whatever the production reader certifies: this calls
    pipeline_helpers_vgeneral_phi_v2.read_band_tables (bounce-fingerprint stamp + MANIFEST
    set equality) and then re-reads exactly the files it
    returned.  A wave that appears in two files for the same s2 must agree to
    `check_tol` relative (the same rule the production reader applies), and the
    identity

        delta_geq3 == Gbar_minus_fv - Gbar_1 - Gbar_2

    is re-verified on the delivered tables (it is how the engine defines the
    production object; a violation means a stale or foreign slice file).

    Returns dict with
        s2   (n_s2,)              ascending slice s2 values
        n    (n_waves,)           ascending waves (2..n_max)
        <key> (n_waves, n_s2)     one dense table per BAND_KEYS entry
        files                     the certified slice file list
        sector_traces_n0_n1       {s2: (Gbar_minus_fv_n0, Gbar_minus_fv_n1)}
                                  -- the raw n = 0, 1 waves of each BASE slice
    """
    manifest = os.path.join(data_dir, f'band_manifest_{tag}.json')
    if not os.path.isfile(manifest):
        raise RuntimeError(f'[ABORT] {os.path.basename(manifest)} not found '
                           f'-- without the MANIFEST guard a foreign same-tag '
                           f'slice file could be summed in silently; the band '
                           f'driver writes it.')
    # the run's own params/fields (from the bounce npz) for the potential_id
    # recompute inside read_band_tables
    with np.load(bounce_npz, allow_pickle=True) as _b:
        _params = np.asarray(_b['params'], float)
        _fields = tuple(str(f) for f in np.asarray(_b['fields']).ravel())
    s2_values, tbl_prod, files = read_band_tables(
        data_dir, tag, expect_bounce_fingerprint=bounce_fingerprint(bounce_npz),
        manifest_path=manifest,
        s2_cap=s2_cap, params=_params, fields=_fields)
    if len(files) == 0:
        raise RuntimeError(f'[ABORT] no residual_band_{tag}_s2*.npz in '
                           f'{data_dir}')

    merged = {k: {} for k in BAND_KEYS}          # merged[key][s2][n] = value
    # the raw n = 0, 1 waves the engine stores in every BASE slice
    # (Gbar_minus_fv_n0/_n1): collected in this same pass, one read
    sector_traces_n0_n1 = {}                       # sector_traces_n0_n1[s2] = (n0, n1)
    for f in files:
        with np.load(f, allow_pickle=True) as band_slice_npz:
            s2 = float(band_slice_npz['s2'])
            nv = np.asarray(band_slice_npz['n_values'], int)
            if 'Gbar_minus_fv_n0' in band_slice_npz.files and 'Gbar_minus_fv_n1' in band_slice_npz.files:
                sector_traces_n0_n1[s2] = (float(band_slice_npz['Gbar_minus_fv_n0']),
                                           float(band_slice_npz['Gbar_minus_fv_n1']))
            for key in BAND_KEYS:
                if key not in band_slice_npz.files:
                    raise RuntimeError(
                        f'[ABORT] slice {os.path.basename(f)} carries no '
                        f'"{key}" array -- not a current-engine slice; '
                        f'regenerate the band or plot a newer tag.')
                traces_at_s2 = merged[key].setdefault(s2, {})
                for n, v in zip(nv, np.asarray(band_slice_npz[key], float)):
                    n, v = int(n), float(v)
                    if n in traces_at_s2 and abs(v - traces_at_s2[n]) > check_tol * max(abs(v),
                                                                     1e-30):
                        raise RuntimeError(
                            f'[ABORT] wave n={n} at s2={s2:g} appears twice '
                            f'with conflicting {key} ({traces_at_s2[n]!r} vs {v!r}) -- '
                            f'a stale extension file; delete it.')
                    traces_at_s2[n] = v

    # the wave axis: the band 2..n_max, which every slice must hold
    n_common = sorted(set.intersection(*(set(merged['delta_geq3'][s2])
                                         for s2 in s2_values)))
    n_common = [n for n in n_common if n <= int(n_max)]
    band_n = np.array(n_common, int)
    if band_n.size == 0:
        raise RuntimeError('[ABORT] the slice files share no common wave.')
    if not np.array_equal(band_n, np.arange(2, int(n_max) + 1)):
        raise RuntimeError(f'[ABORT] the shared wave set {band_n[0]}..{band_n[-1]} is '
                           f'not the band 2..{int(n_max)} ({band_n.size} waves) -- a missing '
                           f'wave would bias every n-sum.')

    out = {'s2': np.asarray(s2_values, float), 'n': band_n, 'files': files,
           'sector_traces_n0_n1': sector_traces_n0_n1}
    for key in BAND_KEYS:
        M = np.array([[merged[key][s2][int(n)] for s2 in s2_values] for n in band_n])
        if not np.all(np.isfinite(M)):
            bad = int(np.sum(~np.isfinite(M)))
            raise RuntimeError(f'[ABORT] {bad} non-finite {key} values in the '
                               f'band tables (fail-loud, no nansum).')
        out[key] = M

    # identity check on the delivered tables (engine definition)
    dev = np.max(np.abs(out['delta_geq3']
                        - (out['Gbar_minus_fv'] - out['Gbar_1']
                           - out['Gbar_2'])))
    scale = max(float(np.max(np.abs(out['delta_geq3']))), 1e-30)
    if dev > check_tol * scale:
        raise RuntimeError(
            f'[ABORT] delta_geq3 != Gbar_minus_fv - Gbar_1 - Gbar_2 on the '
            f'delivered '
            f'tables (max dev {dev:.3e}, scale {scale:.3e}) -- stale slices.')
    # and against the production reader, wave by wave
    dev_prod = max(abs(tbl_prod[s2][int(n)] - out['delta_geq3'][i, j])
                   for i, n in enumerate(band_n) for j, s2 in enumerate(s2_values))
    if dev_prod != 0.0:
        raise RuntimeError(f'[ABORT] delta_geq3 table differs from the '
                           f'production reader by {dev_prod:.3e} (must be 0).')
    print(f'[plot_common] band tables: {len(files)} slice files, '
          f'{len(s2_values)} s2 values [{s2_values[0]:g}, {s2_values[-1]:g}], waves '
          f'n = {band_n[0]}..{band_n[-1]} ({band_n.size}); keys {", ".join(BAND_KEYS)}')
    return out


# --------------------------------------------------------------------------- #
#  The raw n = 0, 1 waves on the band s2 grid                                  #
# --------------------------------------------------------------------------- #
def trace_G_minus_G_fv_n0_n1_on_band_grid(tab):
    """The RAW (un-pole-subtracted) n = 0 and n = 1 waves
    Gbar_minus_fv_n on the band slice s2 grid -- what
    traces_vs_s2_n_common.add_trace_G_minus_G_fv_n0_n1 attaches so the
    core-equation figure sums ALL waves n = 0..n_max.

    ONE source: the keys Gbar_minus_fv_n0/_n1 the band engine stores in every
    BASE slice, already collected by band_tables_multikey (tab['sector_traces_n0_n1']).
    A slice set lacking them is not a current-engine band: abort and
    regenerate -- no side cache, no re-solve.

    Returns (Gbar_minus_fv_n0, Gbar_minus_fv_n1) aligned with tab['s2'].  A NaN entry means that slice
    sits exactly ON a sector pole (det C = 0); the caller guards it.
    """
    s2_values = np.asarray(tab['s2'], float)
    sector_traces_n0_n1 = tab['sector_traces_n0_n1']
    missing = [float(s2) for s2 in s2_values if float(s2) not in sector_traces_n0_n1]
    if missing:
        raise RuntimeError(
            f'[ABORT] {len(missing)} band slice(s) store no Gbar_minus_fv_n0/'
            f'_n1 keys (first at s2 = {missing[0]:g}) -- not a current-engine '
            f'band; regenerate the band with band_adaptive_vgeneral_phi_v2.')
    Gbar_minus_fv_n0_n1 = np.array([sector_traces_n0_n1[float(s2)] for s2 in s2_values])
    print(f'[plot_common] n=0,1 traces read from the slice files '
          f'({len(s2_values)} s2 values).')
    return Gbar_minus_fv_n0_n1[:, 0], Gbar_minus_fv_n0_n1[:, 1]


# --------------------------------------------------------------------------- #
#  Whole-run loader                                                            #
# --------------------------------------------------------------------------- #
def load_run(data_dir, tag, bounce_npz, need=('lnD_ren',)):
    """Load the production npz of a run into one dict.

    `need` names the stage_paths keys that MUST exist (missing -> abort);
    every other stage npz is loaded when present and left out when not.
    The bounce is loaded through the production reader load_background, so the
    returned 'bg' carries R, phi_b, fields, N, params, pot_prime, m2,
    m_min2, mu, U_int, tr_U_int3.  Every stage npz is read through
    pipeline_helpers.load_stage_npz: it must carry THIS bounce's fingerprint.

    Returns dict: bg + one entry per stage key + 'tag', 'data_dir',
    'bounce_npz', 'bounce_fingerprint' (the checksum of the bounce file; on
    disk the stage npz spell that key `bounce_sha`).
    """
    if not bounce_npz:
        raise RuntimeError('[ABORT] --bounce-npz (or $VGENERAL_PHI_BOUNCE_NPZ) '
                           'is required.')
    if not os.path.isfile(bounce_npz):
        raise RuntimeError(f'[ABORT] bounce npz not found: {bounce_npz}')
    paths = stage_paths(data_dir, tag)
    run = {'tag': tag, 'data_dir': data_dir, 'bounce_npz': bounce_npz,
           'bounce_fingerprint': bounce_fingerprint(bounce_npz),
           'paths': paths}
    for key in need:
        if not os.path.isfile(paths[key]):
            raise RuntimeError(f'[ABORT] required stage output missing: '
                               f'{paths[key]}')
    for key, path in paths.items():
        if not os.path.isfile(path):
            continue
        # the bounce npz is the reference itself; every stage npz must carry
        # its fingerprint (load_stage_npz aborts on a missing or wrong stamp)
        run[key] = (np.load(path, allow_pickle=True) if key == 'bounce'
                    else load_stage_npz(path, bounce_npz,
                                        'the stage that writes '
                                        + os.path.basename(path)))
    run['bg'] = load_background(bounce_npz)
    return run


# --------------------------------------------------------------------------- #
#  Reproduce the production scalars from the stored arrays                     #
# --------------------------------------------------------------------------- #
def lnD_from_parts(run):
    """Rebuild lnD_ren in the EXACT production term order of
    assemble_lnD_vgeneral_phi_v2 (its seven terms (T1)-(T7)):

        lnD_ren = Sigma_dang                 (T1)
                + G1_addback_band                         (T2) Born-order-1 add-back on [0, s2_max]
                + sum_{n<=n_max} I_n(s2_max)             (T3)
                + G1_addback_high_s2                         (T4) Born-order-1 add-back beyond s2_max
                + sum_{n<=n_max} T_n(s2_max)             (T5)
                + Zcal                                   (T6)
                - A2_fin/2                               (T7)

    from the stored stage npz.  Returns (lnD, parts_dict); parts['A1_fin'] is
    (T2)+(T4), the tadpole finite part as the add-back of the renormalized
    Green's function, and parts['band_sum'] is (T3)+(T5) = sum_n I_n(inf).
    Used by the plots BOTH to draw the curves and to certify -- at the
    production cutoffs the rebuilt number must equal the stored lnD_ren to the
    last digit.
    """
    d = run['lnD_ren']
    deg1 = float(d['deg1'])
    omega2_00 = float(d['omega2_00'])
    omega2_10 = float(d['omega2_10'])
    s2_max, m2fict = float(d['s2_max']), float(d['m2fict'])
    fict, fict_target = fict_block_terms(m2fict, s2_max,
                                         omega2_00, omega2_10, deg1)
    # ISFINITE FIRST: a NaN deviation is FALSE under every comparison, so the bare
    # `abs(...) > 1e-9` form PASSES a bookkeeping that produced NaN and the
    # figure is then drawn from it, with the NaN surfacing later (or never) and
    # no cause named.
    dev = abs(fict - fict_target)
    if not np.isfinite(dev) or dev > 1e-9:
        raise RuntimeError(f'[ABORT] the m2fict-collapse identity is off by '
                           f'{fict - fict_target:.2e} (limit 1e-9, and a '
                           f'non-finite deviation fails too).')
    Delta_0 = float(d['Delta_0'])
    Delta_1 = float(d['Delta_1'])
    # the negative-mode add-back, once for the one negative mode of this
    # index-1 bounce (assemble_lnD_vgeneral_phi_v2 says why the |.| is the
    # convention and not a dropped sign)
    ln_abs_omega2_00 = float(np.log(abs(omega2_00)))
    Sigma_dang = float(Delta_0 + Delta_1 + fict + ln_abs_omega2_00)
    band_n = np.asarray(d['band_n'], int)
    n_max = int(d['n_max'])
    in_band = band_n <= n_max
    band_sum_s2_max = float(np.asarray(d['I_n_s2_max'], float)[in_band].sum())     # (T3)
    T_sum = float(np.asarray(d['T_n_s2_max'], float)[in_band].sum())               # (T5)
    band_sum = float(np.asarray(d['I_n_inf'], float)[in_band].sum())               # (T3)+(T5)
    Zcal = float(d['Zcal'])                                                            # (T6)
    G1_addback_band = float(d['G1_addback_band'])                                       # (T2)
    G1_addback_beyond = float(d['G1_addback_high_s2'])                                     # (T4)
    A2 = float(d['A2_fin'])                                                            # (T7)
    lnD = (Sigma_dang + G1_addback_band + band_sum_s2_max + G1_addback_beyond
           + T_sum + Zcal - 0.5 * A2)
    parts = {'Delta_0': Delta_0, 'Delta_1': Delta_1,
             'fict_block': float(fict),
             'ln_abs_omega2_00': ln_abs_omega2_00,
             'Sigma_dang': Sigma_dang,
             'G1_addback_band': G1_addback_band, 'G1_addback_high_s2': G1_addback_beyond,
             'band_sum_s2_max': band_sum_s2_max, 'T_sum': T_sum,
             'band_sum': band_sum, 'Zcal': Zcal,
             'A1_fin': G1_addback_band + G1_addback_beyond,
             'minus_half_A2_fin': -0.5 * A2, 'lnD_ren': float(lnD)}
    return float(lnD), parts


def certify_lnD(run, atol=1e-9):
    """Rebuild lnD_ren from the stored parts and check it against the stored
    scalar.  Every figure that shows an absolute lnD_ren calls this first, so a
    figure can never be drawn from an inconsistent npz set."""
    lnD, parts = lnD_from_parts(run)
    stored = float(run['lnD_ren']['lnD_ren'])
    dev = lnD - stored
    if not np.isfinite(dev) or abs(dev) > atol:     # a NaN passes every comparison
        raise RuntimeError(f'[ABORT] rebuilt lnD_ren {lnD:.9f} != stored '
                           f'{stored:.9f} (dev {dev:.2e} > {atol:g}).')
    print(f'[plot_common] lnD_ren rebuilt from parts = {lnD:.9f} '
          f'(stored {stored:.9f}, dev {dev:+.1e})')
    return lnD, parts


def tex_number(x):
    """x for a mathtext label, 6 significant digits, with a power of ten
    written out (3.06177 \\times 10^{-12}) instead of the e-notation of %g."""
    text = f'{x:.6g}'
    if 'e' not in text:
        return text
    mantissa, exponent = text.split('e')
    return rf'{mantissa} \times 10^{{{int(exponent)}}}'


# --------------------------------------------------------------------------- #
#  Layout rules that every figure of this folder obeys (ONE owner each)        #
# --------------------------------------------------------------------------- #
def legend_outside(ax, where='right', fontsize=None, ncol=1, below=-0.17, **kw):
    """THE legend rule of this plot folder: a legend NEVER sits on top of the
    data.  It is placed OUTSIDE the axes rectangle, so no curve can ever be
    hidden by it, whatever the data does.

    where='right'  legend to the RIGHT of the axes, aligned with its top edge.
                   Use it for a single wide panel or for panels STACKED in one
                   column -- there is free canvas to the right of both.
    where='below'  legend BELOW the axes, centred under it.  Use it for panels
                   sitting SIDE BY SIDE in a row: a right-hand legend would be
                   squeezed between two panels, a legend below is not.  Its
                   top edge sits at `below` (axes fraction) under the axes:
                   a short panel needs more to clear its x label.

    The axes keep their full size and the legend hangs off the figure; the
    figure is written with bbox_inches='tight' (save(), below), which grows the
    PNG until the legend fits.  Nothing is clipped and nothing is covered.
    fontsize None = the one legend size (rcParams, set above).
    Returns the matplotlib Legend.
    """
    if where == 'right':
        loc, anchor = 'upper left', (1.02, 1.0)
    elif where == 'below':
        loc, anchor = 'upper center', (0.5, below)
    else:
        raise ValueError(f"legend_outside: where must be 'right' or 'below', "
                         f"got {where!r}")
    return ax.legend(loc=loc, bbox_to_anchor=anchor, borderaxespad=0.0,
                     fontsize=fontsize, ncol=ncol, **kw)


def add_standard_plot_cli(ap):
    """The uniform CLI of every plot script: the three production selectors
    plus --out-dir (default: the figure_folder of the run's potential and
    bounce, filled in by resolve_run).  --tag and --bounce-npz left empty are AUTO-DISCOVERED from
    the data dir by resolve_run (so a plot script runs BARE on a finished run
    dir)."""
    ap.add_argument('--data-dir',
                    default=os.environ.get('VGENERAL_PHI_DATA', '.'),
                    help='data directory holding the stage npz '
                         '(default $VGENERAL_PHI_DATA, else .)')
    ap.add_argument('--tag', default=None,
                    help='run tag threaded into every filename (default: '
                         'auto-discovered from the most recent lnD_ren_*.npz, '
                         'else bounce_data_*.npz, in --data-dir)')
    ap.add_argument('--bounce-npz',
                    default=os.environ.get('VGENERAL_PHI_BOUNCE_NPZ'),
                    help='bounce npz of the run (default '
                         '$VGENERAL_PHI_BOUNCE_NPZ, else auto-discovered as '
                         'bounce_data_<tag>.npz in --data-dir)')
    ap.add_argument('--out-dir', default=None,
                    help='figure folder (default: pictures/pictures_<potential>'
                         '_<bounce>/<group>/ inside the pipeline, see '
                         'figure_folder)')
    return ap


def resolve_run(args, group='computation'):
    """AUTO-DISCOVERY of the run selectors, so plot scripts run BARE with no
    flags and no configuration: they plot THE RUN YOU DID LAST.  Order per
    selector: explicit flag > env > the most recent run found under the data
    dir (searched one level deep, so a data ROOT holding many run directories
    works too); an explicit --bounce-npz without --tag names its own run.
    Every resolution is printed; a different run is selected with
    --data-dir or --tag.  Fills args.tag / args.bounce_npz IN PLACE, and
    args.out_dir, when --out-dir is not given, with figure_folder(group, the
    run's bounce npz): group is 'computation' for the codes of this folder
    and 'bounce_gallery' for the bounce gallery."""
    import glob
    data_dir = args.data_dir
    if args.tag is None and args.bounce_npz:
        # an explicit bounce names its run: the tag of bounce_data_<tag>.npz
        # (the most recent run of the folder could be another one)
        name = os.path.basename(args.bounce_npz)
        if not (name.startswith('bounce_data_') and name.endswith('.npz')):
            raise SystemExit(f'[ABORT] --bounce-npz {name} is not named '
                             f'bounce_data_<tag>.npz -- pass --tag as well.')
        args.tag = name[len('bounce_data_'):-4]
        print(f'[plots] tag:   {args.tag}  (from --bounce-npz)')
    if args.tag is None:
        # Only npz produced by THIS pipeline qualify: every v2 output carries
        # the vgeneral_phi_v2 tag in its name (the tag audit at the end of
        # run_vgeneral_phi_v2_production.sh checks exactly that),
        # so files from other pipelines can never be picked up.  Search the
        # data dir itself, then its subdirectories (the data-root case).
        for kind in ('lnD_ren', 'bounce_data'):
            cands = [c for pat in (os.path.join(data_dir, f'{kind}_*.npz'),
                                   os.path.join(data_dir, '*', f'{kind}_*.npz'))
                     for c in glob.glob(pat)
                     if 'vgeneral_phi_v2' in os.path.basename(c)]
            if cands:
                break
        if not cands:
            raise SystemExit(
                f'[ABORT] no vgeneral_phi_v2 run found in {data_dir} or its '
                f'subdirectories -- run the pipeline first, or pass '
                f'--data-dir <run dir>.')
        # THE RUN YOU DID LAST: newest by file time.  Printed, so the choice
        # is never silent; --data-dir/--tag override it.
        newest = max(cands, key=os.path.getmtime)
        args.data_dir = data_dir = os.path.dirname(newest)
        args.tag = os.path.basename(newest)[len(kind) + 1:-4]
        note = '' if len(cands) == 1 else f'  (most recent of {len(cands)})'
        print(f'[plots] using: {newest}{note}')
        print(f'[plots] tag:   {args.tag}')
    if not args.bounce_npz:
        path = os.path.join(data_dir, f'bounce_data_{args.tag}.npz')
        if not os.path.isfile(path):
            cands = sorted(c for c in
                           glob.glob(os.path.join(data_dir,
                                                  'bounce_data_*.npz'))
                           if 'vgeneral_phi_v2' in os.path.basename(c))
            if len(cands) == 1:
                path = cands[0]
            elif cands:
                listing = '\n'.join('    ' + os.path.basename(c)
                                    for c in cands)
                raise SystemExit(
                    f'[ABORT] ambiguous bounce npz (bounce_data_{args.tag}'
                    f'.npz not found, {len(cands)} candidates in {data_dir}) '
                    f'-- pass --bounce-npz.  Candidates:\n{listing}')
            else:
                raise SystemExit(
                    f'[ABORT] no bounce_data_*.npz in {data_dir} (looked for '
                    f'bounce_data_{args.tag}.npz) -- pass --bounce-npz.')
        args.bounce_npz = path
        print(f'[plots] using: {args.bounce_npz}')
    args.out_dir = args.out_dir or figure_folder(group, args.bounce_npz)
    return args


def figure_folder(group, bounce_npz=None):
    """THE folder a figure is saved in -- one rule for every plot code:

        pictures/pictures_<potential>_<bounce>/<group>/   a figure of one run
        pictures/<group>/                                 bounce_npz=None: the
                                                          thin-wall test bed
                                                          (its own one-field
                                                          model, no run)

    SIDE EFFECT, when bounce_npz is given: every call creates the folder
    pictures/pictures_<potential>_<bounce>/ and rewrites the note
    potential_and_bounce.txt in it, which names the fields, the bounce and
    every non-zero coupling (6 significant digits).  With bounce_npz=None
    it only returns the path.

    pictures/ is the pipeline's own folder (PICTURES).  <group> is the plot
    family: the name of the plot code's folder without 'plot_'
    (computation, add_on_greensfct, computation_aside, test_greensequation,
    test_thin_wall), and bounce_gallery for the bounce gallery.
    <potential> and <bounce> are read from the run's bounce npz; <bounce> is
    the bounce name (e.g. M2_TV0).  <potential> is the
    fields and every non-zero coupling written in front of its monomial, as
    in the polynomial, c the constant term (x-y_1.2c_-2x2_... for
    V = 1.2 - 2 x^2 + ...), each coupling to 6 significant digits: two
    potentials that differ only beyond the 6th digit get the same label and
    share a folder.  A label longer than 60 characters is replaced by
    <N>fields_<n>couplings_<8 hex digits of the sha256 of the exact
    couplings>, which does tell such potentials apart.

    One folder per potential and bounce, not per run: two runs of the same
    potential and bounce (e.g. a rerun of one model with other cutoffs)
    share one folder, and a re-plot of one overwrites the other's figures
    of the same name; --out-dir keeps them apart.  Every plot code with
    --out-dir calls this as  args.out_dir or figure_folder(...):  the
    given --out-dir, else this folder."""
    if bounce_npz is None:
        return os.path.join(PICTURES, group)
    with np.load(bounce_npz, allow_pickle=True) as bounce_npz_data:
        fields = tuple(str(f) for f in np.asarray(bounce_npz_data['fields']).ravel())
        params = np.asarray(bounce_npz_data['params'], float)
        bounce_pair_name = bounce_name(bounce_npz_data)
    # the non-zero couplings (name, value) in the pipeline's monomial order
    couplings = [(name, value) for name, value in zip(monomial_order(fields), params)
                 if value != 0.0]
    potential = '-'.join(fields) + '_' + '_'.join(f'{value:g}{name}'
                                                  for name, value in couplings)
    if len(potential) > 60:
        digest = hashlib.sha256(params.tobytes()).hexdigest()[:8]
        potential = f'{len(fields)}fields_{len(couplings)}couplings_{digest}'
    run_folder = os.path.join(PICTURES, f'pictures_{potential}_{bounce_pair_name}')
    os.makedirs(run_folder, exist_ok=True)
    with open(os.path.join(run_folder, 'potential_and_bounce.txt'), 'w') as note:
        note.write('\n'.join(['fields:  ' + ', '.join(fields),
                              'bounce:  ' + bounce_pair_name,
                              'potential (non-zero couplings):']
                             + [f'    {name} = {value:g}' for name, value in couplings])
                   + '\n')
    return os.path.join(run_folder, group)


def thick_titles(fig):
    """Write every title of the figure (the figure title and each panel's
    title) thicker: a thin outline in the title's own colour, which thickens
    the words and the formulas alike (a bold font would leave every formula
    of a title thin)."""
    for title in [ax.title for ax in fig.axes] + [getattr(fig, '_suptitle', None)]:
        if title is not None and title.get_text():
            title.set_path_effects([patheffects.withStroke(
                linewidth=0.6, foreground=title.get_color())])


def save(fig, out_dir, name, dpi=150):
    """Write one figure into out_dir and print its path -- the one writer of
    a figure; every title is written thicker first (thick_titles).

    bbox_inches='tight' is NOT cosmetic here: legend_outside() deliberately
    places the legend outside the axes, and this is what grows the PNG until
    it fits instead of clipping it.
    """
    thick_titles(fig)
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, name)
    fig.savefig(path, dpi=dpi, bbox_inches='tight')
    print(f'[OK] wrote {path}')
    return path
