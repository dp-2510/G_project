#!/usr/bin/env python3
"""gbar_sector_common_vgeneral_phi_v2.py -- the machinery the two sector
pictures share: the unrenormalized objects of one sector npz (with the
production identities re-verified on the stored arrays) and the renderer that
draws them.  Not a plot code itself: plot_Gbar_n0 and plot_Gbar_n1 each own a
main() and one picture, and import these two from here.

Everything below is the SAME for n = 0 and n = 1; only the sector npz and the
picture name differ.
"""
import numpy as np

from plot_common_vgeneral_phi_v2 import POINT_SIZE, legend_outside


# one colour per role + neutral ink
COLOR_GBAR_BOUNCE, COLOR_P_N, COLOR_GBAR_BOUNCE_MINUS_P_N = '#2a78d6', '#eb6834', '#1baf7a'
COLOR_POLE_MARK = '#52514e'
#: the two curves without data points (the dangerous-mode contribution and
#: the difference) are drawn a bit thicker than the data curve, so they stay
#: readable on the grid
THICK_LINE_WIDTH = 2.2



def unrenormalized_sector_facts(sector_npz, n):
    """The unrenormalized objects of one sector npz, with the production
    identities re-verified on the stored arrays (fail-loud)."""
    s2 = np.asarray(sector_npz['s2_sector'], float)
    bounce = (np.asarray(sector_npz['Gbar_minus_fv'], float)
              + np.asarray(sector_npz['Gbar_fv'], float))
    deg = float(sector_npz['deg'])
    omega2 = float(sector_npz['omega2_00']) if n == 0 else float(sector_npz['omega2_10'])
    pole = deg / (omega2 + s2)
    sub = bounce - pole
    for name, a in (('Gbar_bounce', bounce), ('P_n', pole)):
        if not np.all(np.isfinite(a)):
            raise RuntimeError(f'[ABORT] non-finite {name} in sector n={n} '
                               f'(fail-loud, no nansum).')
    # the stored unrenormalized keys must belong to the certified
    # production wave ...
    dev = np.max(np.abs(np.asarray(sector_npz['delta_geq3'], float)
                        - (np.asarray(sector_npz['Gbar_minus_fv'], float)
                           - np.asarray(sector_npz['Gbar_1'], float)
                           - np.asarray(sector_npz['Gbar_2'], float))))
    if dev > 1e-10 * max(float(np.max(np.abs(sector_npz['delta_geq3']))), 1.0):
        raise RuntimeError(f'[ABORT] sector n={n}: delta_geq3 != Gbar_minus_fv '
                           f'- Gbar_1 - Gbar_2 (max dev {dev:.3e}) -- the '
                           f'stored unrenormalized keys do not belong to the '
                           f'production wave.')
    # ... and P_n here must be EXACTLY the pole term production subtracted
    dev = np.max(np.abs(np.asarray(sector_npz['delta_sub'], float)
                        - (np.asarray(sector_npz['delta_geq3'], float) - pole)))
    if dev > 1e-12 * max(float(np.max(np.abs(sector_npz['delta_sub']))), 1.0):
        raise RuntimeError(f'[ABORT] sector n={n}: delta_sub != delta_geq3 - '
                           f'deg/(omega2+s2) (max dev {dev:.3e}).')
    return dict(n=n, s2=s2, Gbar_bounce=bounce, P_n=pole,
                Gbar_bounce_minus_P_n=sub, deg=deg,
                omega2=omega2, s2_c=-omega2,    # the det C = 0 crossing
                nodes=int(len(s2)),
                resid=float(sector_npz['crossing_resid']),
                s2_max=float(sector_npz['s2_max']),
                mode=('negative mode' if n == 0 else 'zero mode'),
                # the numbering omega^2_{n,r} of the labels: wave n first
                sym=(r'\omega^2_{0,0}' if n == 0 else r'\omega^2_{1,0}'))


def _gap_at_pole(s2, y, s2_c):
    """s2 and y with a NaN inserted at the pole s2_c, so that no line joins
    the two sides of the pole (+inf on one side, -inf on the other)."""
    k = int(np.searchsorted(s2, s2_c))
    return np.insert(s2, k, s2_c), np.insert(np.asarray(y, float), k, np.nan)


def draw_sector_unren_curves(ax, sector_facts):
    """The three unrenormalized curves of one sector against s^2 (log axis),
    in the notation of the legend: Gbar_{n}(s^2) with its data points (the s^2
    nodes), the dangerous-mode contribution Pbar_{n}(s^2) = deg/(omega^2 + s^2)
    and their difference.

    The y axis is LINEAR for both sectors, with the signed values, so the
    subtraction reads directly off the picture (the gap between Gbar_n and
    Pbar_n IS the difference), and it is cut at 4 times the largest
    |difference| -- the parts beyond leave the picture:
      * n = 0: the pole lies ON the plotted range (the negative mode, at
        s^2 = -omega^2_{0,0} > 0); Gbar_n and Pbar_n change sign across it,
        so their negative parts lie below zero and the pole is a jump from
        -inf to +inf (no line joins its two sides; the axis is cut
        symmetrically);
      * n = 1: the pole lies at s^2 ~ 0, below the grid floor (the zero
        modes; nothing is drawn for it); every value is positive, so the
        axis starts just below zero and the ~4/s^2 rise of Gbar_n and Pbar_n
        towards small s^2 leaves it at the top.
    No title: the plot code sets it."""
    n = sector_facts['n']
    s2, s2_c = sector_facts['s2'], sector_facts['s2_c']
    G_n = r'\bar G_{n=%d}(s^2)' % n
    P_n = r'\bar P_{n=%d}(s^2)' % n
    # the difference: the wave n without its dangerous-mode contribution
    G_n_projected = r'\bar G^{\prime}_{n=%d}(s^2)' % n
    # Pbar_n is a formula, not data: drawn on a fine log grid, so the sparse
    # stretches of the sector grid leave no corners in it
    s2_fine = np.geomspace(s2[0], s2[-1], 2000)
    P_fine = sector_facts['deg'] / (sector_facts['omega2'] + s2_fine)
    # (s2, values, colour, label, line style, crosses the pole?)  Gbar_n with
    # its data points: edgeless, on TOP of the other two curves where they run
    # together; points and lines in the one set of sizes (plot_common)
    curves = (
        (s2, sector_facts['Gbar_bounce'], COLOR_GBAR_BOUNCE, rf'${G_n}$',
         dict(marker='o', ms=POINT_SIZE, mew=0, zorder=4), True),
        (s2_fine, P_fine, COLOR_P_N,
         rf'${P_n} = {sector_facts["deg"]:g}/({sector_facts["sym"]} + s^2)$',
         dict(lw=THICK_LINE_WIDTH, zorder=3), True),
        (s2, sector_facts['Gbar_bounce_minus_P_n'], COLOR_GBAR_BOUNCE_MINUS_P_N,
         rf'${G_n_projected} = {G_n} - {P_n}$',
         dict(lw=THICK_LINE_WIDTH, zorder=3), False))
    pole_on_axis = s2_c > float(s2[0])
    for x, y, color, label, style, crosses_pole in curves:
        x_plot, y_plot = (_gap_at_pole(x, y, s2_c) if crosses_pole and pole_on_axis
                          else (x, y))
        ax.plot(x_plot, y_plot, color=color, label=label, **style)
    if pole_on_axis:
        ax.axvline(s2_c, color=COLOR_POLE_MARK, ls='-.', zorder=2,
                   label=(f'{sector_facts["mode"]} pole at $s^2_z = '
                          f'-{sector_facts["sym"]}$'))
    # cut at 4 times the subtracted curve: the pole (n = 0) and the ~4/s^2
    # rise (n = 1) are infinite; symmetric where values are negative, else
    # from just below zero
    lim = 4.0 * float(np.max(np.abs(sector_facts['Gbar_bounce_minus_P_n'])))
    negative = any(float(np.min(y)) < 0.0 for _, y, *_ in curves)
    ax.set_ylim(-lim if negative else -0.05 * lim, lim)
    ax.set_xscale('log')
    ax.grid(True, which='both')     # the one grid style of plot_common
    ax.set_xlabel(r'$s^2$')
    # the legend BELOW the axes: it costs one short band, the axes keep the page
    legend_outside(ax, 'below')
