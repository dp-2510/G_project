#!/usr/bin/env python3
"""plot_hessian_components_vgeneral_phi_v2.py -- ONE picture: ALL components of the primed
Hessian along the bounce, with the wall of every bounce component marked.

WHAT THE PICTURE SHOWS
----------------------
The FULL Hessian with the bounce inserted is, in the primed frame,

    [H_b(r)]_ij = [H'(phi'_b(r))]_ij ,  i, j = 1..N ,  symmetric,

computed by computation_aside/hessian_components_wall_vgeneral_phi_v2.py.  ALL
components are drawn -- diagonal and off-diagonal together, which is what
makes this figure the complete statement about the Hessian along the bounce:

  * the DIAGONAL components are the curvature along each primed axis.  Each
    one must relax onto its own false-vacuum channel mass m_i^2 far outside
    the wall, since H_b(inf) = Lambda = diag(m_1^2, ..., m_N^2).

  * the OFF-DIAGONAL components are the field MIXING.  They vanish in the
    false vacuum and in between are the eigenvalue splitting times the
    rotation ANGLE of the eigenframe: for two fields
    b(r) = (Delta lambda / 2) sin(2 mixing_angle_1), so b is large where a
    large splitting meets a rotated eigenframe; the splitting peaks where the
    bounce passes near the saddle point.

The false-vacuum value of each component (m_i^2 on the diagonal, zero off it,
legend [H'(phi^v)]_ij = delta_ij m_i^2) is drawn as a dotted horizontal line in
the curve's own colour, so the relaxation can be read off directly.

THE WALL POSITIONS are the point of the figure: one vertical line per bounce
COMPONENT, at the radius where that component turns over fastest,
r_wall,i = argmax_r |d phi'_i / dr|, labelled 'steepest slope of phi'_i'
(no number; no line for a component that never moves: its stored r_wall,i
is NaN).

NOTATION of the axis and the legend: the y axis holds the
entries of the hessian matrix [H'(phibar^b(r))]_ij, the legend names each
entry [H'(phi^b_1(r), .., phi^b_N(r))]_ij.  No title: the bounce is named in
the file name.
The claim they test is that the off-diagonal Hessian peaks AT the wall -- the
mixing is a wall effect -- and keeping one line per component shows a field
that switches on ahead of or behind the others instead of averaging it away.

LAYOUT, and why: for TWO fields the three independent components are
a(r), b(r), c(r) of H_b = [[a, b], [b, c]], and they belong on ONE axes,
where they are directly
comparable (that comparison is the whole content of the two-field story); ONE
field has a single component, drawn on one axes the same way.  For three
fields or more there are N(N+1)/2 components of very different
size, so each gets its own panel in a triangular grid; H_b is symmetric,
lower triangle would only repeat the upper one and is left empty.

reads:  hessian_components_wall_<tag>.npz -- R, H_b,
        r_W_per_field, r_W, m2, N
OUTPUT  <out-dir>/plot_hessian_components_vgeneral_phi_v2_<bounce>.png
"""
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

sys.dont_write_bytecode = True
_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

from plot_aside_common_vgeneral_phi_v2 import (bounce_label, legend_outside,
                               parse_cli, read_aside_npz, save)
from plot_common_vgeneral_phi_v2 import DOTTED_LINE_WIDTH

#: colours of the per-component wall lines (one per bounce component)
WALL_COLORS = ('crimson', 'darkgreen', 'purple', 'saddlebrown', 'teal')

#: the legend entry of the dotted false-vacuum lines
FV_VALUE_LABEL = r"$[H'(\phi^v)]_{ij} = \delta_{ij}\, m_i^2$"


def independent_components_with_fv(m2):
    """The N(N+1)/2 independent entries of the symmetric Hessian: for each,
    the index pair, the component label (H' with the two primed channel
    numbers, counted from 1, on the bounce components phi^b_1(r) ..
    phi^b_N(r)), and the false-vacuum value the component relaxes onto
    (m_i^2 on the diagonal, zero off it).
    The (i, j) order, i <= j row by row, is the same as in the component
    table printed by computation_aside/hessian_components_wall_vgeneral_phi_v2.py."""
    N = len(m2)
    bounce = ', '.join(rf'\phi^b_{{{k + 1}}}(r)' for k in range(N))
    # channel numbers = indices + 1 (phi'_1 is index 0)
    return [(i, j, rf"$[H'({bounce})]_{{{i + 1}{j + 1}}}$",
             float(m2[i]) if i == j else 0.0)
            for i in range(N) for j in range(i, N)]


def draw_component(ax, R, component, fv_value, label):
    """One component curve plus the dotted line marking the false-vacuum value
    it must relax onto far outside the wall."""
    line, = ax.plot(R, component, label=label)
    ax.axhline(fv_value, color=line.get_color(), ls=':', lw=DOTTED_LINE_WIDTH)
    ax.grid(True)
    return line


def draw_walls(ax, r_W_per_field, label_them):
    """The wall of every bounce component, as a vertical line per component.
    Labelled only where the figure's single legend is built, so the same lines
    can be repeated in every panel without repeating the legend."""
    for i, r_w in enumerate(r_W_per_field):
        if np.isnan(r_w):            # this component never moves: no wall
            continue
        ax.axvline(r_w, color=WALL_COLORS[i % len(WALL_COLORS)], ls='--',
                   label=(rf"steepest slope of $\phi'_{{{i + 1}}}$"
                          if label_them else '_nolegend_'))


def figure_two_fields(R, H_b, components, r_W_per_field):
    """One or two fields: every independent component in ONE axes -- a, b, c
    of H_b = [[a, b], [b, c]] for two fields, the single entry for
    one field."""
    # a picture <= 8 in wide with its legend (sizes: plot_common)
    fig, ax = plt.subplots(figsize=(4.8, 3.8))
    for i, j, label, fv_value in components:
        draw_component(ax, R, H_b[:, i, j], fv_value, label)
    ax.plot([], [], color='0.4', ls=':', lw=DOTTED_LINE_WIDTH, label=FV_VALUE_LABEL)
    draw_walls(ax, r_W_per_field, label_them=True)
    ax.set_xlabel('$r$')
    ax.set_ylabel(r"Entries of the hessian $[H'(\overline{\phi}^{\,b}(r))]_{ij}$")
    legend_outside(ax)
    return fig


def figure_component_grid(R, H_b, components, r_W_per_field, N):
    """Three fields or more: one panel per
    independent component, placed at its own (i, j) slot of the matrix, so the
    grid is read like H_b itself; the lower triangle stays empty because
    H_b is symmetric.  Every panel carries its own axes and wall lines
    -- components differ by orders of magnitude, which is why they are not on
    one axes -- and ONE legend beside the grid labels the wall lines for all
    of them (the panel titles name the components)."""
    fig, axes = plt.subplots(N, N, figsize=(3.8 * N, 2.8 * N), squeeze=False)
    for ax in axes.ravel():
        ax.set_visible(False)
    for i, j, label, fv_value in components:
        ax = axes[i, j]
        ax.set_visible(True)
        first_panel = (i == 0 and j == 0)
        draw_component(ax, R, H_b[:, i, j], fv_value, '_nolegend_')
        draw_walls(ax, r_W_per_field, label_them=first_panel)
        ax.set_title(label)
        ax.set_xlabel('$r$')
        if i == j:
            ax.set_ylabel(r"Entry of the hessian $[H'(\overline{\phi}^{\,b}(r))]_{ij}$")
        if first_panel:
            ax.plot([], [], color='0.4', ls=':', lw=DOTTED_LINE_WIDTH, label=FV_VALUE_LABEL)
    fig.suptitle("The hessian is symmetric: the empty lower triangle repeats "
                 "the upper one")
    fig.tight_layout()
    # one legend for the whole grid, placed beside it so no panel is covered
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper left', bbox_to_anchor=(1.0, 0.98),
               borderaxespad=0.0, fontsize=9)
    return fig


def main():
    args = parse_cli('aside figure: all components of H_b, the full '
                    'Hessian with the bounce inserted, with the walls marked')
    data = read_aside_npz(args, 'hessian_components_wall')
    R = data['R']
    H_b = data['H_b']
    r_W_per_field = np.asarray(data['r_W_per_field'], float)
    N = int(data['N'])
    components = independent_components_with_fv(data['m2'])

    # the physics reading: where every component does its work, against the
    # wall of every bounce component
    for i, r_w in enumerate(r_W_per_field):
        wall = 'none (the bounce does not move it)' if np.isnan(r_w) else f'r = {r_w:.4f}'
        print(f"[plot] wall of component phi'_{i + 1}: {wall}")
    print(f'[plot] full profile (norm of the slope vector): r_W = '
          f'{float(data["r_W"]):.4f}')
    for i, j, _label, fv_value in components:
        component = H_b[:, i, j]
        i_max_departure = int(np.argmax(np.abs(component - fv_value)))
        print(f"[plot] H_b_{i + 1}{j + 1}: false-vacuum value "
              f"{fv_value:+.6f}, largest departure {component[i_max_departure]:+.6f} at "
              f"r = {R[i_max_departure]:.4f}")

    bounce_pair_name = bounce_label(args)
    if N <= 2:
        fig = figure_two_fields(R, H_b, components, r_W_per_field)
    else:
        fig = figure_component_grid(R, H_b, components, r_W_per_field, N)
    save(fig, args.out_dir, f'plot_hessian_components_vgeneral_phi_v2_{bounce_pair_name}.png')


if __name__ == '__main__':
    main()
