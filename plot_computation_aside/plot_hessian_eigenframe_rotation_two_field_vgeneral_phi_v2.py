#!/usr/bin/env python3
"""plot_hessian_eigenframe_rotation_two_field_vgeneral_phi_v2.py -- ONE picture,
two panels: the rotation of the eigenframe of the primed Hessian along the
bounce of a TWO-FIELD run, the angle theta'(r) that diagonalizes H'(phi'_b(r)).

WHAT THE PICTURE SHOWS
----------------------
At the false vacuum the primed frame diagonalizes the Hessian (that is how
it was built: the rotation angle theta^v = 1/2 atan2(2 H_12, H_11 - H_22) of
the original-frame Hessian).  Along the bounce the primed Hessian
H'(phi'_b(r)) is no longer diagonal, and the same formula gives the extra
rotation that diagonalizes it at radius r
(plot_aside_common.hessian_rotation_two_field):

    theta'(r) = 1/2 atan2( 2 H'_12(r), H'_11(r) - H'_22(r) ),
    e'^(1)(r) = (cos theta', sin theta'),   e'^(2)(r) = (-sin theta', cos theta'),

the local principal curvature axes of the landscape at the bounce point
(the primed channels numbered with phi'_1 the HEAVIEST).  As with theta^v,
e'^(1) belongs to the larger eigenvalue, so far outside the wall, where
H' -> diag(m'_1^2, m'_2^2), theta' -> 0 and the eigenvectors lie on the
primed axes.

  TOP     the two eigenvectors as unit arrows in the (phi'_1, phi'_2) plane,
          one pair at each of N_ARROWS radii evenly spaced in r, in the box
          (r, phi'_1, phi'_2) -- with the phi'_1 axis pointing toward the
          viewer, so that the arrows of the two eigenvectors do not cross;
  BOTTOM  theta'(r) in degrees, in (-90, 90] and continuous as long as the
          two eigenvalues do not cross, with the radius of the steepest
          bounce slope it is compared against.

Across the wall the arrows turn and theta' leaves 0 -- that rotation is the
field mixing as long as the two eigenvalues stay apart.  Where they cross
(e.g. a potential without mixing, H'_12 = 0), theta' jumps by 90 degrees
and e'^(1) swaps axes: a relabeling of the larger eigenvalue, not mixing.
A run with N != 2 aborts here.

reads:  hessian_eigenproblem_<tag>.npz -- R, H_b, eigenvalues, r_W, N
OUTPUT  <out-dir>/plot_hessian_eigenframe_rotation_two_field_vgeneral_phi_v2_<bounce>.png
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

from plot_aside_common_vgeneral_phi_v2 import (bounce_label,
                                               hessian_rotation_two_field,
                                               legend_outside, parse_cli,
                                               read_aside_npz, save)

DEGREE = 180.0 / np.pi

#: arrows drawn per eigenvector, at radii spaced evenly in r over the bounce grid
N_ARROWS = 36


def draw_eigenvector_arrows(ax, R, eigenvectors):
    """TOP panel, a 3-D axes: at N_ARROWS radii one unit arrow per
    eigenvector, rooted on the r axis.  The radii are evenly spaced in r (the
    bounce grid itself is dense at the wall, so grid indices would bunch the
    arrows there): for each, the first grid point at or beyond it
    (np.searchsorted).  e'^(1) is blue and e'^(2) orange."""
    sample = np.unique(np.searchsorted(R, np.linspace(R[0], R[-1], N_ARROWS)))
    sample = np.minimum(sample, len(R) - 1)
    r = R[sample]
    labels = (r"$e'^{(1)}(r) = (\cos\theta', \sin\theta')$",
              r"$e'^{(2)}(r) = (-\sin\theta', \cos\theta')$")
    for alpha in range(2):
        color = f'C{alpha}'
        ax.quiver(r, np.zeros_like(r), np.zeros_like(r), np.zeros_like(r),
                  eigenvectors[sample, 0, alpha], eigenvectors[sample, 1, alpha],
                  color=color, arrow_length_ratio=0.12)
        ax.plot([], [], color=color, label=labels[alpha])
    ax.set_xlabel('$r$')
    ax.set_ylabel(r"$\phi'_1$")
    ax.set_zlabel(r"$\phi'_2$")
    ax.set_ylim(-1.0, 1.0)
    ax.set_zlim(-1.0, 1.0)
    # the phi'_1 axis points toward the viewer: the e'^(1) arrows come out of
    # the box and the e'^(2) arrows point up, so the two never cross (the
    # default view, turned by 180 degrees with r turned back to grow to the
    # right as in the lower panel)
    ax.invert_yaxis()
    legend_outside(ax)


def draw_rotation_angle(ax, R, theta, r_W):
    """BOTTOM panel: theta'(r) in degrees and the radius of the steepest
    bounce slope."""
    ax.plot(R, theta * DEGREE,
            label=r"$\theta'(r) = \frac{1}{2}\,\mathrm{atan2}(2 [H']_{12},\,"
                  r" [H']_{11} - [H']_{22})$,"
                  "\n" r"$H' = H'(\overline{\phi}^{\,b}(r))$")
    ax.axvline(r_W, color='crimson', ls='--',
               label=rf'$r_w = {r_W:.2f}$')
    ax.set_xlabel('$r$')
    ax.set_ylabel(r"Rotation angle $\theta'(r)$ (degree)")
    ax.set_ylim(-95.0, 95.0)
    ax.grid(True)
    legend_outside(ax)


def main():
    args = parse_cli('aside figure: the rotation of the primed Hessian '
                     'eigenframe along the bounce: the '
                     'eigenvectors as arrows and the rotation angle (N = 2)')
    data = read_aside_npz(args, 'hessian_eigenproblem')
    if int(data['N']) != 2:
        raise SystemExit(
            f'[ABORT] the eigenframe rotation angle exists for N = 2 only '
            f'(this run has N = {int(data["N"])}).')
    R = data['R']
    theta, eigenvectors, _eigenvalues = hessian_rotation_two_field(
        data['H_b'], data['eigenvalues'])
    r_W = float(data['r_W'])
    bounce_pair_name = bounce_label(args)

    # the physics reading: where the eigenvectors start and end, the largest
    # rotation against the wall; a jump of theta' is reported: by 180 degrees
    # a wrap at +-90 degrees (the same axis), by 90 degrees a crossing of the
    # two eigenvalues (a relabeling of the larger one, not mixing)
    for alpha in range(2):
        print(f"[plot] e'^({alpha + 1}) at r = 0 "
              f"({eigenvectors[0, 0, alpha]:+.4f}, {eigenvectors[0, 1, alpha]:+.4f}), "
              f"at r_max ({eigenvectors[-1, 0, alpha]:+.4f}, {eigenvectors[-1, 1, alpha]:+.4f})")
    i_peak = int(np.argmax(np.abs(theta)))
    print(f"[plot] theta' largest at r = {R[i_peak]:.4f}: {theta[i_peak] * DEGREE:+.3f} deg;"
          f"  far value {theta[-1] * DEGREE:+.3f} deg;  steepest bounce slope at r_W = {r_W:.4f}")
    for i in np.where(np.abs(np.diff(theta)) > np.pi / 4)[0]:
        jump = abs(theta[i + 1] - theta[i]) * DEGREE
        print(f"[plot] theta' jumps by {jump:.0f} deg at r = {R[i]:.4f}: "
              f"by 180 a wrap at +-90 deg (the same axis), by 90 a "
              f"crossing of the two eigenvalues (e'^(1) swaps axes: a "
              f"relabeling of the larger eigenvalue, not mixing)")

    fig = plt.figure(figsize=(7.7, 9.0))     # a picture <= 8 in wide
    panels = fig.add_gridspec(2, 1, height_ratios=(2.0, 1.0), hspace=0.05)
    draw_eigenvector_arrows(fig.add_subplot(panels[0], projection='3d'),
                            R, eigenvectors)
    # the lower panel narrower than the box above (60 % of the width), centred
    lower = panels[1].subgridspec(1, 3, width_ratios=(0.2, 0.6, 0.2))
    draw_rotation_angle(fig.add_subplot(lower[1]), R, theta, r_W)
    save(fig, args.out_dir,
         f'plot_hessian_eigenframe_rotation_two_field_vgeneral_phi_v2_{bounce_pair_name}.png')


if __name__ == '__main__':
    main()
