#!/usr/bin/env python3
"""plot_hessian_eigenvalues_vgeneral_phi_v2.py -- ONE picture: the eigenvalues of the primed
Hessian along the bounce.

WHAT THE PICTURE SHOWS
----------------------
The primed frame is the false-vacuum eigenframe, phi = phi_F + L phi', in
which the Hessian reads H'(phi') = L^T H(phi_F + L phi') L (L orthogonal;
PrimedPotential.H of computation/potential_vgeneral_phi_v2.py).  Along
the bounce profile phi'_b(r) the symmetric eigenproblem

    H'(phi'_b(r)) V^alpha(r) = lambda^alpha(r) V^alpha(r)

was solved at every radius by computation_aside/hessian_eigenproblem_vgeneral_phi_v2.py.
This figure draws its eigenvalues: the eigenvalue curves lambda^alpha(r).

The two physics statements the curves carry:

  * FAR OUTSIDE THE WALL the bounce has relaxed into the false vacuum,
    phi'_b -> 0, so H' -> Lambda = diag(m_1^2, ..., m_N^2) and every
    eigenvalue curve settles on its own channel mass^2.  Those false-vacuum
    values are drawn as horizontal reference lines; each curve meeting its
    own line IS that statement, made visible.

  * NEAR THE WALL the lowest eigenvalue curve dips and goes negative: the
    bounce passes close to a SADDLE POINT of the primed potential, where the
    Hessian has mixed signature.  That negative curvature is what the
    eigenvector rotation of the companion figures feeds on.

The radius r_W of the steepest slope of the bounce is marked as a dashed
vertical line, as in the eigenframe-rotation figure, so the dip of the
lowest curve can be read against the wall.

NOTATION: the y axis reads "Eigenvalues of [H'(phibar^b(r))]_ij", the
curves lambda'_k(r), the reference lines m'_k^2 = value, and the primed
channels are numbered with phi'_1 the HEAVIEST (e^(1) belongs to the larger
eigenvalue).  For TWO fields the eigenvalues follow from the rotation angle
theta' of the eigenframe along the bounce
(plot_aside_common.hessian_rotation_two_field):
lambda'_1 belongs to e'^(1) = (cos theta', sin theta') and is the larger
eigenvalue at every radius.  For N != 2 the stored continuous eigenvalue
curves are drawn, heaviest first.  No title: the bounce is named in the
file name.

reads:  hessian_eigenproblem_<tag>.npz -- R, eigenvalues, m2, r_W, N and,
        for N = 2, H_b
OUTPUT  <out-dir>/plot_hessian_eigenvalues_vgeneral_phi_v2_<bounce>.png
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
from plot_common_vgeneral_phi_v2 import DOTTED_LINE_WIDTH, LINE_WIDTH

#: line styles of the false-vacuum reference lines, one per channel
FV_STYLES = ('--', ':', '-.', (0, (3, 1, 1, 1)))


def draw_eigenvalues(R, eigenvalues, m2, r_W):
    """One axes: the eigenvalues lambda'_k(r) of the primed Hessian, the
    false-vacuum values m'_k^2 they must relax onto far outside the wall and
    the radius r_W of the steepest bounce slope.  No title: the bounce is
    named in the file name."""
    N = eigenvalues.shape[1]
    # a picture <= 8 in wide with its legend (sizes: plot_common)
    fig, ax = plt.subplots(figsize=(5.6, 3.9))
    for k in range(N):
        ax.plot(R, eigenvalues[:, k],
                label=rf"$\lambda'_{{{k + 1}}}(r)$")
    for k in range(N):
        ax.axhline(m2[k], color='k',
                   ls=FV_STYLES[k % len(FV_STYLES)],
                   # a dotted line thicker (plot_common.DOTTED_LINE_WIDTH)
                   lw=(DOTTED_LINE_WIDTH if FV_STYLES[k % len(FV_STYLES)] == ':'
                       else LINE_WIDTH),
                   label=rf"$m_{{{k + 1}}}'^{{\,2}} = {m2[k]:.4f}$")
    ax.axvline(r_W, color='crimson', ls='--',
               label=rf'$r_w = {r_W:.2f}$')
    ax.set_xlabel('$r$')
    ax.set_ylabel(r"Eigenvalues of $[H'(\overline{\phi}^{\,b}(r))]_{ij}$")
    ax.grid(True)
    legend_outside(ax)
    return fig


def main():
    args = parse_cli('aside figure: eigenvalues of the primed Hessian along '
                    'the bounce')
    data = read_aside_npz(args, 'hessian_eigenproblem')
    R = data['R']
    eigenvalues = data['eigenvalues']
    m2 = data['m2']
    r_W = float(data['r_W'])
    bounce_pair_name = bounce_label(args)
    if int(data['N']) == 2:
        # two fields: lambda'_1 and lambda'_2 belong to e'^(1) and e'^(2) of
        # the rotated eigenframe along the bounce (lambda'_1 the larger one);
        # otherwise the stored continuous eigenvalue curves
        _theta, _eigenvectors, eigenvalues = hessian_rotation_two_field(
            data['H_b'], eigenvalues)

    # the two readings the picture makes: where each eigenvalue curve ends
    # up, and how deep the lowest one dips (the saddle the bounce passes)
    for k in range(eigenvalues.shape[1]):
        print(f'[plot] eigenvalue curve {k + 1}: lambda -> '
              f'{eigenvalues[-1, k]:+.6f} at '
              f'r_max  (false-vacuum m^2 = {m2[k]:+.6f}),  minimum '
              f'{eigenvalues[:, k].min():+.6f} at r = '
              f'{R[int(np.argmin(eigenvalues[:, k]))]:.4f}')
    print(f'[plot] steepest bounce slope at r_W = {r_W:.4f}')

    fig = draw_eigenvalues(R, eigenvalues, m2, r_W)
    save(fig, args.out_dir, f'plot_hessian_eigenvalues_vgeneral_phi_v2_{bounce_pair_name}.png')


if __name__ == '__main__':
    main()
