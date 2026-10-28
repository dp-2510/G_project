#!/usr/bin/env python3
"""plot_aside_common_vgeneral_phi_v2.py -- the small shared read/write layer of the aside
figures.  It is NOT a plot code: it draws nothing and writes no picture.

The three plot codes in this folder READ ONLY the npz written by
computation_aside/ (hessian_eigenproblem_vgeneral_phi_v2.py and
hessian_components_wall_vgeneral_phi_v2.py).  Nothing is re-solved here: no bounce,
no Hessian -- every number drawn, except the two-field rotation theta',
e'^(alpha), lambda' that hessian_rotation_two_field below forms in
closed form from the stored Hessian (checked against the stored
eigenvalues), was computed and stored by
the aside computation stage.  This file holds the two things all three codes
would otherwise repeat, and RE-EXPORTS the production plot layer's owners of
the figure conventions (plot_computation/plot_common_vgeneral_phi_v2:
bounce_label, legend_outside, save):

    parse_cli(description)  the run selectors --data-dir / --tag /
                            --bounce-npz and the figure folder --out-dir
                            (default: the run's computation_aside folder of
                            plot_common.figure_folder)
    read_aside_npz(...)     one aside npz, fail-loud, provenance printed

THE FIGURE-NAMING RULE (bounce_label): a picture is named after the PLOT
CODE that drew it plus the BOUNCE it shows -- the vacuum pair, read from the
bounce npz key pair_name (e.g. 'M2_TV0') -- never the full run tag and
never a date, so the same figure of the same bounce keeps one name across
runs.
"""
import argparse
import os
import sys

import numpy as np

sys.dont_write_bytecode = True
_HERE = os.path.dirname(os.path.abspath(__file__))
for _p in ('plot_computation', 'computation'):
    _d = os.path.normpath(os.path.join(_HERE, '..', _p))
    if _d not in sys.path:
        sys.path.insert(0, _d)
from plot_common_vgeneral_phi_v2 import (bounce_label, figure_folder,
                                         legend_outside, save)


def parse_cli(description):
    """The uniform command line of the three aside figures: which run to read
    (--data-dir, --tag), which bounce names the picture and its folder
    (--bounce-npz), and where the picture goes (--out-dir; default the run's
    folder pictures/pictures_<potential>_<bounce>/computation_aside/)."""
    ap = argparse.ArgumentParser(description=description)
    ap.add_argument('--data-dir',
                    default=os.environ.get('VGENERAL_PHI_DATA', '.'),
                    help='data directory holding the computation_aside npz '
                         '(default $VGENERAL_PHI_DATA)')
    ap.add_argument('--tag', default='vgeneral_phi_v2',
                    help='run tag threaded into the npz file names')
    ap.add_argument('--bounce-npz',
                    default=os.environ.get('VGENERAL_PHI_BOUNCE_NPZ'),
                    help='bounce npz of the run; it names the picture and '
                         'its folder (default $VGENERAL_PHI_BOUNCE_NPZ, else '
                         'bounce_data_<tag>.npz in --data-dir)')
    ap.add_argument('--out-dir', default=None,
                    help='figure folder (default: pictures/pictures_'
                         '<potential>_<bounce>/computation_aside/ inside the '
                         'pipeline)')
    args = ap.parse_args()
    if args.bounce_npz is None:
        args.bounce_npz = os.path.join(args.data_dir,
                                       f'bounce_data_{args.tag}.npz')
    if not os.path.isfile(args.bounce_npz):
        raise SystemExit(f'[ABORT] bounce npz {args.bounce_npz} not found '
                         f'-- pass --tag <run tag> or --bounce-npz <file>')
    args.out_dir = args.out_dir or figure_folder('computation_aside', args.bounce_npz)
    return args


def read_aside_npz(args, stem):
    """Read <stem>_<tag>.npz written by computation_aside/, as a plain dict.

    Fail-loud: a missing file names the computation code that produces it
    instead of drawing an empty figure.  The provenance the computation stage
    stamped in (which bounce file, which potential) is printed, so a picture
    can always be traced back to the bounce it was made from."""
    path = os.path.join(args.data_dir, f'{stem}_{args.tag}.npz')
    if not os.path.isfile(path):
        raise SystemExit(
            f'[ABORT] {os.path.basename(path)} not found in {args.data_dir} '
            f'-- run the matching computation_aside code first (this figure '
            f'only reads its npz, it computes nothing itself).')
    with np.load(path, allow_pickle=True) as npz:
        data = {key: npz[key] for key in npz.files}
    print(f'[plot] read {path}')
    print(f'[plot] bounce fingerprint {str(data["bounce_sha"])[:12]}...,  '
          f'potential {data["potential_id"]},  code {data["code_version"]}')
    return data


def hessian_rotation_two_field(H_b, eigenvalues):
    """The eigenframe rotation of the primed Hessian along the bounce of
    a TWO-FIELD run (phi'_1 the heavier channel): the rotation angle
    theta'(r) of the local eigenframe of H'(phi'_b(r)) against the primed
    axes,

        theta'(r) = 1/2 atan2( 2 H'_12, H'_11 - H'_22 ),

    its eigenvectors e'^(1) = (cos theta', sin theta'), e'^(2) = (-sin theta',
    cos theta') and the eigenvalues lambda'_1 = e'^(1).H'.e'^(1), lambda'_2 =
    e'^(2).H'.e'^(2), the diagonal of L^T H' L.  Here atan2 puts
    theta' in (-90, 90] degrees and makes e'^(1) the eigenvector of the
    LARGER eigenvalue at every radius, lambda'_1 >= lambda'_2; far outside
    the wall H' -> diag(m'_1^2, m'_2^2) with m'_1 > m'_2, so theta' -> 0.
    theta' is continuous as long as the two eigenvalues do not cross (it
    could only wrap at +-90 degrees, where H'_12 changes sign while
    H'_11 < H'_22).
    `eigenvalues` are the stored eigenvalue curves (from eigh): the two
    lambda' must be the same two numbers at every radius (fail-loud check of
    the formula).
    Returns theta' (nR,) in rad, eigenvectors (nR, 2, 2) with column alpha =
    e'^(alpha), eigenvalues (nR, 2) = (lambda'_1, lambda'_2)."""
    a, b, c = H_b[:, 0, 0], H_b[:, 0, 1], H_b[:, 1, 1]
    theta = 0.5 * np.arctan2(2.0 * b, a - c)
    e1 = np.stack([np.cos(theta), np.sin(theta)], axis=1)
    e2 = np.stack([-np.sin(theta), np.cos(theta)], axis=1)
    lam = np.stack([np.einsum('ri,rij,rj->r', e1, H_b, e1),
                    np.einsum('ri,rij,rj->r', e2, H_b, e2)], axis=1)
    mismatch = np.max(np.abs(np.sort(lam, axis=1) - np.sort(eigenvalues, axis=1)))
    if mismatch > 1e-9 * (1.0 + np.max(np.abs(eigenvalues))):
        raise SystemExit(f'[ABORT] the theta\' decomposition of H\' disagrees '
                         f'with the stored eigenvalues by {mismatch:.2e}.')
    return theta, np.stack([e1, e2], axis=2), lam
