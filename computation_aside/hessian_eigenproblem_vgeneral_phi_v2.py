#!/usr/bin/env python3
"""hessian_eigenproblem_vgeneral_phi_v2.py -- the eigenvalues and eigenvectors
of the primed Hessian H_b(r) along the bounce, and the mixing angles built
from them (README tree: Level 2, side branch computation_aside/; no driver
stage; nothing in lnD_ren or G_ren reads what it writes).

READS    --bounce-npz: the run's bounce_data_<tag>.npz (Level 1), through
         load_background (R, phi'_b, m2, U_int); H_b = U_int + diag(m2)
         comes from hessian_along_bounce of hessian_components_wall_vgeneral_phi_v2.
WRITES   hessian_eigenproblem_<tag>.npz in --data-dir, read by the plot codes
         plot_hessian_eigenvalues and
         plot_hessian_eigenframe_rotation_two_field of plot_computation_aside/:
    R                 (nR,)      the bounce grid
    phi_b             (nR,N)     the primed bounce profile phi'_b(r) (the
                                 bounce npz key of the same object)
    H_b               (nR,N,N)   the FULL Hessian with the bounce inserted
    eigenvalues       (nR,N)     lambda^alpha(r), continuous eigenvalue curves
    eigenvectors      (nR,N,N)   column alpha = V^alpha(r), sign-continuous
    tilt_angles       (nR,N)     tilt_k(r) in rad, N-general
    tilt_angle_max    (nR,)      max_k tilt_k(r), the largest tilt over the
                                 eigenvectors at each radius: the one-number
                                 mixing measure
    mixing_angle_1    (nR,)      signed angle of V^1, rad  ] N = 2 only,
    mixing_angle_2    (nR,)      signed angle of V^2, rad  ] absent otherwise
    r_W               ()         THE wall radius (pipeline_helpers'
                                 wall_radius: radius of steepest bounce slope)
    m2                (N,)       false-vacuum channel masses^2 = diag(Lambda)
    N                            the number of fields
    bounce_sha, potential_id, code_version    provenance
    The rotation angle theta' that the plot code
    plot_hessian_eigenframe_rotation_two_field draws is mixing_angle_1,
    modulo 180 deg, as long as the two eigenvalues do not cross: both are
    measured from the primed axis e^(1) of the heavier channel, and
    eigenvalue curve 1 is the larger eigenvalue.

HOW (the steps of main)
  (1) the bounce the user selected, through load_background;
  (2) H_b(r_i) = H'(phi'_b(r_i)) at every radius r_i of the bounce grid;
  (3) the eigenproblem at every radius (numpy's symmetric solver eigh), then
      the labels of the eigenvalue curves, from r_max inward
      (order_eigenvalue_curves_continuously);
  (4) the sign of every eigenvector made continuous in r
      (fix_eigenvector_signs), and a warning wherever two eigenvalue curves
      touch (check_eigenvalue_curves_are_separated);
  (5) the angles below, and the wall radius r_W they are compared with.

WHAT IT COMPUTES
----------------
At every radius r_i of the bounce grid, the eigenproblem of the full Hessian
with the bounce inserted, H'_ij(r_i) = [H_b(r_i)]_ij (primed frame:
phi = phi_F + L phi', L the orthogonal matrix whose columns e^(1)..e^(N) are
the false-vacuum Hessian eigenvectors, heaviest channel first,
H'(phi') = L^T H(phi_F + L phi') L):

    H_b(r_i) V^alpha(r_i) = lambda^alpha(r_i) V^alpha(r_i),
    alpha = 1..N,   V^alpha(r_i) orthonormal.

  * lambda^alpha(r) is the local curvature of the potential along the
    eigen-direction V^alpha(r).  Far outside the wall phi'_b -> 0, so
    H_b -> Lambda = diag(m_1^2, ..., m_N^2) and eigenvalue curve alpha
    settles on the channel mass^2 m_alpha^2.  Across the wall the bounce
    crosses the barrier between the two vacua: there the potential curves
    downward along the bounce's path, so the lowest eigenvalue goes NEGATIVE.

  * V^alpha(r) is the local eigen-direction.  The off-diagonal entries of
    H_b rotate it away from the primed axis e^(alpha); that rotation is
    the mixing of the fields.

THE ANGLES
----------
An eigenvector is an AXIS (V^alpha and -V^alpha solve the same problem): the
tilt ignores its sign, the signed angles use the sign-continuous
eigenvectors of step (4).

  tilt_angles (any N)   tilt_k(r) = arccos |<V^k(r), e^(k)>|  in [0, pi/2],
                        the tilt of V^k away from its OWN primed axis
                        (eigenvalue curve k is labelled on axis k at r_max,
                        order_eigenvalue_curves_continuously).  For distinct
                        channel masses every tilt_k tends to zero far
                        outside the wall; for two fields both tilts are
                        equal.

  mixing_angle_1,       mixing_angle_alpha(r) = atan2(V_2^alpha(r), V_1^alpha(r)),
  mixing_angle_2 (N = 2)  V_i the component along the primed axis e^(i): the
                        signed angle of V^alpha measured from e^(1),
                        unwrapped across the +-pi cut of atan2.  It is the
                        angle that diagonalizes the 2 x 2 Hessian,
                        tan(2 mixing_angle_alpha) = 2 H'_12 / (H'_11 - H'_22),
                        and mixing_angle_2 - mixing_angle_1 = pi/2 modulo pi
                        (the eigenvectors are orthonormal).
"""
import argparse
import os
import sys

import numpy as np
from scipy.optimize import linear_sum_assignment

sys.dont_write_bytecode = True
_HERE = os.path.dirname(os.path.abspath(__file__))
_COMP = os.path.join(os.path.dirname(_HERE), 'computation')
for _p in (_COMP, _HERE):
    if _p not in sys.path:
        sys.path.insert(0, _p)

# Step 2 -- the primed Hessian along the bounce: imported, hessian_along_bounce
from hessian_components_wall_vgeneral_phi_v2 import hessian_along_bounce
from pipeline_helpers_vgeneral_phi_v2 import (add_standard_cli,
                                              atomic_savez,
                                              load_background,
                                              provenance_stamp,
                                              wall_radius)


# --------------------------------------------------------------------------- #
#  Step 3 -- the labelling of the eigenvalue curves                           #
# --------------------------------------------------------------------------- #
def order_eigenvalue_curves_continuously(eigenvalues, eigenvectors):
    """Label the eigenvalue curves from the false-vacuum end inward.  At
    every radius r_i, from r_max down to r = 0, the columns are relabelled so
    that eigenvalue curve alpha continues eigenvalue curve alpha of the
    previous radius: the one-to-one pairing with the largest SUM of overlaps
    |<V^alpha(previous), V^beta(r_i)>| (scipy's linear_sum_assignment;
    absolute value: the signs are fixed afterwards).  The eigensolver sorts
    each radius independently, so where two eigenvalue curves cross their
    columns come back exchanged, and this pairing undoes that.  The
    'previous radius' of r_max is the set of primed axes e^(alpha): there
    the bounce is on its tail, H_b = Lambda + U_int(r_max), and eigenvalue
    curve alpha starts as the eigenvalue curve closest to channel alpha.
    For distinct channel masses that is a definite eigenvector lying
    (almost) on the axis; for equal masses the eigenvectors of that subspace
    are a matter of basis, and the pairing is a naming convention.  (The
    eigensolver's ascending order is not used for this: it runs opposite to
    the channel numbering, heaviest first, and for equal or nearly equal
    channel masses the small tail U_int(r_max) can reorder it.)"""
    eigenvalues = np.array(eigenvalues, float, copy=True)
    eigenvectors = np.array(eigenvectors, float, copy=True)
    V_prev = np.eye(eigenvectors.shape[1])      # the primed axes e^(alpha), as columns
    for i in range(len(eigenvectors) - 1, -1, -1):
        # overlap[a, b] = |<V^a(previous radius), V^b(r_i)>|, the sum over
        # the component index k of V_prev[k, a] * eigenvectors[i][k, b]
        overlap = np.abs(np.einsum('ka,kb->ab', V_prev, eigenvectors[i]))
        prev, new = linear_sum_assignment(-overlap)   # best total pairing
        # linear_sum_assignment returns the chosen pairs as (prev[t], new[t]);
        # sorting them by the PREVIOUS eigenvalue curve index gives order[a] =
        # the column at this radius that continues eigenvalue curve a.
        order = new[np.argsort(prev)]
        eigenvectors[i] = eigenvectors[i][:, order]
        eigenvalues[i] = eigenvalues[i][order]
        V_prev = eigenvectors[i]
    return eigenvalues, eigenvectors


# --------------------------------------------------------------------------- #
#  Step 4 -- sign continuity                                                  #
# --------------------------------------------------------------------------- #
def fix_eigenvector_signs(eigenvectors):
    """Make every eigenvector continuous in r: V and -V solve the same
    eigenproblem, so sweep outward from r = 0 and flip V^alpha(r_i) whenever
    <V^alpha(r_{i-1}), V^alpha(r_i)> < 0; the sign at the first grid point is
    the anchor of the sweep."""
    eigenvectors = np.array(eigenvectors, float, copy=True)
    for i in range(1, len(eigenvectors)):
        points_backwards = np.sum(eigenvectors[i - 1] * eigenvectors[i], axis=0) < 0.0
        eigenvectors[i][:, points_backwards] *= -1.0
    return eigenvectors


def check_eigenvalue_curves_are_separated(R, eigenvalues, rel_tol=1e-6):
    """Warn if two eigenvalue curves touch anywhere on the grid: at a
    crossing the eigenvectors of the degenerate pair are not unique, so every
    angle built from them is meaningless there.  The gap is judged against
    the model's own spectral spread, never against an absolute number."""
    if eigenvalues.shape[1] < 2:
        return float('inf'), float('nan')            # one field: no gap
    # the smallest separation of ANY two eigenvalue curves at each radius:
    # the eigenvalue curves are labelled by continuity, not by size, so the
    # gaps are taken between neighbours of the SORTED eigenvalues
    gap_of_r = np.diff(np.sort(eigenvalues, axis=1), axis=1).min(axis=1)
    spread = float(eigenvalues.max() - eigenvalues.min())
    i = int(np.argmin(gap_of_r))
    min_gap, r_min_gap = float(gap_of_r[i]), float(R[i])
    if min_gap < rel_tol * spread:
        print(f'[WARN DEGENERACY] eigenvalue curves approach to '
              f'{min_gap:.3e} at r = {r_min_gap:.4f}, i.e. below '
              f'{rel_tol:g} of the spectral spread {spread:.3e}: the '
              f'eigenvectors are not unique there, so the mixing angles '
              f'near that radius are not meaningful.')
    return min_gap, r_min_gap


# --------------------------------------------------------------------------- #
#  Step 5 -- the mixing angles                                                #
# --------------------------------------------------------------------------- #
def tilt_angles_to_primed_axes(eigenvectors):
    """tilt_k(r) = arccos |<V^k(r), e^(k)>|, shape (nR, N), range [0, pi/2]:
    the tilt of V^k away from its own primed axis e^(k), i.e. the k-th
    component of the k-th eigenvector (absolute value: V and -V are the
    same axis).  Eigenvalue curve k is labelled on e^(k) at r_max, where
    H_b = Lambda up to the bounce's tail
    (order_eigenvalue_curves_continuously), so for distinct channel masses
    every tilt_k tends to zero far out (for equal masses the pairing is a
    convention and the far-out tilt depends on the basis)."""
    overlap = np.abs(np.diagonal(eigenvectors, axis1=1, axis2=2))   # <V^k, e^(k)>
    return np.arccos(np.clip(overlap, 0.0, 1.0))


def mixing_angles_two_field(eigenvectors):
    """The signed, continuous mixing angles mixing_angle_1(r),
    mixing_angle_2(r) for N = 2: mixing_angle_alpha = atan2(V_2^alpha,
    V_1^alpha) (numpy arctan2; V_i = the component along the primed axis
    e^(i), so the angle is measured from e^(1)) of the sign-continuous
    eigenvectors, with the 2 pi jumps at the atan2 cut removed (numpy
    unwrap); anchored at the eigensolver's value at r = 0."""
    angles = []
    for alpha in (0, 1):
        raw = np.arctan2(eigenvectors[:, 1, alpha], eigenvectors[:, 0, alpha])
        angles.append(np.unwrap(raw))
    return angles[0], angles[1]


# --------------------------------------------------------------------------- #
#  Driver                                                                     #
# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(
        description='aside: eigenvalue problem of the primed Hessian along '
                    'the bounce (eigenvalues, eigenvectors, mixing angles)')
    add_standard_cli(ap)
    args = ap.parse_args()

    # (1) the bounce the user selected: grid R, primed profile phi_b, fields
    background = load_background(args.bounce_npz)
    R = background['R']
    phi_b = background['phi_b']
    fields = background['fields']
    N = background['N']

    # (2) the primed Hessian along that bounce
    H_b = hessian_along_bounce(background)

    # (3) the eigenproblem H_b V^alpha = lambda^alpha V^alpha at every
    #     radius: H_b is real symmetric, so the symmetric solver eigh
    #     returns real eigenvalues in ASCENDING order and orthonormal
    #     eigenvectors (column alpha of eigenvectors[i] is V^alpha(r_i));
    #     then the labelling of the eigenvalue curves (from the false-vacuum
    #     end inward, where eigenvalue curve k is channel k)
    eigenvalues, eigenvectors = np.linalg.eigh(H_b)
    eigenvalues, eigenvectors = order_eigenvalue_curves_continuously(
        eigenvalues, eigenvectors)

    # (4) the sign of every eigenvector, made continuous in r; then the
    #     warning wherever two eigenvalue curves touch (there the
    #     eigenvectors, and every angle built from them, are not unique)
    eigenvectors = fix_eigenvector_signs(eigenvectors)
    min_gap, r_min_gap = check_eigenvalue_curves_are_separated(R, eigenvalues)

    # (5) the tilt angles (any N), their maximum over the eigenvectors, and
    #     the wall radius they are compared with; for N = 2 the signed
    #     mixing angles follow below
    tilt_angles = tilt_angles_to_primed_axes(eigenvectors)
    tilt_angle_max = tilt_angles.max(axis=1)
    r_W = wall_radius(R, phi_b)

    saved = dict(R=R, phi_b=phi_b, H_b=H_b,
                 eigenvalues=eigenvalues, eigenvectors=eigenvectors,
                 tilt_angles=tilt_angles,
                 tilt_angle_max=tilt_angle_max,
                 r_W=r_W, m2=np.asarray(background['m2'], float),
                 N=N)
    if N == 2:
        mixing_angle_1, mixing_angle_2 = mixing_angles_two_field(eigenvectors)
        saved['mixing_angle_1'] = mixing_angle_1
        saved['mixing_angle_2'] = mixing_angle_2
    else:
        # the signed planar angle is a two-field object; other N carry
        # the N-general content in tilt_angles
        print(f'[EIGENPROBLEM] signed planar mixing angle exists only for '
              f'N = 2 -- skipped for N = {N} (tilt_angles carry the '
              f'N-general content)')

    out = os.path.join(args.data_dir, f'hessian_eigenproblem_{args.tag}.npz')
    atomic_savez(out, **saved,
                 **provenance_stamp(args.bounce_npz, background['m2'],
                                    background['params'], background['fields']))

    # ---- human summary ---------------------------------------------------
    print(f'[EIGENPROBLEM] N = {N} fields {tuple(fields)},  '
          f'{len(R)} grid points,  r in [{R[0]:.4f}, {R[-1]:.4f}]')
    print(f'[EIGENPROBLEM] wall radius r_W (steepest bounce slope) = '
          f'{r_W:.4f}')
    print('[EIGENPROBLEM] eigenvalues of H_b along the bounce:')
    for k in range(N):
        print(f'    eigenvalue curve {k + 1}: lambda in '
              f'[{eigenvalues[:, k].min():+.6f}, {eigenvalues[:, k].max():+.6f}]'
              f'   -> {eigenvalues[-1, k]:+.6f} at r_max   '
              f'(false-vacuum m^2 = {background["m2"][k]:+.6f})')
    if N > 1:
        print(f'[EIGENPROBLEM] smallest separation of two eigenvalue curves '
              f'{min_gap:.6f} at r = {r_min_gap:.4f}')
    print('[EIGENPROBLEM] tilt_angles: tilt of the eigenvector of each '
          'eigenvalue curve off its primed axis:')
    for k in range(N):
        i_max = int(np.argmax(tilt_angles[:, k]))
        print(f'    eigenvalue curve {k + 1}: tilt in '
              f'[{np.degrees(tilt_angles[:, k].min()):6.3f}, '
              f'{np.degrees(tilt_angles[:, k].max()):6.3f}] deg,   peak at '
              f'r = {R[i_max]:.4f}')
    i_max = int(np.argmax(tilt_angle_max))
    print(f'[EIGENPROBLEM] tilt_angle_max: largest tilt over '
          f'all eigenvectors: '
          f'{np.degrees(tilt_angle_max[i_max]):.3f} deg '
          f'at r = {R[i_max]:.4f}   '
          f'(wall at r_W = {r_W:.4f})')
    if N > 1:
        # the eigenvalue splitting Delta_lambda = lambda_max - lambda_min;
        # for two fields it sets the size of the off-diagonal component,
        # |H'_12| = (Delta_lambda / 2) |sin(2 mixing_angle_1)| (hessian_components_wall)
        Delta_lambda = eigenvalues.max(axis=1) - eigenvalues.min(axis=1)
        i_pk = int(np.argmax(Delta_lambda))
        print(f'[EIGENPROBLEM] eigenvalue splitting Delta_lambda = '
              f'lambda_max - lambda_min peaks at {Delta_lambda[i_pk]:.6f} '
              f'(r = {R[i_pk]:.4f}; lowest eigenvalue there '
              f'{eigenvalues[i_pk].min():+.6f})')
    if N == 2:
        # mixing_angle_1 where the two diagonal components come closest (at
        # H'_11 = H'_22 with H'_12 != 0 the eigenvectors are (1, +-1)/sqrt 2:
        # 45 deg modulo 90, maximal mixing) and at the radius of the lowest
        # eigenvalue
        H11_minus_H22 = H_b[:, 0, 0] - H_b[:, 1, 1]
        i_eq = int(np.argmin(np.abs(H11_minus_H22)))
        i_low = int(np.argmin(eigenvalues.min(axis=1)))
        print(f'[EIGENPROBLEM] mixing_angle_1 = {np.degrees(mixing_angle_1[i_eq]):+.3f} '
              f"deg where |H'_11 - H'_22| is smallest (r = {R[i_eq]:.4f}, "
              f"|H'_11 - H'_22| = {abs(H11_minus_H22[i_eq]):.3g}; at H'_11 = H'_22 "
              f"with H'_12 != 0 it is 45 deg modulo 90 = maximal mixing), and "
              f'{np.degrees(mixing_angle_1[i_low]):+.3f} deg at the lowest '
              f'eigenvalue (r = {R[i_low]:.4f})')
        print(f'[EIGENPROBLEM] signed mixing angles: mixing_angle_1 in '
              f'[{np.degrees(mixing_angle_1.min()):.3f}, {np.degrees(mixing_angle_1.max()):.3f}] '
              f'deg,  mixing_angle_2 in [{np.degrees(mixing_angle_2.min()):.3f}, '
              f'{np.degrees(mixing_angle_2.max()):.3f}] deg')
    print(f'[OK] wrote {out}')


if __name__ == '__main__':
    main()
