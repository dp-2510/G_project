#!/usr/bin/env python3
"""hessian_components_wall_vgeneral_phi_v2.py -- the components of the primed
Hessian H_b(r) along the bounce, and the wall radius of every primed
bounce component (README tree: Level 2, side branch computation_aside/; no
driver stage; nothing in lnD_ren or G_ren reads what it writes).

READS    --bounce-npz: the run's bounce_data_<tag>.npz (Level 1), through
         load_background: the bounce grid R, the primed profile phi'_b, the
         channel masses m2 and the interaction Hessian
         U_int = H'(phi'_b) - diag(m2).
WRITES   hessian_components_wall_<tag>.npz in --data-dir, read by
         plot_computation_aside/plot_hessian_components_vgeneral_phi_v2.py
         (hessian_eigenproblem_vgeneral_phi_v2.py imports the function
         hessian_along_bounce from this file):
    R               (nR,)     the bounce grid
    phi_b           (nR,N)    the primed bounce profile phi'_b(r) (the bounce
                              npz key of the same object)
    dphi_b_dr       (nR,N)    its radial derivative, component by component
    H_b             (nR,N,N)  the FULL Hessian with the bounce inserted
    r_W_per_field   (N,)      the wall radius r_W,i of every primed component
                              phi'_i; NaN: phi'_i never moves
    wall_slope      (N,)      |d phi'_i / dr| at that radius
    r_W             ()        THE wall radius of the bounce (pipeline_helpers'
                              wall_radius, from the norm of the slope vector)
    m2              (N,)      false-vacuum channel masses^2 = diag(Lambda)
    N                         the number of fields
    bounce_sha, potential_id, code_version    provenance

HOW
(a) THE HESSIAN COMPONENTS (hessian_along_bounce).  In the primed frame the
    fields are written as phi = phi_F + L phi', with L the orthogonal matrix
    of false-vacuum Hessian eigenvectors, so the Hessian reads
    H'(phi') = L^T H(phi_F + L phi') L and is diagonal at the false vacuum,
    H'(0) = Lambda = diag(m_1^2, ..., m_N^2).  On the bounce grid the code
    forms all components of the FULL Hessian with the bounce inserted,

        H'_ij(r) = [H_b(r)]_ij = [H'(phi'_b(r))]_ij
                 = [diag(m2) + U_int(r)]_ij ,

    from the U_int built by load_background.  The DIAGONAL components are the
    curvature along each primed axis; each relaxes onto its own m_i^2 far
    outside the wall.  The OFF-DIAGONAL components couple the field
    components: they vanish at the false vacuum, where the primed axes are
    the eigen-directions, and for two fields
    H'_12(r) = (Delta lambda / 2) sin(2 mixing_angle_1), with Delta lambda =
    lambda^1 - lambda^2 the eigenvalue splitting and mixing_angle_1 the
    rotation of the eigenframe away from the primed axes (both from
    hessian_eigenproblem_vgeneral_phi_v2): H'_12 is large where a large
    splitting meets a rotated eigenframe.

(b) THE WALL POSITION OF EVERY COMPONENT (profile_slopes,
    wall_position_per_component).  The wall is where the bounce profile
    turns over fastest.  Taken component by component, with d phi'_i / dr
    from second-order differences on the bounce grid,

        r_W,i = argmax_r | d phi'_i / dr | ,

    one radius per PRIMED component: phi'_i is the bounce along the i-th
    false-vacuum mass direction e^(i) (column i of L, mass^2 m_i^2), which
    is the declared field i only when column i of L is that field's unit
    vector, so the printed tables name the channels by index and m_i^2.  A
    component that turns over earlier or later than the others stays
    visible instead of being averaged away.  The pipeline's single wall
    radius r_W, from the norm of the full profile slope, is printed
    alongside, and the component table gives the radius where each Hessian
    component departs most from its false-vacuum value, to read against
    these walls.
"""
import argparse
import os
import sys

import numpy as np

sys.dont_write_bytecode = True
_HERE = os.path.dirname(os.path.abspath(__file__))
_COMP = os.path.join(os.path.dirname(_HERE), 'computation')
for _p in (_COMP, _HERE):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from pipeline_helpers_vgeneral_phi_v2 import (ROUNDOFF_REL,
                                              add_standard_cli,
                                              atomic_savez,
                                              load_background,
                                              provenance_stamp,
                                              require_finite,
                                              wall_radius)


# --------------------------------------------------------------------------- #
#  (a) the primed Hessian along the bounce                                    #
# --------------------------------------------------------------------------- #
def hessian_along_bounce(background):
    """H_b(r) = H'(phi'_b(r)) = diag(m2) + U_int(r) on the bounce grid,
    shape (nR, N, N), from the m2 and U_int that load_background built."""
    H_b = (np.asarray(background['U_int'], float)
           + np.diag(np.asarray(background['m2'], float)))
    require_finite(H_b, '[ABORT] primed Hessian along the bounce')
    return H_b


# --------------------------------------------------------------------------- #
#  (b) the wall position of every bounce component                            #
# --------------------------------------------------------------------------- #
def profile_slopes(R, phi_b):
    """d phi'_i / dr on the bounce grid, shape (nR, N): numpy's gradient with
    the grid R passed in, i.e. second-order differences weighted for the
    unequal spacing of the bounce mesh (the bounce solver refines it where
    the wall is)."""
    dphi_b_dr = np.gradient(phi_b, R, axis=0)
    require_finite(dphi_b_dr,
                   '[ABORT] radial derivative of the bounce profile')
    return dphi_b_dr


def wall_position_per_component(R, dphi_b_dr):
    """The wall radii r_W,i = argmax_r |d phi'_i / dr|, one per primed
    component, and the slope |d phi'_i / dr| reached there.  A component
    that never moves (its largest slope is at most ROUNDOFF_REL of the
    largest slope of all components: round-off) has no wall: NaN."""
    steepness = np.abs(dphi_b_dr)               # (nR, N)
    i_wall = np.argmax(steepness, axis=0)       # one grid index per component
    r_W_per_field = np.asarray(R, float)[i_wall]
    wall_slope = steepness[i_wall, np.arange(steepness.shape[1])]
    # no wall for a component that never moves
    r_W_per_field[wall_slope <= ROUNDOFF_REL * wall_slope.max()] = np.nan
    return r_W_per_field, wall_slope


def main():
    ap = argparse.ArgumentParser(
        description='aside: components of the primed Hessian along the '
                    'bounce and the wall position of every bounce component')
    add_standard_cli(ap)
    args = ap.parse_args()

    # the bounce the user selected: grid R, primed profile phi_b, masses m2
    background = load_background(args.bounce_npz)
    R = background['R']
    phi_b = background['phi_b']
    fields = background['fields']
    N = background['N']
    m2 = np.asarray(background['m2'], float)

    # (a) all components of the primed Hessian along that bounce
    H_b = hessian_along_bounce(background)

    # (b) the wall of each bounce component, plus the pipeline's own wall
    dphi_b_dr = profile_slopes(R, phi_b)
    r_W_per_field, wall_slope = wall_position_per_component(R, dphi_b_dr)
    r_W = wall_radius(R, phi_b)

    out = os.path.join(args.data_dir,
                       f'hessian_components_wall_{args.tag}.npz')
    atomic_savez(out,
                 R=R, phi_b=phi_b, dphi_b_dr=dphi_b_dr,
                 H_b=H_b,
                 r_W_per_field=r_W_per_field,
                 wall_slope=wall_slope,
                 r_W=r_W, m2=m2,
                 N=N,
                 **provenance_stamp(args.bounce_npz, m2,
                                    background['params'], background['fields']))

    # ---- human summary ---------------------------------------------------
    print(f'[HESSIAN-WALL] N = {N} fields {tuple(fields)},  {len(R)} grid '
          f'points,  r in [{R[0]:.4f}, {R[-1]:.4f}]')
    print('[HESSIAN-WALL] wall position of every primed bounce component '
          '(channel i = false-vacuum mass direction e^(i), channel 1 the '
          'heaviest):')
    print('    channel |    m_i^2    |   wall r   |  |dphi\'_i/dr| there')
    print('    --------+-------------+------------+-------------------')
    for i in range(N):
        print(f"    {i + 1:>7} | {m2[i]:+11.6f} | {r_W_per_field[i]:10.4f} | "
              f"{wall_slope[i]:16.6f}")
    print(f'    full profile (norm of the slope vector): r_W = '
          f'{r_W:.4f}')
    # The Hessian components, entries i <= j: the false-vacuum value (m_i^2
    # on the diagonal, zero off it), the value at r_max (close to it: the
    # grid ends on the bounce's tail) and the largest departure from it,
    # with its radius (to read against the walls above).
    print('[HESSIAN-WALL] components of H_b along the bounce '
          '(false-vacuum value, largest departure from it):')
    print('     component   |  FV value  |  at r_max  |   extremum   |  at r')
    print('    -------------+------------+------------+--------------+------')
    for i in range(N):
        for j in range(i, N):        # H_b is symmetric: entries i <= j
            component = H_b[:, i, j]
            fv_value = m2[i] if i == j else 0.0
            i_max_departure = int(np.argmax(np.abs(component - fv_value)))
            label = f"H'_{i + 1}{j + 1}"        # primed channels i, j
            print(f'    {label:>12} | {fv_value:+10.6f} | {component[-1]:+10.6f} | '
                  f'{component[i_max_departure]:+12.6f} | {R[i_max_departure]:7.4f}')
    print(f'[OK] wrote {out}')


if __name__ == '__main__':
    main()
