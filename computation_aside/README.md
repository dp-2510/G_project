# computation_aside -- the field-space Hessian along the bounce

This folder describes how the field-space Hessian changes along the bounce.  In the primed frame
the fields are written as `phi = phi_F + L phi'`, with `L` the orthogonal matrix of eigenvectors
of the Hessian at the false vacuum; there the primed Hessian is diagonal,
`H'(0) = diag(m_1^2, ..., m_N^2)`.  With the bounce inserted, the full Hessian

    H_b(r) = H'(phi_b(r)) = diag(m2) + U_int(r)        (U_int: the interaction Hessian)

is in general neither diagonal nor constant in `r`.  It is the potential term of the fluctuation
operator `M_b = -Box + H'(phi_b(r))`; its off-diagonal entries couple the field components.  The
first code computes its components and the wall radius of every primed component of the bounce;
the second its eigenvalues `lambda^alpha(r)` and eigenvectors `V^alpha(r)`.  The eigenvalues are
the local curvature of the potential seen by the bounce; the rotation of the eigenvectors away
from the primed axes is the mixing of the fields.  Far outside the wall every eigenvalue settles
on a channel mass `m_alpha^2`; for distinct masses every eigenvector settles on its own primed
axis, for equal masses the eigenvectors of that subspace are a matter of basis and the pairing
of the eigenvalue curves is a naming convention.

Nothing computed here enters `lnD_ren` or `G_ren`.

## Where the folder hangs in the workflow

Both codes read only the bounce file of the main computation (`run_vgeneral_phi_v2_production.sh`),
`bounce_data_<tag>.npz` (stage 0, level 1), in the run's data directory.  From it they rebuild the
radial grid `R`, the bounce `phi_b` in the primed frame, the masses `m2` and the insertion `U_int`.
They write their results into the same data directory; no other computation reads them.

## The tree

```
LEVEL 1   bounce_data_<tag>.npz                            stage 0 of the main computation
              |
              v
LEVEL 2   needs only the bounce -- two independent branches
  (h) hessian_components_wall_vgeneral_phi_v2.py   H_b(r), wall radius of every component
        -> hessian_components_wall_<tag>.npz
  (i) hessian_eigenproblem_vgeneral_phi_v2.py      lambda^alpha(r), V^alpha(r), mixing angles
        -> hessian_eigenproblem_<tag>.npz
```

The letters are those of the main README's workflow tree, where these two codes are 2(h) and
2(i).  Code 2(i) builds `H_b` with the function `hessian_along_bounce` of code 2(h); it does not
read the output file of 2(h).

## Running the codes

Run the codes from the pipeline folder, once the bounce exists.  Both take the options of every stage
of the main computation:

* `--bounce-npz`: the bounce file (default `$VGENERAL_PHI_BOUNCE_NPZ`);
* `--data-dir`: the run's data directory (default `$VGENERAL_PHI_DATA`, else the current folder); the output is written there;
* `--tag`: the run tag used in the output file name (default `vgeneral_phi_v2`); pass the tag of the
  run.

## The codes

**`hessian_components_wall_vgeneral_phi_v2.py`** -- every component of `H_b(r)` along the bounce,
and the wall radius of every primed component of the bounce.
* It computes `H_b = diag(m2) + U_int` on the bounce grid (`hessian_along_bounce`).  The
  diagonal entries are the curvature along each primed axis; the off-diagonal entries couple the
  field components.
* It computes the wall radius of every primed component separately,
  `r_W,i = argmax_r |d phi'_b,i/dr|`, from the radial derivative of the bounce, and the pipeline's
  single wall radius `r_W` (the radius of the steepest slope of the whole bounce) for reference.
  Component `i` is the bounce along the false-vacuum mass direction `e^(i)` (column `i` of `L`,
  mass `m_i^2`), not the declared field of the same position.
* Output: `hessian_components_wall_<tag>.npz` with `R`, `phi_b` (the bounce in the primed
  frame), `dphi_b_dr`, `H_b`, `r_W_per_field` (the wall radius `r_W,i` of every primed
  component; NaN for a component the bounce does not move, i.e. whose largest slope is at most
  `ROUNDOFF_REL` of the largest slope of all components), `wall_slope` (the slope
  `|d phi'_b,i/dr|` at `r_W,i`), `r_W`, `m2`, `N`.

```
python3 computation_aside/hessian_components_wall_vgeneral_phi_v2.py \
        --bounce-npz results/bounce_data_<tag>.npz --data-dir results --tag <tag>
```

**`hessian_eigenproblem_vgeneral_phi_v2.py`** -- the eigenvalues `lambda^alpha(r)` and the eigenvectors
`V^alpha(r)` of `H_b(r)` at every radius of the bounce grid, and the mixing angles built from
them.
* It solves the symmetric eigenvalue problem `H_b(r) V^alpha(r) = lambda^alpha(r) V^alpha(r)`
  at every radius.  The eigenvalue curves are labelled continuously from the false-vacuum end
  inward, where eigenvalue curve `k` is channel `k`, and the sign of every eigenvector is made
  continuous in `r`, sweeping outward from `r = 0`; a radius where two eigenvalue curves touch
  is reported (the eigenvectors are not unique there).
* It computes the tilt of every eigenvector away from its own primed axis `e^(k)`,
  `tilt_k(r) = arccos |<V^k(r), e^(k)>|`, for any number of fields, and for two fields the signed
  continuous angles of the two eigenvectors measured from the primed axis `e^(1)`,
  `mixing_angle_alpha = atan2(V_2^alpha, V_1^alpha)` (`V_i` = the component along `e^(i)`); the two
  eigenvectors are orthonormal, so the angles differ by `pi/2` modulo `pi` at every radius.
* Output: `hessian_eigenproblem_<tag>.npz` with `R`, `phi_b`, `H_b`, `eigenvalues`,
  `eigenvectors` (column `alpha` is `V^alpha`), `tilt_angles`, `tilt_angle_max` (the largest
  tilt over the eigenvectors at each radius), `r_W`, `m2`, `N`, and for two fields
  `mixing_angle_1`, `mixing_angle_2`.

```
python3 computation_aside/hessian_eigenproblem_vgeneral_phi_v2.py \
        --bounce-npz results/bounce_data_<tag>.npz --data-dir results --tag <tag>
```

## Plot codes (`plot_computation_aside/`)

The plot codes draw the figures of the Hessian along the bounce.  Each reads one
output file of this folder, found by the same `--data-dir` and `--tag`:

```
python3 plot_computation_aside/<plot code>.py --data-dir results --tag <tag>
```

* `plot_hessian_components_vgeneral_phi_v2.py` -- every component of `H_b(r)`
  against `r`, with its false-vacuum values and the wall radius `r_W,i` of every bounce component as
  a dashed vertical line; from `hessian_components_wall_<tag>.npz`.
* `plot_hessian_eigenvalues_vgeneral_phi_v2.py` -- the eigenvalues `lambda^alpha(r)`
  against `r`, with the false-vacuum masses `m_i^2` as reference lines and the radius `r_W`
  of the steepest bounce slope as a dashed vertical line; from
  `hessian_eigenproblem_<tag>.npz`.
* `plot_hessian_eigenframe_rotation_two_field_vgeneral_phi_v2.py` -- two fields only:
  the eigenvectors as arrows along the bounce and their rotation angle `theta'(r)` against the
  primed axes; same file.

`plot_aside_common_vgeneral_phi_v2.py` is the shared helper of these codes and draws nothing; its
`hessian_rotation_two_field` forms, for two fields, the rotation angle `theta'(r)`, the eigenvectors
`e'^(alpha)(r)` and the eigenvalues `lambda'_alpha(r)` in closed form from the stored `H_b`
(checked against the stored eigenvalues).  The figures number the primed channels as the npz files
do, heaviest first (`phi'_1` the heaviest).  The figures are written as PNG files to the folder of
the run's potential and bounce, `pictures/pictures_<potential>_<bounce>/computation_aside/` inside
the pipeline (`--out-dir` chooses another folder), named `<plot code>_<bounce name>.png`.
