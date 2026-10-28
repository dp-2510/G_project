# plot_computation -- the figures of a run

The codes of this folder draw the objects of one finished run of the main computation.
Every code reads the stored stage npz of the run: every figure draws stored arrays or
numbers rebuilt from them with the routines of `computation/`, and stops when a rebuilt
number does not reproduce the stored one.  One code draws one picture (the bounce code the
gallery of every bounce of the run, two pictures per bounce, one for N = 1 and N > 3), and
the picture carries the name of the code.

| code | picture |
|---|---|
| `plot_bounce_gallery_vgeneral_phi_v2.py` | `..._<bounce>_potential.png` of the run's bounce |
| `plot_bounce_gallery_vgeneral_phi_v2.py` | `..._<bounce>_profiles.png` of the run's bounce |
| `plot_det_C_vgeneral_phi_v2.py` | `..._<bounce>.png` |
| `plot_Gbar_n0_vgeneral_phi_v2.py` | `..._<bounce>.png` |
| `plot_Gbar_n1_vgeneral_phi_v2.py` | `..._<bounce>.png` |
| `plot_core_eq_fictitious_vs_s2_vgeneral_phi_v2.py` | `..._<bounce>.png` |
| `plot_lnD_ren_<part>_vs_s2_vgeneral_phi_v2.py` | `..._<bounce>.png` |
| `plot_lnD_ren_<part>_vs_n_vgeneral_phi_v2.py` | `..._<bounce>.png` |
| `plot_lnD_ren_high_n_fit_vgeneral_phi_v2.py` | `..._<bounce>.png` |

(`...` = the code's name without `.py`; `<part>` = `numerical_part`, `analytical_part`,
`numerical_plus_analytical_part`, in this order.)

## Running the codes

One code on a finished run:

```
python3 plot_computation/<plot code>.py --data-dir <DATA_DIR>
```

`--data-dir` defaults to `$VGENERAL_PHI_DATA`, else the current folder.
`--tag` and `--bounce-npz` are found by `plot_common.resolve_run` when they are not given:
the tag comes from the most recent `lnD_ren_*.npz` (else `bounce_data_*.npz`) in the
data directory or one level below it, and the bounce file is `bounce_data_<tag>.npz`,
else the only `bounce_data_*.npz` of that folder.  The choice is printed.  The code
stops when it finds no run or no bounce file, and when several bounce files leave the
choice open (it then lists them).  `--out-dir` is the figure folder (default: the
folder of the run's potential and bounce, see *Picture folders* below).

All codes of a run, in this order (bounce, det C, sectors, the unrenormalised
projected core equation, the final result in parts):

```
plot_computation/run_diagnostic_plots_vgeneral_phi_v2.sh <DATA_DIR> [TAG] [BOUNCE_NPZ] [OUT_DIR]
```

---

## The rules every figure follows

- **Picture names.**  Every picture is `<plot code>_<bounce>.png`: the code's own name
  and the bounce name (the key `pair_name` of the bounce npz, e.g. `M2_TV0`, through
  `plot_common.bounce_label`); the bounce code writes `<plot code>_<NAME>_potential.png`
  and `<plot code>_<NAME>_profiles.png` for every bounce NAME of the run.  A picture name
  never carries the run tag or a date; the input files keep the full tag.
- **Picture folders.**  Every picture of a run is saved in
  `pictures/pictures_<potential>_<bounce>/<group>/` inside the pipeline
  (`plot_common.figure_folder`, the default of `--out-dir`), named from the run's bounce
  npz: `<potential>` is the fields and every non-zero coupling written in front of its
  monomial, to 6 significant digits (`x-y_1.2c_-2x2_-2y2_1xy2_1x4_0.15x2y2_1y4` for the
  example run), replaced by `<N>fields_<n>couplings_<8 hex digits>` when it is longer
  than 60 characters; `<bounce>` is the bounce name; the group is `computation/`, or
  `bounce_gallery/` for the bounce gallery.  `pictures/pictures_<potential>_<bounce>/`
  also holds the note `potential_and_bounce.txt` (the fields, the bounce and every non-zero
  coupling).  Two runs of the same potential and bounce share this folder.
- **One code, one picture.**  Shared computations live in the helper modules
  `plot_common_vgeneral_phi_v2.py`, `traces_vs_s2_n_common_vgeneral_phi_v2.py`,
  `gbar_sector_common_vgeneral_phi_v2.py` and `lnD_ren_parts_common_vgeneral_phi_v2.py`;
  a plot code is its command line, a call into the helper and the drawing of its one
  object.
- **One set of sizes.**  `plot_common` sets the grid and the sizes of every figure once
  (matplotlib rcParams): an opaque light-grey grid (0.72 grey, 1.3 pt, minor lines alike),
  tick labels 11 pt, legends 12 pt, axis labels and every title 13 pt (titles written
  thicker by `plot_common.thick_titles`, a thin outline that thickens words and formulas
  alike, called by `save` and the add-on's `savefig`), every line 1.5 pt (in the two sector
  figures the two curves without data points are thicker, `gbar_sector_common.THICK_LINE_WIDTH`
  = 2.2 pt, so they stay readable on the grid) and every dotted line `DOTTED_LINE_WIDTH` =
  2.2 pt (dots read lighter), every data point `POINT_SIZE` = 4 pt.  Analytical contributions
  are drawn as lines, only the pipeline's numerical data carry points.
  Every picture is meant to be printed at one scale, `\includegraphics[scale=0.69]`, so a
  picture stays at most 8 in wide and these sizes print alike in every figure; plot codes
  set no sizes of their own.
- **Legends.**  A legend never covers the data: every legend sits outside the axes, below
  each panel (`plot_common.legend_outside`, or the legend of `part_figure` in the six
  final-result figures), and to the right of the axes for the profiles of the bounce
  gallery; `save()` writes with `bbox_inches='tight'`, so the picture grows to fit it.  A
  legend entry carries the curve's name.  No figure writes a footer.
- **Axes.**  Every axis is linear or logarithmic; no figure has a symlog or signed-log
  axis (det C is drawn as |det C|, and its y label says so).  An axis label names the
  quantity and its symbol (`s^2`, `Lambda^2`, `n_c`, `r`, `ln D'`); the two sector
  figures carry no y label (their legend names the three curves) and cut their y range
  at 4 times the subtracted curve.  Titles name the object.  The cutoff of the sums over
  the partial waves is written n_c (upright), the index of a wave n (italic).

---

## The figures, in the order they are drawn

### The bounce gallery: `plot_bounce_gallery_vgeneral_phi_v2.py`

**Shows** every solved bounce of the run -- the gallery of all bounces of the potential, from
which a bounce is chosen by its name (main README, Step 2) -- two pictures per bounce
solution (one for N = 1 and N > 3): in the original fieldspace frame the potential with the
bounce path (red) from the true vacuum to the false vacuum, V drawn without
its constant term c (N = 2: the surface V next to the landscape -V, the inverted
potential as a contour map, both with the path; N = 3: the path in field
space, V as its colour; no panel titles) and the bounce profiles phi'_i(r) of the bounce
from the true vacuum to the false vacuum, in the primed fieldspace frame, with the
tree-level action S_E (for N = 1 with the two vacua as levels; for N > 3 the profiles
alone).  The vacua carry the labels of the selection menu, one label per vacuum and the
same in every picture: TV0 the deepest minimum (the highest hilltop of -V), M1, M2, ... the
next ones upwards, FV0 the highest (the lowest hilltop of -V).  Every bounce carries the
name of the two vacua it connects, `<falseLabel>_<trueLabel>` (e.g. `M2_TV0`), in the file
names of its pictures: the name that chooses the bounce of a run.  A bounce that fails the census (SPECIAL) is dashed and marked.  The drawing is
`plot_bounce` of `computation/bounce_gallery_vgeneral_phi_v2.py`, the same renderer the
selection gallery of stage 0 (`<data-dir>/gallery/`) uses.
**Reads** the run's own bounce file (found by `resolve_run`), which names the figure
folder and the run's bounce and gives the potential and the vacuum-search box
(`--seed-range` overrides it), and every `bounce_data_*.npz` in `--gen-dir` (default
`<data-dir>/bounce_gen`, else `<data-dir>` itself) of that potential.  A data dir whose
selection was refused has no run bounce file: its pictures are the stage-0 gallery
`<data-dir>/gallery/`.
**Writes** `plot_bounce_gallery_vgeneral_phi_v2_<NAME>_potential.png` and
`plot_bounce_gallery_vgeneral_phi_v2_<NAME>_profiles.png` for every bounce NAME in
`pictures/pictures_<potential>_<bounce>/bounce_gallery/`, and the table of all bounces
(name, S_E, census) on the terminal.

### det C: `plot_det_C_vgeneral_phi_v2.py`

**Shows** |det C| against s2 at the scan points for the wave n = 0 (top panel, "Negative
mode eigenvalue search") and the wave n = 1 (bottom panel, "Zero mode eigenvalue search"),
both on the same linear s2 window, from -hi/9 to hi (hi = the top of the stored n = 0
scan, rounded up to a multiple of 0.5 x 10^k with 10^k <= top < 10^(k+1)), with the
located roots marked: det C(n, s2) = 0 at
s2 = -omega^2 means that omega^2 is an eigenvalue of the wave (the negative mode at
n = 0, the translation zero mode at n = 1).  C is the matching matrix of the engine's two
branches, the constant Wronskian evaluated at r_max from the value and the slope of the
regular branch,
`C[a,b] = delta_ab + h^-[a,b](r_max) + r_max I_nu K_nu(kappa_a r_max) h^-'[a,b](r_max)`.
**Reads** `eig_n0_<TAG>.npz` and `eig_n1_<TAG>.npz`: the scan arrays `s2_scan`,
`det_scan`, `lambda_scan`, the roots `omega2_00` / `omega2_10`, the located crossings
`omega2_00_census` / `omega2_10_census` and `crossing_resid`.  No h-ODE is solved.
Checks: `det_scan` equals the product of the eigenvalue curves in `lambda_scan` (skipped,
with a printed note, when a complex eigenvalue pair is present, key `max_abs_imag`); every
stored root lies inside its scan window.
**Writes** `plot_det_C_vgeneral_phi_v2_<bounce>.png`.

### The sectors: `plot_Gbar_n0_vgeneral_phi_v2.py`, `plot_Gbar_n1_vgeneral_phi_v2.py`

**Show** the unrenormalised wave n = 0 (n = 1) as three curves against s2: the bounce
trace `Gbar_bounce = Gbar_minus_fv + Gbar_fv` (box-dependent: its radial integral runs up
to r_max; drawn with its data points), the pole of the mode, `P_0 = 1/(omega2_00 + s2)`
for the negative mode (deg = 1) and `P_1 = 4/(omega2_10 + s2)` for the four translation
zero modes (deg = 4), and the difference `Gbar_bounce - P_n`, the trace with the mode
subtracted.  The legend uses symbols such as Gbar_{n=0}, Pbar_{n=0}, omega^2_{1,0}.
Both have a linear y axis with the signed values, cut at 4 times the subtracted curve
(n = 0: the traces change sign at the pole; n = 1: the ~4/s2 rise leaves the axis at the top).
omega2_10 ~ 0, so the n = 1 pole lies below the grid floor and is not drawn.
The facts, the checks and the drawing of both codes come from
`gbar_sector_common_vgeneral_phi_v2.py`.
**Read** `sector_n0_<TAG>.npz` / `sector_n1_<TAG>.npz` (`s2_sector`, `Gbar_minus_fv`,
`Gbar_fv`, `Gbar_1`, `Gbar_2`, `delta_geq3`, `delta_sub`, `deg`, `omega2_00` /
`omega2_10`) and the bounce npz for r_max.  Checks: `delta_geq3 = Gbar_minus_fv - Gbar_1
- Gbar_2` and `delta_sub = delta_geq3 - P_n` on the stored arrays, so the drawn pole is
the one the sector stage subtracted.
**Write** `plot_Gbar_n0_vgeneral_phi_v2_<bounce>.png` and
`plot_Gbar_n1_vgeneral_phi_v2_<bounce>.png`.

### The unrenormalised projected core equation: `plot_core_eq_fictitious_vs_s2_vgeneral_phi_v2.py`

**Shows** the projected core equation against the upper limit Lambda^2 <= s2_max of its s2
integral: the left side Jbar''(Lambda^2) - Jbar''(0) (y label) is the integral from 0 (the
production quadrature `I_n_from_delta`, sign flipped to + Int) of the integrand summed over
n = 0..n_max -- [G - G_FV] with the two pole terms subtracted and the five fictitious modes
added -- plus 5 ln(m_f^2); points at the stored s2, the title with the cutoffs n_max and
s2_max, the legend the right side alone.
**Reads** the residual-band slices (manifest- and fingerprint-checked by `plot_common`),
summed by the helper `traces_vs_s2_n_common_vgeneral_phi_v2.py` into the n = 2..n_max
curve, and the raw n = 0, 1 traces (the slice keys `Gbar_minus_fv_n0/_n1`) with the pole
terms `1/(omega2_00 + s2)` and `4/(omega2_10 + s2)` subtracted at the stored band nodes
(nothing interpolated) and the five fictitious modes `1/(m_f^2 + s2)` (n = 0) and
`4/(m_f^2 + s2)` (n = 1) added (one fictitious mode replaces one removed mode).  The
numbers come from `band_integrals_<TAG>.npz` (n_max, s2_max), `eig_n0_<TAG>.npz` and
`eig_n1_<TAG>.npz` (`omega2_00`, `omega2_10`) and `lnD_ren_<TAG>.npz` (m_f^2 = `m2fict`,
the 4 = `deg1`), so the figure needs the finished run.  Nothing is refitted.
Checks: the Born split
`trace_delta_geq3 = trace_Gbar_minus_fv - trace_Gbar_1 - trace_Gbar_2` on the summed
curves.
**Writes** `plot_core_eq_fictitious_vs_s2_vgeneral_phi_v2_<bounce>.png`.

### The final result in parts: `plot_lnD_ren_<part>_vs_<cutoff>_vgeneral_phi_v2.py`

Six codes, one picture each: `<part>` = `numerical_part`, `analytical_part`,
`numerical_plus_analytical_part`, `<cutoff>` = `s2` (Lambda^2, the cutoff of every
s^2-integral) or `n` (n_c, the cutoff of every sum over the partial waves).
**Shows** one part of the final result in its numerical form,
`ln D' = Sigma_dang + sum_{n=2}^{n_max} I_n(s2_max) + sum_{n=2}^{n_max} T_n(s2_max) + Zcal
+ ln D|_{k=1} + ln D|_{k=2}` (Sigma_dang = (T1), I_n(s2_max) = (T3), T_n(s2_max) = (T5)
`T_n_s2_max`, Zcal = (T6) `Zcal`, ln D|_{k=1} = (T2)+(T4), ln D|_{k=2} = (T7); the
coefficients a, b, c of the tail fit I_n(s^2 = inf) ~ a/(n+1)^3 + b/(n+1)^5 + c/(n+1)^7 are
`fit_a`, `fit_b`, `fit_c`), against one cutoff, a single panel titled
"<Part> vs. Cutoff Lambda^2 / n_c", all six drawn the same way (the same legends and axis
texts, only the variable differs; against n_c, which counts partial waves, every curve is
points at the integers -- a sum over whole waves has no value in between -- against Lambda^2
the computed part is its data points joined by a line and every other curve a line):
*Numerical part* = what the pipeline computes inside its cutoffs (blue: the dangerous waves
n = 0, 1 and the waves n = 2..n_max integrated up to the cutoff, plus the two mu-dependent
terms, closed forms without a cutoff), the cutoff (dotted) and the blue curve evaluated there
(big dot); against Lambda^2 the dangerous part is Sigma_dang integrated up to Lambda^2 (the
two pole-subtracted dangerous waves by the production sector rule,
`running_sector_integral`, plus the five fictitious modes at m_f^2 = `m2fict`), against n_c
wave by wave (n = 0 the negative-mode wave, n = 1 plus the zero-mode wave, each with its
share 1 and 4 of the fictitious modes); both blue curves end on the same number.
*Analytical part* = what lies beyond the cutoffs, each part followed by its complete value in
one legend column: what Sigma_dang adds beyond s2_max (purple: its pole and fictitious terms,
closed form), sum_n T_n (green, from s2_max on: each wave's high-s2 tail spread over s2 by the
uniform high-s2 shape) and Zcal (orange); against Lambda^2 the first two run towards their
complete values (dotted) while Zcal, a sum over whole s^2-integrals, does not depend on
Lambda^2 and is its dotted constant alone; against n_c the sums of the dangerous waves and of
T_n END at n = 1 and at n_max, where they reach their complete values (big dots), and Zcal
runs on as the cumulative sum of the fit towards its complete value (dotted); the cutoff
(dotted) carries its value in the legend, s2_max or n_max, as in the other two parts.
*Numerical + analytical part* = the computed part and, on top of it, the computed part plus
everything the cutoffs miss (black, `numerical_plus_missed_curves`), running onto ln D'
(red), from s2_max/10 on (Lambda^2) and from n_max/2 on (n_c).
The dangerous waves carry the notation of Sigma_dang split at s2_max: -Int_0^{s2_max} ds^2
sigma_dang(s^2) - 5 ln(m_f^2) + ln|omega^2_{0,0}| in the computed part,
-Int_{s2_max}^inf ds^2 sigma_dang(s^2) beyond it (the same labels in all six figures:
`DANG_UP_TO`, `DANG_BEYOND`, `NUMERICAL_PART_LABEL`, `MISSED`).  The numbers and the
drawing live in `lnD_ren_parts_common_vgeneral_phi_v2.py` (`contributions`, the one owner
of the numbers; `numerical_plus_missed_curves`; `cutoff_view`, `draw_numerical`,
`draw_analytical`, `draw_total`, `part_figure`; the legend below the panel).
**Reads** `lnD_ren_<TAG>.npz` (its `certify_lnD` runs first), `sector_n0/_n1_<TAG>.npz`
(`s2_sector`, `delta_sub`), the band tables, `band_integrals_<TAG>.npz` and the bounce
file (R, tr U_int^3 and m_min^2, for the high-s^2 shapes of the two tails).  Checks: the
two dangerous waves add up to
Sigma_dang; the running sector integrals end on the stored Delta's; the running band sum
ends on sum_n I_n(s2_max); I_n(s2_max) + T_n(s2_max) = I_n(s^2 = inf); the stored
T_n(s2_max) are today's completion of the band tables; the end values against both cutoffs
reproduce the stored lnD_ren (1e-9).
**Writes** `plot_lnD_ren_<part>_vs_<cutoff>_vgeneral_phi_v2_<bounce>.png`.

### The high-n fit: `plot_lnD_ren_high_n_fit_vgeneral_phi_v2.py`

**Shows** the fit behind the high-n tail Zcal = (T6): (n+1)^3 I_n(s^2 = inf) of the computed
waves n = 2..n_max as points (I_n(s^2 = inf) = I_n(s2_max) + T_n(s2_max), `I_n_inf`), the fit
a/(n+1)^3 + b/(n+1)^5 + c/(n+1)^7 (`fit_a`, `fit_b`, `fit_c`, made on the waves `fit_lo`..n_max)
scaled by (n+1)^3 -- a + b/(n+1)^2 + c/(n+1)^4 -- as a line on its window and continued to
3 n_max, the window's start (dotted) and end n_max (dashed), and the two n -> infinity limits
as horizontal lines: the fitted a (dotted) and the analytic leading coefficient `A3` =
(1/16) Int dr r^5 tr U_int^3 (dashed), whose agreement validates the fit.  One panel titled
"High n-fitfunction", the legend below it with every number.
**Reads** `lnD_ren_<TAG>.npz` (`band_n`, `I_n_inf`, `fit_a/b/c`, `fit_lo`, `A3`, `n_max`,
`Zcal`); nothing is refitted.  Checks: finite arrays; the stored Zcal is the zeta sum
a zeta(3, n_max+2) + b zeta(5, n_max+2) + c zeta(7, n_max+2) of the stored coefficients
(1e-9); the fit window lies inside the computed waves.
**Writes** `plot_lnD_ren_high_n_fit_vgeneral_phi_v2_<bounce>.png`.

---

## The shared read layer

`plot_common_vgeneral_phi_v2.py` draws nothing.  It merges the band slices for all
per-wave keys, loads a whole run in one call (`load_run`), rebuilds Sigma_dang and
lnD_ren from the stored parts (`certify_lnD`, term by term as the seven terms (T1)-(T7)
of the assembly), provides the n = 0, 1 traces on the band grid
(`trace_G_minus_G_fv_n0_n1_on_band_grid`), the common command line (`add_standard_plot_cli`)
with the run lookup (`resolve_run`), the figure folder of a run (`figure_folder`, the
default of `--out-dir`), the saver (`save`, `bbox_inches='tight'`), the naming helper
`bounce_label` (the bounce name from the bounce npz `pair_name`, through
`pipeline_helpers.bounce_name`) and the layout rule `legend_outside`.  Its checks: the
MANIFEST and bounce-fingerprint check of the band
slices, the per-wave identity `delta_geq3 = Gbar_minus_fv - Gbar_1 - Gbar_2` on the
delivered tables, no non-finite value anywhere, and the exact reproduction of every
stored scalar.
