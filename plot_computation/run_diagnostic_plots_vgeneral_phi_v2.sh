#!/usr/bin/env bash
# run_diagnostic_plots_vgeneral_phi_v2.sh -- render the figures of
# plot_computation/ of one vgeneral_phi_v2 run, ONE PLOT CODE PER PICTURE, in
# downflow order: bounce gallery -> det C -> the two dangerous sectors -> the
# unrenormalized projected core equation -> the final result in its three
# parts against each cutoff -> the high-n fit behind its tail Zcal.  The Hessian figures (plot_computation_aside/)
# and the Green's-function figures
# (plot_add_on_greensfct/) are drawn from their own folders' files (README,
# Step 3).  No figure before the final result calls plot_common.certify_lnD:
# det C reads eig_n0/eig_n1, the two sector figures sector_n0/sector_n1, and
# the core-equation figure the band slices, band_integrals, eig_n0/eig_n1
# and lnD_ren (m2fict, deg1), so it needs the finished run.
#
# USAGE
#   plot_computation/run_diagnostic_plots_vgeneral_phi_v2.sh <DATA_DIR> [TAG] [BOUNCE_NPZ] [OUT_DIR]
#
#   DATA_DIR    directory holding the stage npz of the run
#   TAG         run tag (default: AUTO-DISCOVERED from the most recent
#               lnD_ren_*.npz, else bounce_data_*.npz, in DATA_DIR)
#   BOUNCE_NPZ  default <DATA_DIR>/bounce_data_<TAG>.npz (auto-discovered)
#   OUT_DIR     default: every plot code's own folder of this run,
#               pictures/pictures_<potential>_<bounce>/<group>/ inside the
#               pipeline (plot_common.figure_folder)
#
#   Bare invocation on a finished run dir:
#       plot_computation/run_diagnostic_plots_vgeneral_phi_v2.sh <DATA_DIR>
#
# PNG NAMES: every picture is <plot-code-name>_<bounce>.png, where <bounce>
# is the bounce name from the bounce npz (pair_name, e.g. M2_TV0) -- NEVER
# the run tag or a date (the input npz lookup keeps the full TAG).
#
# Every figure script reads ONLY stored production npz and aborts if a
# re-derived scalar does not reproduce the stored one.  This script draws
# exactly the figures listed below.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
PY="${PY:-python3}"                 # override with PY=<interpreter> if wanted

DATA_DIR="${1:?usage: $0 <DATA_DIR> [TAG] [BOUNCE_NPZ] [OUT_DIR]}"
TAG="${2:-}"                        # empty = auto-discovered by the scripts
BOUNCE="${3:-}"                     # empty = auto-discovered by the scripts
OUT="${4:-}"                        # empty = the run's folder in pictures/

echo "[plots] python  : $PY"
echo "[plots] data-dir: $DATA_DIR"
echo "[plots] tag     : ${TAG:-(auto-discover)}"
echo "[plots] bounce  : ${BOUNCE:-(auto-discover)}"
echo "[plots] out-dir : ${OUT:-(each figure in its run folder below $ROOT/pictures)}"
COMMON=(--data-dir "$DATA_DIR")
[ -n "$TAG" ]    && COMMON+=(--tag "$TAG")
[ -n "$BOUNCE" ] && COMMON+=(--bounce-npz "$BOUNCE")
[ -n "$OUT" ]    && COMMON+=(--out-dir "$OUT")

# stamp BEFORE the first figure: the final listing shows exactly the PNGs
# this driver run (re)wrote, whatever bounce name they resolved to
STAMP="$(mktemp "${TMPDIR:-/tmp}/vg2_plots_stamp.XXXXXX")"
trap 'rm -f "$STAMP"' EXIT

echo; echo "===== plot_bounce_gallery_vgeneral_phi_v2: the pictures of every bounce solution ====="
# the bounce gallery: every computed bounce solution of the run, census-marked
# -- the potential with the path (the surface V next to the landscape -V for
# N = 2, the field-space path for N = 3) and the primed profiles with the
# tree-level action S_E
"$PY" "$HERE/plot_bounce_gallery_vgeneral_phi_v2.py" "${COMMON[@]}"

echo; echo "===== plot_det_C_vgeneral_phi_v2: det C zero crossings (n = 0, 1) ====="
# det C (the constant Wronskian, built at r_max) vs s^2
# for sectors n = 0, 1: HOW the zero-crossing condition saw the negative mode
# and the zero mode
"$PY" "$HERE/plot_det_C_vgeneral_phi_v2.py" "${COMMON[@]}"

echo; echo "===== plot_Gbar_n0_vgeneral_phi_v2: sector n = 0, unrenormalized ====="
# the unsubtracted bounce-alone trace, the negative-mode
# pole contribution P_0 = 1/(omega2_00 + s^2), and the negative-mode-subtracted
# trace vs s^2
"$PY" "$HERE/plot_Gbar_n0_vgeneral_phi_v2.py" "${COMMON[@]}"

echo; echo "===== plot_Gbar_n1_vgeneral_phi_v2: sector n = 1, unrenormalized ====="
# the same three curves with the four translation zero
# modes: P_1 = 4/(omega2_10 + s^2) and the zero-mode-subtracted trace vs s^2
"$PY" "$HERE/plot_Gbar_n1_vgeneral_phi_v2.py" "${COMMON[@]}"

echo; echo "===== plot_core_eq_fictitious_vs_s2_vgeneral_phi_v2: the unrenormalized projected core equation ====="
# the left-hand side of the unrenormalized projected core
# equation against the upper limit Lambda^2 <= s2_max of its s^2-integral,
# from the certified band slices, band_integrals, eig_n0/eig_n1 and lnD_ren
# (m2fict, deg1)
"$PY" "$HERE/plot_core_eq_fictitious_vs_s2_vgeneral_phi_v2.py" "${COMMON[@]}"

# the final result in its three parts (numerical,
# analytical, both summed), one picture per part, first against Lambda^2,
# then against n_c
for cutoff in s2 n; do
    for part in numerical_part analytical_part numerical_plus_analytical_part; do
        echo; echo "===== plot_lnD_ren_${part}_vs_${cutoff}_vgeneral_phi_v2 ====="
        "$PY" "$HERE/plot_lnD_ren_${part}_vs_${cutoff}_vgeneral_phi_v2.py" "${COMMON[@]}"
    done
done

echo; echo "===== plot_lnD_ren_high_n_fit_vgeneral_phi_v2: the fit behind the high-n tail ====="
# (n+1)^3 I_n(s^2 = inf) of the computed waves against n, the fit a/(n+1)^3
# + b/(n+1)^5 + c/(n+1)^7 scaled the same way, and the two n -> infinity
# limits, the fitted a and the analytic A3
"$PY" "$HERE/plot_lnD_ren_high_n_fit_vgeneral_phi_v2.py" "${COMMON[@]}"

echo; echo "[plots] DONE"
# list exactly the PNGs THIS run (re)wrote (newer than the start stamp); the
# trailing / also follows a pictures/ that is a link to a folder elsewhere
PICS="${OUT:-$ROOT/pictures}"
if [ -d "$PICS" ]; then
  find "$PICS/" -name 'plot_*.png' -newer "$STAMP" | sort | sed 's/^/    /'
fi
