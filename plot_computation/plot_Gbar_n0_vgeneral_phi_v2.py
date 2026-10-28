#!/usr/bin/env python3
"""plot_Gbar_n0_vgeneral_phi_v2.py -- sector n = 0,
UNRENORMALIZED -- the traced Green's function of the wave, the negative-mode
contribution, and the negative-mode-subtracted trace vs s^2.

THE THREE CURVES (all read from sector_n0_<tag>.npz -- nothing recomputed)

  Gbar_{n=0}(s^2)    the unsubtracted, unrenormalized trace of the wave n = 0,
                     the n = 0 term of the sum over all waves (deg-weighted,
                     the radial integral up to r_max), added here from the two
                     stored traces: Gbar_{n=0} = Gbar_minus_fv + Gbar_fv.  It
                     carries the bound-state pole of the negative mode.
                     BOX-DEPENDENT: the trace alone grows with the radial
                     box.  Drawn with its data points (the s^2 nodes).

  Pbar_{n=0}(s^2) = 1 / (omega^2_{0,0} + s^2)    the NEGATIVE-MODE
                     contribution: the one bound state of the n = 0 sector
                     enters the trace as a simple pole at s^2 = -omega^2_{0,0}
                     > 0 with residue deg = (n+1)^2 = 1 (Coleman's single
                     negative mode; omega2_00 from the
                     eigenvalue_neg_vgeneral_phi_v2 det C crossing) -- the
                     term subtracted when the mode is projected out.

  Gbar_{n=0} - Pbar_{n=0}   the negative-mode-SUBTRACTED unrenormalized
                     trace: smooth through the pole -- the mode has been
                     removed.

THE PICTURE  title "Negative mode contribution"; a log s^2 axis and a LINEAR
y axis with the signed values: Gbar_{n=0} and Pbar_{n=0} change sign across
the pole, so their negative parts lie below zero and the pole is a jump from
-inf to +inf (marked by a dash-dotted line; the axis is cut symmetrically at
4 times the subtracted curve); Pbar_{n=0} is drawn on a fine s^2 grid (it
is a formula); x label s^2, no y label and no footer.

GATES (fail-loud)  the stored production identities are re-verified on the
delivered arrays before anything is drawn:
    delta_geq3 == Gbar_minus_fv - Gbar_1 - Gbar_2     (keys belong to the wave)
    delta_sub == delta_geq3 - P_0                     (same pole term as production)
so the Pbar_{n=0} drawn here is EXACTLY the pole the pipeline subtracted.

OUTPUT   plot_Gbar_n0_vgeneral_phi_v2_<bounce>.png   in --out-dir (default
the run's pictures/pictures_<potential>_<bounce>/computation/); <bounce> is
the bounce name from the bounce npz (pair_name, e.g. M2_TV0) -- figure
names never carry the run tag/date.
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.dont_write_bytecode = True
_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

from gbar_sector_common_vgeneral_phi_v2 import (draw_sector_unren_curves,
                                                unrenormalized_sector_facts)
from plot_common_vgeneral_phi_v2 import (add_standard_plot_cli, bounce_label,
                                         load_run, resolve_run, save)


def main():
    ap = argparse.ArgumentParser(
        description='n=0 unrenormalized bounce-alone '
                    'trace, negative-mode contribution, and its subtraction')
    add_standard_plot_cli(ap)
    args = ap.parse_args()
    resolve_run(args)      # bare invocation: auto-discover tag + bounce npz

    run = load_run(args.data_dir, args.tag, args.bounce_npz,
                   need=('sector_n0',))
    sector_facts = unrenormalized_sector_facts(run['sector_n0'], 0)
    r_max = float(run['bg']['R'][-1])
    print(f'[plot_Gbar_n0] n=0: deg = {sector_facts["deg"]:g}, omega2_00 = '
          f'{sector_facts["omega2"]:+.10e} -> pole at s^2 = {sector_facts["s2_c"]:.10g}; '
          f'{sector_facts["nodes"]} nodes on [{sector_facts["s2"][0]:g}, {sector_facts["s2"][-1]:g}]; '
          f'r_max = {r_max:.6g}')

    # the legend hangs BELOW the axes (legend_outside); bbox_inches='tight'
    # in save() grows the PNG downwards to fit it
    fig, ax = plt.subplots(figsize=(7.0, 4.4))     # a picture <= 8 in wide
    draw_sector_unren_curves(ax, sector_facts)
    ax.set_title('Negative mode contribution')
    save(fig, args.out_dir,
         f'plot_Gbar_n0_vgeneral_phi_v2_{bounce_label(args)}.png')
    plt.close(fig)


if __name__ == '__main__':
    main()
