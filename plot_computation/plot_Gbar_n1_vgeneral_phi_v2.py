#!/usr/bin/env python3
"""plot_Gbar_n1_vgeneral_phi_v2.py -- sector n = 1,
UNRENORMALIZED -- the traced Green's function of the wave, the zero-mode
contribution, and the zero-mode-subtracted trace vs s^2.

Same three objects as plot_Gbar_n0_vgeneral_phi_v2, read from
sector_n1_<tag>.npz:

  Gbar_{n=1}(s^2)    the unsubtracted, unrenormalized trace of the wave n = 1,
                     the n = 1 term of the sum over all waves (deg-weighted;
                     box-dependent), drawn with its data points.
  Pbar_{n=1}(s^2) = 4 / (omega^2_{1,0} + s^2)   the ZERO-MODE contribution:
                     the four translation zero modes (deg = (n+1)^2 = 4) enter
                     the trace as a simple pole at s^2 = -omega^2_{1,0} ~ 0
                     (omega2_10 from the eigenvalue_zero_vgeneral_phi_v2 det C
                     crossing; essentially zero for a true bounce, so the pole
                     sits below the grid floor, off the log axis to the left:
                     nothing is drawn for it) -- the term subtracted when the
                     modes are projected out.
  Gbar_{n=1} - Pbar_{n=1}   the zero-mode-SUBTRACTED unrenormalized trace.

THE PICTURE  title "Zero mode contribution"; a log s^2 axis and a LINEAR y
axis, as for n = 0: every value is positive, so the axis starts just below
zero and is cut at 4 times the subtracted curve -- the ~4/s^2 rise of
Gbar_{n=1} and Pbar_{n=1} towards small s^2 leaves it at the top, and the gap
between them is the subtracted curve; Pbar_{n=1} on a fine s^2 grid (it is a
formula); x label s^2, no y label and no footer.

The sector facts, the fail-loud production-identity gates and the drawing are
IMPORTED from gbar_sector_common_vgeneral_phi_v2 (one owner, shared with
plot_Gbar_n0_vgeneral_phi_v2; only the sector differs).

OUTPUT   plot_Gbar_n1_vgeneral_phi_v2_<bounce>.png   in --out-dir (default
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
        description='n=1 unrenormalized bounce-alone '
                    'trace, zero-mode contribution, and its subtraction')
    add_standard_plot_cli(ap)
    args = ap.parse_args()
    resolve_run(args)      # bare invocation: auto-discover tag + bounce npz

    run = load_run(args.data_dir, args.tag, args.bounce_npz,
                   need=('sector_n1',))
    f = unrenormalized_sector_facts(run['sector_n1'], 1)
    r_max = float(run['bg']['R'][-1])
    print(f'[plot_Gbar_n1] n=1: deg = {f["deg"]:g}, omega2_10 = '
          f'{f["omega2"]:+.10e} -> pole at s^2 = {f["s2_c"]:.10g}; '
          f'{f["nodes"]} nodes on [{f["s2"][0]:g}, {f["s2"][-1]:g}]; '
          f'r_max = {r_max:.6g}')

    # same geometry as plot_Gbar_n0 (the legend hangs BELOW the axes)
    fig, ax = plt.subplots(figsize=(7.0, 4.4))     # a picture <= 8 in wide
    draw_sector_unren_curves(ax, f)
    ax.set_title('Zero mode contribution')
    save(fig, args.out_dir,
         f'plot_Gbar_n1_vgeneral_phi_v2_{bounce_label(args)}.png')
    plt.close(fig)


if __name__ == '__main__':
    main()
