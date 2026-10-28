#!/usr/bin/env python3
"""plot_lnD_ren_analytical_part_vs_s2_vgeneral_phi_v2.py -- the analytical part of the
final result ln D' against Lambda^2, the cutoff of every s^2-integral:
    what lies beyond the cutoffs -- the high-s^2 tail, the high-n tail and
    the dangerous waves beyond s2_max -- each running towards its
    complete value (dotted).

Drawn by lnD_ren_parts_common.part_figure (the six figures of the final
result -- numerical part, analytical part and both summed, against Lambda^2
and against n_c -- share their panels, legends and axis texts); the
contributions are computed ONCE, by lnD_ren_parts_common.contributions (with
all its gates).

Output: <out-dir>/plot_lnD_ren_analytical_part_vs_s2_vgeneral_phi_v2_<bounce>.png
(<bounce> = the bounce name from the bounce npz pair_name, e.g. M2_TV0;
figure names never carry the run tag/date).
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

from lnD_ren_parts_common_vgeneral_phi_v2 import contributions, part_figure
from plot_common_vgeneral_phi_v2 import (add_standard_plot_cli, bounce_label,
                                         resolve_run, save)


def main():
    ap = argparse.ArgumentParser(
        description='vgeneral_phi_v2: the analytical part of the final result ln D\' '
                    'against Lambda^2')
    add_standard_plot_cli(ap)
    args = ap.parse_args()
    resolve_run(args)      # bare invocation: auto-discover tag + bounce npz
    fig = part_figure(contributions(args), 's2', 'analytical')
    save(fig, args.out_dir,
         f'plot_lnD_ren_analytical_part_vs_s2_vgeneral_phi_v2_{bounce_label(args)}.png')
    plt.close(fig)


if __name__ == '__main__':
    main()
