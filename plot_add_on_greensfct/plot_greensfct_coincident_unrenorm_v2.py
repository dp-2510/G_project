#!/usr/bin/env python3
"""plot_greensfct_coincident_unrenorm_v2.py -- every component of the
UNrenormalised [G(x, x, s^2)]_ij against r, at the full cutoff n_max AND at
n_half = n_max/2 in the same panel, so the regulator (n_max) dependence --
i.e. the divergence -- is visible.  The value of s^2 stands in the figure
title alone ("Unrenormalized coincident x-space Green's function at s^2 =
..."), each panel is titled [G(x, x)]_ij, the shared legend sits right
below the panel grid:

    python3 plot_greensfct_coincident_unrenorm_v2.py --s2 0.5

reads: greensfct_coincident_unrenorm_v2_<tag>_s2<v>.npz (written by
       add_on_greensfct/greensfct_coincident_unrenorm_v2.py) -- r_grid,
       G_unren, G_unren_half, n_max, n_half, s2
"""
import argparse

import numpy as np
import matplotlib.pyplot as plt
from plot_toolbox_v2 import (add_plot_cli, load_npz,
                             savefig, s2_filetag, trim_edge, resolve_tag,
                             bounce_label)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    add_plot_cli(ap)
    args = ap.parse_args()
    resolve_tag(args)
    coincident_unrenorm_npz = load_npz(args.data_dir,
                 f'greensfct_coincident_unrenorm_v2_{args.tag}'
                 f'_s2{s2_filetag(args.s2)}.npz')
    r_grid = coincident_unrenorm_npz['r_grid']
    G_unren = coincident_unrenorm_npz['G_unren']
    G_unren_half = coincident_unrenorm_npz['G_unren_half']

    k = trim_edge(r_grid)                # drop the r->0 blow-up edge
    r_grid, G_unren, G_unren_half = r_grid[k:], G_unren[k:], G_unren_half[k:]
    N = G_unren.shape[-1]

    # a picture <= 8 in wide (sizes: plot_common, the one set)
    fig, axes = plt.subplots(N, N, figsize=(3.7 * N, 2.9 * N), sharex=True,
                             gridspec_kw={'hspace': 0.38})
    axes = np.atleast_2d(axes)
    for i in range(N):
        for j in range(N):
            ax = axes[i, j]
            ax.plot(r_grid, G_unren[:, i, j],
                    label=rf'$G_{{\rm unren}}$: waves $n=0$ to ${int(coincident_unrenorm_npz["n_max"])}$')
            ax.plot(r_grid, G_unren_half[:, i, j], ls='--',
                    label=rf'$G_{{\rm unren}}$: waves $n=0$ to '
                          rf'${int(coincident_unrenorm_npz["n_half"])}$')
            # channel numbers = indices + 1 (phi'_1 is index 0)
            ax.set_title(rf'$[G(x,x)]_{{{i + 1}{j + 1}}}$')
            ax.set_xlabel(r'$r$'); ax.grid(True)
    fig.suptitle(r"Unrenormalized coincident $x$-space Green's function at "
                 rf'$s^2={float(coincident_unrenorm_npz["s2"]):g}$', y=1.0062)
    # ONE shared legend OUTSIDE, right below the panel grid (every panel
    # draws the same curves; the legend names the two wave ranges summed)
    hh, ll = axes[0, 0].get_legend_handles_labels()
    fig.legend(hh, ll, loc='upper center', bbox_to_anchor=(0.5, 0.02),
               ncol=len(ll))
    savefig(fig, args, f'plot_greensfct_coincident_unrenorm_v2_{bounce_label(args)}'
                 f'_s2{s2_filetag(args.s2)}.png')


if __name__ == '__main__':
    main()
