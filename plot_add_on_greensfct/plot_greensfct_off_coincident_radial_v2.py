#!/usr/bin/env python3
"""plot_greensfct_off_coincident_radial_v2.py -- every component of
[G_n(r, r', s^2)]_ij against r and r' as a 3-D surface ('..._3d_surface.png',
all components on ONE scale with one colour bar: the 99th percentile of the
diagonal components, at least the largest off-diagonal value), from the npz
written by add_on_greensfct/greensfct_off_coincident_radial_v2.py:

    python3 plot_greensfct_off_coincident_radial_v2.py --s2 0.5 --n 3

reads: greensfct_off_coincident_radial_v2_<tag>_n<k>_s2<v>.npz -- r_grid,
       G_n, n, s2
"""
import argparse

import numpy as np

from plot_toolbox_v2 import (add_plot_cli, load_npz,
                             savefig, s2_filetag, trim_edge,
                             resolve_tag, surface_grid, bounce_label)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    add_plot_cli(ap, need_n=True)
    args = ap.parse_args()
    resolve_tag(args)
    radial_npz = load_npz(args.data_dir,
                 f'greensfct_off_coincident_radial_v2_{args.tag}_n{args.n}'
                 f'_s2{s2_filetag(args.s2)}.npz')
    # r_grid (nr,), G_n (nr, nr, N, N)
    r_grid, G_n = radial_npz['r_grid'], radial_npz['G_n']
    k = trim_edge(r_grid)                     # drop the r->0 blow-up edge
    r_grid, G_n = r_grid[k:], G_n[k:, k:]

    # this object is ONE partial wave, not a wave sum: the suptitle names the
    # wave n and s^2 (the stored matrix is the full resolvent of the wave)
    sup = (rf"Radial off-coincident Green's function at $n={int(radial_npz['n'])}$, "
           rf'$s^2={float(radial_npz["s2"]):g}$')

    # the panel title: the component alone
    title = r"$[G_n(r,r')]_{{ij}}$"

    # a 3-D surface grid (shared helper), all panels on ONE scale with one
    # colour bar: the 99th percentile of the diagonal components, so that
    # these are cut flat on at most 1% of the grid points (their growth
    # along r = r' towards r = 0), and at least the largest off-diagonal
    # value, so that the off-diagonal components are never cut
    N = G_n.shape[-1]
    diagonal = max(np.nanpercentile(np.abs(G_n[:, :, i, i]), 99)
                   for i in range(N))
    offdiagonal = max((np.nanmax(np.abs(G_n[:, :, i, j]))
                       for i in range(N) for j in range(N) if i != j),
                      default=0.0)
    fig = surface_grid(r_grid, r_grid, G_n, suptitle=sup, title_fmt=title,
                       limit=max(diagonal, offdiagonal))
    savefig(fig, args, f'plot_greensfct_off_coincident_radial_v2_'
                 f'{bounce_label(args)}_n{args.n}'
                 f'_s2{s2_filetag(args.s2)}_3d_surface.png')


if __name__ == '__main__':
    main()
