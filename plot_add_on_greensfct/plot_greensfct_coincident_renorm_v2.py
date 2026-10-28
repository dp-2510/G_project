#!/usr/bin/env python3
"""plot_greensfct_coincident_renorm_v2.py -- every component of the RENORMALISED
coincident Green's function against r, with its three pieces drawn separately
so their assembly is visible, in the notation of the figure:

    [G_{dim-reg, MSbar}(x, x, s^2)]_ij
        = 1/(2 pi^2) sum_{n=0}^{n_max} (n+1)^2 [G^{(k>=2)}_n(r, r, s^2)]_ij
          + [G^{(0)}_{dim-reg, MSbar}(s^2, mu)]_ij
          + [G^{(1)}_{dim-reg, MSbar}(r, s^2, mu)]_ij

(the numerical remainder = Dyson orders >= 2 summed over the waves n = 0 to
n_max, and the two analytic add-backs of Dyson orders 0 and 1).  The values
of s^2 and mu^2 stand in the figure title alone ("Renormalized coincident
x-space Green's function at s^2 = ..., mu^2 = ..."); each panel is titled
[G_{dim-reg, MSbar}(x, x)]_ij, the shared legend (the four objects of the
formula above, without the values) sits right below the panel grid:

    python3 plot_greensfct_coincident_renorm_v2.py --s2 0.5

reads: greensfct_coincident_renorm_v2_<tag>_s2<v>.npz (written by
       add_on_greensfct/greensfct_coincident_renorm_v2.py) -- r_grid, G_ren,
       G0_fv_ren_addback, G1_ren_addback, G_remainder, n_max, s2, mu
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
    coincident_renorm_npz = load_npz(args.data_dir,
                 f'greensfct_coincident_renorm_v2_{args.tag}'
                 f'_s2{s2_filetag(args.s2)}.npz')
    r_grid = coincident_renorm_npz['r_grid']
    G_ren, G0_fv_ren_addback, G1_ren_addback, G_remainder = (
        coincident_renorm_npz[key] for key in
        ('G_ren', 'G0_fv_ren_addback', 'G1_ren_addback', 'G_remainder'))
    k = trim_edge(r_grid)                # drop the r->0 blow-up edge
    r_grid, G_ren, G1_ren_addback, G_remainder = (
        r_grid[k:], G_ren[k:], G1_ren_addback[k:], G_remainder[k:])
    N = G_ren.shape[-1]

    # the symbols of the docstring formula: the scheme subscript dim-reg over
    # MSbar (a LaTeX \substack; mathtext has none, and a \genfrac subscript
    # rides up to superscript height, \underset keeps dim-reg in the
    # subscript position); the values of s^2 and mu go in the figure title
    scheme = r'\underset{\overline{MS}}{\text{dim-reg}}'
    n_max = int(coincident_renorm_npz['n_max'])

    # the rows a bit further apart than the default: the two-line scheme
    # subscript of the panel titles would touch the r label above them
    # a picture <= 8 in wide (sizes: plot_common, the one set)
    fig, axes = plt.subplots(N, N, figsize=(3.7 * N, 2.9 * N), sharex=True,
                             gridspec_kw={'hspace': 0.38})
    axes = np.atleast_2d(axes)
    for i in range(N):
        for j in range(N):
            ax = axes[i, j]
            ax.plot(r_grid, G_ren[:, i, j], color='k',
                    label=rf'$[G_{{{scheme}}}(x,x)]_{{ij}}$')
            ax.plot(r_grid, G_remainder[:, i, j],
                    color='tab:blue',
                    label=rf'$\frac{{1}}{{2\pi^2}}\sum_{{n=0}}^{{{n_max}}}'
                          rf'(n+1)^2\,[G^{{(k\geq 2)}}_n(r,r)]_{{ij}}$')
            ax.axhline(G0_fv_ren_addback[i, j], color='tab:green',
                       label=rf'$[G^{{(0)}}_{{{scheme}}}]_{{ij}}$')
            ax.plot(r_grid, G1_ren_addback[:, i, j],
                    color='tab:orange',
                    label=rf'$[G^{{(1)}}_{{{scheme}}}(r)]_{{ij}}$')
            # channel numbers = indices + 1 (phi'_1 is index 0); aligned at
            # its bottom, so the MSbar under dim-reg stays above the frame
            ax.set_title(rf'$[G_{{{scheme}}}(x,x)]_{{{i + 1}{j + 1}}}$',
                         va='bottom')
            ax.set_xlabel(r'$r$'); ax.grid(True)
    fig.suptitle(r"Renormalized coincident $x$-space Green's function at "
                 rf'$s^2={float(coincident_renorm_npz["s2"]):g}$, '
                 rf'$\mu^2={float(coincident_renorm_npz["mu"]) ** 2:.6g}$', y=1.0372)
    # ONE shared legend OUTSIDE, right below the panel grid (every panel draws
    # the same four objects), read row by row in the order of the docstring
    # formula: the legend fills column by column, so the handles are
    # interleaved
    hh, ll = axes[0, 0].get_legend_handles_labels()
    order = (0, 2, 1, 3)
    fig.legend([hh[o] for o in order], [ll[o] for o in order],
               loc='upper center', bbox_to_anchor=(0.5, 0.02), ncol=2)
    savefig(fig, args, f'plot_greensfct_coincident_renorm_v2_{bounce_label(args)}'
                 f'_s2{s2_filetag(args.s2)}.png')


if __name__ == '__main__':
    main()
