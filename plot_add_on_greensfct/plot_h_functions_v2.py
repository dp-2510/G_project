#!/usr/bin/env python3
"""plot_h_functions_v2.py -- ONE figure: the h-matrix of BOTH branches at ONE
(n, s^2), all N x N components against r.  h is the deviation of the exact
radial solutions from the free Bessel carriers: the exact mode functions are

    f^-_{beta,i}(r) = B^-_i(r) [delta_{i beta} + h^-_{i beta}(r)],   B^-_i = I_nu(kappa_i r)/r
    f^+_{beta,i}(r) = B^+_i(r) [delta_{i beta} + h^+_{i beta}(r)],   B^+_i = K_nu(kappa_i r)/r

(first index the field component i, second the
basis function beta), with h^-(0) = 0 grown OUTWARD (regular branch) and
h^+(inf) = 0 grown INWARD (decaying branch).  h = 0 everywhere would mean
the bounce's Hessian insertion U_int(r) does nothing at this (n, s^2); the
off-diagonal components are pure field mixing.

NOTATION of the figure: the component (i, beta) of the branch +- is the
h-basisfunction h^{+-,alpha}_{n,i}(r) with alpha = beta, channel 1 the
HEAVIEST.  Title "h^{+-,alpha}_{n,i}(r)-basisfunctions at n, s^2", one
panel per branch, the legend below.

The two h-matrices are the ones the compute code
add_on_greensfct/greensfct_off_coincident_radial_v2.py solved (the ONE
solve_branches call of build_G_n) and stored next to its G_n; this code
only reads and draws them.

    python3 plot_h_functions_v2.py --s2 0.5 --n 3

(The free Bessel carriers these h-matrices dress have their own plot code,
plot_fv_free_modes_v2.py -- same CLI, same npz.)

reads: greensfct_off_coincident_radial_v2_<tag>_n<k>_s2<v>.npz -- r_grid,
       h_minus, h_plus, kappa_i, n, s2
"""
import argparse

import numpy as np
import matplotlib.pyplot as plt
from plot_toolbox_v2 import (add_plot_cli, load_npz, savefig,
                             s2_filetag, trim_edge, resolve_tag, bounce_label)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    add_plot_cli(ap, need_n=True)
    args = ap.parse_args()
    resolve_tag(args)
    h_functions_npz = load_npz(args.data_dir,
                 f'greensfct_off_coincident_radial_v2_{args.tag}_n{args.n}'
                 f'_s2{s2_filetag(args.s2)}.npz')
    r_grid = h_functions_npz['r_grid']
    h_minus = h_functions_npz['h_minus']   # (nr, N, N) each
    h_plus = h_functions_npz['h_plus']
    kappa_i = h_functions_npz['kappa_i']
    n = int(h_functions_npz['n'])
    N = len(kappa_i)
    print('[CHANNELS] kappa_i = sqrt(m_i^2 + s^2): '
          + ',  '.join(f'kappa_{i + 1} = {k:.6f}' for i, k in enumerate(kappa_i)))

    # drop the first ~4% of the grid: next to the engine's start radius the
    # inward branch inherits the K_nu blow-up and would drown the y-scale
    k0 = trim_edge(r_grid)
    r_grid, h_minus, h_plus = r_grid[k0:], h_minus[k0:], h_plus[k0:]

    # ---- THE figure: all N*N components h_ij(r) of BOTH branches ----------
    # left panel the regular (outward) branch h^-, right panel the decaying
    # (inward) branch h^+; one colour per component (i, j), shared panels
    # a picture <= 8 in wide (sizes: plot_common, the one set)
    fig, (axm, axp) = plt.subplots(1, 2, figsize=(7.6, 3.6), sharex=True)
    fig.subplots_adjust(wspace=0.55)     # the two panels well apart
    colors = plt.cm.tab10(np.arange(N * N) % 10)
    for i in range(N):                             # row i = field component
        for alpha in range(N):                     # column = basis function
            c = colors[i * N + alpha]
            # channel numbers = indices + 1 (phi'_1 is index 0), the legend
            # name h^{+-,alpha}_{n,i}(r) of the component (i, alpha)
            axm.plot(r_grid, h_minus[:, i, alpha], color=c,
                     label=rf'$h^{{-,{alpha + 1}}}_{{{n},{i + 1}}}(r)$')
            axp.plot(r_grid, h_plus[:, i, alpha], ls='--', color=c,
                     label=rf'$h^{{+,{alpha + 1}}}_{{{n},{i + 1}}}(r)$')
    axm.set_ylabel(r'$h^{-}(r)$-basisfunctions')
    axp.set_ylabel(r'$h^{+}(r)$-basisfunctions')
    for ax in (axm, axp):
        ax.set_xlabel(r'$r$')
        ax.grid(True)
    fig.suptitle(rf'$h^{{\pm,\alpha}}_{{n,i}}(r)$-basisfunctions at $n={n}$, '
                 rf'$s^2={float(h_functions_npz["s2"]):g}$', y=1.005)
    # ONE shared legend OUTSIDE, below the axes (both panels' curves): the h^-
    # row above the h^+ row (the legend fills column by column, so the two
    # branches are interleaved)
    hm, lm = axm.get_legend_handles_labels()
    hp, lp = axp.get_legend_handles_labels()
    fig.legend([h for pair in zip(hm, hp) for h in pair],
               [l for pair in zip(lm, lp) for l in pair], loc='upper center',
               bbox_to_anchor=(0.5, -0.02), ncol=N * N)
    savefig(fig, args, f'plot_h_functions_v2_{bounce_label(args)}'
                 f'_n{args.n}_s2{s2_filetag(args.s2)}.png')


if __name__ == '__main__':
    main()
