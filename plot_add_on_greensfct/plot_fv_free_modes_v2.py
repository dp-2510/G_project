#!/usr/bin/env python3
"""plot_fv_free_modes_v2.py -- ONE figure: the free-channel carriers of ONE
(n, s^2), the Bessel solutions the h-matrices dress, in the notation of
the legend:

    B^-_{n+1}(r) = I_{n+1}(kappa_i r)/r   (regular at r = 0)
    B^+_{n+1}(r) = K_{n+1}(kappa_i r)/r   (decaying at r = inf)

kappa_i = sqrt(m_i^2 + s^2), one pair per channel i (the channel enters
through kappa_i only; channel 1 is the HEAVIEST).  These are the
FALSE-VACUUM mode functions; the exact radial solutions of the bounce
problem are
f^-_{beta,i}(r) = B^-_i(r) [delta_{i beta} + h^-_{i beta}(r)] (and the same
with +; first index the field component, second the basis function), and
the h-matrices have their own plot code, plot_h_functions_v2.py.

DRAWN AS THE STORED LOGS.  The compute code stores the carriers as their
natural logarithms, log_B_minus_fv = ln[I_nu(kappa_i r)/r] and
log_B_plus_fv = ln[K_nu(kappa_i r)/r], because K_nu(kappa r)/r overflows float64
at small r for large nu.  This figure draws exactly those stored arrays,
ln(B^pm(r)) against r on a linear axis (one curve per channel i and
branch): the same picture as B^pm on a log axis, with nothing
exponentiated.  Title "Free fluctuation modes at n, s^2"; the legend sits
outside the axes, at the top right.

    python3 plot_fv_free_modes_v2.py --s2 0.5 --n 3

reads: greensfct_off_coincident_radial_v2_<tag>_n<k>_s2<v>.npz (written by
       add_on_greensfct/greensfct_off_coincident_radial_v2.py) -- r_grid,
       log_B_minus_fv, log_B_plus_fv, kappa_i, n, s2
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
    free_modes_npz = load_npz(args.data_dir,
                 f'greensfct_off_coincident_radial_v2_{args.tag}_n{args.n}'
                 f'_s2{s2_filetag(args.s2)}.npz')
    r_grid = free_modes_npz['r_grid']
    log_B_minus_fv = free_modes_npz['log_B_minus_fv']   # (nr, N)
    log_B_plus_fv = free_modes_npz['log_B_plus_fv']
    kappa_i = free_modes_npz['kappa_i']
    N = len(kappa_i)
    nu = int(free_modes_npz['n']) + 1
    print('[CHANNELS] kappa_i = sqrt(m_i^2 + s^2): '
          + ',  '.join(f'kappa_{i + 1} = {k:.6f}' for i, k in enumerate(kappa_i)))

    # drop the first ~4% of the grid (the K_nu blow-up next to the engine's
    # start radius would drown the y-scale), same trim as every plot here
    k0 = trim_edge(r_grid)
    r_grid, log_B_minus_fv, log_B_plus_fv = r_grid[k0:], log_B_minus_fv[k0:], log_B_plus_fv[k0:]

    # ---- THE figure: one curve per channel and branch ---------------------
    fig, ax = plt.subplots(figsize=(5.6, 4.0))     # a picture <= 8 in wide
    ch_colors = plt.cm.tab10(np.arange(N) % 10)
    for i in range(N):                                   # one colour per
        c = ch_colors[i]                                 # channel
        # channel number = index + 1 (phi'_1 is index 0); the legend writes
        # B^pm_{n+1}, the channel enters through kappa_i
        ax.plot(r_grid, log_B_minus_fv[:, i], color=c,
                label=rf'$B^-_{{{nu}}}(r) = I_{{{nu}}}(\kappa_{{{i + 1}}} r)/r$')
        ax.plot(r_grid, log_B_plus_fv[:, i], ls='--', color=c,
                label=rf'$B^+_{{{nu}}}(r) = K_{{{nu}}}(\kappa_{{{i + 1}}} r)/r$')
    ax.set_xlabel(r'$r$')
    ax.set_ylabel(r'$\ln(B^{\pm}(r))$')
    ax.grid(True)
    ax.set_title(rf'Free fluctuation modes at $n={int(free_modes_npz["n"])}$, '
                 rf'$s^2={float(free_modes_npz["s2"]):g}$')
    # legend OUTSIDE the axes, to the right at the top
    ax.legend(loc='upper left', bbox_to_anchor=(1.02, 1.0), borderaxespad=0.0)
    savefig(fig, args, f'plot_fv_free_modes_v2_{bounce_label(args)}'
                 f'_n{args.n}_s2{s2_filetag(args.s2)}.png')


if __name__ == '__main__':
    main()
