#!/usr/bin/env python3
"""plot_core_eq_fictitious_vs_s2_vgeneral_phi_v2.py -- ONE picture: the
right-hand side of the core equation with the dangerous modes subtracted and
the five fictitious modes added, against the upper limit Lambda^2 of its s^2
integral, up to s2_max:

    Int_0^{Lambda^2} ds^2 ( Gbar(s^2) - Gbar^v(s^2) - 1/(omega^2_{0,0} + s^2)
                            - 4/(omega^2_{1,0} + s^2) + 5/(m_f^2 + s^2) )
    + 5 ln(m_f^2)

THE INTEGRAND sums ALL waves n = 0..n_max: the raw n = 0, 1 traces of the
BASE band slices, with their pole terms subtracted at the stored s^2 nodes,
plus the band waves n = 2..n_max (traces_vs_s2_n_common.add_core_eq_fictitious);
m_f^2 is the fictitious mass^2 of the assembly (lnD_ren npz, m2fict; every
m_f^2 > 0 gives the same total once the s^2 integral runs to infinity).
THE INTEGRAL is the production quadrature (pipeline_quadrature.I_n_from_delta,
with its sign flipped back to + Int): log-Simpson over the stored nodes and
the exact strip [0, s2_min] below the first one, applied to the truncated
grid at every Lambda^2 -- so the first point sits at the third stored s^2.
It has no plateau: at fixed n_max the integrand decays slowly, and the n-sum
itself grows ~ ln n_max -- the UV divergence the renormalization removes.

THE PICTURE  title "Unrenormalized projected core-equation with cutoffs:
n_max = ..., s2_max = ..."; the running integral at the stored s^2 nodes as
points joined by a thin line; log Lambda^2 axis, linear values; y label the
left side of the core equation, Jbar''(Lambda^2) - Jbar''(0); x label
Lambda^2; the legend below the axes holds the right side alone; no footer.

The stored-curve read + certification and the integrand are IMPORTED from
traces_vs_s2_n_common_vgeneral_phi_v2 (one owner; data = the certified band
slices, nothing refitted).

Output: <out-dir>/plot_core_eq_fictitious_vs_s2_vgeneral_phi_v2_<bounce>.png
"""
import os
import sys

import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.dont_write_bytecode = True
_HERE = os.path.dirname(os.path.abspath(__file__))
_COMP = os.path.normpath(os.path.join(_HERE, '..', 'computation'))
for _p in (_HERE, _COMP):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from pipeline_quadrature_vgeneral_phi_v2 import I_n_from_delta
from plot_common_vgeneral_phi_v2 import POINT_SIZE, legend_outside, save
from traces_vs_s2_n_common_vgeneral_phi_v2 import (COLOR_TRUNCATED,
                                                   add_core_eq_fictitious,
                                                   add_trace_G_minus_G_fv_n0_n1,
                                                   load_stored_traces,
                                                   make_cli, picture_name)

CODE_STEM = 'plot_core_eq_fictitious_vs_s2_vgeneral_phi_v2'


def main():
    args = make_cli('vgeneral_phi_v2: the projected core equation -- the s^2 '
                    'integral of its integrand (summed over n = 0..n_max) '
                    'from 0 to Lambda^2 <= s2_max, plus 5 ln(m_f^2)')
    traces = load_stored_traces(args)
    add_trace_G_minus_G_fv_n0_n1(traces, args)      # the raw n = 0, 1 traces + poles
    add_core_eq_fictitious(traces)

    s2, integrand = traces['s2'], traces['core_eq_fict_all']
    deg1, m2_f = traces['deg1'], traces['m2_f']
    n_fict = 1.0 + deg1                              # one per removed mode
    # the running integral: the production quadrature on the grid truncated at
    # every stored Lambda^2 (it needs three nodes, so it starts at the third);
    # I_n_from_delta returns -Int, the core equation has +Int
    Lambda2 = s2[2:]
    core_eq = np.array([-float(I_n_from_delta(s2[:j + 1], integrand[:j + 1]))
                        for j in range(2, len(s2))]) + n_fict * np.log(m2_f)
    if not np.all(np.isfinite(core_eq)):
        raise RuntimeError('[ABORT] non-finite running integral of the core '
                           'equation (fail-loud).')
    print(f'[plot_core_eq_fictitious_vs_s2] core equation at Lambda^2 = s2_max = '
          f'{Lambda2[-1]:g}: {core_eq[-1]:.6f}  (of which {n_fict:g} ln(m_f^2) = '
          f'{n_fict * np.log(m2_f):.6f}); n = 0..{traces["n_max"]}')

    # the legend: the right side of the core equation alone
    label = (r'$\int_0^{\Lambda^2} ds^2 \left( \overline{G}(s^2) - '
             r'\overline{G}^{\,v}(s^2) - \frac{1}{\omega^2_{0,0}+s^2}'
             rf' - \frac{{{deg1:g}}}{{\omega^2_{{1,0}}+s^2}}'
             rf' + \frac{{{n_fict:g}}}{{m_f^2+s^2}} \right)'
             rf' + {n_fict:g}\,\ln(m_f^2)$')
    s2_max = float(traces['run']['band_integrals']['s2_max'])
    fig, ax = plt.subplots(figsize=(7.0, 4.4))     # a picture <= 8 in wide
    ax.plot(Lambda2, core_eq, '-', color=COLOR_TRUNCATED, marker='o',
            ms=POINT_SIZE, mew=0, label=label)
    ax.set_xscale('log')
    ax.grid(True, which='both')     # the one grid style of plot_common
    ax.set_xlabel(r'$\Lambda^2$')
    # the left side of the core equation
    ax.set_ylabel(r"$\overline{J}''(\Lambda^2) - \overline{J}''(0)$")
    ax.set_title(rf'Unrenormalized projected core-equation with cutoffs: '
                 rf'$n_{{max}} = {traces["n_max"]}$, $s^2_{{max}} = {s2_max:g}$')
    legend_outside(ax, 'below')
    save(fig, args.out_dir, picture_name(CODE_STEM, args))
    plt.close(fig)


if __name__ == '__main__':
    main()
