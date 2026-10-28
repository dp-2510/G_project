#!/usr/bin/env python3
"""plot_toolbox_v2.py -- the shared toolbox of the plot codes: locate the
npz written by the matching add_on_greensfct code (same CLI), and save every
figure into the run's folder pictures/pictures_<potential>_<bounce>/
add_on_greensfct/ inside the pipeline (plot_common.figure_folder)."""
import glob
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.ticker import MaxNLocator

_HERE = os.path.dirname(os.path.abspath(__file__))
_ADDON = os.path.join(os.path.dirname(_HERE), 'add_on_greensfct')
_COMP = os.path.join(os.path.dirname(_HERE), 'computation')
_PLOTS = os.path.join(os.path.dirname(_HERE), 'plot_computation')
for _p in (_COMP, _ADDON, _PLOTS):
    if _p not in sys.path:
        sys.path.insert(0, _p)
# s2_filetag / bounce_name: the pipeline's two filename pieces
# (pipeline_helpers_vgeneral_phi_v2), re-exported here so the plot codes
# import from this toolbox only
from pipeline_helpers_vgeneral_phi_v2 import bounce_name, s2_filetag
from plot_common_vgeneral_phi_v2 import figure_folder, thick_titles


def add_plot_cli(parser, need_n=False):
    parser.add_argument('--s2', type=float, required=True,
                        help='the s^2 of the npz to plot (same as the '
                             'compute code)')
    if need_n:
        parser.add_argument('--n', type=int, required=True,
                            help='the wave n of the npz to plot')
    parser.add_argument('--data-dir',
                        default=os.environ.get('VGENERAL_PHI_DATA', '.'),
                        help='where the compute code wrote its npz')
    parser.add_argument('--tag', default=None,
                        help='run tag in the npz filenames (default: derived '
                             'from the unique bounce_data_<tag>.npz in '
                             '--data-dir)')
    parser.add_argument('--out-dir', default=None,
                        help='figure folder (default: pictures/pictures_'
                             '<potential>_<bounce>/add_on_greensfct/ inside '
                             'the pipeline)')
    return parser


def resolve_tag(args):
    """Fill a missing --tag from the data dir: the unique bounce_data_<tag>
    .npz wins; ambiguity aborts listing the candidates (never guessed)."""
    if getattr(args, 'tag', None) is not None:
        return args
    hits = sorted(glob.glob(os.path.join(args.data_dir, 'bounce_data_*.npz')))
    tags = sorted({re.match(r'bounce_data_(.+)\.npz$',
                            os.path.basename(h)).group(1) for h in hits})
    if len(tags) == 1:
        args.tag = tags[0]
        print(f'[TAG] auto-discovered run tag: {args.tag}')
        return args
    raise SystemExit('[ABORT] --tag not given and it cannot be derived '
                     f'uniquely from {args.data_dir!r} (candidates: {tags}) '
                     '-- pass --tag.')


def closest_stored_s2(data_dir, stem_prefix, s2):
    """The stored s^2 tags of every npz `<stem_prefix>_s2<tag>*.npz` in
    data_dir and the one closest to `s2` -- what the abort of a missing
    s^2-tagged file reports, so the user re-runs with a value that exists.
    Returns (closest or None, sorted list of stored values)."""
    stored = []
    for path in glob.glob(os.path.join(data_dir, f'{stem_prefix}_s2*.npz')):
        tag = os.path.basename(path).split('_s2', 1)[1].split('.npz')[0]
        try:
            stored.append(float(tag.replace('p', '.')))
        except ValueError:
            continue
    stored = sorted(set(stored))
    if not stored:
        return None, stored
    return min(stored, key=lambda v: abs(v - s2)), stored


def load_npz(data_dir, stem):
    path = os.path.join(data_dir, stem)
    if not os.path.isfile(path):
        hint = ''
        if '_s2' in stem:
            prefix = stem.split('_s2', 1)[0]
            want = stem.split('_s2', 1)[1].split('.npz')[0]
            try:
                near, stored = closest_stored_s2(data_dir, prefix,
                                                 float(want.replace('p', '.')))
            except ValueError:
                near, stored = None, []
            hint = (f' s2-value not available; closest stored s2 = {near:g} '
                    f'(stored: {", ".join(f"{v:g}" for v in stored)}).'
                    if near is not None else
                    ' s2-value not available; no stored s2 of this code yet.')
        raise RuntimeError(f'[ABORT] {path} not found.{hint} Run the matching '
                           f'add_on_greensfct compute code first.')
    print(f'[READ] {path}')
    return np.load(path, allow_pickle=True)


def savefig(fig, args, name):
    """Write one figure into --out-dir, by default the run's folder
    pictures/pictures_<potential>_<bounce>/add_on_greensfct/ named by the
    run's bounce npz (plot_common.figure_folder), and close it."""
    out_dir = args.out_dir or figure_folder('add_on_greensfct',
                                            run_bounce_npz(args))
    thick_titles(fig)                    # every title thicker (plot_common)
    os.makedirs(out_dir, exist_ok=True)
    out = os.path.join(out_dir, name)
    fig.savefig(out, dpi=160, bbox_inches='tight')
    plt.close(fig)
    print(f'[WROTE] {out}')


def surface_grid(r, rp, M, suptitle, limit,
                 title_fmt=r'$[G(r,r\prime)]_{{ij}}$'):
    """An N x N grid of surface panels (mpl_toolkits Axes3D,
    ax.plot_surface) of the matrix-valued M[a, b, i, j] over x = r (first
    axis of M) and y = r' (second axis), one panel per component (i, j),
    every grid point drawn.  All panels share ONE scale: the height, the z
    axis and the colour run up to `limit` in every panel (a component
    larger than it is cut flat there, exactly where its colour saturates),
    from 0 when no entry of M is negative, else from -limit; ONE colour bar
    right of the grid, in the middle, reads off the heights (the z axes
    carry no ticks).  The colours are RdBu_r centred on 0 (white = 0).  NaN
    entries (e.g. a masked diagonal) leave a gap in the surface.

    title_fmt is a template for the panel titles: the marker {ij} is
    replaced by the two channel numbers of the panel's component.  Returns
    the figure."""
    N = M.shape[-1]
    if not limit > 0:                  # an identically zero matrix: drawn
        limit = 1.0                    # flat and white
    lower = -limit if np.nanmin(M) < 0 else 0.0
    cmap = plt.get_cmap('RdBu_r')
    norm = mcolors.Normalize(-limit, limit)
    X, Y = np.meshgrid(r, rp, indexing='ij')                  # X=r, Y=r'
    # the layout in fractions of the figure, the same for every N: panel
    # rectangles 0.94/N wide on a column pitch 0.92/N and 0.88/N high on a
    # row pitch 1/N, each 3-D box zoomed by 1.1 to fill its rectangle; the
    # colour bar in the column right of the grid; <= 8 in wide at N = 2
    # (sizes: plot_common, the one set)
    fig = plt.figure(figsize=(3.725 * N, 3.525 * N))
    for i in range(N):
        for j in range(N):
            ax = fig.add_axes([j * 0.92 / N, (N - 1 - i) / N,
                               0.94 / N, 0.88 / N], projection='3d')
            # height AND colour are the component value, both clipped to
            # the shared scale; NaN heights are simply not drawn
            z = np.clip(M[:, :, i, j], lower, limit)
            ax.plot_surface(X, Y, z, facecolors=cmap(norm(z)), rstride=1,
                            cstride=1, linewidth=0, antialiased=False,
                            shade=False)
            ax.set_zlim(lower, limit)
            ax.set_zticks([])
            ax.set_xlabel(r'$r$')
            ax.set_ylabel(r"$r'$")
            ax.set_title(title_fmt.replace('{ij}', f'{i + 1}{j + 1}'),
                         y=0.97)
            ax.set_box_aspect(None, zoom=1.1)
            # no background: a neighbour never hides this panel's labels
            ax.patch.set_visible(False)
    top = 1.0 - 1.0 / N + 0.88 / N             # top edge of the panel grid
    cax = fig.add_axes([0.945, top / 2 - 0.36 / N, 0.022, 0.72 / N])
    bar = fig.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=norm), cax=cax)
    bar.ax.set_ylim(lower, limit)              # the range the panels use
    bar.locator = MaxNLocator(nbins=4)
    bar.update_ticks()
    # the figure title 0.19 in above the panel grid (any N); anchored at its
    # BOTTOM, so a title of several lines grows upwards
    fig.suptitle(suptitle, y=top + 0.054 / N, va='bottom')
    # the saved picture is cut to its tight bounding box, which ignores the
    # labels of 3-D axes: an invisible mark below the figure's bottom edge
    # keeps the bottom row's r labels inside the picture
    fig.text(0.5, -0.02 / N, '.', alpha=0.0)
    return fig


def trim_edge(r, frac=0.04):
    """Index from which to plot: the first ~4% of the grid sits in the
    K_nu(kappa r)/r blow-up (~ r^-(nu+1) as r -> 0) next to the engine's
    start radius r_start = R[1]/2
    (greens_function_vgeneral_phi_v2), where the partial-wave remainder is far from converged and would drown
    every panel's y-scale.  The stored npz always keeps the full grid."""
    return max(3, int(frac * len(r)))


def run_bounce_npz(args):
    """The run's bounce npz, bounce_data_<tag>.npz in --data-dir: it names
    every figure (bounce_label) and the figure folder (savefig)."""
    return os.path.join(args.data_dir, f'bounce_data_{args.tag}.npz')


def bounce_label(args):
    """The BOUNCE NAME for figure filenames (e.g. 'M2_TV0'): the run's
    bounce npz (run_bounce_npz), read through pipeline_helpers.bounce_name
    (key pair_name).  Figure names carry the bounce name + the (n, s^2)
    tags -- enough to tell two parameter choices apart -- never the run
    date.  Fail-loud: no bounce npz, or one without pair_name, aborts."""
    p = run_bounce_npz(args)
    if not os.path.isfile(p):
        raise SystemExit(f'[ABORT] bounce npz {p} not found -- the figure '
                         f'name needs its pair_name.')
    with np.load(p, allow_pickle=True) as bounce_npz_data:
        return bounce_name(bounce_npz_data)
