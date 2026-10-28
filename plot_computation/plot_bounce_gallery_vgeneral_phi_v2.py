#!/usr/bin/env python3
"""plot_bounce_gallery_vgeneral_phi_v2.py -- the bounce PLOT stage: render
ALL computed bounce solutions of a run, not only the production pick -- the
gallery of every bounce of the potential, the same pictures as the selection
gallery of stage 0 from which the bounce of a run is chosen by its name.

ONE RENDERER: the pictures are drawn by plot_bounce of
computation/bounce_gallery_vgeneral_phi_v2.py, the same renderer the stage-0
selection gallery uses; this code rebuilds the vacuum catalogue, resolves
the labels and calls plot_bounce for every bounce file of the run, with
its own picture names and folder.

WHAT IT DRAWS, per bounce npz (two PNGs; one for N=1 and N>3): the
PROFILES picture -- the bounce profiles phi'(r) of the bounce from the true
vacuum to the false vacuum, in the primed fieldspace frame: the channels
phi'_i(r) (phi'_1 the heaviest false-vacuum mass channel) against the
Euclidean radius r, no title (the bounce NAME is in the file name),
the legend outside at the top right holding only the tree-level action S_E
and one colour per channel (N=1: the two vacua are dotted levels,
labelled at the right end, not legend entries) -- and the POTENTIAL
picture, in the original fieldspace frame: N=2 the 3-D surface V with the
bounce path (red) from the true vacuum to the false vacuum, next to the
landscape -V (the inverted potential, a contour map) with the same path;
N=3 the field-space curve (V
as color), V without its constant term c; no panel titles.  The potential
axes are the declared fields, named phi_1, phi_2, ... in their declared
order.  The census eigenvalues are printed on the terminal of every run;
SPECIAL bounces are dashed + [SPECIAL]-marked (never enter lnD).

LABELS (bounce_vgeneral_phi_v2 is the one owner of the label functions AND
of catalogue_indices; names are <falseLabel>_<trueLabel> from the global
vacuum-catalogue labels TV0/M<i>/FV0): resolved by COORDINATES against the
current catalogue, never by the stored indices.

IN   the run's own bounce npz (--bounce-npz, auto-discovered in --data-dir
     by plot_common.resolve_run): it names the figure folder and gives the
     potential (params, fields) and the vacuum-search box (seed_range, unless
     --seed-range is given).  Then every bounce_data_*.npz in --gen-dir
     (default <data-dir>/bounce_gen, falling back to <data-dir> itself);
     files of a different potential in that dir are reported and skipped.
     Schema: bounce_vgeneral_phi_v2.bounce_npz_payload (every file carries
     the census; a file written before the current census schema, without
     its census keys, is drawn like an accepted bounce -- nothing marks it
     in the picture -- and listed as OLDER SCHEMA in the terminal table).
     A data dir whose selection was refused holds no run bounce npz: its
     pictures are the stage-0 gallery <data-dir>/gallery/.
OUT  TWO pictures per bounce solution,
     plot_bounce_gallery_vgeneral_phi_v2_<NAME>_profiles.png and
     plot_bounce_gallery_vgeneral_phi_v2_<NAME>_potential.png (N = 2, 3)
     (NAME = the bounce name, e.g. M2_TV0 -- figure names never carry the
     run tag/date), into --out-dir (default: pictures/pictures_<potential>_
     <bounce>/bounce_gallery/ inside the pipeline,
     plot_common.figure_folder), plus the terminal summary table (name, S_E,
     census verdict) by S_E.
"""
import argparse
import glob
import os
import sys

import numpy as np

sys.dont_write_bytecode = True
os.environ.setdefault('MPLBACKEND', 'Agg')
_HERE = os.path.dirname(os.path.abspath(__file__))
_COMP = os.path.normpath(os.path.join(_HERE, '..', 'computation'))
if _COMP not in sys.path:
    sys.path.insert(0, _COMP)

from bounce_vgeneral_phi_v2 import (
    catalogue_indices,
    minima_catalogue,
    same_potential,
)
from bounce_gallery_vgeneral_phi_v2 import plot_bounce
from plot_common_vgeneral_phi_v2 import add_standard_plot_cli, resolve_run


def load_bounces(gen_dir):
    """All bounce_data_*.npz of gen_dir as plain dicts, sorted by filename."""
    files = sorted(glob.glob(os.path.join(gen_dir, 'bounce_data_*.npz')))
    out = []
    for f in files:
        with np.load(f, allow_pickle=True) as bounce_npz_data:
            b = {k: bounce_npz_data[k] for k in bounce_npz_data.files}
        b['_file'] = f
        out.append(b)
    return out


def main():
    ap = argparse.ArgumentParser(
        description='vgeneral_phi_v2 bounce PLOT stage: render ALL computed '
                    'bounce solutions (profiles + surface/field-space path, '
                    'census-marked)')
    add_standard_plot_cli(ap)
    ap.add_argument('--gen-dir', default=None,
                    help='folder of bounce_data_*.npz (default: '
                         '<data-dir>/bounce_gen when it exists, else '
                         '<data-dir> itself)')
    ap.add_argument('--seed-range', default=None, metavar='LO,HI',
                    help='vacuum-search seed box of the SOLVE (default: '
                         'the box bounce_vgeneral_phi_v2 stored in the bounce '
                         'npz; the catalogue labels must come from the same '
                         'box); write a negative LO as --seed-range=LO,HI, '
                         'else it is read as an option')
    args = ap.parse_args()
    # the run's bounce npz and, from it, the folder of the gallery figures
    resolve_run(args, 'bounce_gallery')

    seed_range = None
    if args.seed_range is not None:
        seed_range = tuple(float(v) for v in args.seed_range.split(','))

    gen_dir = args.gen_dir
    if gen_dir is None:
        cand = os.path.join(args.data_dir, 'bounce_gen')
        gen_dir = cand if os.path.isdir(cand) else args.data_dir
    bounces = load_bounces(gen_dir)
    if not bounces:
        raise SystemExit(f'[ABORT GALLERY-PLOT] no bounce_data_*.npz in '
                         f'{gen_dir} -- run bounce_vgeneral_phi_v2 first.')

    # the run's potential comes FROM the run's bounce npz, the file that
    # names the figure folder (no coupling re-entry); bounce files of another
    # potential are skipped below.
    with np.load(args.bounce_npz, allow_pickle=True) as run_bounce:
        params = np.asarray(run_bounce['params'], float)
        fields = tuple(str(v) for v in np.asarray(run_bounce['fields']).ravel())
        stored_box = np.asarray(run_bounce.get('seed_range', ()), float)
    # no --seed-range given: the box the solve stored (an empty or absent
    # seed_range = the box derived from the couplings)
    if seed_range is None and stored_box.size == 2:
        seed_range = tuple(stored_box)
    minima = minima_catalogue(params, fields, seed_range=seed_range)
    print(f'[gallery-plot] {len(bounces)} bounce file(s) in {gen_dir}; '
          f'fields {fields}, {len(minima)} minima; figures -> {args.out_dir}')

    gallery_rows = []
    for b in bounces:
        fname = os.path.basename(b['_file'])
        if not same_potential(b, params, fields):
            print(f'  [gallery-plot] SKIPPING {fname}: solved for a DIFFERENT '
                  f'potential than the run bounce {args.bounce_npz}')
            continue
        iF, iT = catalogue_indices(minima,
                                   np.asarray(b['false_vac'], float),
                                   np.asarray(b['true_vac'], float), fname)
        # the pictures carry THIS plot code's name plus the bounce name
        # alone -- no run tag and no date ever appears in a figure name
        gallery_rows.append(plot_bounce(
            b, minima, iF, iT, fields, params, args.out_dir,
            png_prefix='plot_bounce_gallery_vgeneral_phi_v2_'))

    if not gallery_rows:
        raise SystemExit(f'[ABORT GALLERY-PLOT] no bounce file in {gen_dir} '
                         f'belongs to the potential of the run bounce '
                         f'{args.bounce_npz} -- select the run with '
                         f'--data-dir/--tag/--bounce-npz.')
    gallery_rows.sort(key=lambda entry: entry[1])   # by tree-level action S_E, the smallest first
    print('\nALL computed bounce solutions (by tree-level action):')
    for name, S_E, verdict, _pngs in gallery_rows:
        print(f'  {name:<12s}  S_E = {S_E:12.4f}   census: {verdict}')


if __name__ == '__main__':
    main()
