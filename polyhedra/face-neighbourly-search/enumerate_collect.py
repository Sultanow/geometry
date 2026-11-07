"""Merge the cube results: number of cubes done, labelled solutions, isomorphism classes of candidate maps with their
statistics (face sizes, pairs sharing two or three edges); cross-check against the main case (K9 embeddings with three
quadrilaterals, facenb/enumerate.py), whose candidates must all be contained. Writes candidates.json and summary.md."""
import argparse, glob, json, os, collections
from facenb.enumerate2 import canonical_map
ap = argparse.ArgumentParser(); ap.add_argument('--root', required=True); ap.add_argument('--crosscheck', action='store_true')
ap.add_argument('--cubes', type=int, default=256)
a = ap.parse_args()
files = glob.glob(os.path.join(a.root, 'cube_*.json')); allc = set(); lab = 0
for f in files:
    r = json.load(open(f)); lab += r['labelled']; allc |= {tuple(tuple(tuple(t) for t in g) for g in F) for F in r['classes']}
def stats(F):
    E = {}
    for p, f in enumerate(F):
        for i in range(len(f)): E.setdefault(frozenset((f[i], f[(i + 1) % len(f)])), set()).add(p)
    share = collections.Counter(frozenset(ps) for ps in E.values())
    return tuple(sorted(len(f) for f in F)), tuple(sorted(v for v in share.values() if v > 1))
lines = [f'# Candidate enumeration: {len(files)}/{a.cubes} cubes, {lab} labelled solutions, {len(allc)} classes', '']
if len(files) < a.cubes: lines.append('**incomplete** - not all cubes are done yet\n')
mult = collections.Counter(stats(F)[1] for F in allc)
lines += ['| multiplicities of the pairs sharing more than one edge | classes |', '|---|---|'] + [f'| {list(k)} | {v} |' for k, v in sorted(mult.items())]
if a.crosscheck:
    from facenb.enumerate import cubes as qcubes, enumerate_cubes, expand
    found, _, _ = enumerate_cubes(qcubes()); main = set()
    for F in found:
        for m in expand([tuple(f) for f in F]): main.add(canonical_map([[tuple(t) for t in g] for g in m]))
    lines += ['', f'cross-check: main case (three quadrilaterals) gives {len(main)} classes; contained in the complete enumeration: '
              f'{len(main & allc)}/{len(main)}' + (' OK' if main <= allc else ' **MISMATCH**')]
json.dump([[list(map(list, g)) for g in F] for F in sorted(allc)], open(os.path.join(a.root, 'candidates.json'), 'w'))
open(os.path.join(a.root, 'summary.md'), 'w').write('\n'.join(lines) + '\n'); print('\n'.join(lines))
