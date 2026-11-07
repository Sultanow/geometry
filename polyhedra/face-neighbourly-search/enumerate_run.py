"""Enumerate cubes of the complete candidate enumeration (facenb/enumerate2.py) and store one JSON per cube:
python enumerate_run.py --cubes 0,1,2 --out DIR [--n 9 --X 3]"""
import argparse, json, os, time
from facenb.enumerate2 import Enumerator2, add_cube_switches, canonical_map
ap = argparse.ArgumentParser(); ap.add_argument('--cubes', required=True); ap.add_argument('--out', required=True)
ap.add_argument('--n', type=int, default=9); ap.add_argument('--X', type=int, default=3)
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
E = Enumerator2(n=a.n, X=a.X); sw = add_cube_switches(E)
for cid in [int(x) for x in a.cubes.split(',')]:
    path = os.path.join(a.out, f'cube_{cid}.json')
    if os.path.exists(path): continue
    assum = [sw[i] if not (cid >> i) & 1 else -sw[i] for i in range(len(sw))]; found = set(); lab = 0; t0 = time.time()
    for ch in E.solutions(assumptions=assum):
        lab += 1; found.add(canonical_map(E.dual(ch)))
    json.dump(dict(cube=cid, n=a.n, X=a.X, labelled=lab, seconds=round(time.time() - t0, 1),
                   classes=[[list(map(list, f)) for f in F] for F in sorted(found)]), open(path + '.tmp', 'w'))
    os.replace(path + '.tmp', path); print(f'cube {cid}: {lab} labelled, {len(found)} classes, {time.time() - t0:.0f} s', flush=True)
