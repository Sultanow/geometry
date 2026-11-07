"""Collect the hits of all workers: statistics per family, combinatorial types (WL hash, then exact isomorphism), the known
types (tetrahedron, Szilassi, Mizhaev, Roest-Vigh) and, for every NEW type, an integer realisation verified in exact rational
arithmetic. Writes summary.md, types.json and, per new type, new_type_<k>.json into --out."""
import argparse, glob, json, os, numpy as np
from facenb.arrsat import solve, describe
from facenb.combtype import same_type
from facenb.exact import verify_exact
from facenb.known import name_of
from facenb.integerize import integerize
ap = argparse.ArgumentParser(); ap.add_argument('--root', required=True); ap.add_argument('--out', required=True)
ap.add_argument('--max-scale', type=int, default=400)
ap.add_argument('--no-exact', action='store_true', help='quick status: skip integer realisation and exact verification of new types')
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
stats, hits = {}, []
for f in glob.glob(os.path.join(a.root, '*', 'stats_*.json')):
    s = json.load(open(f)); t = stats.setdefault(s['family'], dict(tasks=0, samples=0, usable=0, degenerate=0, timeouts=0, solutions=0, float_rejected=0, cpu_h=0))
    t['tasks'] += 1; t['cpu_h'] += s['seconds'] / 3600
    for k in ('samples', 'usable', 'degenerate', 'timeouts', 'solutions', 'float_rejected'): t[k] += s[k]
for f in glob.glob(os.path.join(a.root, '*', 'hits_*.jsonl')):
    hits += [json.loads(l) for l in open(f) if l.strip()]
types = []                                                   # list of dict(rep=hit, members=[...])
for h in hits:
    fc = [[tuple(t) for t in f] for f in h['faces']]
    for T in types:
        if T['wl'] == h['wl'] and same_type(T['faces'], fc): T['members'].append(h); break
    else: types.append(dict(wl=h['wl'], faces=fc, rep=h, members=[h]))
lines = ['# Face-neighbourly search: summary', '', '| family | tasks | CPU h | samples | usable | degenerate | SAT timeouts | solutions | rejected in floats |', '|---|---|---|---|---|---|---|---|---|']
for fam, t in sorted(stats.items()):
    lines.append(f"| {fam} | {t['tasks']} | {t['cpu_h']:.1f} | {t['samples']} | {t['usable']} | {t['degenerate']} | {t['timeouts']} | {t['solutions']} | {t['float_rejected']} |")
lines += ['', '| type | faces | V | E | genus | face sizes | pairs sharing two | known as | hits (families) |', '|---|---|---|---|---|---|---|---|---|']
out_types = []; new_k = 0
for k, T in enumerate(sorted(types, key=lambda T: (T['rep']['F'], -len(T['members'])))):
    h = T['rep']; name = name_of(T['faces']); fams = {}
    for m in T['members']: fams[m['family']] = fams.get(m['family'], 0) + 1
    rec = dict(type=k, wl=T['wl'], known=name, F=h['F'], V=h['V'], E=h['E'], genus=h['genus'], face_sizes=h['face_sizes'],
               pairs_sharing_two=h['pairs_sharing_two'], families=fams, example=dict(planes=h['planes'], faces=h['faces']))
    if name is None:
        new_k += 1; res = 'skipped' if a.no_exact else integerize(h['planes'], T['faces'], a.max_scale)
        if res == 'skipped': rec['exact'] = 'not yet checked (status run)'; res = False
        if res is False: pass
        elif res is None: rec['exact'] = 'no integer realisation found up to the scale limit (float-verified only)'
        else:
            Q, fcs, r, K = res
            rec['exact'] = 'VERIFIED in exact rational arithmetic'
            cert = dict(planes=Q, faces=[[list(t) for t in f] for f in fcs], vertices={str(list(t)): [str(c) for c in x] for t, x in r['vertices'].items()},
                        report={kk: str(v) for kk, v in r.items() if kk != 'vertices'})
            json.dump(cert, open(os.path.join(a.out, f'new_type_{new_k}.json'), 'w'), indent=1)
    out_types.append(rec)
    lines.append(f"| {k} | {h['F']} | {h['V']} | {h['E']} | {h['genus']:.0f} | {h['face_sizes']} | {h['pairs_sharing_two']} | {name or '**NEW** — ' + rec['exact']} | "
                 + ', '.join(f'{f}: {c}' for f, c in sorted(fams.items())) + ' |')
json.dump(out_types, open(os.path.join(a.out, 'types.json'), 'w'), indent=1)
open(os.path.join(a.out, 'summary.md'), 'w').write('\n'.join(lines) + '\n'); print('\n'.join(lines))
