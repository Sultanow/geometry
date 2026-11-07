"""One search worker (one core): sample plane arrangements of a family until the time budget is spent; for every arrangement
find up to --max-solutions face-neighbourly polyhedra with the SAT construction, check each in floating point (Roest-Vigh
interval test), classify it (Weisfeiler-Lehman hash of the incidence structure + exact isomorphism against the known types)
and log it. Output, in --out: hits_<seed>.jsonl (one line per new type in this task, at most --keep per type) and
stats_<seed>.json (counts, refreshed every 5 minutes)."""
import argparse, json, os, time, zlib, numpy as np, networkx as nx
from facenb.arrsat import solve, describe
from facenb.families import FAMILIES, perturb
from facenb.verify2 import verify2
from facenb.combtype import incidence
from facenb.known import name_of
ap = argparse.ArgumentParser()
ap.add_argument('--family', required=True); ap.add_argument('--seed', type=int, required=True)
ap.add_argument('--seconds', type=float, default=3600); ap.add_argument('--out', required=True)
ap.add_argument('--max-solutions', type=int, default=4); ap.add_argument('--time-limit', type=float, default=60)
ap.add_argument('--keep', type=int, default=3)
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
rng = np.random.default_rng([a.seed, zlib.crc32(a.family.encode())])   # reproducible per (family, seed)
sampler = FAMILIES[a.family]
st = dict(family=a.family, seed=a.seed, samples=0, usable=0, degenerate=0, timeouts=0, solutions=0, float_rejected=0, by_type={}, seconds=0)
kept = {}; t0 = time.time(); t_last = t0
hits = open(os.path.join(a.out, f'hits_{a.seed}.jsonl'), 'a')
def flush():
    st['seconds'] = round(time.time() - t0)
    with open(os.path.join(a.out, f'stats_{a.seed}.json'), 'w') as f: json.dump(st, f)
while time.time() - t0 < a.seconds:
    P, par = sampler(rng); st['samples'] += 1
    try:
        A, sols = solve(P, max_solutions=a.max_solutions, max_rounds=400, time_limit=a.time_limit)
    except ValueError:
        try:
            P = perturb(P, rng, 1e-4); A, sols = solve(P, max_solutions=a.max_solutions, max_rounds=400, time_limit=a.time_limit)
        except ValueError:
            st['degenerate'] += 1; continue
    st['usable'] += 1; st['timeouts'] += int(A.get('timeout', False))
    for ch in sols:
        st['solutions'] += 1
        d = describe(A, ch)
        r = verify2(P.tolist(), d['faces'], d['coords'])
        if not r['simple'] or r['pairs_bad']:
            st['float_rejected'] += 1
            if st['float_rejected'] <= 20:
                rec = dict(planes=np.asarray(P).tolist(), faces=[[list(t) for t in f] for f in d['faces']], simple=r['simple'], bad_pairs=len(r['pairs_bad']))
                with open(os.path.join(a.out, f'rejected_{a.seed}.jsonl'), 'a') as fr: print(json.dumps(rec), file=fr)
            continue
        key = nx.weisfeiler_lehman_graph_hash(incidence(d['faces']), node_attr='kind', iterations=4)
        st['by_type'][key] = st['by_type'].get(key, 0) + 1
        if kept.get(key, 0) < a.keep:
            kept[key] = kept.get(key, 0) + 1
            rec = dict(family=a.family, seed=a.seed, wl=key, known=name_of(d['faces']), V=d['V'], E=d['E'], F=d['F'], genus=d['genus'],
                       face_sizes=sorted(len(f) for f in d['faces']), pairs_sharing_two=sum(v == 2 for v in d['share'].values()),
                       planes=np.asarray(P).tolist(), params=par, faces=[[list(t) for t in f] for f in d['faces']])
            hits.write(json.dumps(rec) + '\n'); hits.flush()
    if time.time() - t_last > 300: flush(); t_last = time.time()
flush(); hits.close()
