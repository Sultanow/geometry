"""One guided-search worker (one core): repeated stochastic hill climbs on the parameters of a family, uphill in the MaxSAT
score (the largest number of pairs of faces that can share an edge at once). Full scores are face-neighbourly polyhedra and
go through the same checks and output format as search.py (hits_<seed>.jsonl); stats_<seed>.json is compatible with
collect.py and adds the climb statistics (best score, histogram of final scores) for progress.py."""
import argparse, json, os, time, zlib, numpy as np, networkx as nx
from facenb.arrsat import build, describe
from facenb.guided import score_timed, score_anytime, score_genus, score_faces
from facenb.gfamilies import GFAMILIES
from facenb.verify2 import verify2
from facenb.combtype import incidence
from facenb.known import name_of
ap = argparse.ArgumentParser()
ap.add_argument('--family', required=True); ap.add_argument('--seed', type=int, required=True)
ap.add_argument('--seconds', type=float, default=3600); ap.add_argument('--out', required=True)
ap.add_argument('--eval-timeout', type=float, default=120); ap.add_argument('--iters', type=int, default=200)
ap.add_argument('--patience', type=int, default=30); ap.add_argument('--step', type=float, default=0.08)
ap.add_argument('--scorer', choices=['anytime', 'exact', 'genus', 'faces'], default='anytime')
ap.add_argument('--hop', type=float, default=0.5, help='probability to restart near the best point of this task instead of a random one')
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
rng = np.random.default_rng([a.seed, zlib.crc32(('guided' + a.family).encode())])
lo, hi, init, make = GFAMILIES[a.family]
st = dict(family=a.family, seed=a.seed, mode='guided', samples=0, usable=0, degenerate=0, timeouts=0, solutions=0, float_rejected=0,
          by_type={}, seconds=0, climbs=0, best=-10 ** 9, nmax=None, final_hist={})
hits = open(os.path.join(a.out, f'hits_{a.seed}.jsonl'), 'a'); t0 = time.time(); t_last = t0; kept = {}
def flush():
    st['seconds'] = round(time.time() - t0)
    with open(os.path.join(a.out, f'stats_{a.seed}.json'), 'w') as f: json.dump(st, f)
def evaluate(x, ctx):
    P = make(x, ctx); st['samples'] += 1
    if a.scorer == 'faces':
        try: s, ch, info = score_faces(P, a.eval_timeout)
        except AssertionError: s, ch, info = -1000, None, {}
        if not info.get('proven', False): st['timeouts'] += 1
        if s <= -1000: st['degenerate'] += 1
        else: st['usable'] += 1
        return P, s, ch
    if a.scorer == 'genus':
        try: s, ch, info = score_genus(P, a.eval_timeout)
        except AssertionError: s, ch, info = -999, None, {}
        if s <= -997: s = -100 if s != -999 else -1                            # no surface with all planes / degenerate
        if not info.get('proven', False): st['timeouts'] += 1
    elif a.scorer == 'anytime':
        try: s, ch, proven = score_anytime(P, a.eval_timeout)
        except AssertionError: s, ch, proven = -3, None, False                 # encoding inconsistency: never trust, skip
        if s >= 0 and not proven: st['timeouts'] += 1                          # best found so far, optimum not proven
    else:
        s, ch = score_timed(P, a.eval_timeout)
        if s == -2: st['timeouts'] += 1
    if s == -1: st['degenerate'] += 1
    elif s >= 0: st['usable'] += 1
    return P, s, ch
def record_hit(P, ch, x, ctx):
    A = build(P); d = describe(A, ch); st['solutions'] += 1
    r = verify2(P.tolist(), d['faces'], d['coords'])
    if not r['simple'] or r['pairs_bad']: st['float_rejected'] += 1; return
    key = nx.weisfeiler_lehman_graph_hash(incidence(d['faces']), node_attr='kind', iterations=4)
    st['by_type'][key] = st['by_type'].get(key, 0) + 1
    if kept.get(key, 0) < 3:
        kept[key] = kept.get(key, 0) + 1
        rec = dict(family=a.family, seed=a.seed, wl=key, known=name_of(d['faces']), V=d['V'], E=d['E'], F=d['F'], genus=d['genus'],
                   face_sizes=sorted(len(f) for f in d['faces']), pairs_sharing_two=sum(v == 2 for v in d['share'].values()),
                   planes=np.asarray(P).tolist(), params=dict(x=list(map(float, x)), ctx={k: v for k, v in ctx.items() if k not in ('noise', 'lo', 'hi')}),
                   faces=[[list(t) for t in f] for f in d['faces']])
        hits.write(json.dumps(rec) + '\n'); hits.flush()
best_start = None                                                             # (x, ctx) of the best climb end in this task
while time.time() - t0 < a.seconds:
    if best_start is not None and rng.random() < a.hop:
        x0, ctx = best_start; L0, H0 = ctx.get('lo', lo), ctx.get('hi', hi); x = np.clip(x0 + 4 * a.step * (H0 - L0) * rng.normal(size=len(x0)), L0, H0)   # hop near the best
    else:
        x, ctx = init(rng)
    L, H = ctx.get('lo', lo), ctx.get('hi', hi)
    P, s, ch = evaluate(x, ctx)
    if s < 0 and a.scorer not in ('genus', 'faces'): continue
    if s <= -100: continue
    n = len(P); nmax = n * (n - 1) // 2; st['nmax'] = nmax; st['climbs'] += 1
    step, stuck = a.step, 0
    for it in range(a.iters):
        if s >= nmax or time.time() - t0 > a.seconds: break
        y = np.clip(x + step * (H - L) * rng.normal(size=len(x)), L, H)
        Py, t, chy = evaluate(y, ctx)
        if t > s: x, s, P, ch, stuck = y, t, Py, chy, 0
        elif t == s: x, P, ch, stuck = y, Py, chy, stuck + 1
        else: stuck += 1
        if stuck > a.patience:
            step *= 0.5; stuck = 0
            if step < a.step / 16: break
    st['final_hist'][str(s)] = st['final_hist'].get(str(s), 0) + 1
    if s >= st['best'] and (s > 0 or a.scorer == 'faces'): best_start = (x.copy(), ctx)
    if s > st['best']:
        st['best'] = s
        json.dump(dict(score=s, nmax=nmax, planes=np.asarray(P).tolist(), x=list(map(float, x))), open(os.path.join(a.out, f'best_{a.seed}.json'), 'w'))
    if s >= nmax and ch is not None: record_hit(P, ch, x, ctx)
    if time.time() - t_last > 300: flush(); t_last = time.time()
flush(); hits.close()
