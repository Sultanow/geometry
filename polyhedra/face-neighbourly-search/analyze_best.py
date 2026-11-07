"""Analyse the best configurations of a guided run (best_<seed>.json with score >= --min-score): recompute the surface with
a longer limit, describe it (V, E, genus, face sizes, pairs sharing two edges), list the MISSING pairs (involving the added
plane, index >= 8, or two of the original eight), check it independently (Roest-Vigh interval test) and exactly (rounded
integer planes, rational arithmetic). Writes analysis.jsonl and summary.md into --out."""
import argparse, glob, json, os, itertools, numpy as np
from collections import Counter
from multiprocessing import Pool
from facenb.arrsat import build, describe
from facenb.guided import score_genus, score_faces
from facenb.verify2 import verify2
from facenb.exact import verify_exact
ap = argparse.ArgumentParser()
ap.add_argument('--root', required=True); ap.add_argument('--out', required=True); ap.add_argument('--min-score', type=int, default=34)
ap.add_argument('--timeout', type=float, default=120); ap.add_argument('--scorer', choices=['genus', 'faces'], default='genus'); ap.add_argument('--workers', type=int, default=8); ap.add_argument('--max-scale', type=int, default=300)
a = ap.parse_args()
def exact_check(P, faces):
    P = np.array(P, float); P = P / np.abs(P).max(axis=1, keepdims=True)
    for K in range(8, a.max_scale + 1, 4):
        Q = np.rint(K * P).astype(int)
        if any(not q[:3].any() for q in Q): continue
        try: r = verify_exact(Q.tolist(), faces)
        except Exception: continue
        if all(r[k] for k in ('no_parallel', 'general_position', 'simple', 'closed', 'orientable', 'manifold_links', 'embedded')):
            return dict(scale=K, planes=Q.tolist(), genus=str(r['genus']), V=r['V'], E=r['E'])
    return None
def work(path):
    b = json.load(open(path)); P = b['planes']; n = len(P)
    s, ch, info = (score_faces if a.scorer == 'faces' else score_genus)(P, a.timeout)
    rec = dict(file=os.path.relpath(path, a.root), stored=b['score'], rescored=s, proven=info.get('proven'))
    if ch is None or s < 0: return rec
    A = build(P); d = describe(A, ch); share = {tuple(sorted(k)): v for k, v in d['share'].items()}
    missing = [pq for pq in itertools.combinations(range(n), 2) if share.get(pq, 0) == 0]
    r = verify2(P, d['faces'], d['coords'])
    rec.update(V=d['V'], E=d['E'], genus=d['genus'], face_sizes=[len(f) for f in d['faces']], pairs=n * (n - 1) // 2 - len(missing),
               missing=missing, missing_kind=['new-old' if max(pq) >= 8 else 'old-old' for pq in missing],
               doubled=sum(v == 2 for v in share.values()), float_ok=bool(r['simple'] and not r['pairs_bad']),
               exact=exact_check(P, d['faces']), planes=P, faces=[[list(t) for t in f] for f in d['faces']])
    return rec
if __name__ == '__main__':
    os.makedirs(a.out, exist_ok=True)
    files = [f for f in glob.glob(os.path.join(a.root, '*', 'best_*.json')) if json.load(open(f))['score'] >= a.min_score]
    print(f'{len(files)} configurations with score >= {a.min_score}')
    with Pool(a.workers) as pool: recs = pool.map(work, files)
    with open(os.path.join(a.out, 'analysis.jsonl'), 'w') as f:
        for r in recs: f.write(json.dumps(r) + '\n')
    ok = [r for r in recs if 'missing' in r]
    lines = [f'# Analysis of {len(files)} configurations (score >= {a.min_score})', '',
             f'- recomputed: {len(ok)}; independently verified (floats): {sum(r["float_ok"] for r in ok)}; exactly verified (integer planes): {sum(r["exact"] is not None for r in ok)}',
             f'- genus: {sorted({r["genus"] for r in ok})}; pairs: {sorted({r["pairs"] for r in ok})}', '',
             '| kinds of the missing pairs | count |', '|---|---|']
    for k, v in Counter(tuple(sorted(r['missing_kind'])) for r in ok).most_common(): lines.append(f'| {", ".join(k) or "none"} | {v} |')
    lines += ['', '| file | pairs | V | E | genus | face sizes | doubled | missing | float | exact |', '|---|---|---|---|---|---|---|---|---|---|']
    for r in sorted(ok, key=lambda r: -r['pairs']):
        lines.append(f"| {r['file']} | {r['pairs']} | {r['V']} | {r['E']} | {r['genus']:.0f} | {r['face_sizes']} | {r['doubled']} | {r['missing']} | {r['float_ok']} | {'K=' + str(r['exact']['scale']) if r['exact'] else '-'} |")
    open(os.path.join(a.out, 'summary.md'), 'w').write('\n'.join(lines) + '\n'); print('\n'.join(lines))
