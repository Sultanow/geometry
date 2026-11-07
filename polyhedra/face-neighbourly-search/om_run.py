"""Oriented-matroid SAT test of candidate maps: python om_run.py --candidates FILE --ids 0,1,... --out DIR [--limit SECONDS]
Per candidate k: cnf_k.cnf.gz (the full formula in DIMACS, for independent re-solving), result_k.json
(SAT / UNSAT / undecided, seconds, and for SAT the chirotope signs)."""
import argparse, gzip, json, os, time, multiprocessing as mp
from facenb.omsat import OMSat, add_map, add_reflex, reorientation_breaking
from pysat.solvers import Solver
ap = argparse.ArgumentParser(); ap.add_argument('--candidates', required=True); ap.add_argument('--ids', required=True)
ap.add_argument('--out', required=True); ap.add_argument('--limit', type=float, default=20 * 3600)
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
C = json.load(open(a.candidates))
def work(k, conn):
    faces = [[tuple(t) for t in f] for f in C[k]]; O = OMSat(len(faces)); O.grassmann_pluecker(); add_map(O, faces); add_reflex(O, faces); chosen = reorientation_breaking(O)
    with gzip.open(os.path.join(a.out, f'cnf_{k}.cnf.gz'), 'wt') as f:
        f.write(f'c candidate {k}; variables 1..{len(O.var)} are chi of the sorted 4-subsets in lexicographic order of {{0..{O.n}}} ({O.n} = plane at infinity)\n')
        f.write(f'p cnf {O.top} {len(O.cl)}\n'); f.writelines(' '.join(map(str, c)) + ' 0\n' for c in O.cl)
    s = Solver(name='cadical153', bootstrap_with=O.cl); t0 = time.time(); ok = s.solve(); dt = time.time() - t0
    out = dict(candidate=k, encoding='v2 (simple faces, no piercing, at most one reflex corner per vertex)', result='SAT' if ok else 'UNSAT', seconds=round(dt, 1), variables=O.top, clauses=len(O.cl), breaking=[list(S) for S in chosen])
    if ok:
        m = set(l for l in s.get_model() if l > 0); out['chi'] = {','.join(map(str, S)): (1 if v in m else -1) for S, v in O.var.items()}
    s.delete(); conn.send(out)
for k in [int(x) for x in a.ids.split(',')]:
    path = os.path.join(a.out, f'result_{k}.json')
    if os.path.exists(path): continue
    ctx = mp.get_context('fork'); r, w = ctx.Pipe(False); pr = ctx.Process(target=work, args=(k, w)); t0 = time.time(); pr.start(); w.close()
    out = r.recv() if r.poll(a.limit) else dict(candidate=k, result='undecided', seconds=round(time.time() - t0, 1), limit=a.limit)
    if pr.is_alive(): pr.kill()
    pr.join(); json.dump(out, open(path, 'w')); print(f"candidate {k}: {out['result']} ({out['seconds']} s)", flush=True)
