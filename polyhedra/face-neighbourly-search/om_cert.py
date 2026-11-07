"""Certified decision of candidate maps: python om_cert.py --candidates FILE --ids 0,1 --out DIR --drat-trim PATH
Per candidate: v2 encoding (simple faces, no piercing, at most one reflex corner per vertex, reorientation breaking) solved
with a DRAT proof that drat-trim checks independently; the alternative encoding (other ray lines) solved as a cross-check.
Writes cert_<k>.json and the formula and the proof as cnf_<k>.cnf.gz, proof_<k>.drat.gz."""
import argparse, gzip, json, os, shutil, subprocess, time
from pysat.solvers import Solver
from facenb.omsat import OMSat, add_map, add_reflex, add_map_alt, add_reflex_alt, reorientation_breaking
ap = argparse.ArgumentParser(); ap.add_argument('--candidates', required=True); ap.add_argument('--ids', required=True)
ap.add_argument('--out', required=True); ap.add_argument('--drat-trim', required=True); ap.add_argument('--tmp', default=os.environ.get('TMPDIR', '/tmp'))
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True); C = json.load(open(a.candidates))
def build(faces, alt=False):
    O = OMSat(len(faces)); O.grassmann_pluecker()
    (add_map_alt if alt else add_map)(O, faces); (add_reflex_alt if alt else add_reflex)(O, faces); reorientation_breaking(O); return O
for k in [int(x) for x in a.ids.split(',')]:
    path = os.path.join(a.out, f'cert_{k}.json')
    if os.path.exists(path): continue
    faces = [[tuple(t) for t in f] for f in C[k]]; rec = dict(candidate=k)
    O = build(faces); cnf = os.path.join(a.tmp, f'cnf_{k}.cnf'); prf = os.path.join(a.tmp, f'proof_{k}.drat')
    with open(cnf, 'w') as f: f.write(f'p cnf {O.top} {len(O.cl)}\n'); f.writelines(' '.join(map(str, c)) + ' 0\n' for c in O.cl)
    t0 = time.time(); s = Solver(name='cadical153', bootstrap_with=O.cl, with_proof=True); ok = s.solve(); rec['seconds'] = round(time.time() - t0, 1)
    rec['result'] = 'SAT' if ok else 'UNSAT'
    if not ok:
        with open(prf, 'w') as f: f.writelines(l + '\n' for l in s.get_proof())
        t1 = time.time(); r = subprocess.run([a.drat_trim, cnf, prf, '-t', '20000'], capture_output=True, text=True)
        rec['drat_trim'] = 'VERIFIED' if 's VERIFIED' in r.stdout else ('NOT VERIFIED: ' + r.stdout[-300:]); rec['check_seconds'] = round(time.time() - t1, 1)
        for src, dst in ((cnf, f'cnf_{k}.cnf.gz'), (prf, f'proof_{k}.drat.gz')):
            with open(src, 'rb') as fi, gzip.open(os.path.join(a.out, dst), 'wb') as fo: shutil.copyfileobj(fi, fo)
            os.remove(src)
    else:
        m = set(l for l in s.get_model() if l > 0); rec['chi'] = {','.join(map(str, S)): (1 if v in m else -1) for S, v in O.var.items()}
        with open(cnf, 'rb') as fi, gzip.open(os.path.join(a.out, f'cnf_{k}.cnf.gz'), 'wb') as fo: shutil.copyfileobj(fi, fo)
        os.remove(cnf)
    s.delete()
    O2 = build(faces, alt=True); s2 = Solver(name='cadical153', bootstrap_with=O2.cl); t2 = time.time(); rec['alt_result'] = 'SAT' if s2.solve() else 'UNSAT'
    rec['alt_seconds'] = round(time.time() - t2, 1); s2.delete()
    json.dump(rec, open(path, 'w')); print(json.dumps({x: rec[x] for x in rec if x != 'chi'}), flush=True)
