"""Controls for the oriented-matroid encoding. Prints OM CONTROLS PASSED or stops at the first failure.
Soundness (real valid polyhedra satisfy the encoding, also after symmetry breaking) is what makes UNSAT a proof."""
import sys, numpy as np
sys.path.insert(0, '.')
from facenb.omsat import OMSat, om_solve3 as om_solve2, chirotope_of, reorientation_breaking
from facenb.arrsat import solve, describe, build
from facenb.known import ROESTVIGH_PLANES, MIZHAEV_VERTICES, MIZHAEV_FACES
from facenb.families import FAMILIES
from facenb.gfamilies import GFAMILIES
from facenb.guided import score_genus
from facenb.verify2 import verify2
from facenb.combtype import same_type
from pysat.solvers import Solver
def check(cond, msg):
    print(('ok    ' if cond else 'FAIL  ') + msg)
    if not cond: sys.exit(1)
rng = np.random.default_rng(0)
O = OMSat(9); O.grassmann_pluecker(); good = bad = 0
for _ in range(3):
    chi = chirotope_of([[*rng.normal(size=3), rng.normal()] for _ in range(9)])
    s = Solver(name='cadical153', bootstrap_with=O.cl); good += s.solve(assumptions=[O.var[S] if v > 0 else -O.var[S] for S, v in chi.items()]); s.delete()
    s = Solver(name='cadical153', bootstrap_with=O.cl); bad += not s.solve(assumptions=[O.var[S] if rng.random() < .5 else -O.var[S] for S in O.var]); s.delete()
check(good == 3 and bad == 3, 'Grassmann-Pluecker: real chirotopes accepted (3/3), random sign vectors rejected (3/3)')
def valid(P, faces):
    coords = {t: np.linalg.solve(np.array(P)[list(t), :3], np.array(P)[list(t), 3]) for g in faces for t in g}; r = verify2(np.array(P).tolist(), faces, coords)
    return r['simple'] and not r['pairs_bad']
def reorient(chi, chosen, n):
    for mask in range(1 << n):
        f = [(mask >> e) & 1 for e in range(n)] + [0]; new = {S: v * (-1) ** sum(f[e] for e in S) for S, v in chi.items()}
        if all(new[S] > 0 for S in chosen): return new
cases = []
A, s = solve(ROESTVIGH_PLANES, max_solutions=1); cases.append(('Roest-Vigh', np.array(ROESTVIGH_PLANES, float), describe(A, s[0])['faces']))
V = np.array(MIZHAEV_VERTICES, dtype=np.int64); pm = []
for f in MIZHAEV_FACES:
    nrm = np.cross(V[f[1]] - V[f[0]], V[f[2]] - V[f[0]]); pm.append([*map(int, nrm), int(nrm @ V[f[0]])])
A, s = solve(pm, max_solutions=1); cases.append(('Mizhaev', np.array(pm, float), describe(A, s[0])['faces']))
while True:
    P, _ = FAMILIES['C2_7'](rng); P = P + 1e-3 * rng.normal(size=P.shape)
    try: A, s = solve(P, max_solutions=1)
    except ValueError: continue
    if s and all(v is not None for v in chirotope_of(P.tolist()).values()): cases.append(('Szilassi', P, describe(A, s[0])['faces'])); break
k = 0
while k < 2:                                                      # nine-plane embedded polyhedra (genus 3) from the growth start
    lo, hi, init, make = GFAMILIES['GROW_R9']; x, ctx = init(rng); P = make(x, ctx) + 1e-4 * rng.normal(size=(9, 4))
    sc, ch, info = score_genus(P, 30)
    if ch is None or sc < 0: continue
    faces = describe(build(P), ch)['faces']
    if valid(P, faces) and all(v is not None for v in chirotope_of(P.tolist()).values()): cases.append((f'nine planes #{k + 1}', P, faces)); k += 1
for name, P, faces in cases:
    n = len(faces); ok_free, O, _ = om_solve2(faces); chi = chirotope_of(P.tolist())
    chi2 = reorient(chi, reorientation_breaking(OMSat(n)), n)
    ok_real, _, _ = om_solve2(faces, assumptions=[O.var[S] if v > 0 else -O.var[S] for S, v in chi2.items()])
    rej = tested = 0
    for t in range(6):
        Q = P + 0.25 * np.abs(P).max() * rng.normal(size=P.shape)
        if valid(Q, faces): continue
        chq = chirotope_of(Q.tolist())
        if any(v is None for v in chq.values()): continue
        tested += 1; okq, _, _ = om_solve2(faces, assumptions=[O.var[S] if v > 0 else -O.var[S] for S, v in chq.items()], breaking=False); rej += not okq
    check(ok_free and ok_real and rej == tested and tested > 0, f'{name}: SAT; its real chirotope (reoriented) satisfies the encoding (incl. reflex condition) with symmetry breaking; invalid geometries rejected {rej}/{tested}')
print('OM CONTROLS PASSED')
