"""Controls that must pass before any search result is trusted. Prints ALL CONTROLS PASSED or stops at the first failure."""
import sys, numpy as np, networkx as nx
from fractions import Fraction
sys.path.insert(0, '.')
from facenb.arrsat import solve, describe
from facenb.combtype import same_type, incidence
from facenb.exact import verify_exact
from facenb.verify2 import verify2
from facenb.known import MIZHAEV_VERTICES, MIZHAEV_FACES, ROESTVIGH_PLANES, known_types
def check(cond, msg):
    print(('ok    ' if cond else 'FAIL  ') + msg)
    if not cond: sys.exit(1)
V = np.array(MIZHAEV_VERTICES, dtype=np.int64); pm = []
for f in MIZHAEV_FACES:
    nrm = np.cross(V[f[1]] - V[f[0]], V[f[2]] - V[f[0]]); nrm //= np.gcd.reduce(np.abs(nrm)); pm.append([*map(int, nrm), int(nrm @ V[f[0]])])
A, s = solve(pm, max_solutions=10); check(len(s) == 1, "Mizhaev's planes carry exactly one face-neighbourly polyhedron")
d = describe(A, s[0]); check(same_type(d['faces'], MIZHAEV_FACES), 'and it is Mizhaev\'s polyhedron')
r = verify_exact(pm, d['faces']); check(sorted(tuple(int(c) for c in x) for x in r['vertices'].values()) == sorted(map(tuple, V.tolist())), 'exact vertices = Mizhaev\'s integer vertices')
check(r['embedded'] and r['simple'] and r['genus'] == 3 and r['pairs_sharing_two'] == 8, 'Mizhaev passes the exact verification')
A, s = solve(ROESTVIGH_PLANES, max_solutions=10); check(len(s) == 1, "Roest-Vigh's planes carry exactly one")
r = verify_exact(ROESTVIGH_PLANES, describe(A, s[0])['faces'])
check(r['embedded'] and r['volume'] == Fraction(351105133, 18170460), 'Roest-Vigh passes, with the exact volume of their paper')
K = known_types()
autos = lambda fc: sum(1 for _ in nx.algorithms.isomorphism.GraphMatcher(incidence(fc), incidence(fc), node_match=lambda x, y: x['kind'] == y['kind']).isomorphisms_iter())
check(autos(K['Mizhaev']) == 4 and autos(K['Roest-Vigh']) == 8, 'automorphism group orders 4 and 8 as stated by Roest-Vigh')
check(not same_type(K['Mizhaev'], K['Roest-Vigh']), 'the two known types are distinguished')
flip = [f[::-1] if k % 2 else f for k, f in enumerate(MIZHAEV_FACES)]
check(same_type(flip, MIZHAEV_FACES), 'type comparison does not depend on face orientation')
P = np.array(pm, float); rng = np.random.default_rng(3); bad = 0
for _ in range(10):
    Q = P.copy(); Q[:, :3] += 0.25 * np.abs(Q[:, :3]).max() * rng.normal(size=(8, 3)) / 5
    coords = {t: np.linalg.solve(Q[list(t), :3], Q[list(t), 3]) for f in d['faces'] for t in f}
    rr = verify2(Q.tolist(), d['faces'], coords); bad += int((not rr['simple']) or bool(rr['pairs_bad']))
check(bad == 10, 'wrong geometries with the right combinatorics are rejected (10/10)')
print('ALL CONTROLS PASSED')
