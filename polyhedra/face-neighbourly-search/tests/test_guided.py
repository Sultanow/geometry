"""Controls for the guided search. Prints GUIDED CONTROLS PASSED or stops at the first failure."""
import sys, numpy as np
sys.path.insert(0, '.')
from facenb.guided import score_conn, score_timed
from facenb.gfamilies import GFAMILIES
from facenb.arrsat import describe
from facenb.combtype import same_type
from facenb.known import ROESTVIGH_PLANES, MIZHAEV_VERTICES, MIZHAEV_FACES, known_types
def check(cond, msg):
    print(('ok    ' if cond else 'FAIL  ') + msg)
    if not cond: sys.exit(1)
V = np.array(MIZHAEV_VERTICES, dtype=np.int64); pm = []
for f in MIZHAEV_FACES:
    nrm = np.cross(V[f[1]] - V[f[0]], V[f[2]] - V[f[0]]); nrm //= np.gcd.reduce(np.abs(nrm)); pm.append([*map(int, nrm), int(nrm @ V[f[0]])])
s, (A, ch) = score_conn(pm); check(s == 28 and same_type(describe(A, ch)['faces'], MIZHAEV_FACES), "score 28/28 on Mizhaev's planes, and the optimum is Mizhaev's polyhedron")
s, (A, ch) = score_conn(ROESTVIGH_PLANES); check(s == 28 and same_type(describe(A, ch)['faces'], known_types()['Roest-Vigh']), "score 28/28 on Roest-Vigh's planes, optimum is theirs")
P = np.array(pm, float); P[:, 3] += 0.3 * np.abs(P[:, 3]).max() * np.array([1, -1, 1, -1, 1, -1, 1, -1]); s, _ = score_conn(P)
check(0 <= s < 28, f'a distorted arrangement scores below 28 (got {s})')
check(score_timed(ROESTVIGH_PLANES, 30)[0] == 28 and score_timed(ROESTVIGH_PLANES, 0.005)[0] == -2, 'evaluation in a child process: result and timeout')
rng = np.random.default_rng(0)
for name, (lo, hi, init, make) in GFAMILIES.items():
    if name.startswith('SEED'): continue                     # needs the analysis file of a previous run (checked by submit_G.sh)
    x, ctx = init(rng); P = make(x, ctx); n = int(__import__('re').findall(r'\d+', name)[-1])
    L, H = ctx.get('lo', lo), ctx.get('hi', hi); check(len(P) == n and np.all(x >= L) and np.all(x <= H), f'family {name}: {n} planes, parameters inside their box')

from facenb.guided import score_genus, target_genus
check([target_genus(n) for n in (7, 8, 9, 10, 11, 12)] == [1, 2, 3, 4, 5, 6], 'target genus 1..6 for 7..12 faces')
s, _, info = score_genus(ROESTVIGH_PLANES, 60); check(s == 28 and info['stage'] == 2, f'genus-aware score 28/28 on Roest-Vigh ({info})')
lo_, hi_, init_, make_ = GFAMILIES['GROW_R9']; x, ctx = init_(np.random.default_rng(4)); P = make_(x, ctx)
check(len(P) == 9 and np.all(x >= ctx['lo']) and np.all(x <= ctx['hi']), 'growth family: nine planes, start inside its window')

from facenb.guided import score_faces
s, _, info = score_faces(ROESTVIGH_PLANES, 60); check(s == 28 and info['stage'] == 3, f'three-stage score 28/28 on Roest-Vigh ({info})')
print('GUIDED CONTROLS PASSED')
