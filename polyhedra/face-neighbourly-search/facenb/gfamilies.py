"""Parametrised families for the guided search. Each family: (lo, hi, init, make) with a continuous parameter vector x in
[lo, hi] and per-climb discrete choices ctx (which axis, which single plane, a fixed tiny perturbation where the symmetry
forces four planes through a point). make(x, ctx) -> planes keeps the symmetry while x moves."""
import numpy as np
from .families import nvec, orbit, S4, D2, C2
PI = np.pi
def _noise(rng, n, on): return 1e-4 * rng.normal(size=(n, 4)) if on else None
def _apply(P, ctx):
    P = np.array(P, float)
    return P + ctx['noise'] if ctx.get('noise') is not None else P
SINGLES = [[0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0]]
def _axis_pair(psi, axis):
    pair = [[0, np.cos(psi), np.sin(psi), 0], [0, np.cos(psi), -np.sin(psi), 0]]; perm = [[0, 1, 2], [1, 2, 0], [2, 0, 1]][axis]
    return [[q[perm[0]], q[perm[1]], q[perm[2]], 0] for q in pair]
def S4_family(n):
    lo = [0.05, 0.05, 0, -3]; hi = [PI - 0.05, PI - 0.05, 2 * PI, 3]
    if n in (10, 11): lo, hi = lo + [0], hi + [PI / 2]
    def make(x, ctx):
        P = orbit(nvec(x[0], 0), 1., S4) + orbit(nvec(x[1], x[2]), x[3], S4)
        if n in (10, 11): P += [[np.cos(x[4]), np.sin(x[4]), 0, 0], [-np.sin(x[4]), np.cos(x[4]), 0, 0]]
        if n in (9, 11): P += [[0, 0, 1, 0]]
        return _apply(P, ctx)
    def init(rng):
        lo_, hi_ = np.array(lo, float), np.array(hi, float)
        return lo_ + (hi_ - lo_) * rng.random(len(lo)), dict(noise=_noise(rng, n, n in (10, 11)))
    return np.array(lo, float), np.array(hi, float), init, make
def D2_family(n):
    lo = [0.05, 0, 0.05, 0, -3]; hi = [PI - 0.05, 2 * PI, PI - 0.05, 2 * PI, 3]
    if n in (10, 11): lo, hi = lo + [0.05], hi + [PI / 2 - 0.05]
    def make(x, ctx):
        P = orbit(nvec(x[0], x[1]), 1., D2) + orbit(nvec(x[2], x[3]), x[4], D2)
        if n in (10, 11): P += _axis_pair(x[5], ctx['axis'])
        if n in (9, 11): P += [SINGLES[ctx['single']]]
        return _apply(P, ctx)
    def init(rng):
        lo_, hi_ = np.array(lo, float), np.array(hi, float)
        return lo_ + (hi_ - lo_) * rng.random(len(lo)), dict(axis=int(rng.integers(3)), single=int(rng.integers(3)), noise=_noise(rng, n, n in (10, 11)))
    return np.array(lo, float), np.array(hi, float), init, make
def C2_family(n):
    m = n // 2; fixed = n % 2 == 1
    lo = ([-1.] if fixed else []) + [0.05, 0, -2] * m; hi = ([1.] if fixed else []) + [PI - 0.05, 2 * PI, 2] * m
    def make(x, ctx):
        P = []; k = 0
        if fixed:
            P.append([0, 0, 1, x[0]] if ctx['kind'] == 'perp' else [np.cos(PI * (x[0] + 1) / 2), np.sin(PI * (x[0] + 1) / 2), 0, 0]); k = 1
        for j in range(m): P += orbit(nvec(x[k + 3 * j], x[k + 3 * j + 1]), x[k + 3 * j + 2], C2)
        return np.array(P, float)
    def init(rng):
        lo_, hi_ = np.array(lo, float), np.array(hi, float)
        return lo_ + (hi_ - lo_) * rng.random(len(lo)), dict(kind='perp' if rng.random() < 0.5 else 'contain')
    return np.array(lo, float), np.array(hi, float), init, make
GFAMILIES = {f'S4_{n}': S4_family(n) for n in (8, 9, 10, 11)}
GFAMILIES.update({f'D2_{n}': D2_family(n) for n in (8, 9, 10, 11)})
GFAMILIES.update({f'C2_{n}': C2_family(n) for n in (7, 8, 9, 10, 11)})

def ASYM_family(n):
    """n planes without any symmetry: x = (polar angle, azimuth, offset) per plane"""
    lo = [0.05, 0, -1] * n; hi = [PI - 0.05, 2 * PI, 1] * n
    def make(x, ctx): return np.array([[*nvec(x[3 * i], x[3 * i + 1]), x[3 * i + 2]] for i in range(n)])
    def init(rng):
        lo_, hi_ = np.array(lo, float), np.array(hi, float); return lo_ + (hi_ - lo_) * rng.random(len(lo)), {}
    return np.array(lo, float), np.array(hi, float), init, make
GFAMILIES.update({f'ASYM_{n}': ASYM_family(n) for n in (8, 9, 10, 11)})

# --- growth from the known eight-faced polyhedra: add a plane that truncates a vertex (genus stays 3, a ninth face appears) ---
def _plane_params(n, d):
    n = np.asarray(n, float); k = np.linalg.norm(n); n, d = n / k, d / k
    return [np.arccos(np.clip(n[2], -1, 1)), np.arctan2(n[1], n[0]) % (2 * PI), d]
def _seed(name):
    """(planes scaled to unit size, vertex coordinates, neighbour lists) of a known eight-faced polyhedron"""
    from .known import MIZHAEV_VERTICES, MIZHAEV_FACES, ROESTVIGH_PLANES
    from .arrsat import solve, describe
    if name == 'Mizhaev':
        X = np.array(MIZHAEV_VERTICES, float); faces = MIZHAEV_FACES; pts = {i: X[i] for i in range(len(X))}
    else:
        A, s = solve(ROESTVIGH_PLANES, max_solutions=1); d = describe(A, s[0]); faces = d['faces']; pts = d['coords']
    keys = list(pts); scale = max(np.abs(pts[k]).max() for k in keys)
    P = []
    for f in faces:
        a, b, c = (pts[f[0]], pts[f[1]], pts[f[2]]); nrm = np.cross(b - a, c - a); P.append(_plane_params(nrm, nrm @ a / scale))
    nb = {k: set() for k in keys}
    for f in faces:
        for i in range(len(f)): nb[f[i]].add(f[i - 1]); nb[f[i]].add(f[(i + 1) % len(f)])
    return np.array(P, float), {k: pts[k] / scale for k in keys}, nb
def GROW_family(seed_name, n_extra=1, window=(0.35, 0.35, 0.35)):
    base, pts, nb = _seed(seed_name); keys = list(pts); n0 = len(base)
    def make(x, ctx): return np.array([[*nvec(x[3 * i], x[3 * i + 1]), x[3 * i + 2]] for i in range(n0 + n_extra)])
    def init(rng):
        x = list(base.reshape(-1)); ctx = dict(cut=[])
        for _ in range(n_extra):
            v = keys[int(rng.integers(len(keys)))]; w = list(nb[v]); t = rng.uniform(0.1, 0.4)
            q = [pts[v] + t * (pts[u] - pts[v]) for u in w[:3]]; nrm = np.cross(q[1] - q[0], q[2] - q[0])
            x += _plane_params(nrm, nrm @ q[0]); ctx['cut'].append(str(v))
        x = np.array(x, float); win = np.tile(window, n0 + n_extra)
        ctx['lo'], ctx['hi'] = x - win, x + win
        return x, ctx
    lo = np.full(3 * (n0 + n_extra), -10.0); hi = np.full(3 * (n0 + n_extra), 10.0)
    return lo, hi, init, make
GFAMILIES.update(GROW_M9=GROW_family('Mizhaev'), GROW_R9=GROW_family('Roest-Vigh'))

# --- second generation: start from certified near-misses (analysis.jsonl of a previous run, e.g. 34/36 with nine faces) ---
def SEED_family(window=(0.12, 0.12, 0.12)):
    """seeds: records with 'planes' and 'pairs' >= FN_SEED_MIN (default 34) in FN_SEEDS (default $FN_OUT/analysis_F/analysis.jsonl)"""
    import os, json
    cache = {}
    def seeds():
        if 'list' not in cache:
            path = os.environ.get('FN_SEEDS') or os.path.join(os.environ.get('FN_OUT', '.'), 'analysis_F', 'analysis.jsonl')
            mn = int(os.environ.get('FN_SEED_MIN', '34'))
            cache['list'] = [r for r in (json.loads(l) for l in open(path) if l.strip()) if r.get('pairs', -1) >= mn and r.get('planes')]
            if not cache['list']: raise RuntimeError(f'no seeds with >= {mn} pairs in {path}')
        return cache['list']
    def make(x, ctx): return np.array([[*nvec(x[3 * i], x[3 * i + 1]), x[3 * i + 2]] for i in range(len(x) // 3)])
    def init(rng):
        S = seeds(); r = S[int(rng.integers(len(S)))]
        x = np.array([v for pl in r['planes'] for v in _plane_params(pl[:3], pl[3])], float); win = np.tile(window, len(x) // 3)
        return x, dict(seed=r['file'], lo=x - win, hi=x + win)
    return np.full(27, -10.0), np.full(27, 10.0), init, make
GFAMILIES['SEED_9'] = SEED_family()
