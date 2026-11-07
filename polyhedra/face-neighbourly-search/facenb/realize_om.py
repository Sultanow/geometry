"""Realisation guided by a chirotope: minimise  W * sign_violation(x) + penalty(x)  where sign_violation sums, over all
4-subsets S of the n planes and the plane at infinity, max(0, margin - sigma_S * det(h_S) / prod|h|) for the target signs
sigma_S (a SAT model of facenb/omsat.py), and penalty is the exact embedding penalty of realize_fast (0 = valid)."""
import itertools, numpy as np
from .realize_fast import make_penalty
def sign_terms(n):
    S = list(itertools.combinations(range(n + 1), 4)); return S, np.array(S)
def make_objective(faces, chi, W=20.0, margin=0.02):
    n = len(faces); pen, M = make_penalty(faces); S, Sa = sign_terms(n)
    sig = np.array([chi[s] for s in S], float)
    def H(x):
        P = np.asarray(x, float).reshape(n, 4); Hm = np.zeros((n + 1, 4)); Hm[:n, :3] = P[:, :3]; Hm[:n, 3] = -P[:, 3]; Hm[n] = [0, 0, 0, 1]
        return Hm / np.linalg.norm(Hm, axis=1, keepdims=True)
    def violation(x):
        d = np.linalg.det(H(x)[Sa]); return float(np.sum(np.maximum(0.0, margin - sig * d)))
    def f(x): return W * violation(x) + pen(x)
    return f, violation, pen, M
def realise_om(faces, chi, rng, time_budget=120, maxfevals=30000, sigma=0.5, W=20.0):
    import cma, time
    f, viol, pen, M = make_objective(faces, chi, W=W); n = len(faces); t0 = time.time(); best = (np.inf, None); runs = 0
    while time.time() - t0 < time_budget:
        x0 = np.concatenate([np.concatenate([rng.normal(size=3), [rng.uniform(-1, 1)]]) for _ in range(n)])
        es = cma.CMAEvolutionStrategy(x0, sigma, {'maxfevals': maxfevals, 'verbose': -9, 'seed': int(rng.integers(1 << 30)), 'tolfun': 1e-14})
        hit = None
        while not es.stop() and time.time() - t0 < time_budget:
            xs = es.ask(); fs = [f(z) for z in xs]; es.tell(xs, fs); k = int(np.argmin(fs))
            if pen(xs[k]) == 0: hit = np.asarray(xs[k]); break
        runs += 1
        xb = hit if hit is not None else np.asarray(es.result.xbest); p = pen(xb)
        if p < best[0]: best = (p, xb, viol(xb))
        if p == 0: break
    return best, runs
