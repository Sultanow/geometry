"""Numba version of realize.penalty (identical definition, compiled): use make_penalty(faces) -> f(x)."""
import numpy as np
from numba import njit
from .realize import prepare
@njit(cache=True)
def _pen(x, trip, Fidx, Flen, Eij, Epq, n, NBidx):
    P = x.reshape(n, 4); N = np.empty((n, 3)); D = np.empty(n)
    for p in range(n):
        nn = np.sqrt(P[p, 0] ** 2 + P[p, 1] ** 2 + P[p, 2] ** 2)
        if nn < 1e-6: return 1e6
        N[p, 0] = P[p, 0] / nn; N[p, 1] = P[p, 1] / nn; N[p, 2] = P[p, 2] / nn; D[p] = P[p, 3] / nn
    V = trip.shape[0]; X = np.empty((V, 3)); bad = 0
    for k in range(V):
        a, b, c = trip[k, 0], trip[k, 1], trip[k, 2]
        m00, m01, m02 = N[a, 0], N[a, 1], N[a, 2]; m10, m11, m12 = N[b, 0], N[b, 1], N[b, 2]; m20, m21, m22 = N[c, 0], N[c, 1], N[c, 2]
        det = m00 * (m11 * m22 - m12 * m21) - m01 * (m10 * m22 - m12 * m20) + m02 * (m10 * m21 - m11 * m20)
        if abs(det) < 1e-4: bad += 1; continue
        r0, r1, r2 = D[a], D[b], D[c]
        X[k, 0] = (r0 * (m11 * m22 - m12 * m21) - m01 * (r1 * m22 - m12 * r2) + m02 * (r1 * m21 - m11 * r2)) / det
        X[k, 1] = (m00 * (r1 * m22 - m12 * r2) - r0 * (m10 * m22 - m12 * m20) + m02 * (m10 * r2 - r1 * m20)) / det
        X[k, 2] = (m00 * (m11 * r2 - r1 * m21) - m01 * (m10 * r2 - r1 * m20) + r0 * (m10 * m21 - m11 * m20)) / det
    if bad > 0: return 1e5 + 1e3 * bad
    cx, cy, cz = 0.0, 0.0, 0.0
    for k in range(V): cx += X[k, 0]; cy += X[k, 1]; cz += X[k, 2]
    cx /= V; cy /= V; cz /= V
    R = np.empty(V)
    for k in range(V): R[k] = np.sqrt((X[k, 0] - cx) ** 2 + (X[k, 1] - cy) ** 2 + (X[k, 2] - cz) ** 2)
    scale = np.median(R) + 1e-9; pen = 0.0
    for k in range(V):
        if R[k] / scale > 20: pen += R[k] / scale - 20
    maxl = Fidx.shape[1]; Q = np.empty((n, maxl, 2)); U = np.empty((n, 3)); W = np.empty((n, 3))
    for p in range(n):                                   # frames and 2D polygons
        nx, ny, nz = N[p, 0], N[p, 1], N[p, 2]
        if abs(nx) < 0.9: ax, ay, az = 1.0, 0.0, 0.0
        else: ax, ay, az = 0.0, 1.0, 0.0
        ux, uy, uz = ny * az - nz * ay, nz * ax - nx * az, nx * ay - ny * ax; un = np.sqrt(ux * ux + uy * uy + uz * uz); ux /= un; uy /= un; uz /= un
        wx, wy, wz = ny * uz - nz * uy, nz * ux - nx * uz, nx * uy - ny * ux
        U[p, 0], U[p, 1], U[p, 2] = ux, uy, uz; W[p, 0], W[p, 1], W[p, 2] = wx, wy, wz
        for i in range(Flen[p]):
            k = Fidx[p, i]; Q[p, i, 0] = X[k, 0] * ux + X[k, 1] * uy + X[k, 2] * uz; Q[p, i, 1] = X[k, 0] * wx + X[k, 1] * wy + X[k, 2] * wz
    for p in range(n):                                   # self-intersections
        m = Flen[p]
        for i in range(m):
            a0x, a0y = Q[p, i, 0], Q[p, i, 1]; a1x, a1y = Q[p, (i + 1) % m, 0], Q[p, (i + 1) % m, 1]
            for j in range(i + 2, m):
                if (j + 1) % m == i: continue
                b0x, b0y = Q[p, j, 0], Q[p, j, 1]; b1x, b1y = Q[p, (j + 1) % m, 0], Q[p, (j + 1) % m, 1]
                if NBidx[p, i] == NBidx[p, j]:
                    dx, dy = a1x - a0x, a1y - a0y; dl = np.sqrt(dx * dx + dy * dy)
                    if dl < 1e-30: dl = 1e-30
                    dx /= dl; dy /= dl
                    t1 = (b0x - a0x) * dx + (b0y - a0y) * dy; t2 = (b1x - a0x) * dx + (b1y - a0y) * dy
                    lo2, hi2 = min(t1, t2), max(t1, t2); ov = min(dl, hi2) - max(0.0, lo2)
                    if ov > 1e-9 * scale: pen += 1 + ov / scale
                    continue
                o1 = (a1x - a0x) * (b0y - a0y) - (a1y - a0y) * (b0x - a0x); o2 = (a1x - a0x) * (b1y - a0y) - (a1y - a0y) * (b1x - a0x)
                o3 = (b1x - b0x) * (a0y - b0y) - (b1y - b0y) * (a0x - b0x); o4 = (b1x - b0x) * (a1y - b0y) - (b1y - b0y) * (a1x - b0x)
                mo = min(min(abs(o1), abs(o2)), min(abs(o3), abs(o4)))
                if o1 * o2 < 0 and o3 * o4 < 0 and mo > 1e-12 * scale ** 2: pen += 1 + mo / scale ** 2
    for e in range(Eij.shape[0]):                        # edge-face piercing
        i, j = Eij[e, 0], Eij[e, 1]
        for f in range(n):
            if f == Epq[e, 0] or f == Epq[e, 1]: continue
            if trip[i, 0] == f or trip[i, 1] == f or trip[i, 2] == f or trip[j, 0] == f or trip[j, 1] == f or trip[j, 2] == f: continue
            si = N[f, 0] * X[i, 0] + N[f, 1] * X[i, 1] + N[f, 2] * X[i, 2] - D[f]; sj = N[f, 0] * X[j, 0] + N[f, 1] * X[j, 1] + N[f, 2] * X[j, 2] - D[f]
            if si * sj >= 0: continue
            t = si / (si - sj); yx3 = X[i, 0] + t * (X[j, 0] - X[i, 0]); yy3 = X[i, 1] + t * (X[j, 1] - X[i, 1]); yz3 = X[i, 2] + t * (X[j, 2] - X[i, 2])
            yx = yx3 * U[f, 0] + yy3 * U[f, 1] + yz3 * U[f, 2]; yy = yx3 * W[f, 0] + yy3 * W[f, 1] + yz3 * W[f, 2]
            m = Flen[f]; inside = False; dmin = 1e300
            for k in range(m):
                a0x, a0y = Q[f, k, 0], Q[f, k, 1]; a1x, a1y = Q[f, (k + 1) % m, 0], Q[f, (k + 1) % m, 1]
                if (a0y > yy) != (a1y > yy):
                    if yx < a0x + (yy - a0y) * (a1x - a0x) / (a1y - a0y): inside = not inside
                abx, aby = a1x - a0x, a1y - a0y; L2 = abx * abx + aby * aby
                tt = ((yx - a0x) * abx + (yy - a0y) * aby) / L2 if L2 > 1e-30 else 0.0
                tt = min(1.0, max(0.0, tt)); dx, dy = yx - (a0x + tt * abx), yy - (a0y + tt * aby); d = np.sqrt(dx * dx + dy * dy)
                if d < dmin: dmin = d
            if inside: pen += 1 + dmin / scale
    return pen
def make_penalty(faces):
    M = prepare(faces); n = M['n']; maxl = max(len(f) for f in M['F'])
    Fidx = np.zeros((n, maxl), np.int64); Flen = np.array([len(f) for f in M['F']], np.int64)
    for p, f in enumerate(M['F']): Fidx[p, :len(f)] = f
    Eij = np.array([e for e, _ in M['E']], np.int64); Epq = np.array([ps for _, ps in M['E']], np.int64); trip = M['trip'].astype(np.int64)
    NBidx = np.full((n, maxl), -1, np.int64)
    for p, nb in enumerate(M['NB']): NBidx[p, :len(nb)] = nb
    return (lambda x: _pen(np.ascontiguousarray(x, dtype=np.float64), trip, Fidx, Flen, Eij, Epq, n, NBidx)), M
def realise_fast(faces, rng, restarts=1000, maxfevals=20000, sigma=0.6, time_budget=None, verbose=False):
    import cma, time
    f, M = make_penalty(faces); n = M['n']; best = (np.inf, None); t0 = time.time(); runs = 0
    for r in range(restarts):
        if time_budget and time.time() - t0 > time_budget: break
        x0 = np.concatenate([np.concatenate([rng.normal(size=3), [rng.uniform(-1, 1)]]) for _ in range(n)])
        es = cma.CMAEvolutionStrategy(x0, sigma, {'maxfevals': maxfevals, 'verbose': -9, 'seed': int(rng.integers(1 << 30)), 'tolfun': 1e-12})
        hit = None
        while not es.stop():
            xs = es.ask(); fs = [f(z) for z in xs]; es.tell(xs, fs)
            k = int(np.argmin(fs))
            if fs[k] == 0: hit = np.asarray(xs[k]); break
            if time_budget and time.time() - t0 > time_budget: break
        runs += 1; fb = 0.0 if hit is not None else f(es.result.xbest)
        xb = hit if hit is not None else np.asarray(es.result.xbest)
        if verbose: print(f'  restart {r}: best {fb:.4f}')
        if fb < best[0]: best = (fb, xb)
        if fb == 0: break
    return best, M, runs
