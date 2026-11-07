"""Geometric realisation of a FIXED combinatorial type. The map is given as face cycles of vertex triples (every vertex is the
intersection of the planes of the three faces meeting there). Unknowns: the n face planes (a, b, c, d), 4n numbers.
penalty(x) >= 0 adds up
  - edge-face piercing: for every edge and every face not containing it, if the edge crosses the face's plane inside the face
    polygon, the distance of the crossing point to the polygon boundary (+ a constant)
  - self-intersection of faces: for crossing non-adjacent edges of the same face, the smallest orientation depth (+ a constant)
  - degeneracy: nearly parallel planes / far away vertices
and is 0 exactly when the planes realise the map as an embedded polyhedron with simple faces (checked again afterwards by the
Roest-Vigh interval test and exactly). Minimised with CMA-ES from random starts."""
import numpy as np
def prepare(faces):
    n = len(faces); trip = sorted({tuple(sorted(t)) for f in faces for t in f}); ix = {t: i for i, t in enumerate(trip)}
    F = [[ix[tuple(sorted(t))] for t in f] for f in faces]
    edges = {}
    for p, f in enumerate(F):
        for i in range(len(f)): edges.setdefault(frozenset((f[i], f[(i + 1) % len(f)])), set()).add(p)
    E = [(tuple(e), tuple(sorted(ps))) for e, ps in edges.items()]
    NB = []                                                            # neighbour face across each side of each face
    for p, f in enumerate(faces):
        NB.append([(set(tuple(sorted(f[i]))) & set(tuple(sorted(f[(i + 1) % len(f)]))) - {p}).pop() for i in range(len(f))])
    return dict(n=n, trip=np.array(trip), F=F, E=E, NB=NB)
def frame(nrm):
    n = nrm / np.linalg.norm(nrm); u = np.cross(n, [1, 0, 0] if abs(n[0]) < 0.9 else [0, 1, 0]); u /= np.linalg.norm(u); return n, u, np.cross(n, u)
def cross2(a, b): return a[..., 0] * b[..., 1] - a[..., 1] * b[..., 0]
def seg_dist(P, A, B):
    AB = B - A; t = np.clip(((P - A) @ AB) / max(AB @ AB, 1e-30), 0, 1); return np.linalg.norm(P - (A + t * AB))
def penalty(x, M, details=False):
    n = M['n']; P = x.reshape(n, 4); N = P[:, :3]; D = P[:, 3]
    nn = np.linalg.norm(N, axis=1)
    if np.any(nn < 1e-6): return 1e6
    N = N / nn[:, None]; D = D / nn
    T = M['trip']; A = N[T]; b = D[T]
    det = np.linalg.det(A)
    if np.any(np.abs(det) < 1e-4): return 1e5 + 1e3 * np.sum(np.abs(det) < 1e-4)
    X = np.linalg.solve(A, b[..., None])[..., 0]
    scale = np.median(np.linalg.norm(X - X.mean(0), axis=1)) + 1e-9
    pen = 0.0; parts = dict(far=0.0, simple=0.0, pierce=0.0)
    far = np.linalg.norm(X - X.mean(0), axis=1) / scale; parts['far'] = float(np.sum(np.maximum(far - 20, 0)))
    polys = []
    for p, f in enumerate(M['F']):                                     # 2D polygons in their planes
        nrm, u, v = frame(N[p]); Q = np.stack([X[f] @ u, X[f] @ v], 1); polys.append((Q, u, v))
        m = len(f)
        for i in range(m):                                             # self-intersection
            a0, a1 = Q[i], Q[(i + 1) % m]
            for j in range(i + 2, m):
                if (j + 1) % m == i: continue
                b0, b1 = Q[j], Q[(j + 1) % m]
                if M['NB'][p][i] == M['NB'][p][j]:                     # two sides on the same line (doubled pair): no overlap
                    dvec = a1 - a0; dl = np.sqrt(dvec @ dvec); dvec = dvec / max(dl, 1e-30)
                    s1 = sorted((0.0, dl)); s2 = sorted(((b0 - a0) @ dvec, (b1 - a0) @ dvec))
                    ov = min(s1[1], s2[1]) - max(s1[0], s2[0])
                    if ov > 1e-9 * scale: parts['simple'] += 1 + ov / scale
                    continue
                o1, o2 = cross2(a1 - a0, b0 - a0), cross2(a1 - a0, b1 - a0); o3, o4 = cross2(b1 - b0, a0 - b0), cross2(b1 - b0, a1 - b0)
                tol = 1e-12 * scale ** 2
                if o1 * o2 < 0 and o3 * o4 < 0 and min(abs(o1), abs(o2), abs(o3), abs(o4)) > tol: parts['simple'] += 1 + min(abs(o1), abs(o2), abs(o3), abs(o4)) / scale ** 2
    for (i, j), ps in M['E']:                                          # edge-face piercing
        Xi, Xj = X[i], X[j]
        touch = set(M['trip'][i]) | set(M['trip'][j])                     # faces through an endpoint: no piercing test
        for f in range(n):
            if f in ps or f in touch: continue
            si, sj = N[f] @ Xi - D[f], N[f] @ Xj - D[f]
            if si * sj >= 0: continue
            Y = Xi + si / (si - sj) * (Xj - Xi); Q, u, v = polys[f]; y = np.array([Y @ u, Y @ v])
            inside = False; m = len(Q)
            for k in range(m):
                a0, a1 = Q[k], Q[(k + 1) % m]
                if (a0[1] > y[1]) != (a1[1] > y[1]) and y[0] < a0[0] + (y[1] - a0[1]) * (a1[0] - a0[0]) / (a1[1] - a0[1]): inside = not inside
            if inside:
                dmin = min(seg_dist(y, Q[k], Q[(k + 1) % m]) for k in range(m)); parts['pierce'] += 1 + dmin / scale
    parts['reflex'] = 0.0                                              # at most one reflex corner per vertex
    reflex_at = {}
    for p_, f in enumerate(M['F']):
        Q = polys[p_][0]; m = len(f); area = 0.5 * sum(Q[i, 0] * Q[(i + 1) % m, 1] - Q[(i + 1) % m, 0] * Q[i, 1] for i in range(m))
        for i in range(m):
            tr = cross2(Q[i] - Q[i - 1], Q[(i + 1) % m] - Q[i])
            if tr * area < 0: reflex_at.setdefault(f[i], []).append(abs(tr))
    for v, lst in reflex_at.items():
        if len(lst) >= 2: parts['reflex'] += (len(lst) - 1) * (1 + min(lst) / scale ** 2)
    pen = parts['far'] + parts['simple'] + parts['pierce'] + parts['reflex']
    return (pen, parts, X) if details else pen
def realise(faces, rng, restarts=10, maxfevals=20000, sigma=0.6, verbose=False, time_budget=None):
    import cma, time
    M = prepare(faces); n = M['n']; best = (np.inf, None); t0 = time.time()
    for r in range(restarts):
        if time_budget and time.time() - t0 > time_budget: break
        x0 = np.concatenate([np.concatenate([rng.normal(size=3), [rng.uniform(-1, 1)]]) for _ in range(n)])
        es = cma.CMAEvolutionStrategy(x0, sigma, {'maxfevals': maxfevals, 'verbose': -9, 'seed': int(rng.integers(1 << 30)), 'tolfun': 1e-12})
        while not es.stop():
            xs = es.ask(); fs = [penalty(np.asarray(z), M) for z in xs]; es.tell(xs, fs)
            if min(fs) == 0 or (time_budget and time.time() - t0 > time_budget): break
        xb = es.result.xbest if es.result.xbest is not None else xs[int(np.argmin(fs))]
        fb = penalty(np.asarray(xb), M)
        if min(fs) == 0: xb = np.asarray(xs[int(np.argmin(fs))]); fb = 0.0
        if verbose: print(f'  restart {r}: best penalty {fb:.4f}, evaluations {es.countevals}')
        if fb < best[0]: best = (fb, np.asarray(xb))
        if fb == 0: break
    return best, M
