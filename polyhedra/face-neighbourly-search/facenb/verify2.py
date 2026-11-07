"""Embedding check as in Roest-Vigh: every face is a simple polygon, and for every pair of faces the intersection of the two
closed polygons, computed on the line where their planes meet, equals the union of the edges they share. Polygons are grown
by a hair (eps) so that edges lying exactly on the line are not lost to rounding; lengths are compared in absolute units."""
import numpy as np
from shapely.geometry import Polygon, LineString
def frame(pl):
    n = np.array(pl[:3], float); nn = np.linalg.norm(n); n /= nn; O = n * pl[3] / nn
    u = np.cross(n, [1, 0, 0] if abs(n[0]) < 0.9 else [0, 1, 0]); u /= np.linalg.norm(u); v = np.cross(n, u); return O, u, v
def intervals(poly2, O, u, v, P0, dh, eps):
    a = np.array([(P0 - O) @ u, (P0 - O) @ v]); b = np.array([dh @ u, dh @ v]); b /= np.linalg.norm(b)
    L = 1e3 * (1 + np.abs(np.array(poly2.bounds)).max() + np.linalg.norm(a))
    I = poly2.buffer(eps, join_style=2).intersection(LineString([a - L * b, a + L * b]))
    geoms = [I] if I.geom_type in ('LineString', 'Point') else list(getattr(I, 'geoms', []))
    out = []
    for g in geoms:
        if g.is_empty: continue
        ss = [(np.array(c) - a) @ b for c in g.coords]; out.append((min(ss), max(ss)))
    return sorted(out)
def verify2(planes, faces, coords):
    X = np.array(list(coords.values())); scale = 1 + np.abs(X).max(); eps = 1e-9 * scale; tol = 1e-6 * scale
    res = dict(simple=True, pairs_ok=0, pairs_bad=[]); polys = []
    for p, f in enumerate(faces):
        O, u, v = frame(planes[p]); P2 = Polygon([((coords[t] - O) @ u, (coords[t] - O) @ v) for t in f])
        if not (P2.is_valid and P2.exterior.is_simple): res['simple'] = False
        polys.append((P2, O, u, v))
    E = {}
    for p, f in enumerate(faces):
        for i in range(len(f)): E.setdefault(frozenset((f[i], f[(i + 1) % len(f)])), set()).add(p)
    for p in range(len(faces)):
        for q in range(p + 1, len(faces)):
            np_, nq = np.array(planes[p][:3], float), np.array(planes[q][:3], float); d = np.cross(np_, nq); dh = d / np.linalg.norm(d)
            P0 = np.linalg.solve(np.array([np_, nq, d]), [planes[p][3], planes[q][3], 0.0])
            Ip = intervals(*polys[p], P0, dh, eps); Iq = intervals(*polys[q], P0, dh, eps)
            pos = []
            for a0, a1 in Ip:
                for b0, b1 in Iq:
                    lo, hi = max(a0, b0), min(a1, b1)
                    if hi - lo > tol: pos.append((lo, hi))
            shared = sorted(tuple(sorted(((coords[a] - P0) @ dh, (coords[b] - P0) @ dh))) for e, fs in E.items() if fs == {p, q} for a, b in [tuple(e)])
            ok = len(pos) == len(shared) and all(abs(x[0] - y[0]) < tol and abs(x[1] - y[1]) < tol for x, y in zip(sorted(pos), shared))
            if ok: res['pairs_ok'] += 1
            else: res['pairs_bad'].append((p, q, pos, shared))
    return res
