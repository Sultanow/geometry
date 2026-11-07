"""Exact verification of a polyhedron with integer face planes, in rational arithmetic only (fractions.Fraction).
Input: integer planes [a, b, c, d] (ax + by + cz = d) and, for every face p, its boundary cycle as triples of plane indices
(each vertex is the intersection of the three planes of the faces meeting there). Checks, all exact:
 1 general position: no two planes parallel, every three meet in one point, no four through a common point
 2 every face is a simple polygon (no two non-consecutive edges meet, no consecutive edges overlap); reflex corners counted
 3 closed and orientable: every edge lies in exactly two face cycles, once in each direction after a consistent orientation
 4 vertex links are single cycles (closed 2-manifold); Euler characteristic and genus
 5 embedding: for every pair of faces, the intersection of the two closed polygons, computed exactly on the line where their
   planes meet as a union of closed intervals, equals the union of the edges they share
 6 every two faces share an edge; number of pairs sharing two
 7 signed volume (positive after orienting outward)"""
from fractions import Fraction as Fr
import itertools
def det3(M):
    return (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
            + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))
def solve3(M, r):
    D = det3(M)
    if D == 0: return None
    out = []
    for k in range(3):
        Mk = [list(row) for row in M]
        for i in range(3): Mk[i][k] = r[i]
        out.append(Fr(det3(Mk), D))
    return tuple(out)
def cross(a, b): return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])
def dot(a, b): return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
def sub(a, b): return tuple(x - y for x, y in zip(a, b))
def orient2(a, b, c): return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
def on_seg2(p, a, b):
    return orient2(a, b, p) == 0 and min(a[0], b[0]) <= p[0] <= max(a[0], b[0]) and min(a[1], b[1]) <= p[1] <= max(a[1], b[1])
def seg_meet2(a, b, c, d):
    o1, o2, o3, o4 = orient2(a, b, c), orient2(a, b, d), orient2(c, d, a), orient2(c, d, b)
    if o1 * o2 < 0 and o3 * o4 < 0: return True
    return on_seg2(c, a, b) or on_seg2(d, a, b) or on_seg2(a, c, d) or on_seg2(b, c, d)
def verify_exact(planes, faces):
    planes = [tuple(int(c) for c in p) for p in planes]; n = len(planes); N = [p[:3] for p in planes]; D = [p[3] for p in planes]
    rep = {}
    # 1 general position
    rep['no_parallel'] = all(cross(N[a], N[b]) != (0, 0, 0) for a, b in itertools.combinations(range(n), 2))
    X = {}
    gp = True
    for t in itertools.combinations(range(n), 3):
        x = solve3([N[i] for i in t], [D[i] for i in t])
        if x is None: gp = False; continue
        X[t] = x
        if any(dot(N[s], x) == D[s] for s in range(n) if s not in t): gp = False
    rep['general_position'] = gp
    V = sorted({t for f in faces for t in f}); rep['V'] = len(V)
    for t in V:
        assert all(dot(N[i], X[t]) == D[i] for i in t)
    # 2 simple polygons (in a coordinate projection that is bijective on the plane)
    proj, P2 = [], []
    simple, reflex = True, []
    for p, f in enumerate(faces):
        k = max(range(3), key=lambda i: abs(N[p][i])); keep = [i for i in range(3) if i != k]; proj.append(keep)
        pts = [(X[t][keep[0]], X[t][keep[1]]) for t in f]; P2.append(pts); m = len(pts)
        for i, j in itertools.combinations(range(m), 2):
            if (j - i) % m in (1, m - 1): continue
            if seg_meet2(pts[i], pts[(i + 1) % m], pts[j], pts[(j + 1) % m]): simple = False
        for i in range(m):
            a, b, c = pts[i - 1], pts[i], pts[(i + 1) % m]
            if orient2(a, b, c) == 0 and (min(a[0], c[0]) > b[0] or max(a[0], c[0]) < b[0] or min(a[1], c[1]) > b[1] or max(a[1], c[1]) < b[1] or True):
                if not (on_seg2(b, a, c)): simple = False                  # consecutive edges fold back onto each other
        area2 = sum(pts[i][0] * pts[(i + 1) % m][1] - pts[(i + 1) % m][0] * pts[i][1] for i in range(m))
        reflex.append(sum(1 for i in range(m) if orient2(pts[i - 1], pts[i], pts[(i + 1) % m]) * area2 < 0))
    rep['simple'] = simple; rep['reflex_corners'] = reflex
    # 3 closed and orientable
    owner = {}
    for p, f in enumerate(faces):
        for i in range(len(f)): owner.setdefault(frozenset((f[i], f[(i + 1) % len(f)])), []).append(p)
    rep['E'] = len(owner); rep['closed'] = all(len(v) == 2 for v in owner.values())
    sign = {0: 1}; queue = [0]; orientable = True
    def directed(p, s):
        f = faces[p] if s == 1 else faces[p][::-1]; return {(f[i], f[(i + 1) % len(f)]) for i in range(len(f))}
    while queue:
        p = queue.pop(); dp = directed(p, sign[p])
        for e, fs in owner.items():
            if p not in fs: continue
            q = fs[0] if fs[1] == p else fs[1]; a, b = tuple(e)
            dir_p = (a, b) if (a, b) in dp else (b, a)
            for s in (1, -1):
                if (dir_p[1], dir_p[0]) in directed(q, s):
                    if q in sign and sign[q] != s: orientable = False
                    if q not in sign: sign[q] = s; queue.append(q)
    rep['orientable'] = orientable and len(sign) == n
    # 4 vertex links
    links = True
    for v in V:
        inc = [frozenset(e) for e in owner if v in e]
        fl = {}
        for e in inc:
            for p in owner[e]: fl.setdefault(p, []).append(e)
        G = {}
        for p, es in fl.items():
            if len(es) != 2: links = False; continue
            G.setdefault(es[0], []).append(es[1]); G.setdefault(es[1], []).append(es[0])
        if any(len(x) != 2 for x in G.values()): links = False; continue
        seen, stack = set(), [inc[0]]
        while stack:
            e = stack.pop()
            if e in seen: continue
            seen.add(e); stack += G[e]
        if len(seen) != len(inc): links = False
    rep['manifold_links'] = links
    chi = rep['V'] - rep['E'] + n; rep['euler'] = chi; rep['genus'] = Fr(2 - chi, 2)
    # 5 embedding on the lines of pairs of planes
    def inside_closed(p, x3):
        q = (x3[proj[p][0]], x3[proj[p][1]]); pts = P2[p]; m = len(pts)
        if any(on_seg2(q, pts[i], pts[(i + 1) % m]) for i in range(m)): return True
        c = False
        for i in range(m):
            a, b = pts[i], pts[(i + 1) % m]
            if (a[1] > q[1]) != (b[1] > q[1]):
                xint = a[0] + (q[1] - a[1]) * (b[0] - a[0]) / (b[1] - a[1])
                if xint > q[0]: c = not c
        return c
    def on_line(p, q, d, P0):
        """closed set face_p ∩ L as merged closed intervals in the parameter s = d·x"""
        f = faces[p]; m = len(f); g = lambda t: dot(N[q], X[t]) - D[q]
        bps = set()
        for i in range(m):
            a, b = f[i], f[(i + 1) % m]; ga, gb = g(a), g(b)
            if ga == 0: bps.add(dot(d, X[a]))
            if ga * gb < 0:
                lam = Fr(ga, 1) / (ga - gb); x = tuple(X[a][k] + lam * (X[b][k] - X[a][k]) for k in range(3)); bps.add(dot(d, x))
        bps = sorted(bps); s0 = dot(d, P0); dd = dot(d, d)
        pt = lambda s: tuple(P0[k] + (s - s0) / dd * d[k] for k in range(3))
        iv = []
        for i in range(len(bps) - 1):
            if inside_closed(p, pt((bps[i] + bps[i + 1]) / 2)): iv.append([bps[i], bps[i + 1]])
        for b in bps:
            if not any(lo <= b <= hi for lo, hi in iv) and inside_closed(p, pt(b)): iv.append([b, b])
        iv.sort(); merged = []
        for lo, hi in iv:
            if merged and lo <= merged[-1][1]: merged[-1][1] = max(merged[-1][1], hi)
            else: merged.append([lo, hi])
        return merged
    emb = True; share = {}
    for p, q in itertools.combinations(range(n), 2):
        d = cross(N[p], N[q]); P0 = solve3([N[p], N[q], d], [D[p], D[q], 0])
        Ip, Iq = on_line(p, q, d, P0), on_line(q, p, d, P0)
        inter = []
        for a in Ip:
            for b in Iq:
                lo, hi = max(a[0], b[0]), min(a[1], b[1])
                if lo <= hi: inter.append([lo, hi])
        inter.sort()
        shared = sorted([sorted([dot(d, X[a]), dot(d, X[b])]) for e, fs in owner.items() if set(fs) == {p, q} for a, b in [tuple(e)]])
        share[(p, q)] = len(shared)
        merged = []
        for lo, hi in shared:
            if merged and lo <= merged[-1][1]: merged[-1][1] = max(merged[-1][1], hi)
            else: merged.append([lo, hi])
        if inter != merged: emb = False
    rep['embedded'] = emb
    rep['face_neighbourly'] = all(v >= 1 for v in share.values()); rep['pairs_sharing_two'] = sum(v == 2 for v in share.values())
    # 7 signed volume with the consistent orientation (divergence theorem, fan triangulation of each face)
    vol = Fr(0)
    for p, f in enumerate(faces):
        cyc = f if sign.get(p, 1) == 1 else f[::-1]
        for i in range(1, len(cyc) - 1):
            vol += Fr(dot(X[cyc[0]], cross(X[cyc[i]], X[cyc[i + 1]])), 6)
    rep['volume'] = abs(vol)
    rep['vertices'] = {t: X[t] for t in V}
    return rep
