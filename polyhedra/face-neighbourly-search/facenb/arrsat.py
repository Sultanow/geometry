"""Face-neighbourly polyhedra from a plane arrangement.
Given n planes in general position, every polyhedron whose faces lie in these planes (one face per plane) is a choice of
bounded 2-cells of the arrangement. The surface is closed and embedded iff (i) every 1-cell (segment of a line L_pq) has
0 or 2 chosen incident cells, (ii) no face is pinched at a vertex (two opposite quadrants in one plane), (iii) the link of
every vertex is a single cycle (no two disjoint trihedral corners), and (iv) the cells chosen in each plane form one simple
polygon. Face-neighbourly: every pair of planes shares at least one edge on its line. (i)-(iii) and neighbourliness are
CNF clauses; (iv) is enforced lazily with valid cuts."""
import itertools, numpy as np
from shapely.geometry import LineString, Polygon, box
from shapely.ops import unary_union, polygonize
from pysat.solvers import Cadical153

def build(planes, tol=1e-7):
    P = np.asarray(planes, float); n = len(P); N, D = P[:, :3], P[:, 3]
    for a, b in itertools.combinations(range(n), 2):                           # parallel planes can never share an edge
        if np.linalg.norm(np.cross(N[a], N[b])) < 1e-9 * np.linalg.norm(N[a]) * np.linalg.norm(N[b]): raise ValueError('parallel planes')
    tri = {}
    for t in itertools.combinations(range(n), 3):
        M = N[list(t)]
        if abs(np.linalg.det(M)) < 1e-9 * np.prod(np.linalg.norm(M, axis=1)): continue   # line of two planes parallel to the third: no vertex
        tri[t] = np.linalg.solve(M, D[list(t)])
    for t, x in tri.items():                                     # no fourth plane through a triple point
        for s in range(n):
            if s not in t and abs(N[s] @ x - D[s]) < tol * np.linalg.norm(N[s]) * (1 + np.abs(x).max()): raise ValueError('four planes meet')
    cells, seginc, vinc = [], {}, {}
    span = 10 * (1 + max(np.abs(x).max() for x in tri.values()))
    for p in range(n):
        nn = N[p] / np.linalg.norm(N[p]); O = nn * D[p] / np.linalg.norm(N[p])
        u = np.cross(nn, [1, 0, 0] if abs(nn[0]) < 0.9 else [0, 1, 0]); u /= np.linalg.norm(u); v = np.cross(nn, u)
        to2 = lambda X: np.array([(X - O) @ u, (X - O) @ v])
        pts2 = {t: to2(x) for t, x in tri.items() if p in t}
        lines = []
        for q in range(n):
            if q == p: continue
            a, b, c = N[q] @ u, N[q] @ v, D[q] - N[q] @ O
            nrm = np.hypot(a, b); a, b, c = a / nrm, b / nrm, c / nrm
            base = np.array([a, b]) * c; dirv = np.array([-b, a])
            lines.append(LineString([base - 4 * span * dirv, base + 4 * span * dirv]))
        B = box(-2 * span, -2 * span, 2 * span, 2 * span)
        polys = list(polygonize(unary_union(lines + [B.exterior])))
        keys = list(pts2.keys()); arr = np.array([pts2[k] for k in keys])
        # vertex tolerance from the LOCAL geometry of this plane (a far-away triple point must not coarsen it)
        dmin = min(np.linalg.norm(arr[i] - arr[j]) for i in range(len(arr)) for j in range(i + 1, len(arr))) if len(arr) > 1 else 1.0
        tol_v = min(1e-6 * span, 0.25 * dmin)
        if dmin < 1e-9: raise ValueError('two triple points (almost) coincide')
        for poly in polys:
            if poly.exterior.distance(B.exterior) < 1e-9 * span: continue          # unbounded cell (touches the box)
            co = np.array(poly.exterior.coords)[:-1]
            if len(co) < 3 or Polygon(co).area <= 0: continue                        # only truly empty polygons are skipped
            ts = []
            for c2 in co:
                dist = np.linalg.norm(arr - c2, axis=1); k = int(np.argmin(dist))
                if dist[k] > tol_v: raise ValueError('cell vertex is not a triple point')
                ts.append(keys[k])
            segs, corners = [], {}
            for i in range(len(ts)):
                t1, t2 = ts[i], ts[(i + 1) % len(ts)]
                common = (set(t1) & set(t2)) - {p}
                assert len(common) == 1, (t1, t2)
                q = common.pop(); r1 = (set(t1) - {p, q}).pop(); r2 = (set(t2) - {p, q}).pop()
                segs.append((min(p, q), max(p, q), frozenset((r1, r2))))
            for i in range(len(ts)):
                corners[ts[i]] = (segs[i - 1], segs[i])                            # the two segments at this corner
            cid = len(cells); cells.append(dict(plane=p, tri=ts, segs=segs, corners=corners, poly=Polygon(co), frame=(O, u, v)))
            for s in segs: seginc.setdefault(s, []).append(cid)
            for t, (sa, sb) in corners.items(): vinc.setdefault(t, []).append((cid, sa, sb))
    return dict(n=n, planes=P, tri=tri, cells=cells, seginc=seginc, vinc=vinc)

def clauses(A):
    cells, seginc, vinc, n = A['cells'], A['seginc'], A['vinc'], A['n']
    x = lambda c: c + 1; nxt = [len(cells) + 1]; cl = []
    def newvar():
        nxt[0] += 1; return nxt[0] - 1
    for s, inc in seginc.items():                                                 # (i) 0 or 2 chosen cells per segment
        for k in range(1, len(inc) + 1):
            if k in (0, 2): continue
            for ch in itertools.combinations(inc, k):
                cl.append([-x(c) for c in ch] + [x(c) for c in inc if c not in ch])
    for t, inc in vinc.items():
        byplane = {}
        for (c, sa, sb) in inc: byplane.setdefault(cells[c]['plane'], []).append((c, frozenset((sa, sb))))
        for p, lst in byplane.items():                                            # (ii) no pinch in a plane
            for (c1, n1), (c2, n2) in itertools.combinations(lst, 2):
                if not (n1 & n2):
                    others = [c for c, _ in lst if c not in (c1, c2)]
                    cl.append([-x(c1), -x(c2)] + [x(c) for c in others])
        cyc = []                                                                  # (iii) trihedral 3-cycles at t
        planes = sorted(byplane)
        if len(planes) == 3:
            for a in byplane[planes[0]]:
                for b in byplane[planes[1]]:
                    for c in byplane[planes[2]]:
                        nodes = a[1] | b[1] | c[1]
                        if len(nodes) == 3 and len(a[1] & b[1]) == 1 and len(a[1] & c[1]) == 1 and len(b[1] & c[1]) == 1:
                            cyc.append((frozenset(nodes), (a[0], b[0], c[0])))
        for (n1, c1), (n2, c2) in itertools.combinations(cyc, 2):
            if not (n1 & n2): cl.append([-x(c) for c in c1 + c2])
    edge_of_pair = {}
    for s, inc in seginc.items():                                                 # neighbourliness
        p, q, _ = s
        cp = [c for c in inc if cells[c]['plane'] == p]; cq = [c for c in inc if cells[c]['plane'] == q]
        if not cp or not cq: continue
        e = newvar()
        for side in (cp, cq):
            cl.append([-e] + [x(c) for c in side])
            if len(side) == 2: cl.append([-e, -x(side[0]), -x(side[1])])
        edge_of_pair.setdefault((p, q), []).append(e)
    for p, q in itertools.combinations(range(n), 2):
        es = edge_of_pair.get((p, q), [])
        if not es: return None                                                    # this pair can never share an edge
        cl.append(es)
    return cl

def face_check(A, chosen):
    """returns a list of violated planes with a valid cut each, or [] if every plane holds one simple polygon"""
    cells = A['cells']; cuts = []
    for p in range(A['n']):
        S = [c for c in chosen if cells[c]['plane'] == p]
        U = unary_union([cells[c]['poly'] for c in S]).buffer(0)
        comps = list(U.geoms) if U.geom_type == 'MultiPolygon' else [U]
        bad = len(comps) != 1 or any(len(g.interiors) for g in comps)
        if not bad: continue
        segS = {s for c in S for s in cells[c]['segs']}
        nbr = [c for c, cc in enumerate(cells) if cc['plane'] == p and c not in S and set(cc['segs']) & segS]
        if len(comps) != 1: cuts.append([-(c + 1) for c in S] + [c + 1 for c in nbr])
        else:
            hole = [c for c, cc in enumerate(cells) if cc['plane'] == p and c not in S and U.contains(cc['poly'].representative_point())]
            cuts.append([-(c + 1) for c in S] + [c + 1 for c in hole])
    return cuts

def solve(planes, max_solutions=10, max_rounds=100000, time_limit=None):
    """all face-neighbourly polyhedra on the arrangement (up to max_solutions); time_limit in seconds for the SAT part"""
    import threading, time as _t
    A = build(planes); cl = clauses(A)
    if cl is None: return A, []
    sols = []; s = Cadical153(bootstrap_with=cl); rounds = 0; t_end = None if time_limit is None else _t.time() + time_limit
    while len(sols) < max_solutions and rounds < max_rounds:
        if t_end is not None:
            left = t_end - _t.time()
            if left <= 0: A['timeout'] = True; break
            timer = threading.Timer(left, s.interrupt); timer.start()
            ok = s.solve_limited(expect_interrupt=True); timer.cancel()
            if ok is None: A['timeout'] = True; break
        else:
            ok = s.solve()
        if not ok: break
        rounds += 1
        m = s.get_model(); chosen = [c for c in range(len(A['cells'])) if m[c] > 0]
        cuts = face_check(A, chosen)
        if cuts:
            for c in cuts: s.add_clause(c)
            continue
        sols.append(chosen); s.add_clause([-(c + 1) for c in chosen] + [c + 1 for c in range(len(A['cells'])) if c not in chosen])
    s.delete(); return A, sols

def describe(A, chosen):
    """faces as cycles of 3D corner points; V, E, F, genus; number of shared edges per pair"""
    cells, tri, n = A['cells'], A['tri'], A['n']
    faces = []
    for p in range(n):
        S = [c for c in chosen if cells[c]['plane'] == p]
        U = unary_union([cells[c]['poly'] for c in S]).buffer(0)
        O, u, v = cells[S[0]]['frame']
        co = np.array(U.exterior.coords)[:-1]
        pts = {t: x for t, x in tri.items() if p in t}; keys = list(pts); arr = np.array([[(pts[k] - O) @ u, (pts[k] - O) @ v] for k in keys])
        cyc = [keys[int(np.argmin(np.linalg.norm(arr - c2, axis=1)))] for c2 in co]
        # keep true corners only (a polyhedron vertex lies on two lines of this face's boundary)
        corner = []
        for i, t in enumerate(cyc):
            a, b = set(cyc[i - 1]), set(cyc[(i + 1) % len(cyc)])
            if len((a & set(t)) - {p}) and len((b & set(t)) - {p}) and ((a & set(t)) - {p}) != ((b & set(t)) - {p}): corner.append(t)
        faces.append(corner)
    V = {t for f in faces for t in f}
    edges = {}
    for p, f in enumerate(faces):
        for i in range(len(f)):
            t1, t2 = f[i], f[(i + 1) % len(f)]
            q = ((set(t1) & set(t2)) - {p}).pop(); edges.setdefault(frozenset((t1, t2)), set()).update({p, q})
    E = len(edges); F = n; chi = len(V) - E + F
    share = {}
    for e, pq in edges.items(): share[frozenset(pq)] = share.get(frozenset(pq), 0) + 1
    return dict(faces=faces, V=len(V), E=E, F=F, genus=(2 - chi) / 2, share=share, coords={t: tri[t] for t in V})
