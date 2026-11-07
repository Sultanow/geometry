"""Guided search. score(planes) = the largest number of pairs of planes that share an edge, over all closed embedded surfaces
built from cells of the arrangement with at most one simple polygon per plane (MaxSAT: the surface constraints are hard,
"pair p,q shares an edge" is soft; face simplicity is enforced lazily with the same cuts as in arrsat). A plane arrangement
carries a face-neighbourly polyhedron iff score = n(n-1)/2. climb() moves the parameters of a family uphill in this score."""
import itertools, numpy as np
from pysat.formula import WCNF
from pysat.examples.rc2 import RC2
from .arrsat import build, face_check
def surface_clauses(A):
    """the hard clauses of arrsat.clauses without the neighbourliness requirement, plus pair indicator variables"""
    cells, seginc, vinc, n = A['cells'], A['seginc'], A['vinc'], A['n']
    x = lambda c: c + 1; nxt = [len(cells) + 1]; cl = []
    def newvar():
        nxt[0] += 1; return nxt[0] - 1
    for s, inc in seginc.items():
        for k in range(1, len(inc) + 1):
            if k == 2: continue
            for ch in itertools.combinations(inc, k): cl.append([-x(c) for c in ch] + [x(c) for c in inc if c not in ch])
    for t, inc in vinc.items():
        byplane = {}
        for (c, sa, sb) in inc: byplane.setdefault(cells[c]['plane'], []).append((c, frozenset((sa, sb))))
        for p, lst in byplane.items():
            for (c1, n1), (c2, n2) in itertools.combinations(lst, 2):
                if not (n1 & n2): cl.append([-x(c1), -x(c2)] + [x(c) for c, _ in lst if c not in (c1, c2)])
        cyc = []; planes = sorted(byplane)
        if len(planes) == 3:
            for a in byplane[planes[0]]:
                for b in byplane[planes[1]]:
                    for c in byplane[planes[2]]:
                        nodes = a[1] | b[1] | c[1]
                        if len(nodes) == 3 and len(a[1] & b[1]) == 1 and len(a[1] & c[1]) == 1 and len(b[1] & c[1]) == 1:
                            cyc.append((frozenset(nodes), (a[0], b[0], c[0])))
        for (n1, c1), (n2, c2) in itertools.combinations(cyc, 2):
            if not (n1 & n2): cl.append([-x(c) for c in c1 + c2])
    pair = {}
    for s, inc in seginc.items():
        p, q, _ = s
        cp = [c for c in inc if cells[c]['plane'] == p]; cq = [c for c in inc if cells[c]['plane'] == q]
        if not cp or not cq: continue
        e = newvar()
        for side in (cp, cq):
            cl.append([-e] + [x(c) for c in side])
            if len(side) == 2: cl.append([-e, -x(side[0]), -x(side[1])])
        pair.setdefault((p, q), []).append(e)
    y = {}
    for pq, es in pair.items():
        y[pq] = newvar(); cl.append([-y[pq]] + es)
    return cl, y
def score(planes, max_cuts=200):
    try: A = build(planes)
    except ValueError: return -1, None
    cl, y = surface_clauses(A); w = WCNF()
    for c in cl: w.append(c)
    for v in y.values(): w.append([v], weight=1)
    rc2 = RC2(w); best = None
    for _ in range(max_cuts):
        m = rc2.compute()
        if m is None: break
        chosen = [c for c in range(len(A['cells'])) if m[c] > 0]
        cuts = face_check(A, chosen)
        if not cuts: best = (sum(1 for v in y.values() if m[v - 1] > 0), chosen); break
        for c in cuts: rc2.add_clause(c)
    rc2.delete()
    return (best[0], (A, best[1])) if best else (0, None)
def climb(make, x0, lo, hi, rng, iters=200, step=0.15, target=None):
    """stochastic hill climbing on the parameters x of make(x) -> planes; sideways moves allowed"""
    x = np.array(x0, float); s, _ = score(make(x)); hist = [s]; evals = 1
    for it in range(iters):
        if target is not None and s >= target: break
        y = np.clip(x + step * (hi - lo) * rng.normal(size=len(x)), lo, hi)
        t, _ = score(make(y)); evals += 1
        if t >= s: x, s = y, t
        hist.append(s)
    return x, s, hist, evals

def connectivity_clauses(A, first_var):
    """every plane: the chosen cells form one region connected through segments (root + reachability levels), and every
    unchosen cell is reachable from outside through unchosen cells (no holes). Returns (clauses, next free variable)."""
    cells = A['cells']; nv = [first_var]; cl = []
    def new():
        nv[0] += 1; return nv[0] - 1
    x = lambda c: c + 1
    for p in range(A['n']):
        C = [c for c, cc in enumerate(cells) if cc['plane'] == p]
        if not C: continue
        segs = {c: set(cells[c]['segs']) for c in C}
        adj = {c: [d for d in C if d != c and segs[c] & segs[d]] for c in C}
        outer = [c for c in C if any(len([d for d in A['seginc'][s] if cells[d]['plane'] == p]) < 2 for s in segs[c])]
        L = len(C)
        root = {c: new() for c in C}
        for a_, b_ in itertools.combinations(C, 2): cl.append([-root[a_], -root[b_]])          # at most one root
        for c in C:
            cl.append([-root[c], x(c)]); cl.append([-x(c)] + [root[d] for d in C])               # a chosen cell implies a root
        reach = {c: [new() for _ in range(L + 1)] for c in C}
        for c in C:
            cl.append([-reach[c][0], root[c]])
            for k in range(1, L + 1):
                cl.append([-reach[c][k], x(c)])
                cl.append([-reach[c][k], reach[c][k - 1]] + [reach[d][k - 1] for d in adj[c]])
            cl.append([-x(c), reach[c][L]])
        out = {c: [new() for _ in range(L + 1)] for c in C}
        for c in C:
            cl.append([-out[c][0], -x(c)])
            if c not in outer: cl.append([-out[c][0]])
            for k in range(1, L + 1):
                cl.append([-out[c][k], -x(c)])
                cl.append([-out[c][k], out[c][k - 1]] + [out[d][k - 1] for d in adj[c]])
            cl.append([x(c), out[c][L]])
    return cl, nv[0]
def score_conn(planes):
    """MaxSAT score with the connectivity encoding (no lazy cuts needed); returns (score, (A, chosen))"""
    try: A = build(planes)
    except ValueError: return -1, None
    cl, y = surface_clauses(A)
    top = max([abs(l) for c in cl for l in c] + list(y.values()))
    cc, _ = connectivity_clauses(A, top + 1)
    w = WCNF()
    for c in cl + cc: w.append(c)
    for v in y.values(): w.append([v], weight=1)
    rc2 = RC2(w); m = rc2.compute(); rc2.delete()
    if m is None: return 0, None
    chosen = [c for c in range(len(A['cells'])) if m[c] > 0]
    assert not face_check(A, chosen), 'connectivity encoding let a non-simple face through'
    return sum(1 for v in y.values() if m[v - 1] > 0), (A, chosen)

def score_timed(planes, timeout):
    """score_conn in a forked child process that is killed after `timeout` seconds (RC2 cannot be interrupted).
    Returns (score, chosen cells or None); score -2 on timeout, -1 for a degenerate arrangement."""
    import multiprocessing as mp
    ctx = mp.get_context('fork'); recv, send = ctx.Pipe(duplex=False)
    def child(P, conn):
        try:
            s, res = score_conn(P); conn.send((s, None if res is None else res[1]))
        except Exception as e:
            conn.send((-3, repr(e)))
        conn.close()
    pr = ctx.Process(target=child, args=(planes, send)); pr.start(); send.close()
    if recv.poll(timeout):
        try: out = recv.recv()
        except EOFError: out = (-3, 'child died')
    else:
        out = (-2, None)
    if pr.is_alive(): pr.kill()
    pr.join(); recv.close()
    return out

def score_anytime(planes, time_limit, solver='g4'):
    """Anytime MaxSAT by linear search: find a surface, then require at least one more connected pair, and so on.
    Returns (score, chosen, proven): the best score found within time_limit seconds, the chosen cells of that surface, and
    whether it is proven optimal. -1 for a degenerate arrangement. Interruptible (CaDiCaL), so no child process is needed."""
    import threading, time as _t
    from pysat.card import ITotalizer
    from pysat.solvers import Solver                                   # Glucose: interruptible (CaDiCaL is not, in PySAT)
    t_end = _t.time() + time_limit
    try: A = build(planes)
    except ValueError: return -1, None, False
    cl, y = surface_clauses(A)
    top = max([abs(l) for c in cl for l in c] + list(y.values()))
    cc, top2 = connectivity_clauses(A, top + 1)
    ys = list(y.values())
    if not ys: return 0, None, True
    # at least `need` pairs  <=>  at most n - need pairs MISSING: a totalizer over the negated pair variables,
    # whose outputs encode upper bounds (rhs[j] is forced true once j + 1 inputs are true; assume -rhs[j] for "at most j")
    n = len(ys); tot = ITotalizer(lits=[-v for v in ys], ubound=n, top_id=top2)
    s = Solver(name=solver, bootstrap_with=cl + cc + tot.cnf.clauses)
    best, best_ch, proven, need = 0, None, False, 1
    while need <= n:
        left = t_end - _t.time()
        if left <= 0: break
        timer = threading.Timer(left, s.interrupt); timer.start()
        ok = s.solve_limited(assumptions=[-tot.rhs[n - need]], expect_interrupt=True); timer.cancel()
        if ok is None: break                                              # time is up: keep the best so far
        if not ok: proven = True; break                                   # no surface with `need` pairs: best is optimal
        m = s.get_model(); val = {abs(l): l > 0 for l in m}
        best = sum(1 for v in ys if val.get(v, False)); best_ch = [c for c in range(len(A['cells'])) if val.get(c + 1, False)]
        need = best + 1
    if need > n: proven = True
    s.delete(); tot.delete()
    if best_ch is not None: assert not face_check(A, best_ch), 'connectivity encoding let a non-simple face through'
    return best, best_ch, proven

def target_genus(n):
    """smallest genus g for which a trivalent surface with n faces has enough edges, 3(n - 2 + 2g) >= n(n-1)/2"""
    g = 0
    while 3 * (n - 2 + 2 * g) < n * (n - 1) // 2: g += 1
    return g
def score_genus(planes, time_limit, g_target=None, solver='g4'):
    """Two-stage anytime score aimed at the right kind of surface. Every plane must carry a face. Stage 1 maximises the number
    of surface vertices (corners where three faces meet); a trivalent surface with n faces and V vertices has genus
    1 - (n - V/2)/2, so genus >= g_target means V >= 2(n - 2 + 2 g_target). If that is not reached, the score is the deficit
    V - V_target (< 0). Stage 2 maximises the connected pairs among surfaces with V >= V_target; the score is that number
    (n(n-1)/2 = a face-neighbourly polyhedron). Returns (score, chosen, info)."""
    import threading, time as _t
    from pysat.card import ITotalizer
    from pysat.solvers import Solver
    t_end = _t.time() + time_limit
    try: A = build(planes)
    except ValueError: return -999, None, dict(stage=0)
    n = A['n']; g_t = target_genus(n) if g_target is None else g_target; V_t = 2 * (n - 2 + 2 * g_t)
    cl, y = surface_clauses(A)
    top = max([abs(l) for c in cl for l in c] + list(y.values()))
    cc, nxt = connectivity_clauses(A, top + 1)
    cells = A['cells']
    for p in range(n): cl.append([c + 1 for c, cc_ in enumerate(cells) if cc_['plane'] == p])          # every plane carries a face
    corners = []
    for t, inc in A['vinc'].items():                                                                     # corner indicators
        byplane = {}
        for (c, sa, sb) in inc: byplane.setdefault(cells[c]['plane'], []).append(c)
        if len(byplane) < 3: continue
        z = nxt; nxt += 1; corners.append(z)
        for p, lst in byplane.items(): cl.append([-z] + [c + 1 for c in lst])
    ys = list(y.values())
    totV = ITotalizer(lits=[-z for z in corners], ubound=len(corners), top_id=nxt)
    totP = ITotalizer(lits=[-v for v in ys], ubound=len(ys), top_id=totV.top_id)
    s = Solver(name=solver, bootstrap_with=cl + cc + totV.cnf.clauses + totP.cnf.clauses)
    def run(assum):
        left = t_end - _t.time()
        if left <= 0: return None
        timer = threading.Timer(left, s.interrupt); timer.start(); ok = s.solve_limited(assumptions=assum, expect_interrupt=True); timer.cancel()
        return ok
    def atleast(tot, lits, k):                                    # assumption "at least k of lits"
        return [] if k <= 0 else ([-tot.rhs[len(lits) - k]] if k <= len(lits) else None)
    info = dict(stage=1, V_target=V_t, g_target=g_t, V=None, proven=False)
    ok = run([])
    if not ok: s.delete(); return (-998 if ok is False else -997), None, info                         # no surface uses all planes
    m = s.get_model(); val = {abs(l): l > 0 for l in m}; V = sum(val.get(z, False) for z in corners); chosen = [c for c in range(len(cells)) if val.get(c + 1, False)]
    while V < V_t:                                                                                      # stage 1: more corners
        ok = run(atleast(totV, corners, V + 1))
        if not ok:
            info.update(V=V, proven=ok is False); s.delete(); return V - V_t, chosen, info
        m = s.get_model(); val = {abs(l): l > 0 for l in m}; V = sum(val.get(z, False) for z in corners); chosen = [c for c in range(len(cells)) if val.get(c + 1, False)]
    info.update(stage=2, V=V)
    P = sum(val.get(v, False) for v in ys); base = atleast(totV, corners, V_t)
    while P < len(ys):                                                                                  # stage 2: more pairs
        ok = run(base + atleast(totP, ys, P + 1))
        if not ok: info['proven'] = ok is False; break
        m = s.get_model(); val = {abs(l): l > 0 for l in m}; P = sum(val.get(v, False) for v in ys); chosen = [c for c in range(len(cells)) if val.get(c + 1, False)]
        info['V'] = sum(val.get(z, False) for z in corners)
    else:
        info['proven'] = True
    s.delete()
    assert not face_check(A, chosen), 'connectivity encoding let a non-simple face through'
    return P, chosen, info

def score_faces(planes, time_limit, min_edges=None, solver='g4'):
    """Three-stage anytime score for the last step to a face-neighbourly polyhedron (every face needs >= n-1 edges).
      genus target not reached           -> -200
      genus ok, some face too small      -> sum_p min(corners_p, n-1) - n(n-1)   (in [-n(n-1), -1]; -2 = one hexagon among 9 faces)
      every face has >= n-1 corners       -> number of connected pairs (n(n-1)/2 = face-neighbourly)
    -500: no closed surface uses all planes; -1000: degenerate arrangement. Returns (score, chosen, info)."""
    import threading, time as _t
    from pysat.card import ITotalizer
    from pysat.solvers import Solver
    t_end = _t.time() + time_limit
    try: A = build(planes)
    except ValueError: return -1000, None, dict(stage=0)
    n = A['n']; k_min = n - 1 if min_edges is None else min_edges; V_t = 2 * (n - 2 + 2 * target_genus(n)); cells = A['cells']
    cl, y = surface_clauses(A)
    top = max([abs(l) for c in cl for l in c] + list(y.values()))
    cc, nxt = connectivity_clauses(A, top + 1)
    for p in range(n): cl.append([c + 1 for c, cc_ in enumerate(cells) if cc_['plane'] == p])
    corners, at_plane = [], {p: [] for p in range(n)}
    for t, inc in A['vinc'].items():
        bp = {}
        for (c, sa, sb) in inc: bp.setdefault(cells[c]['plane'], []).append(c)
        if len(bp) < 3: continue
        z = nxt; nxt += 1; corners.append(z)
        for p, lst in bp.items(): cl.append([-z] + [c + 1 for c in lst]); at_plane[p].append(z)
    def new_tot(lits, top_id): return ITotalizer(lits=[-l for l in lits], ubound=len(lits), top_id=top_id)
    totV = new_tot(corners, nxt); top_id = totV.top_id
    totF = {}
    for p in range(n):
        if at_plane[p]: totF[p] = new_tot(at_plane[p], top_id); top_id = totF[p].top_id
    levels = []                                                     # L_{p,k}: plane p has at least k corners, k = 1..k_min
    lvl_cl = []
    for p in range(n):
        m = len(at_plane[p])
        for k in range(1, k_min + 1):
            top_id += 1; L = top_id; levels.append(L)
            lvl_cl.append([-L, -totF[p].rhs[m - k]] if m >= k else [-L])
    totL = new_tot(levels, top_id); top_id = totL.top_id
    ys = list(y.values()); totP = new_tot(ys, top_id)
    s = Solver(name=solver, bootstrap_with=cl + cc + lvl_cl + totV.cnf.clauses + totL.cnf.clauses + totP.cnf.clauses
                                         + [c for t in totF.values() for c in t.cnf.clauses])
    def atleast(tot, lits, k): return [] if k <= 0 else [-tot.rhs[len(lits) - k]]
    def run(assum):
        left = t_end - _t.time()
        if left <= 0: return None
        tm = threading.Timer(left, s.interrupt); tm.start(); ok = s.solve_limited(assumptions=assum, expect_interrupt=True); tm.cancel(); return ok
    def model():
        val = {abs(l): l > 0 for l in s.get_model()}
        return val, [c for c in range(len(cells)) if val.get(c + 1, False)]
    info = dict(stage=1, V_target=V_t, min_edges=k_min, proven=False)
    gV = atleast(totV, corners, V_t)
    ok = run([])
    if not ok: s.delete(); return (-500 if ok is False else -499), None, info
    ok = run(gV)
    if not ok: info['proven'] = ok is False; s.delete(); return -200, None, info
    val, chosen = model(); info['stage'] = 2
    full = gV + [a for p in range(n) for a in atleast(totF[p], at_plane[p], k_min)] if all(p in totF for p in range(n)) else None
    F8 = sum(val.get(L, False) for L in levels)
    if full is not None:
        okf = run(full)
        if okf: val, chosen = model(); F8 = len(levels)
    while F8 < len(levels):                                        # stage 2: fewer missing corners
        ok = run(gV + atleast(totL, levels, F8 + 1))
        if not ok: info['proven'] = ok is False; s.delete(); return F8 - len(levels), chosen, info
        val, chosen = model(); F8 = sum(val.get(L, False) for L in levels)
    info['stage'] = 3; base = full if full is not None else gV + atleast(totL, levels, len(levels))
    P = sum(val.get(v, False) for v in ys)
    while P < len(ys):                                             # stage 3: more connected pairs
        ok = run(base + atleast(totP, ys, P + 1))
        if not ok: info['proven'] = ok is False; break
        val, chosen = model(); P = sum(val.get(v, False) for v in ys)
    else:
        info['proven'] = True
    s.delete()
    assert not face_check(A, chosen), 'connectivity encoding let a non-simple face through'
    return P, chosen, info
