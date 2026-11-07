"""Oriented-matroid SAT test for a FIXED combinatorial type (face cycles of vertex triples, n faces).
Planes h_i = (a, b, c, -d) in R^4 (plane ax + by + cz = d), plus the plane at infinity h_inf = (0, 0, 0, 1) as element n.
Variables: chi(S) = sign det(h_S) for every 4-subset S of {0..n} (uniform: the realisation space of a valid polyhedron is
open, so all determinants may be assumed non-zero). Vertex v = (i, j, k) lies on the positive side of plane l iff
chi(i,j,k,l) * chi(i,j,k,inf) > 0. Constraints:
  - 3-term Grassmann-Pluecker relations (they characterise uniform chirotopes)
  - every face polygon is simple (non-adjacent sides do not meet; collinear sides on the line to a doubly adjacent face
    do not overlap)
  - no edge pierces a face (edge crosses the plane of f, and the crossing point lies inside polygon f by ray casting along
    a line of the arrangement, half-open rule for vertices on the ray)
UNSAT means: the type has NO realisation as an embedded polyhedron with planar faces (necessary conditions violated).
SAT gives a candidate chirotope (realisability is not implied)."""
import itertools
from pysat.solvers import Solver
def perm_parity(seq):
    seq = list(seq); par = 0
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            if seq[i] > seq[j]: par ^= 1
    return par
class OMSat:
    def __init__(self, n):
        self.n = n; self.inf = n; self.m = n + 1; self.top = 0; self.cl = []; self.var = {}; self._xor = {}; self._and = {}
        for S in itertools.combinations(range(self.m), 4): self.top += 1; self.var[S] = self.top
        self.true = self.new(); self.cl.append([self.true])
    def new(self): self.top += 1; return self.top
    def chi(self, *t):
        """literal true <=> det(h_t) > 0 for the ordered tuple t (4 distinct elements)"""
        S = tuple(sorted(t)); v = self.var[S]; return -v if perm_parity(t) else v
    def xor(self, a, b):
        sgn = 1
        if a < 0: a, sgn = -a, -sgn
        if b < 0: b, sgn = -b, -sgn
        if a > b: a, b = b, a
        if a == b: return -self.true if sgn > 0 else self.true
        key = (a, b)
        if key not in self._xor:
            x = self.new(); self._xor[key] = x
            self.cl += [[-x, a, b], [-x, -a, -b], [x, -a, b], [x, a, -b]]
        x = self._xor[key]; return x if sgn > 0 else -x
    def xor_many(self, lits):
        acc = -self.true
        for l in lits: acc = self.xor(acc, l)
        return acc
    def and_(self, *lits):
        lits = tuple(sorted(set(lits)))
        if len(lits) == 1: return lits[0]
        if lits not in self._and:
            x = self.new(); self._and[lits] = x
            for l in lits: self.cl.append([-x, l])
            self.cl.append([x] + [-l for l in lits])
        return self._and[lits]
    def side(self, v, l):
        """literal true <=> vertex v = (i, j, k) lies on the positive side of plane l (l not in v)"""
        i, j, k = v; return -self.xor(self.chi(i, j, k, l), self.chi(i, j, k, self.inf))
    def grassmann_pluecker(self):
        for X in itertools.combinations(range(self.m), 2):
            rest = [e for e in range(self.m) if e not in X]
            for a, b, c, d in itertools.combinations(rest, 4):
                t1 = -self.xor(self.chi(*X, a, b), self.chi(*X, c, d))      # sign of chi(Xab) chi(Xcd)
                t2 = self.xor(self.chi(*X, a, c), self.chi(*X, b, d))       # sign of -chi(Xac) chi(Xbd)
                t3 = -self.xor(self.chi(*X, a, d), self.chi(*X, b, c))      # sign of chi(Xad) chi(Xbc)
                self.cl += [[-t1, -t2, -t3], [t1, t2, t3]]                   # the three terms are not all of one sign
def chirotope_of(planes):
    """real chirotope of planes (rows a, b, c, d) plus the plane at infinity: dict sorted 4-subset -> +1/-1 (None if zero)"""
    import numpy as np
    H = [np.array([p[0], p[1], p[2], -p[3]], float) for p in planes] + [np.array([0, 0, 0, 1.0])]
    out = {}
    for S in itertools.combinations(range(len(H)), 4):
        d = np.linalg.det(np.array([H[i] for i in S])); out[S] = None if abs(d) < 1e-9 else (1 if d > 0 else -1)
    return out

def add_map(O, faces):
    """geometric conditions of the combinatorial type `faces` (face p = cycle of sorted vertex triples containing p)"""
    T = lambda t: tuple(sorted(t))
    F = [[T(t) for t in f] for f in faces]; V = {t for f in F for t in f}
    def opposite(u, w, l): return O.xor(O.side(u, l), O.side(w, l))
    def nb(p, u, w): return (set(u) & set(w) - {p}).pop()
    # 1. every face polygon is simple
    for p, f in enumerate(F):
        m = len(f); sides = [(f[i], f[(i + 1) % m], nb(p, f[i], f[(i + 1) % m])) for i in range(m)]
        for i in range(m):
            for j in range(i + 2, m):
                if (j + 1) % m == i: continue
                (a1, b1, x1), (a2, b2, x2) = sides[i], sides[j]
                if x1 != x2:                                                    # meeting point Z = (p, x1, x2)
                    c1 = O.true if (x2 in a1 or x2 in b1) else opposite(a1, b1, x2)
                    c2 = O.true if (x1 in a2 or x1 in b2) else opposite(a2, b2, x1)
                    O.cl.append([-c1, -c2])
                else:                                                           # collinear sides on the line to x1
                    def within(z, a, b):
                        c = (set(a) - {p, x1}).pop(); d = (set(b) - {p, x1}).pop()
                        return O.and_(-O.xor(O.side(z, c), O.side(b, c)), -O.xor(O.side(z, d), O.side(a, d)))
                    for z, (a, b) in ((a2, (a1, b1)), (b2, (a1, b1)), (a1, (a2, b2)), (b1, (a2, b2))):
                        O.cl.append([-within(z, a, b)])
    # 2. no edge pierces a face
    edges = {}
    for p, f in enumerate(F):
        for i in range(len(f)): edges.setdefault(frozenset((f[i], f[(i + 1) % len(f)])), set()).add(p)
    for e, ps in edges.items():
        u, w = tuple(e); p, q = sorted(ps)
        for f in range(O.n):
            if f in ps or f in u or f in w: continue
            cond = opposite(u, w, f); y = T((p, q, f))
            if y in V: O.cl.append([-cond]); continue                          # a vertex of f inside the edge
            cs = []
            for k in range(len(F[f])):
                al, be = F[f][k], F[f][(k + 1) % len(F[f])]; x = nb(f, al, be)
                if x == p: continue                                             # side on the ray line
                if x == q:                                                      # side through the ray origin
                    if p not in al and p not in be: O.cl.append([-cond, -opposite(al, be, p)])
                    continue
                sa = O.true if p in al else O.side(al, p); sb = O.true if p in be else O.side(be, p)
                z = al if p in al else (be if p in be else T((f, p, x)))
                cs.append(O.and_(O.xor(sa, sb), O.side(z, q)))
            O.cl.append([-cond, -O.xor_many(cs)])
def om_solve(faces, assumptions=(), fix_sign=True, solver='cadical153', time_limit=None):
    n = len(faces); O = OMSat(n); O.grassmann_pluecker(); add_map(O, faces)
    if fix_sign: O.cl.append([O.var[(0, 1, 2, 3)]])
    s = Solver(name=solver, bootstrap_with=O.cl)
    ok = s.solve(assumptions=list(assumptions)); model = s.get_model() if ok else None; s.delete()
    return ok, O, model

def reorientation_breaking(O):
    """Reorienting a plane (h -> -h) changes nothing geometrically and flips every chi containing it; reorienting all m
    elements at once changes no chi. Choose m-1 four-subsets whose incidence vectors have rank m-1 over GF(2) and require
    them positive: every reorientation class then has exactly one representative."""
    import numpy as np
    m = O.m; rows = []; chosen = []
    def rank(R):
        M = np.array(R, dtype=np.uint8) % 2; r = 0
        for c in range(M.shape[1]):
            piv = next((i for i in range(r, M.shape[0]) if M[i, c]), None)
            if piv is None: continue
            M[[r, piv]] = M[[piv, r]]
            for i in range(M.shape[0]):
                if i != r and M[i, c]: M[i] ^= M[r]
            r += 1
        return r
    for S in itertools.combinations(range(m), 4):
        v = [1 if e in S else 0 for e in range(m)]
        if rank(rows + [v]) > len(rows): rows.append(v); chosen.append(S)
        if len(rows) == m - 1: break
    for S in chosen: O.cl.append([O.var[S]])
    return chosen
def om_solve2(faces, assumptions=(), solver='cadical153', breaking=True):
    n = len(faces); O = OMSat(n); O.grassmann_pluecker(); add_map(O, faces)
    if breaking: reorientation_breaking(O)
    s = Solver(name=solver, bootstrap_with=O.cl); ok = s.solve(assumptions=list(assumptions)); model = s.get_model() if ok else None; s.delete()
    return ok, O, model

def add_reflex(O, faces):
    """At a vertex of an embedded polyhedral surface at most ONE of the three face angles is reflex (two reflex faces at a
    vertex would both contain the extension of their common edge). reflex_p(V): points of the line L_pq just beyond V (away
    from the neighbour U on L_pq) lie inside face p; decided by parity along that ray (half-open rule for vertices on L_pq)."""
    T = lambda t: tuple(sorted(t))
    F = [[T(t) for t in f] for f in faces]
    def nb(p, u, w): return (set(u) & set(w) - {p}).pop()
    R = {}
    for p, f in enumerate(F):
        m = len(f)
        for i in range(m):
            V = f[i]; U = f[i - 1]; q = nb(p, U, V); r = nb(p, V, f[(i + 1) % m])     # side U-V on L_pq, side V-W on L_pr
            cs = []
            for k in range(m):
                al, be = f[k], f[(k + 1) % m]; x = nb(p, al, be)
                if x == q or V in (al, be) or x == r: continue
                sa = O.true if q in al else O.side(al, q); sb = O.true if q in be else O.side(be, q)
                z = al if q in al else (be if q in be else T((p, q, x)))
                beyond = O.xor(O.side(z, r), O.side(U, r))
                cs.append(O.and_(O.xor(sa, sb), beyond))
            R[(p, V)] = O.xor_many(cs)
    for V in {t for f in F for t in f}:
        a, b, c = V
        for x, y in ((a, b), (a, c), (b, c)): O.cl.append([-R[(x, V)], -R[(y, V)]])
    return R
_old_om_solve2 = om_solve2
def om_solve3(faces, assumptions=(), solver='cadical153', breaking=True):
    n = len(faces); O = OMSat(n); O.grassmann_pluecker(); add_map(O, faces); add_reflex(O, faces)
    if breaking: reorientation_breaking(O)
    s = Solver(name=solver, bootstrap_with=O.cl); ok = s.solve(assumptions=list(assumptions)); model = s.get_model() if ok else None; s.delete()
    return ok, O, model

# ---------- second, independent formulation of the same geometric conditions (for cross-checking) ----------
def add_map_alt(O, faces):
    """as add_map, but the inside test for the crossing point y = (p,q,f) casts the ray along the OTHER line L_fq (instead of
    L_fp) and in the NEGATIVE direction of plane p; face simplicity is identical (it uses no ray)"""
    T = lambda t: tuple(sorted(t))
    F = [[T(t) for t in f] for f in faces]; V = {t for f in F for t in f}
    def opposite(u, w, l): return O.xor(O.side(u, l), O.side(w, l))
    def nb(p, u, w): return (set(u) & set(w) - {p}).pop()
    for p, f in enumerate(F):                                    # 1. simple faces (same as add_map)
        m = len(f); sides = [(f[i], f[(i + 1) % m], nb(p, f[i], f[(i + 1) % m])) for i in range(m)]
        for i in range(m):
            for j in range(i + 2, m):
                if (j + 1) % m == i: continue
                (a1, b1, x1), (a2, b2, x2) = sides[i], sides[j]
                if x1 != x2:
                    c1 = O.true if (x2 in a1 or x2 in b1) else opposite(a1, b1, x2)
                    c2 = O.true if (x1 in a2 or x1 in b2) else opposite(a2, b2, x1)
                    O.cl.append([-c1, -c2])
                else:
                    def within(z, a, b):
                        c = (set(a) - {p, x1}).pop(); d = (set(b) - {p, x1}).pop()
                        return O.and_(-O.xor(O.side(z, c), O.side(b, c)), -O.xor(O.side(z, d), O.side(a, d)))
                    for z, (a, b) in ((a2, (a1, b1)), (b2, (a1, b1)), (a1, (a2, b2)), (b1, (a2, b2))):
                        O.cl.append([-within(z, a, b)])
    edges = {}
    for p, f in enumerate(F):
        for i in range(len(f)): edges.setdefault(frozenset((f[i], f[(i + 1) % len(f)])), set()).add(p)
    for e, ps in edges.items():                                  # 2. no piercing, ray along L_fq towards the negative side of p
        u, w = tuple(e); p0, q0 = sorted(ps); p, q = q0, p0          # roles swapped
        for f in range(O.n):
            if f in ps or f in u or f in w: continue
            cond = opposite(u, w, f); y = T((p, q, f))
            if y in V: O.cl.append([-cond]); continue
            cs = []
            for k in range(len(F[f])):
                al, be = F[f][k], F[f][(k + 1) % len(F[f])]; x = nb(f, al, be)
                if x == p: continue
                if x == q:
                    if p not in al and p not in be: O.cl.append([-cond, -opposite(al, be, p)])
                    continue
                sa = -O.true if p in al else O.side(al, p); sb = -O.true if p in be else O.side(be, p)   # half-open rule, other side
                z = al if p in al else (be if p in be else T((f, p, x)))
                cs.append(O.and_(O.xor(sa, sb), -O.side(z, q)))                                      # negative ray direction
            O.cl.append([-cond, -O.xor_many(cs)])
def add_reflex_alt(O, faces):
    """as add_reflex, but reflexness at V is read off the extension of the NEXT side (on L_pr, beyond V away from W)"""
    T = lambda t: tuple(sorted(t))
    F = [[T(t) for t in f] for f in faces]
    def nb(p, u, w): return (set(u) & set(w) - {p}).pop()
    R = {}
    for p, f in enumerate(F):
        m = len(f)
        for i in range(m):
            V = f[i]; W = f[(i + 1) % m]; r = nb(p, V, W); q = nb(p, f[i - 1], V)          # side V-W on L_pr, side U-V on L_pq
            cs = []
            for k in range(m):
                al, be = f[k], f[(k + 1) % m]; x = nb(p, al, be)
                if x == r or V in (al, be) or x == q: continue
                sa = O.true if r in al else O.side(al, r); sb = O.true if r in be else O.side(be, r)
                z = al if r in al else (be if r in be else T((p, r, x)))
                cs.append(O.and_(O.xor(sa, sb), O.xor(O.side(z, q), O.side(W, q))))
            R[(p, V)] = O.xor_many(cs)
    for V in {t for f in F for t in f}:
        a, b, c = V
        for x, y in ((a, b), (a, c), (b, c)): O.cl.append([-R[(x, V)], -R[(y, V)]])
    return R
def om_solve_alt(faces, assumptions=(), solver='cadical153', breaking=True):
    n = len(faces); O = OMSat(n); O.grassmann_pluecker(); add_map_alt(O, faces); add_reflex_alt(O, faces)
    if breaking: reorientation_breaking(O)
    s = Solver(name=solver, bootstrap_with=O.cl); ok = s.solve(assumptions=list(assumptions)); model = s.get_model() if ok else None; s.delete()
    return ok, O, model
