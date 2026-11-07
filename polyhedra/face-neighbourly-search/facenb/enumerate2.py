"""Complete enumeration of the combinatorial candidates, all multiplicity patterns at once. A face-neighbourly trivalent
polyhedron with n faces is dual to a triangulation T of a closed orientable surface with n vertices in which every pair of
vertices is joined (edge copies allowed, X extra copies in total) and no triple of vertices spans two triangles (three planes
meet in one point). Encoding: edge copies (u,v,k), k = 0..X (k = 0 always present, k >= 1 optional, ordered); a variable for
every triangle with chosen copies of its three sides; every present copy in exactly two triangles, absent copies in none;
every triple at most once; consistent orientation; every vertex link one cycle (lazy cuts). Symmetry breaking: vertex 0 has
no extra copies (one always exists, X extra copies touch at most 2X vertices) and its link reads 1..n-1."""
import itertools
from pysat.card import CardEnc, EncType
from pysat.solvers import Solver
class Enumerator2:
    def __init__(self, n=9, X=3, solver='cadical153'):
        self.n, self.X = n, X; top = [0]
        def new(): top[0] += 1; return top[0]
        pairs = list(itertools.combinations(range(n), 2))
        self.ex = {(p, k): new() for p in pairs for k in range(1, X + 1)}               # extra copies
        def exists_lit(p, k): return None if k == 0 else self.ex[(p, k)]
        self.T = {}                                                                     # (a,b,c,kab,kbc,kac) -> var
        for a, b, c in itertools.combinations(range(n), 3):
            for kab, kbc, kac in itertools.product(range(X + 1), repeat=3):
                self.T[(a, b, c, kab, kbc, kac)] = new()
        cl = []
        for (a, b, c, kab, kbc, kac), t in self.T.items():                               # a triangle needs its copies
            for p, k in (((a, b), kab), ((b, c), kbc), ((a, c), kac)):
                if k: cl.append([-t, self.ex[(p, k)]])
        for p in pairs:
            for k in range(2, X + 1): cl.append([-self.ex[(p, k)], self.ex[(p, k - 1)]])  # copies in order
        enc = CardEnc.equals(lits=list(self.ex.values()), bound=X, top_id=top[0], encoding=EncType.seqcounter); top[0] = max(top[0], enc.nv); cl += enc.clauses
        uses = {}
        for key, t in self.T.items():
            a, b, c, kab, kbc, kac = key
            for p, k in (((a, b), kab), ((b, c), kbc), ((a, c), kac)): uses.setdefault((p, k), []).append((key, t))
        self.orient = {key: new() for key in self.T}
        for (p, k), lst in uses.items():                                                # copy in exactly two triangles if present
            lits = [t for _, t in lst]
            enc = CardEnc.atmost(lits=lits, bound=2, top_id=top[0], encoding=EncType.seqcounter); top[0] = max(top[0], enc.nv); cl += enc.clauses
            enc = CardEnc.atleast(lits=lits, bound=2, top_id=top[0], encoding=EncType.seqcounter); top[0] = max(top[0], enc.nv)
            g = [] if k == 0 else [-self.ex[(p, k)]]
            cl += [c_ + g for c_ in enc.clauses]
            u, v = p
            for (k1, t1), (k2, t2) in itertools.combinations(lst, 2):                    # opposite traversal
                d1, d2 = self._dir(k1, u, v), self._dir(k2, u, v); o1, o2 = self.orient[k1], self.orient[k2]; gg = [-t1, -t2]
                cl += [gg + [o1, o2], gg + [-o1, -o2]] if d1 == d2 else [gg + [o1, -o2], gg + [-o1, o2]]
        for a, b, c in itertools.combinations(range(n), 3):                              # a triple at most once
            lits = [t for key, t in self.T.items() if key[:3] == (a, b, c)]
            enc = CardEnc.atmost(lits=lits, bound=1, top_id=top[0], encoding=EncType.seqcounter); top[0] = max(top[0], enc.nv); cl += enc.clauses
        for u in range(1, n):                                                            # symmetry breaking at vertex 0
            for k in range(1, X + 1): cl.append([-self.ex[((0, u), k)]])
        ring = list(range(1, n)); L = len(ring); want = {frozenset((ring[i], ring[(i + 1) % L])) for i in range(L)}
        for key, t in self.T.items():
            if key[0] == 0 and frozenset(key[1:3]) not in want: cl.append([-t])
        self.s = Solver(name=solver, bootstrap_with=cl); self.uses = uses
    @staticmethod
    def _dir(key, u, v):
        a, b, c = key[:3]; cyc = (a, b, c); i = cyc.index(u); return 1 if cyc[(i + 1) % 3] == v else -1
    @staticmethod
    def copies(key):
        a, b, c, kab, kbc, kac = key; return [((a, b), kab), ((b, c), kbc), ((a, c), kac)]
    def links(self, chosen):
        """cuts for vertex links that are not one cycle; nodes of the link of v are the edge copies at v"""
        cuts = []
        for v in range(self.n):
            at = [key for key in chosen if v in key[:3]]; adj = {}
            for key in at:
                cs = [c for c in self.copies(key) if v in c[0]]
                adj.setdefault(cs[0], []).append((cs[1], key)); adj.setdefault(cs[1], []).append((cs[0], key))
            seen = set()
            for st in adj:
                if st in seen: continue
                comp, stack = set(), [st]
                while stack:
                    x = stack.pop()
                    if x in seen: continue
                    seen.add(x)
                    for y, key in adj[x]: comp.add(key); stack.append(y)
                if len(comp) < len(at): cuts.append([-self.T[key] for key in comp])
        return cuts
    def dual(self, chosen):
        """the polyhedral map: face v = cycle of the triangles (as vertex triples) around v"""
        faces = []
        for v in range(self.n):
            at = [key for key in chosen if v in key[:3]]; by = {}
            for key in at:
                for c in self.copies(key):
                    if v in c[0]: by.setdefault(c, []).append(key)
            cyc = [at[0]]; prev = None
            while len(cyc) < len(at):
                x = cyc[-1]; cs = [c for c in self.copies(x) if v in c[0]]
                nxt = [y for c in cs for y in by[c] if y != x and y != prev and y not in cyc]
                prev = x; cyc.append(nxt[0])
            faces.append([tuple(key[:3]) for key in cyc])
        return faces
    def solutions(self, assumptions=()):
        while self.s.solve(assumptions=list(assumptions)):
            m = set(l for l in self.s.get_model() if l > 0); chosen = [key for key, t in self.T.items() if t in m]
            cuts = self.links(chosen)
            if cuts:
                for c in cuts: self.s.add_clause(c)
                continue
            yield chosen
            self.s.add_clause([-self.T[key] for key in chosen])
def canonical_map(faces):
    """normal form of a dual map (faces = cycles of vertex triples): relabel so that a face with no doubled edges (8 distinct
    neighbours) is 0 and its neighbours in cyclic order are 1..n-1; normalise every face cycle; take the smallest"""
    n = len(faces); best = None
    for v, f in enumerate(faces):
        nb = []
        for i in range(len(f)):
            common = set(f[i]) & set(f[(i + 1) % len(f)]) - {v}; nb.append(common.pop())
        if len(set(nb)) != n - 1 or len(nb) != n - 1: continue
        for st in range(n - 1):
            for d in (1, -1):
                seq = [nb[(st + d * i) % (n - 1)] for i in range(n - 1)]; lab = {v: 0}; lab.update({w: i + 1 for i, w in enumerate(seq)})
                rel = {}
                for w, g in enumerate(faces):
                    cyc = [tuple(sorted(lab[x] for x in t)) for t in g]
                    rots = [cyc[i:] + cyc[:i] for i in range(len(cyc))]; rots += [r[::-1] for r in rots]
                    rel[lab[w]] = min(tuple(r) for r in rots)
                F = tuple(rel[i] for i in range(n))
                if best is None or F < best: best = F
    return best

def add_cube_switches(E):
    """indicator a_i (i = 1..n-1): the triangle at vertex 0 on the link edge (i, i+1) uses copy 0 of that edge (a_i) or an extra
    copy (not a_i). Returns the list of indicator variables; a cube fixes all of them by assumptions."""
    n = E.n; sw = []; nv = E.s.nof_vars() + 1
    for i in range(1, n):
        j = i % (n - 1) + 1; a, b = sorted((i, j)); ai = nv; nv += 1; sw.append(ai)
        for key, t in E.T.items():
            if key[0] == 0 and set(key[1:3]) == {a, b}:
                k = {frozenset((key[0], key[1])): key[3], frozenset((key[1], key[2])): key[4], frozenset((key[0], key[2])): key[5]}[frozenset((a, b))]
                E.s.add_clause([-ai, -t] if k >= 1 else [ai, -t])
    return sw
def run_cubes(cube_ids, n=9, X=3, time_budget=None):
    import time
    E = Enumerator2(n=n, X=X); sw = add_cube_switches(E); found = set(); labelled = 0; done = 0; t0 = time.time()
    for cid in cube_ids:
        assum = [sw[i] if not (cid >> i) & 1 else -sw[i] for i in range(len(sw))]
        for ch in E.solutions(assumptions=assum):
            labelled += 1; found.add(canonical_map(E.dual(ch)))
        done += 1
        if time_budget and time.time() - t0 > time_budget: break
    return found, labelled, done
