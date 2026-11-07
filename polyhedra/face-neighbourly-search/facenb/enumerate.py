"""Combinatorial candidates for a face-neighbourly polyhedron with nine faces and genus 3 (trivalent: 26 vertices, 39 edges).
Its dual is a triangulation of the genus-3 surface with 9 vertices in which every pair of vertices is joined and exactly three
edges are doubled. Removing one copy of each doubled edge merges two triangles into a quadrilateral, so the candidates are
the genus-3 embeddings of K9 with 20 triangles and 3 quadrilaterals (other face patterns, a pentagon or a hexagon, follow
the same scheme). SAT encoding: a variable per possible face (84 triangles, 378 quadrilaterals); every edge in exactly two
faces; exactly three quadrilaterals; every vertex link is ONE cycle (lazy subcycle cuts); a consistent orientation (encoded).
Symmetry breaking: the link of vertex 0 is the cycle 1-2-...-8."""
import itertools
from pysat.card import CardEnc, EncType
from pysat.solvers import Solver
N = 9
def faces_catalogue(n=None):
    n = N if n is None else n; F = []
    for t in itertools.combinations(range(n), 3): F.append(tuple(t))
    for a, b, c, d in itertools.combinations(range(n), 4):
        for cyc in ((a, b, c, d), (a, b, d, c), (a, c, b, d)): F.append(cyc)
    return F
def edges_of(f): return [frozenset((f[i], f[(i + 1) % len(f)])) for i in range(len(f))]
def link_edge(f, v):
    """the edge of the link of v contributed by face f (the two neighbours of v along f)"""
    i = f.index(v); return frozenset((f[i - 1], f[(i + 1) % len(f)]))
def direction(f, u, v):
    """+1 if the cycle f traverses the edge u->v, -1 if v->u"""
    i = f.index(u); return 1 if f[(i + 1) % len(f)] == v else -1
class Enumerator:
    def __init__(self, n_quads=3, solver='cadical153', n=None):
        self.n = N if n is None else n; self.F = faces_catalogue(self.n); self.var = {f: i + 1 for i, f in enumerate(self.F)}; top = len(self.F)
        self.orient = {}; cl = []
        for f in self.F: top += 1; self.orient[f] = top
        inc = {}
        for f in self.F:
            for e in edges_of(f): inc.setdefault(e, []).append(f)
        for e, fs in inc.items():                                        # every edge in exactly two faces
            enc = CardEnc.equals(lits=[self.var[f] for f in fs], bound=2, top_id=top, encoding=EncType.seqcounter); top = max(top, enc.nv); cl += enc.clauses
            u, v = tuple(e)
            for f1, f2 in itertools.combinations(fs, 2):                 # opposite traversal of a shared edge
                same = direction(f1, u, v) == direction(f2, u, v)
                o1, o2 = self.orient[f1], self.orient[f2]; g = [-self.var[f1], -self.var[f2]]
                cl += [g + [o1, o2], g + [-o1, -o2]] if same else [g + [o1, -o2], g + [-o1, o2]]
        quads = [self.var[f] for f in self.F if len(f) == 4]
        enc = CardEnc.equals(lits=quads, bound=n_quads, top_id=top, encoding=EncType.seqcounter); top = max(top, enc.nv); cl += enc.clauses
        ring = list(range(1, self.n)); L = len(ring)                    # symmetry breaking: link(0) = 1-2-...-(n-1)-1
        want = {frozenset((ring[i], ring[(i + 1) % L])) for i in range(L)}
        for f in self.F:
            if 0 in f and link_edge(f, 0) not in want: cl.append([-self.var[f]])
        cl.append([self.orient[next(f for f in self.F if f == (0, 1, 2))], self.orient[next(f for f in self.F if f == (0, 1, 2))]] if False else [])
        self.s = Solver(name=solver, bootstrap_with=[c for c in cl if c])
        self.inc = inc
    def links_ok(self, chosen):
        """returns list of cuts (one per subcycle) or [] if every link is a single cycle"""
        cuts = []
        for v in range(self.n):
            fs = [f for f in chosen if v in f]; adj = {}
            for f in fs:
                a, b = tuple(link_edge(f, v)); adj.setdefault(a, []).append((b, f)); adj.setdefault(b, []).append((a, f))
            seen = set()
            for start in adj:
                if start in seen: continue
                comp_f, stack = set(), [start]
                while stack:
                    x = stack.pop()
                    if x in seen: continue
                    seen.add(x)
                    for y, f in adj[x]: comp_f.add(f); stack.append(y)
                if len(seen) < len(adj) or len(comp_f) < len(fs): cuts.append([-self.var[f] for f in comp_f])
        return cuts
    def solutions(self, limit=None):
        found = 0
        while self.s.solve():
            m = set(l for l in self.s.get_model() if l > 0)
            chosen = [f for f in self.F if self.var[f] in m]
            cuts = self.links_ok(chosen)
            if cuts:
                for c in cuts: self.s.add_clause(c)
                continue
            yield chosen
            self.s.add_clause([-self.var[f] for f in chosen])           # next solution (faces only; orientation is then forced)
            found += 1
            if limit and found >= limit: return

def link_cycle(faces, v):
    """neighbours of v in the cyclic order of its link (a Hamiltonian cycle on the other eight vertices)"""
    adj = {}
    for f in faces:
        if v in f:
            a, b = tuple(link_edge(f, v)); adj.setdefault(a, []).append(b); adj.setdefault(b, []).append(a)
    start = min(adj); cyc = [start]; prev = None
    while len(cyc) < len(adj):
        x = cyc[-1]; nxt = [y for y in adj[x] if y != prev and (len(cyc) < 2 or y != cyc[-2])]
        y = nxt[0] if nxt else adj[x][0]
        if y == start: break
        prev = x; cyc.append(y)
    return cyc
def canon_face(f):
    if len(f) == 3: return tuple(sorted(f))
    rots = [f[i:] + f[:i] for i in range(len(f))]; rots += [r[::-1] for r in rots]
    return min(tuple(r) for r in rots)
def canonical(faces, n=None, count=False):
    """exact normal form: relabel so that some vertex is 0 and its link reads 1..n-1; take the smallest face list.
    With count=True also return how many relabellings attain it (= order of the automorphism group)"""
    n = N if n is None else n; L = n - 1; best = None; hits = 0
    for v in range(n):
        cyc = link_cycle(faces, v)
        for start in range(L):
            for d in (1, -1):
                seq = [cyc[(start + d * i) % L] for i in range(L)]
                lab = {v: 0}; lab.update({w: i + 1 for i, w in enumerate(seq)})
                F = tuple(sorted(canon_face(tuple(lab[x] for x in f)) for f in faces))
                if best is None or F < best: best, hits = F, 1
                elif F == best: hits += 1
    return (best, hits) if count else best
def cubes():
    """all ways to choose the faces at vertex 0 (link 1..8 fixed): per link edge a triangle or a quadrilateral 0-i-x-(i+1)"""
    edges = [(i, i % 8 + 1) for i in range(1, 9)]
    out = []
    for k in range(4):
        for qs in itertools.combinations(range(8), k):
            choices = [[None] if j not in qs else [x for x in range(1, 9) if x not in edges[j]] for j in range(8)]
            for xs in itertools.product(*choices): out.append(tuple(xs))
    return out
def cube_assumptions(E, cube):
    lits = []
    for j, x in enumerate(cube):
        i, i1 = j + 1, (j + 1) % 8 + 1
        if x is None: f = tuple(sorted((0, i, i1)))
        else: f = next(g for g in E.F if len(g) == 4 and canon_face(g) == canon_face((0, i, x, i1)))
        lits.append(E.var[f])
    return lits
def enumerate_cubes(cube_list, time_budget=None):
    """complete enumeration of the given cubes; returns the set of canonical forms and counts"""
    import time
    E = Enumerator(); t0 = time.time(); found = set(); labelled = 0; done = 0
    for cube in cube_list:
        assum = cube_assumptions(E, cube)
        while E.s.solve(assumptions=assum):
            m = set(l for l in E.s.get_model() if l > 0); chosen = [f for f in E.F if E.var[f] in m]
            cuts = E.links_ok(chosen)
            if cuts:
                for c in cuts: E.s.add_clause(c)
                continue
            labelled += 1; found.add(canonical(chosen)); E.s.add_clause([-E.var[f] for f in chosen])
        done += 1
        if time_budget and time.time() - t0 > time_budget: break
    return found, labelled, done

def expand(faces):
    """all triangulations T of a K9 embedding with quadrilaterals: one diagonal per quadrilateral (the doubled edges). Returns
    the dual polyhedral maps as face cycles of vertex triples (the format of arrsat.describe), only where every triple of faces
    meets at most once (three planes meet in one point)."""
    tris0 = [tuple(sorted(f)) for f in faces if len(f) == 3]; quads = [f for f in faces if len(f) == 4]
    out = []
    for diag in itertools.product((0, 1), repeat=len(quads)):
        tris = [(t, [frozenset((t[0], t[1])), frozenset((t[1], t[2])), frozenset((t[0], t[2]))]) for t in tris0]
        diag_ends = {}
        for qi, (q, dg) in enumerate(zip(quads, diag)):
            a, b, c, d = q if dg == 0 else (q[1], q[2], q[3], q[0])           # diagonal a-c
            eid = ('d', qi); diag_ends[eid] = (a, c)
            tris.append((tuple(sorted((a, b, c))), [frozenset((a, b)), frozenset((b, c)), eid]))
            tris.append((tuple(sorted((a, c, d))), [eid, frozenset((c, d)), frozenset((a, d))]))
        if len({t for t, _ in tris}) < len(tris): continue                        # a triple of faces would meet twice
        dual = []
        for v in range(N):                                                        # face v of the polyhedron: triangles around v
            at = [(t, es) for t, es in tris if v in t]
            out_ok = True
            # adjacency of triangles around v through shared edges incident to v
            inc = {}
            for k, (t, es) in enumerate(at):
                for e in es:
                    if isinstance(e, frozenset):
                        if v in e: inc.setdefault(e, []).append(k)
                    else:
                        qa = diag_ends[e]
                        if v in qa: inc.setdefault(e, []).append(k)
            if any(len(x) != 2 for x in inc.values()): out_ok = False
            if not out_ok: break
            nb = {k: [] for k in range(len(at))}
            for e, (k1, k2) in inc.items(): nb[k1].append(k2); nb[k2].append(k1)
            cyc = [0]; prev = None
            while len(cyc) < len(at):
                nxt = [k for k in nb[cyc[-1]] if k != prev and k not in cyc]
                if not nxt: break
                prev = cyc[-1]; cyc.append(nxt[0])
            if len(cyc) != len(at): out_ok = False; break
            dual.append([at[k][0] for k in cyc])
        else:
            out.append(dual)
    return out
