"""Controls for the combinatorial enumeration. Prints ENUMERATION CONTROLS PASSED or stops at the first failure."""
import sys, random
sys.path.insert(0, '.')
from facenb.enumerate import Enumerator, canonical
from facenb.enumerate2 import Enumerator2, canonical_map
def check(cond, msg):
    print(('ok    ' if cond else 'FAIL  ') + msg)
    if not cond: sys.exit(1)
E = Enumerator(n_quads=0, n=7); S = list(E.solutions()); C = {canonical(s, n=7, count=True) for s in S}
check(len(C) == 1 and list(C)[0][1] == 42, 'K7 on the torus: one class, automorphism group of order 42 (Szilassi\'s dual)')
E2 = Enumerator2(n=7, X=0); S2 = list(E2.solutions()); D = [E2.dual(s) for s in S2]
check(len({canonical_map(d) for d in D}) == 1 and sorted(len(f) for f in D[0]) == [6] * 7, 'complete enumerator, 7 faces: one class, seven hexagons (Szilassi)')
E9 = Enumerator2(n=9, X=3); ch = next(E9.solutions()); d = E9.dual(ch); perm = list(range(9)); random.Random(1).shuffle(perm)
dp = [None] * 9
for v, f in enumerate(d): dp[perm[v]] = [tuple(sorted(perm[x] for x in t)) for t in f]
check(canonical_map(d) == canonical_map(dp), 'normal form of a nine-faced candidate is invariant under relabelling')
check(sum(len(f) for f in d) == 78 and len({t for f in d for t in f}) == 26, 'nine-faced candidate: 26 vertices, 39 edges')
print('ENUMERATION CONTROLS PASSED')
