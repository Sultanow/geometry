"""The known face-neighbourly polyhedra with eight faces, for classification of hits.
Seven faces: every face-neighbourly polyhedron with 7 faces and genus 1 has the combinatorial type of the Szilassi polyhedron
(dual of the unique 7-vertex triangulation of the torus)."""
# Mizhaev (arXiv:2609.17700): integer vertices and nonagonal faces (0-based)
MIZHAEV_VERTICES = [[-72.0, 84.0, 18.0], [-36.0, 112.0, 102.0], [72.0, -84.0, 18.0], [36.0, -112.0, 102.0], [0.0, 300.0, 234.0], [-84.0, 48.0, -18.0], [9.0, 147.0, 207.0], [0.0, -300.0, 234.0], [84.0, -48.0, -18.0], [-112.0, -36.0, -102.0], [48.0, 84.0, 18.0], [-147.0, 9.0, -207.0], [-9.0, -147.0, 207.0], [84.0, 72.0, -18.0], [112.0, 36.0, -102.0], [-48.0, -84.0, 18.0], [147.0, -9.0, -207.0], [-18.0, 126.0, 144.0], [-126.0, -18.0, -144.0], [-300.0, 0.0, -234.0], [-84.0, -72.0, -18.0], [18.0, -126.0, 144.0], [126.0, 18.0, -144.0], [300.0, 0.0, -234.0]]
MIZHAEV_FACES = [[5, 9, 15, 21, 17, 6, 4, 0, 1], [8, 14, 10, 17, 21, 12, 7, 2, 3], [6, 4, 7, 2, 16, 22, 8, 14, 13], [12, 7, 4, 0, 11, 18, 5, 9, 20], [0, 1, 10, 17, 6, 13, 23, 19, 11], [2, 3, 15, 21, 12, 20, 19, 23, 16], [13, 23, 16, 22, 18, 5, 1, 10, 14], [9, 20, 19, 11, 18, 22, 8, 3, 15]]
# Roest-Vigh (arXiv:2609.32998): face planes A1, B1 and their images under the half-turns about the coordinate axes
_A1 = (3, -4, -2, 5); _B1 = (-2, 5, -5, 3)
_Rz = lambda p: (-p[0], -p[1], p[2], p[3]); _Rx = lambda p: (p[0], -p[1], -p[2], p[3]); _Ry = lambda p: (-p[0], p[1], -p[2], p[3])
ROESTVIGH_PLANES = [_A1, _Rz(_A1), _Rx(_A1), _Ry(_A1), _B1, _Rz(_B1), _Rx(_B1), _Ry(_B1)]
_cache = {}
def known_types():
    """name -> face cycles; Roest-Vigh's cycles are computed from its planes with the SAT construction"""
    if not _cache:
        from .arrsat import solve, describe
        _cache['Mizhaev'] = MIZHAEV_FACES
        A, s = solve(ROESTVIGH_PLANES, max_solutions=1); _cache['Roest-Vigh'] = describe(A, s[0])['faces']
    return _cache
def name_of(faces):
    """name of a known type, or None"""
    from .combtype import same_type
    E = len({frozenset((f[i], f[(i + 1) % len(f)])) for f in faces for i in range(len(f))})
    if len(faces) == 7 and E == 21: return 'Szilassi'        # 7 faces, every pair exactly one edge: dual of the 7-vertex torus
    if len(faces) == 4 and E == 6: return 'tetrahedron'
    for name, f in known_types().items():
        if len(f) == len(faces) and same_type(f, faces): return name
    return None
