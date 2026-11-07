"""Families of plane arrangements to sample. Each sampler returns (planes, params) with planes as rows (a, b, c, d) of
ax + by + cz = d. Symmetric families follow the two known eight-faced examples; families whose symmetry forces four planes
through a point are slightly perturbed (eps), so every hit is an honest, possibly asymmetric, polyhedron.
  C2_7   : 1 + 3 pairs under a half-turn            (control: must find the Szilassi type)
  S4_8   : 4 + 4 under Mizhaev's rotoreflection T    D2_8 : 4 + 4 under the three coordinate half-turns
  C2_8   : 4 pairs under a half-turn                  ASYM_8 : 8 random planes
  S4_11  : 4 + 4 + 2 + 1 (perturbed)                  D2_11 : 4 + 4 + 2 + 1 (perturbed)
  C2_11  : 1 + 5 pairs under a half-turn              ASYM_11: 11 random planes"""
import numpy as np
def nvec(a, b): return np.array([np.sin(a) * np.cos(b), np.sin(a) * np.sin(b), np.cos(a)])
def rand_normal(rng): v = rng.normal(size=3); return v / np.linalg.norm(v)
T = np.array([[0, 1, 0], [-1, 0, 0], [0, 0, -1.]])
Rx, Ry, Rz = np.diag([1, -1, -1.]), np.diag([-1, 1, -1.]), np.diag([-1, -1, 1.])
def orbit(n, d, mats): return [[*(M @ n), d] for M in mats]
S4 = [np.linalg.matrix_power(T, k) for k in range(4)]; D2 = [np.eye(3), Rz, Rx, Ry]; C2 = [np.eye(3), Rz]
def perturb(P, rng, eps):
    P = np.array(P, float); P[:, :3] += eps * rng.normal(size=(len(P), 3)); P[:, 3] += eps * rng.normal(size=len(P)); return P
def fixed_c2(rng):
    """a plane invariant under the half-turn about z: perpendicular to the axis, or containing it"""
    if rng.random() < 0.5: return [0, 0, 1., rng.uniform(-1, 1)], 'perp'
    phi = rng.uniform(0, np.pi); return [np.cos(phi), np.sin(phi), 0, 0.], 'contain'
def s_C2_7(rng):
    f, kind = fixed_c2(rng); P = [f]
    for _ in range(3): P += orbit(rand_normal(rng), rng.uniform(-2, 2), C2)
    return np.array(P), dict(fixed=kind)
def s_S4_8(rng):
    p = [rng.uniform(0.05, np.pi - 0.05), rng.uniform(0.05, np.pi - 0.05), rng.uniform(0, 2 * np.pi), rng.uniform(-3, 3)]
    return np.array(orbit(nvec(p[0], 0), 1., S4) + orbit(nvec(p[1], p[2]), p[3], S4)), dict(p=p)
def s_D2_8(rng):
    return np.array(orbit(rand_normal(rng), 1., D2) + orbit(rand_normal(rng), rng.uniform(-3, 3), D2)), {}
def s_C2_8(rng):
    P = []
    for _ in range(4): P += orbit(rand_normal(rng), rng.uniform(-2, 2), C2)
    return np.array(P), {}
def s_ASYM(n):
    def f(rng): return np.array([[*rand_normal(rng), rng.uniform(-1, 1)] for _ in range(n)]), {}
    return f
def s_S4_11(rng, eps=1e-4):
    P, par = s_S4_8(rng); phi = rng.uniform(0, np.pi / 2)
    P = list(P) + [[np.cos(phi), np.sin(phi), 0, 0], [-np.sin(phi), np.cos(phi), 0, 0], [0, 0, 1, 0]]
    return perturb(P, rng, eps), par
def s_D2_11(rng, eps=1e-4):
    P, _ = s_D2_8(rng); psi = rng.uniform(0.05, np.pi / 2 - 0.05); axis = int(rng.integers(3))
    pair = [[0, np.cos(psi), np.sin(psi), 0], [0, np.cos(psi), -np.sin(psi), 0]]                 # planes through the x-axis
    perm = [[0, 1, 2], [1, 2, 0], [2, 0, 1]][axis]; pair = [[q[perm[0]], q[perm[1]], q[perm[2]], 0] for q in pair]
    single = [[0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0]][int(rng.integers(3))]
    return perturb(list(P) + pair + [single], rng, eps), dict(axis=axis)
def s_C2_11(rng):
    f, kind = fixed_c2(rng); P = [f]
    for _ in range(5): P += orbit(rand_normal(rng), rng.uniform(-2, 2), C2)
    return np.array(P), dict(fixed=kind)
FAMILIES = dict(C2_7=s_C2_7, S4_8=s_S4_8, D2_8=s_D2_8, C2_8=s_C2_8, ASYM_8=s_ASYM(8),
                S4_11=s_S4_11, D2_11=s_D2_11, C2_11=s_C2_11, ASYM_11=s_ASYM(11))


# --- ladder towards eleven faces: nine and ten faces (extensions of the known 4 + 4 examples) ---
SINGLES = [[0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0]]
def _pair_through_axis_D2(rng):
    psi = rng.uniform(0.05, np.pi / 2 - 0.05); axis = int(rng.integers(3))
    pair = [[0, np.cos(psi), np.sin(psi), 0], [0, np.cos(psi), -np.sin(psi), 0]]
    perm = [[0, 1, 2], [1, 2, 0], [2, 0, 1]][axis]
    return [[q[perm[0]], q[perm[1]], q[perm[2]], 0] for q in pair], axis
def s_S4_9(rng):
    P, par = s_S4_8(rng); return np.array(list(P) + [[0, 0, 1., 0]]), par
def s_D2_9(rng):
    P, _ = s_D2_8(rng); k = int(rng.integers(3)); return np.array(list(P) + [SINGLES[k]], float), dict(single=k)
def s_C2_9(rng):
    f, kind = fixed_c2(rng); P = [f]
    for _ in range(4): P += orbit(rand_normal(rng), rng.uniform(-2, 2), C2)
    return np.array(P), dict(fixed=kind)
def s_S4_10(rng, eps=1e-4):
    P, par = s_S4_8(rng); phi = rng.uniform(0, np.pi / 2)
    return perturb(list(P) + [[np.cos(phi), np.sin(phi), 0, 0], [-np.sin(phi), np.cos(phi), 0, 0]], rng, eps), par
def s_D2_10(rng, eps=1e-4):
    P, _ = s_D2_8(rng); pair, axis = _pair_through_axis_D2(rng); return perturb(list(P) + pair, rng, eps), dict(axis=axis)
def s_C2_10(rng):
    P = []
    for _ in range(5): P += orbit(rand_normal(rng), rng.uniform(-2, 2), C2)
    return np.array(P), {}
FAMILIES.update(S4_9=s_S4_9, D2_9=s_D2_9, C2_9=s_C2_9, S4_10=s_S4_10, D2_10=s_D2_10, C2_10=s_C2_10)
