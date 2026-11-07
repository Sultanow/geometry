"""Integer face planes with the same combinatorial type, verified in exact rational arithmetic."""
import numpy as np
from .arrsat import solve, describe
from .combtype import same_type
from .exact import verify_exact
def integerize(planes, faces, max_scale):
    P = np.array(planes, float); P = P / np.abs(P).max(axis=1, keepdims=True)
    for K in range(4, max_scale + 1, 2):
        Q = np.rint(K * P).astype(int)
        if any(not q[:3].any() for q in Q): continue
        try: A, sols = solve(Q.tolist(), max_solutions=8, max_rounds=400, time_limit=60)
        except ValueError: continue
        for ch in sols:
            d = describe(A, ch)
            if not same_type(d['faces'], faces): continue
            r = verify_exact(Q.tolist(), d['faces'])
            ok = all(r[k] for k in ('no_parallel', 'general_position', 'simple', 'closed', 'orientable', 'manifold_links', 'embedded', 'face_neighbourly'))
            if ok: return Q.tolist(), d['faces'], r, K
    return None
