"""Independent numerical spot-check of the simplex identities underlying
the codex/dimension6 exact filters (d6_theory_filters.md):
  for a centered regular unit K7 in R^6 and arbitrary points x,y:
    (a) x = -sum_i u_i(x) q_i          with u_i(x)=||x-q_i||^2-1
    (b) ||u(x)||^2 = 1 + c_x^2/7       with c_x = 1 + sum_i u_i(x)
    (c) M_xy := ||x-y||^2 - 1 = c_x c_y/7 - u(x).u(y)
    (d) points with c_x=0 lie on the sphere ||x||^2=3/7, and the lift
        p -> (sqrt2 p, 1/sqrt7) sends unit distances to orthogonality.
Random x,y over many trials; tolerances at 1e-9 (double precision).
"""
import numpy as np

rng = np.random.default_rng(7)

# centered regular unit simplex q_1..q_7 in R^6
# build from the standard simplex in R^7 projected to sum-zero hyperplane
E = np.eye(7)
C = E - np.ones((7, 7)) / 7
# rows of C live in the sum-zero hyperplane; pairwise distances equal
# |e_i - e_j| = sqrt(2); rescale to unit edge
Q7 = C / np.sqrt(2)
# orthonormal basis of the hyperplane -> coordinates in R^6
B = np.linalg.svd(C)[0][:, :6]
Q = Q7 @ B          # 7 x 6, wait: project each row
Q = (C @ B) / np.sqrt(2)
d = np.linalg.norm(Q[0] - Q[1])
assert abs(d - 1) < 1e-12, d
assert np.allclose(Q.sum(axis=0), 0, atol=1e-12)
assert np.allclose(np.linalg.norm(Q, axis=1) ** 2, 3 / 7, atol=1e-12)

def u(x):
    return np.array([np.linalg.norm(x - Q[i]) ** 2 - 1 for i in range(7)])

ok = True
for _ in range(2000):
    x = rng.normal(size=6) * rng.uniform(0.2, 2.0)
    y = rng.normal(size=6) * rng.uniform(0.2, 2.0)
    ux, uy = u(x), u(y)
    cx, cy = 1 + ux.sum(), 1 + uy.sum()
    if not np.allclose(x, -ux @ Q, atol=1e-9): ok = False; print("(a) FAIL")
    if abs(ux @ ux - (1 + cx * cx / 7)) > 1e-9: ok = False; print("(b) FAIL")
    Mxy = np.linalg.norm(x - y) ** 2 - 1
    if abs(Mxy - (cx * cy / 7 - ux @ uy)) > 1e-9: ok = False; print("(c) FAIL")
print("identities (a),(b),(c): OK" if ok else "FAILURES above")

# (d): sample points with c_x = 0 (i.e. |x|^2=3/7) and unit distance
okd = True
for _ in range(2000):
    x = rng.normal(size=6); x *= np.sqrt(3 / 7) / np.linalg.norm(x)
    assert abs((1 + u(x).sum())) < 1e-9
    # random second point at unit distance on the same sphere (if exists)
    y = rng.normal(size=6); y *= np.sqrt(3 / 7) / np.linalg.norm(y)
    lx = np.append(np.sqrt(2) * x, 1 / np.sqrt(7))
    ly = np.append(np.sqrt(2) * y, 1 / np.sqrt(7))
    dxy = np.linalg.norm(x - y) ** 2
    ip = lx @ ly
    # ip should equal (2 <x,y> + 1/7) and unit distance <=> ip == 0
    if abs(ip - (2 * x @ y + 1 / 7)) > 1e-9: okd = False; print("(d) lift FAIL")
    if abs(dxy - 1) < 1e-12 and abs(ip) > 1e-9: okd = False; print("(d) ortho FAIL")
print("identity (d): OK" if okd else "FAILURES above")
