"""Exact verification that f(10) >= 26 (construction of 2026-08-06,
DIMENSIONS_1_TO_10.md): the even-parity half-cube (16 points, squared
radius 5/8, in R^5) orthogonally joined with the Petersen equiangular
set (10 points, squared radius 3/8, in R^5).

Everything is verified in exact rational arithmetic at the Gram level:
  (1) the 26x26 Gram matrix G of the union is PSD of rank <= 10
      (block-diagonal; each block is checked PSD of rank 5 by exact
      LDL^T / pivoted elimination over Q);
  (2) all pairwise squared distances G_ii + G_jj - 2 G_ij are computed
      exactly; among any three points, some pair has squared distance 1;
  (3) all 26 points are distinct (no zero squared distance).
PSD + rank <= 10 guarantees a realization by 26 points in R^10 whose
unit-distance pattern matches G exactly.
"""
from fractions import Fraction as Fr
from itertools import combinations, product

def gram_halfcube():
    """16 even-parity sign vectors scaled by 1/(2*sqrt(2)): inner
    products are rational: <x,y> = (5 - 2*hamming)/8."""
    V = [v for v in product((1, -1), repeat=5)
         if v.count(-1) % 2 == 0]
    assert len(V) == 16
    def ip(u, v):
        return Fr(sum(a * b for a, b in zip(u, v)), 8)
    return [[ip(u, v) for v in V] for u in V], V

def gram_petersen():
    """(3/8) * B with B = (4I + 2A - J)/3 on the Kneser(5,2) Petersen."""
    verts = list(combinations(range(5), 2))
    assert len(verts) == 10
    def adj(a, b):
        return len(set(a) & set(b)) == 0
    G = []
    for a in verts:
        row = []
        for b in verts:
            if a == b:
                B = Fr(1)
            elif adj(a, b):
                B = Fr(1, 3)
            else:
                B = Fr(-1, 3)
            row.append(Fr(3, 8) * B)
        G.append(row)
    return G, verts

def psd_rank(G):
    """exact symmetric pivoted elimination: returns rank if PSD, else None"""
    n = len(G)
    M = [row[:] for row in G]
    rank = 0
    used = [False] * n
    for _ in range(n):
        # pick a positive diagonal pivot
        piv = None
        for i in range(n):
            if not used[i] and M[i][i] > 0:
                piv = i; break
        if piv is None:
            # all remaining diagonal entries must be zero, and then all
            # remaining off-diagonals must be zero for PSD
            for i in range(n):
                if used[i]: continue
                if M[i][i] != 0: return None
                for j in range(n):
                    if not used[j] and M[i][j] != 0: return None
            return rank
        used[piv] = True
        rank += 1
        d = M[piv][piv]
        for i in range(n):
            if used[i]: continue
            f = M[i][piv] / d
            for j in range(n):
                if used[j]: continue
                M[i][j] -= f * M[piv][j]
        for j in range(n):
            M[piv][j] = M[j][piv] = Fr(0)
    return rank

def main():
    GX, VX = gram_halfcube()
    GY, VY = gram_petersen()
    rx = psd_rank(GX)
    ry = psd_rank(GY)
    assert rx is not None and rx <= 5, f"half-cube Gram not PSD rank<=5: {rx}"
    assert ry is not None and ry <= 5, f"Petersen Gram not PSD rank<=5: {ry}"
    print(f"half-cube Gram: PSD, rank {rx}; Petersen Gram: PSD, rank {ry}")
    # radii
    assert all(GX[i][i] == Fr(5, 8) for i in range(16))
    assert all(GY[i][i] == Fr(3, 8) for i in range(10))
    # union Gram: block diagonal (orthogonal join)
    n = 26
    G = [[Fr(0)] * n for _ in range(n)]
    for i in range(16):
        for j in range(16):
            G[i][j] = GX[i][j]
    for i in range(10):
        for j in range(10):
            G[16 + i][16 + j] = GY[i][j]
    # squared distances
    def d2(i, j):
        return G[i][i] + G[j][j] - 2 * G[i][j]
    for i, j in combinations(range(n), 2):
        assert d2(i, j) > 0, f"points {i},{j} coincide"
    # cross pairs all unit
    for i in range(16):
        for j in range(16, 26):
            assert d2(i, j) == 1
    # almost-equidistance over all triples
    bad = 0
    for a, b, c in combinations(range(n), 3):
        if not (d2(a, b) == 1 or d2(a, c) == 1 or d2(b, c) == 1):
            bad += 1
    assert bad == 0, f"{bad} triples without a unit pair"
    print("union: 26 distinct points, rank <= 10 PSD Gram, every triple "
          "contains a unit pair.")
    print("==> f(10) >= 26, verified in exact rational arithmetic.")

if __name__ == "__main__":
    main()
