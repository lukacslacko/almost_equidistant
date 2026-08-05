"""Exact verification that f(5) >= 16 (Larman-Rogers construction).

S = (1/sqrt(8)) * { v in {+-1}^5 : v has an odd number of +1 coordinates }.

We verify, in exact rational arithmetic on the UNSCALED vectors (squared
distances scale by 1/8):
  (1) |S| = 16 and all points are distinct;
  (2) any two distinct points differ in exactly 2 or 4 coordinates, giving
      squared distance 8 or 16, i.e. exactly 1 or 2 after scaling;
  (3) among ANY three points, some two are at squared distance 8 (unit
      distance after scaling)  — the almost-equidistant property;
  (4) the graph of NON-unit pairs (Hamming distance 4) is triangle-free and
      5-regular on 16 vertices (the Clebsch graph), consistent with (3).

Everything is integer arithmetic; no floating point is involved.
"""
from itertools import combinations, product

def main():
    V = [v for v in product((1, -1), repeat=5) if sum(1 for x in v if x == 1) % 2 == 1]
    assert len(V) == 16 and len(set(V)) == 16, "need 16 distinct points"

    def sqdist(u, v):  # squared Euclidean distance of unscaled vectors
        return sum((a - b) ** 2 for a, b in zip(u, v))

    def hamming(u, v):
        return sum(1 for a, b in zip(u, v) if a != b)

    # (2) two-distance set: squared distances only 8 (unit) or 16
    for u, v in combinations(V, 2):
        d2, h = sqdist(u, v), hamming(u, v)
        assert (d2, h) in ((8, 2), (16, 4)), (u, v, d2, h)

    # (3) almost-equidistant: any 3 points contain a pair at squared dist 8
    bad = 0
    for a, b, c in combinations(V, 3):
        if not (sqdist(a, b) == 8 or sqdist(a, c) == 8 or sqdist(b, c) == 8):
            bad += 1
    assert bad == 0, f"{bad} triples without a unit pair"

    # (4) the non-unit graph is 5-regular and triangle-free (Clebsch graph)
    idx = {v: i for i, v in enumerate(V)}
    nonunit = {i: set() for i in range(16)}
    for u, v in combinations(V, 2):
        if sqdist(u, v) == 16:
            nonunit[idx[u]].add(idx[v]); nonunit[idx[v]].add(idx[u])
    assert all(len(s) == 5 for s in nonunit.values()), "must be 5-regular"
    for i in range(16):
        for j in nonunit[i]:
            assert not (nonunit[i] & nonunit[j]), "triangle in non-unit graph"

    print("OK: 16 points, exact; any 3 contain a unit pair; non-unit graph is "
          "the 5-regular triangle-free Clebsch graph.")
    print("==> f(5) >= 16 (Larman-Rogers), verified in integer arithmetic.")

if __name__ == "__main__":
    main()
