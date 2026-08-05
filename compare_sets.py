"""Cross-validate the two independent generation routes at n=13:
(A) in-house enumerator (enumaeq5.c): ALL abstract almost-equidistant
    graphs in R^5 on 13 vertices (7623096, matching BPSSV Table 3),
    filtered here to the minimal ones (expect 242, Table 2);
(B) triangleramsey + filter_mtf: minimal graphs directly (242).
The two 242-element sets must agree up to isomorphism, element for element.
"""
import struct, sys
sys.setrecursionlimit(10000)
from enumerate_aeq import isomorphic  # d=4 module's iso is generic

def is_minimal(adj, n):
    full = (1 << n) - 1
    for u in range(n):
        for v in range(u + 1, n):
            if (adj[u] >> v) & 1:
                if (adj[u] | adj[v] | (1 << u) | (1 << v)) == full:
                    return False
    return True

def invariant(adj, n):
    from itertools import combinations
    tri = []
    for v in range(n):
        nb = [i for i in range(n) if (adj[v] >> i) & 1]
        t = sum(1 for i, j in combinations(nb, 2) if (adj[i] >> j) & 1)
        nd = tuple(sorted(bin(adj[u]).count("1") for u in nb))
        tri.append((bin(adj[v]).count("1"), t, nd))
    return tuple(sorted(tri))

def main():
    n = 13
    # (A) minimal graphs from the in-house full enumeration
    A = []
    with open(f"aeq_d5_n{n}.bin", "rb") as f:
        hdr = struct.unpack("<II", f.read(8))
        assert hdr[0] == n
        for _ in range(hdr[1]):
            adj = list(struct.unpack(f"<{n}I", f.read(4 * n)))
            if is_minimal(adj, n):
                A.append(adj)
    # (B) triangleramsey route
    B = []
    with open(f"aeq_d5_n{n}_tr.txt") as f:
        for line in f:
            parts = [int(x) for x in line.split()]
            assert parts[0] == n
            B.append(parts[1:1 + n])
    print(f"A (in-house minimal): {len(A)}, B (triangleramsey): {len(B)}")
    assert len(A) == len(B) == 242
    # match up to isomorphism via invariant buckets
    from collections import defaultdict
    bucketsB = defaultdict(list)
    for g in B:
        bucketsB[invariant(g, n)].append(g)
    unmatched = 0
    usedB = set()
    for g in A:
        found = False
        for k, h in enumerate(bucketsB.get(invariant(g, n), [])):
            key = (invariant(g, n), k)
            if key in usedB:
                continue
            if isomorphic(g, h, n):
                usedB.add(key); found = True; break
        if not found:
            unmatched += 1
    print(f"unmatched: {unmatched}")
    assert unmatched == 0
    print("OK: the 242 minimal graphs on 13 vertices agree (up to iso) "
          "between the two independent generation routes.")

if __name__ == "__main__":
    main()
