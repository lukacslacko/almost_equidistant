"""Validation controls for the d=5 certified kernel (ckernel5.c).

KILL controls (must be certified non-realizable in R^5):
  C1. K_7 — seven pairwise-unit points cannot exist in R^5 (at most d+1=6).
  C2. K_1 + K_{2,2,2,2,2} (apex joined to every vertex of the 5-dimensional
      cross-polytope graph, 11 vertices): the cross-polytope's unit
      realization is unique up to congruence with all vertices at distance
      1/sqrt(2) from its centre, and a common unit neighbour of all ten
      vertices would have to be that centre — contradiction.

SURVIVE controls (realizable; engine must NOT kill, reports survivors):
  S1. K_6 — the unit 5-simplex itself (seed only, trivially realizable).
  S2. K_{2,2,2,2,2} — the unit cross-polytope in R^5 (exists; the engine
      must isolate its two mirror images, like K_{2,2,2,2} did for d=4).
  S3. The 16-vertex half-cube graph (complement of the Clebsch graph) —
      the Larman-Rogers almost-equidistant set realizes it.
"""
import math, time
from cdriver5 import decide5, gen_orders

TWO_PI = 2 * math.pi

def K(n):
    return [((1 << n) - 1) & ~(1 << i) for i in range(n)]

def cross_polytope(k=5):
    """K_{2,2,...,2} on 2k vertices; i and i+k are the non-adjacent pairs."""
    n = 2 * k
    adj = [0] * n
    for i in range(n):
        for j in range(n):
            if i != j and j != (i + k) % n and i != (j + k) % n:
                adj[i] |= 1 << j
    return adj

def apex_cross():
    """K_1 + cross-polytope: vertex 10 adjacent to all of 0..9."""
    base = cross_polytope(5)
    adj = [m | (1 << 10) for m in base] + [(1 << 10) - 1]
    return adj

def halfcube16():
    """Unit-distance graph of the Larman-Rogers set: vertices = odd-weight
    +-1 vectors in R^5, adjacent iff Hamming distance 2 (squared dist 1
    after scaling); complement is the Clebsch graph."""
    from itertools import product
    V = [v for v in product((1, -1), repeat=5)
         if sum(1 for x in v if x == 1) % 2 == 1]
    n = len(V)
    adj = [0] * n
    for i in range(n):
        for j in range(n):
            if i != j and sum(1 for a, b in zip(V[i], V[j]) if a != b) == 2:
                adj[i] |= 1 << j
    return adj

def run(name, adj, n, expect, **kw):
    t = time.time()
    st, nodes, unres = decide5(adj, n, **kw)
    ok = (st == expect)
    print(f"{name}: {st} (expect {expect}) nodes={nodes} unresolved={unres} "
          f"t={time.time()-t:.1f}s  {'PASS' if ok else 'FAIL'}", flush=True)
    return ok

def main():
    ok = True
    ok &= run("C1 K7            ", K(7), 7, "KILLED")
    ok &= run("S1 K6            ", K(6), 6, "SURVIVORS")
    ok &= run("S2 cross-polytope", cross_polytope(), 10, "SURVIVORS",
              theta_min=TWO_PI / (1 << 12), max_nodes=100_000_000)
    ok &= run("C2 apex+cross    ", apex_cross(), 11, "KILLED",
              max_nodes=200_000_000)
    hc = halfcube16()
    print("halfcube orders:", [(s, o) for s, o in gen_orders(hc, 16, kmax=2)][:1])
    ok &= run("S3 half-cube 16  ", hc, 16, "SURVIVORS",
              theta_min=TWO_PI / (1 << 10), max_nodes=400_000_000)
    print("controls5:", "PASS" if ok else "FAIL")
    return 0 if ok else 1

if __name__ == "__main__":
    raise SystemExit(main())
