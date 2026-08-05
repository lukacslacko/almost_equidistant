"""Independent brute-force verification of the K_{1,3,3,3} discrepancy at
n=11 (d=6 class): reads the complete pool of 105071 complement-triangle-
free graphs (tf_n11.bin), and for each graph decides, from the definition
with itertools (no bitmask tricks):
  - has K_8 (subgraph)?
  - has K_{1,3,3,3} as a subgraph: an apex a and three disjoint triples
    T1,T2,T3 (not containing a) with all 27 cross edges present and a
    adjacent to all nine vertices (edges INSIDE triples unrestricted)?
Prints the class count and cross-checks against the fast C checkers.
"""
import struct
from itertools import combinations

def brute_k8(adj, n):
    # clique of size 8 via simple recursion on candidate sets
    def ext(clq_last, cand, size):
        if size == 8: return True
        for v in range(clq_last + 1, n):
            if (cand >> v) & 1:
                if ext(v, cand & adj[v], size + 1): return True
        return False
    return ext(-1, (1 << n) - 1, 0)

def brute_k1333(adj, n):
    verts = range(n)
    for a in verts:
        Na = [v for v in verts if (adj[a] >> v) & 1]
        if len(Na) < 9: continue
        for nine in combinations(Na, 9):
            nine_set = set(nine)
            # partition nine into 3 triples with all cross edges
            rest = list(nine)
            first = rest[0]
            for t1o in combinations(rest[1:], 2):
                T1 = (first,) + t1o
                rem = [v for v in rest if v not in T1]
                f2 = rem[0]
                for t2o in combinations(rem[1:], 2):
                    T2 = (f2,) + t2o
                    T3 = tuple(v for v in rem if v not in T2)
                    ok = True
                    for x in T1:
                        for y in T2 + T3:
                            if not (adj[x] >> y) & 1: ok = False; break
                        if not ok: break
                    if ok:
                        for x in T2:
                            for y in T3:
                                if not (adj[x] >> y) & 1: ok = False; break
                            if not ok: break
                    if ok:
                        return True
    return False

def main():
    with open("tf_n11.bin", "rb") as f:
        n, cnt = struct.unpack("<II", f.read(8))
        assert n == 11
        kept = k8c = k13c = 0
        for i in range(cnt):
            adj = list(struct.unpack(f"<{n}I", f.read(4 * n)))
            if brute_k8(adj, n):
                k8c += 1
            elif brute_k1333(adj, n):
                k13c += 1
            else:
                kept += 1
            if (i + 1) % 20000 == 0:
                print(f"  {i+1}/{cnt} ... kept={kept} K8={k8c} K1333={k13c}",
                      flush=True)
    print(f"BRUTE FORCE: total={cnt} kept={kept} K8={k8c} K1333={k13c}")
    print("fast C checkers said:      kept=103194 K8=1738 K1333=139")
    print("BPSSV Table 3 (d=6, n=11): 103333")

if __name__ == "__main__":
    main()
