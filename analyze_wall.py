"""Structural analysis of the hard tail of level 17.

Hypothesis: the resistant graphs admit DEGENERATE edge-exact
configurations obtained by identifying one non-adjacent vertex pair
(u,v), collapsing the 17-vertex graph onto a realizable 16-vertex graph
-- most naturally a subgraph of the half-cube unit-distance graph (the
extremal 16-point almost-equidistant set). Around such a collapsed
configuration the search sees a near-solution channel in which the
coincidence rules are the only kill mechanism (the d=4 campaign met the
same phenomenon; cf. "multiple points" in Gyorey's thesis).

For each remaining graph and each non-edge (u,v): form G/(u=v) (16
vertices, adjacency = union) and test whether it embeds into the
half-cube graph as a spanning subgraph (backtracking monomorphism,
16 -> 16 vertices, edges to edges).

Usage: python analyze_wall.py idx1 idx2 ...   (or no args: default list)
"""
import sys
from reproduce5 import load_level
from controls5 import halfcube16

HC = halfcube16()
N = 16

def embeds_into_halfcube(adj16):
    """is there an injective map V(G)->V(HC) sending edges to edges?"""
    degs = [bin(a).count("1") for a in adj16]
    order = sorted(range(N), key=lambda v: -degs[v])
    mapping = [-1] * N
    used = [False] * 16
    def bt(i):
        if i == N: return True
        v = order[i]
        for w in range(16):
            if used[w]: continue
            ok = True
            for j in range(i):
                u = order[j]
                if (adj16[v] >> u) & 1 and not ((HC[w] >> mapping[u]) & 1):
                    ok = False; break
            if ok:
                mapping[v] = w; used[w] = True
                if bt(i + 1): return True
                used[w] = False; mapping[v] = -1
        return False
    return bt(0)

def collapses(gi, graphs):
    adj = graphs[gi]
    n = 17
    out = []
    for u in range(n):
        for v in range(u + 1, n):
            if (adj[u] >> v) & 1: continue
            # identify v into u; relabel to 0..15
            keep = [x for x in range(n) if x != v]
            pos = {x: i for i, x in enumerate(keep)}
            m = [0] * N
            for x in keep:
                nb = adj[x]
                if x == u: nb |= adj[v]
                for y in keep:
                    if y == x: continue
                    if (nb >> y) & 1 or (x == u and (adj[v] >> y) & 1):
                        m[pos[x]] |= 1 << pos[y]
            # symmetrize (v's edges now point at u)
            for a in range(N):
                for b in range(N):
                    if (m[a] >> b) & 1: m[b] |= 1 << a
            if embeds_into_halfcube(m):
                out.append((u, v))
    return out

def main():
    graphs = load_level(17)
    idxs = [int(a) for a in sys.argv[1:]] or [11402, 11403, 11420, 7152, 5892, 9384]
    for gi in idxs:
        cc = collapses(gi, graphs)
        print(f"graph {gi}: {len(cc)} collapsing non-edges onto half-cube"
              f"{('  e.g. ' + str(cc[:4])) if cc else ''}", flush=True)

if __name__ == "__main__":
    main()
