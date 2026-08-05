"""Profile the 38 pending level-17 graphs: extract them to a data file
and build a combinatorial + numerical summary table."""
import json, glob
from itertools import combinations
from reproduce5 import load_level

PENDING = [3013, 4190, 4816, 5604, 5892, 7152, 7904, 9384, 9518, 9621,
           10429, 10440, 11377, 11378, 11381, 11384, 11385, 11395, 11396,
           11402, 11403, 11405, 11586, 11602, 11605, 11608, 11677, 11770,
           12313, 12316, 12372, 12519, 12520, 12560, 12572, 12590, 12594,
           12595]

def cliques(adj, n, k):
    res = []
    def ext(clq, cand, start):
        if len(clq) == k: res.append(tuple(clq)); return
        for v in range(start, n):
            if (cand >> v) & 1: ext(clq + [v], cand & adj[v], v + 1)
    ext([], (1 << n) - 1, 0)
    return res

def profile(adj, n=17):
    degs = sorted(bin(a).count("1") for a in adj)
    e = sum(degs) // 2
    k6 = cliques(adj, n, 6)
    k5 = cliques(adj, n, 5)
    # complement (the maximal triangle-free graph)
    full = (1 << n) - 1
    comp = [full & ~adj[i] & ~(1 << i) for i in range(n)]
    cdegs = sorted(bin(c).count("1") for c in comp)
    # non-edges with large common neighbourhoods (the coincidence-prone pairs)
    bigcn = 0
    maxcn = 0
    for u in range(n):
        for v in range(u + 1, n):
            if (adj[u] >> v) & 1: continue
            cn = bin(adj[u] & adj[v]).count("1")
            maxcn = max(maxcn, cn)
            if cn >= 10: bigcn += 1
    return dict(edges=e, degs=degs, nK6=len(k6), nK5=len(k5),
                cdegs=cdegs, bigcn=bigcn, maxcn=maxcn)

def main():
    graphs = load_level(17)
    # sweep residuals
    sweep = {}
    for f in glob.glob("sweep/n17*.out"):
        for line in open(f):
            p = line.split()
            sweep[int(p[0])] = (float(p[2]), float(p[3]), float(p[4]))
    with open("aeq_d5_n17_pending38.txt", "w") as out:
        for gi in PENDING:
            out.write(f"# graph {gi}\n17 " +
                      " ".join(map(str, graphs[gi])) + "\n")
    rows = []
    for gi in PENDING:
        p = profile(graphs[gi])
        bf, bfd, mind = sweep.get(gi, (None, None, None))
        rows.append((gi, p, bf, bfd))
        dmin, dmax = p["degs"][0], p["degs"][-1]
        print(f"{gi:6d} e={p['edges']:3d} deg[{dmin}-{dmax}] "
              f"K6={p['nK6']:2d} K5={p['nK5']:3d} "
              f"cn>=10:{p['bigcn']:2d} maxcn={p['maxcn']:2d} "
              f"bestF={bf:.3f} bestFdistinct={bfd:.3f}")
    # family clustering by complement degree sequence
    from collections import Counter
    fam = Counter(tuple(profile(graphs[gi])["cdegs"]) for gi in PENDING)
    print("\ncomplement degree-sequence families:")
    for k, v in fam.most_common():
        print(f"  {v} graphs with complement degs {list(k)}")

if __name__ == "__main__":
    main()
