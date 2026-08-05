"""Enumerate all abstract almost-equidistant graphs in R^4 (BPSSV Section 3):
graphs G with:
  (1) no independent set of size 3  (complement triangle-free)
  (2) no K6
  (3) no K_{1,3,3} (complete tripartite, parts 1,3,3) as subgraph
Expected counts (BPSSV Table 3, d=4): n=4..17:
  7, 14, 37, 97, 316, 934, 2362, 2814, 944, 59, 4, 1, 1, 0
Minimal ones (Table 2): n=13: 12, n=14: 3, n=15: 1, n=16: 1.
Graphs stored as tuples of adjacency bitmasks.
"""
import sys, json, time
from itertools import combinations

def popcount(x): return bin(x).count("1")

def bits(x):
    out = []
    i = 0
    while x:
        if x & 1: out.append(i)
        x >>= 1; i += 1
    return out

# ---------- full (non-incremental) checks, used for validation ----------
def has_indep3(adj, n):
    for a in range(n):
        for b in range(a + 1, n):
            if not (adj[a] >> b) & 1:
                # need c nonadjacent to both
                mask = ~(adj[a] | adj[b]) & ((1 << n) - 1) & ~(1 << a) & ~(1 << b)
                if mask: return True
    return False

def cliques_of_size(adj, n, k):
    """all k-cliques as bitmasks"""
    res = []
    def ext(clq, cand, size):
        if size == k:
            res.append(clq); return
        c = cand
        while c:
            v = (c & -c).bit_length() - 1
            c &= c - 1
            ext(clq | (1 << v), cand & adj[v] & ~((1 << (v + 1)) - 1), size + 1)
    ext(0, (1 << n) - 1, 0)
    return res

def has_k6(adj, n):
    return len(cliques_of_size(adj, n, 6)) > 0

def has_k133(adj, n):
    full = (1 << n) - 1
    for a in range(n):
        na = adj[a]
        nb = bits(na)
        for X in combinations(nb, 3):
            common = na & adj[X[0]] & adj[X[1]] & adj[X[2]]
            common &= ~((1 << X[0]) | (1 << X[1]) | (1 << X[2]))
            if popcount(common) >= 3:
                return True
    return False

def is_abstract(adj, n):
    return not has_indep3(adj, n) and not has_k6(adj, n) and not has_k133(adj, n)

# ---------- canonical-ish invariant + exact isomorphism ----------
def invariant(adj, n):
    tri = []
    for v in range(n):
        nb = bits(adj[v])
        t = sum(1 for i, j in combinations(nb, 2) if (adj[i] >> j) & 1)
        nd = tuple(sorted(popcount(adj[u]) for u in nb))
        tri.append((popcount(adj[v]), t, nd))
    return tuple(sorted(tri))

def isomorphic(adj1, adj2, n):
    d1 = [popcount(a) for a in adj1]; d2 = [popcount(a) for a in adj2]
    if sorted(d1) != sorted(d2): return False
    # backtracking: map vertices of g1 (ordered by rarity of degree) to g2
    order = sorted(range(n), key=lambda v: (d1[v],))
    mapping = [-1] * n; used = [False] * n
    def bt(idx):
        if idx == n: return True
        v = order[idx]
        for w in range(n):
            if used[w] or d2[w] != d1[v]: continue
            ok = True
            for u in order[:idx]:
                if ((adj1[v] >> u) & 1) != ((adj2[w] >> mapping[u]) & 1):
                    ok = False; break
            if ok:
                mapping[v] = w; used[w] = True
                if bt(idx + 1): return True
                mapping[v] = -1; used[w] = False
        return False
    return bt(0)

# ---------- incremental extension ----------
def extend_all(level_graphs, n):
    """level_graphs: list of adj tuples on n vertices -> all abstract graphs on n+1"""
    m = n + 1
    seen = {}  # invariant -> list of adj tuples
    total_cand = 0
    for adj in level_graphs:
        degs = [popcount(a) for a in adj]
        cl5 = cliques_of_size(adj, n, 5) if m >= 6 else []
        full_old = (1 << n) - 1
        mindeg_bound = min(degs) + 1 if n else n
        lo = max(0, m - 6)
        hi = min(n, mindeg_bound)
        for k in range(lo, hi + 1):
            for Sv in combinations(range(n), k):
                S = 0
                for s in Sv: S |= (1 << s)
                # min-degree of new vertex
                okd = True
                for u in range(n):
                    du = degs[u] + (1 if (S >> u) & 1 else 0)
                    if du < k: okd = False; break
                if not okd: continue
                # alpha <= 2: non-neighbours of v pairwise adjacent
                nonN = full_old & ~S
                okA = True
                nn = bits(nonN)
                for i in range(len(nn)):
                    for j in range(i + 1, len(nn)):
                        if not (adj[nn[i]] >> nn[j]) & 1: okA = False; break
                    if not okA: break
                if not okA: continue
                # K6: 5-clique of parent inside S
                if any((c & S) == c for c in cl5): continue
                # K_{1,3,3} involving v:
                bad = False
                Sb = list(Sv)
                # (a) v as singleton: X triple in S, |N(x1)&N(x2)&N(x3)&S \ X| >= 3
                for X in combinations(Sb, 3):
                    common = S & adj[X[0]] & adj[X[1]] & adj[X[2]]
                    common &= ~((1 << X[0]) | (1 << X[1]) | (1 << X[2]))
                    if popcount(common) >= 3: bad = True; break
                # (b) v in a triple part: a in S; x2,x3 in N(a); Y subset of
                #     N(a)&S&N(x2)&N(x3) minus {x2,x3}, size >= 3
                if not bad:
                    for a in Sb:
                        na = adj[a]
                        nbrs = bits(na)
                        for x2, x3 in combinations(nbrs, 2):
                            comm = na & S & adj[x2] & adj[x3]
                            comm &= ~((1 << x2) | (1 << x3))
                            if popcount(comm) >= 3: bad = True; break
                        if bad: break
                if bad: continue
                # build extended graph
                newadj = list(adj) + [S]
                for u in Sv: newadj[u] |= (1 << n)
                cand = tuple(newadj)
                total_cand += 1
                inv = invariant(cand, m)
                bucket = seen.setdefault(inv, [])
                if not any(isomorphic(cand, g, m) for g in bucket):
                    bucket.append(cand)
    out = [g for b in seen.values() for g in b]
    return out, total_cand

def is_minimal(adj, n):
    """minimal: removing any edge creates an independent 3-set"""
    full = (1 << n) - 1
    for u in range(n):
        for v in range(u + 1, n):
            if (adj[u] >> v) & 1:
                # removable iff N(u)|N(v)|{u,v} covers V
                if (adj[u] | adj[v] | (1 << u) | (1 << v)) == full:
                    return False   # some edge removable -> not minimal
    return True

EXPECTED = {4: 7, 5: 14, 6: 37, 7: 97, 8: 316, 9: 934, 10: 2362, 11: 2814,
            12: 944, 13: 59, 14: 4, 15: 1, 16: 1, 17: 0}

def main():
    t0 = time.time()
    level = [(0,)]  # single vertex, no edges (mask 0)
    n = 1
    results = {}
    while n < 17 and level:
        level, ncand = extend_all(level, n)
        n += 1
        # validation with full checks
        for g in level:
            assert is_abstract(g, n), f"incremental/full mismatch at n={n}"
        nmin = sum(1 for g in level if is_minimal(g, n))
        exp = EXPECTED.get(n, "?")
        print(f"n={n:2d}: {len(level):6d} abstract graphs (expected {exp}), "
              f"{nmin} minimal, candidates {ncand}, t={time.time()-t0:.1f}s", flush=True)
        results[n] = [[g[i] for i in range(n)] for g in level]
        if n >= 12:
            with open(f"aeq_d4_n{n}.json", "w") as f:
                json.dump({"n": n, "graphs": results[n],
                           "minimal": [i for i, g in enumerate(level) if is_minimal(g, n)]}, f)
    print("done", time.time() - t0)

if __name__ == "__main__":
    main()
