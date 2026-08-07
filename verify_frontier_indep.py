"""Independent re-implementation of the extension/frontier pipeline
(ext6.c + escalate_frontier.py) for cross-checking the 2026-08-07 result
   644 residue -> 117,290 distinct level-20 candidates -> 34 true frontier
   -> 3,602 level-21 candidates -> 0 survive.

Independence: this file shares NO code with ext6.c.  Subgraph tests are
re-derived (different loop structures), isomorphism is decided by profile
bucketing + explicit backtracking matching (no canonical labeling at all),
and the all-maximal-completions enumeration uses its own branching.

Phases (run all by default):
  --control14   extensions of aeq_d6_n14.txt: expect 3,969 classes AND
                iso-bijection with aeq_d6_n15.txt
  --control15   extensions of aeq_d6_n15.txt: expect 18,917 classes
  --run         the real thing: residue644.txt -> level 20 -> level 21
"""
import sys
from itertools import combinations
from multiprocessing import Pool

# ---------- basic graph utilities (independent implementations) ----------

def parse(fn):
    out = []
    for line in open(fn):
        p = line.split()
        if not p:
            continue
        n = int(p[0])
        out.append((n, tuple(int(x) for x in p[1:1 + n])))
    return out

def complement(n, a):
    full = (1 << n) - 1
    return tuple(full & ~a[i] & ~(1 << i) for i in range(n))

def triangle_free(n, h):
    for u in range(n):
        m = h[u]
        while m:
            v = m & -m
            vi = v.bit_length() - 1
            m ^= v
            if vi > u and h[u] & h[vi] & ~((1 << (u + 1)) - 1) & ~0:
                pass
            if vi > u and (h[u] & h[vi]):
                return False
    return True

def maximal_tf(n, h):
    for u in range(n):
        for v in range(u + 1, n):
            if not (h[u] >> v) & 1 and not (h[u] & h[v]):
                return False
    return True

def clique_number_ge(n, a, k):
    """Does a (as masks) contain K_k?  Greedy pivoting backtracker."""
    def bk(cand, size):
        if size == k:
            return True
        while cand:
            v = (cand & -cand).bit_length() - 1
            cand &= cand - 1
            if bk(cand & a[v], size + 1):
                return True
        return False
    return bk((1 << n) - 1, 0)

def has_k1333(n, a):
    """K_{1,3,3,3} subgraph, apex-first enumeration (independent of the
    triple-pair-first loop in filter_mtf6.c)."""
    for apex in range(n):
        nb = a[apex]
        cnt = bin(nb).count("1")
        if cnt < 9:
            continue
        verts = [v for v in range(n) if (nb >> v) & 1]
        for t1 in combinations(verts, 3):
            m1 = (1 << t1[0]) | (1 << t1[1]) | (1 << t1[2])
            c1 = a[t1[0]] & a[t1[1]] & a[t1[2]] & nb & ~m1
            if bin(c1).count("1") < 6:
                continue
            v2 = [v for v in verts if (c1 >> v) & 1]
            for t2 in combinations(v2, 3):
                m2 = (1 << t2[0]) | (1 << t2[1]) | (1 << t2[2])
                c2 = a[t2[0]] & a[t2[1]] & a[t2[2]] & c1 & ~m2
                if bin(c2).count("1") < 3:
                    continue
                v3 = [v for v in v2 if (c2 >> v) & 1]
                for t3 in combinations(v3, 3):
                    return True
    return False

# --------------------- isomorphism via backtracking ---------------------

def vprofile(n, a):
    """Per-vertex invariant: (degree, triangles, sorted neighbour degrees)."""
    deg = [bin(a[v]).count("1") for v in range(n)]
    tri = []
    for v in range(n):
        t = 0
        m = a[v]
        while m:
            u = (m & -m).bit_length() - 1
            m &= m - 1
            t += bin(a[v] & a[u]).count("1")
        tri.append(t // 2)
    prof = []
    for v in range(n):
        nbd = tuple(sorted(deg[u] for u in range(n) if (a[v] >> u) & 1))
        prof.append((deg[v], tri[v], nbd))
    return prof

def gprofile(n, a):
    return (n, tuple(sorted(vprofile(n, a))))

def isomorphic(n, a, b):
    pa, pb = vprofile(n, a), vprofile(n, b)
    if sorted(pa) != sorted(pb):
        return False
    order = sorted(range(n), key=lambda v: (pa[v], -bin(a[v]).count("1")))
    cand0 = [[u for u in range(n) if pb[u] == pa[v]] for v in range(n)]
    mapping = [-1] * n
    used = [False] * n

    def bt(i):
        if i == n:
            return True
        v = order[i]
        for u in cand0[v]:
            if used[u]:
                continue
            ok = True
            for j in range(i):
                w = order[j]
                if ((a[v] >> w) & 1) != ((b[u] >> mapping[w]) & 1):
                    ok = False
                    break
            if ok:
                mapping[v] = u
                used[u] = True
                if bt(i + 1):
                    return True
                used[u] = False
                mapping[v] = -1
        return False
    return bt(0)

class IsoSet:
    """Set of graphs up to isomorphism (profile buckets + pairwise tests)."""
    def __init__(self):
        self.buckets = {}
        self.count = 0
    def add(self, n, a):
        """Returns True if new."""
        key = gprofile(n, a)
        lst = self.buckets.setdefault(key, [])
        for b in lst:
            if isomorphic(n, a, b):
                return False
        lst.append(a)
        self.count += 1
        return True
    def contains(self, n, a):
        key = gprofile(n, a)
        for b in self.buckets.get(key, []):
            if isomorphic(n, a, b):
                return True
        return False

# --------------------------- extension step ---------------------------

def extensions_of(args):
    """All valid one-vertex extensions of R (independent derivation used
    only for N-completeness: N ranges over ALL subsets of size 1..7; each
    built C is checked from scratch)."""
    n, R = args
    full = (1 << n) - 1
    H0 = complement(n, R)
    out = []
    for k in range(1, 8):
        for combo in combinations(range(n), k):
            N = 0
            for v in combo:
                N |= 1 << v
            ok = True
            for u in range(n):
                if not (N >> u) & 1 and not (H0[u] & N):
                    ok = False
                    break
            if not ok:
                continue
            C = list(R)
            for u in combo:
                add = H0[u] & N
                for v in combo:
                    if v != u and (add >> v) & 1:
                        C[u] |= 1 << v
            for i in range(n):
                if not (N >> i) & 1:
                    C[i] |= 1 << n
            C.append(full & ~N)
            nn = n + 1
            Ct = tuple(C)
            H = complement(nn, Ct)
            if not triangle_free(nn, H):
                raise AssertionError("TF anomaly")
            if not maximal_tf(nn, H):
                continue
            if clique_number_ge(nn, Ct, 8):
                continue
            if has_k1333(nn, Ct):
                continue
            out.append(Ct)
    return out

# ----------------------- completion enumeration -----------------------

def free_pairs(n, h):
    out = []
    for u in range(n):
        for v in range(u + 1, n):
            if not (h[u] >> v) & 1 and not (h[u] & h[v]):
                out.append((u, v))
    return out

def greedy_complete_lastfirst(n, h):
    """A maximal TF completion — deliberately completes the LAST free pair
    first (the C tool completes the first), so the two implementations
    probe different completions where the choice matters."""
    h = list(h)
    while True:
        fp = free_pairs(n, h)
        if not fp:
            return tuple(h)
        u, v = fp[-1]
        h[u] |= 1 << v
        h[v] |= 1 << u

def all_completions(n, h, cap=100000):
    """All maximal TF supergraphs; own branching: process free pairs in
    order; for the first one, either include it, or exclude it by adding a
    blocking cherry through each possible centre."""
    results = []
    seen = set()
    stack = [tuple(h)]
    nodes = 0
    while stack:
        cur = stack.pop()
        nodes += 1
        if nodes > cap:
            raise RuntimeError("completion cap")
        fp = free_pairs(n, cur)
        if not fp:
            if cur not in seen:
                seen.add(cur)
                results.append(cur)
            continue
        u, v = fp[0]
        h2 = list(cur)
        h2[u] |= 1 << v
        h2[v] |= 1 << u
        stack.append(tuple(h2))
        for z in range(n):
            if z == u or z == v:
                continue
            h3 = list(cur)
            good = True
            for (x, y) in ((u, z), (v, z)):
                if (h3[x] >> y) & 1:
                    continue
                if h3[x] & h3[y]:
                    good = False
                    break
                h3[x] |= 1 << y
                h3[y] |= 1 << x
            if good:
                stack.append(tuple(h3))
    return results

# ------------------------------ pipeline ------------------------------

def dedup(graphs, n):
    s = IsoSet()
    for a in graphs:
        s.add(n, a)
    return s

def frontier_step(parents, np, workers):
    """parents: IsoSet reps list at size np. Returns (distinct IsoSet,
    survivors list) at size np+1 after greedy + full-escalation tests."""
    reps = []
    for lst in parents.buckets.values():
        reps.extend(lst)
    with Pool(workers) as pool:
        chunks = pool.map(extensions_of, [(np, R) for R in reps])
    allext = [c for ch in chunks for c in ch]
    nn = np + 1
    dis = dedup(allext, nn)
    print(f"  level {nn}: {len(allext)} with multiplicity, "
          f"{dis.count} distinct", flush=True)
    survivors = []
    nkill = 0
    for lst in dis.buckets.values():
        for C in lst:
            dead = False
            for u in range(nn):
                keep = [i for i in range(nn) if i != u]
                Cd = []
                for i in keep:
                    m = 0
                    for jp, j in enumerate(keep):
                        if (C[i] >> j) & 1:
                            m |= 1 << jp
                    Cd.append(m)
                h = complement(nn - 1, tuple(Cd))
                hm = greedy_complete_lastfirst(nn - 1, h)
                c3 = complement(nn - 1, hm)
                if not parents.contains(nn - 1, c3):
                    dead = True
                    break
            if dead:
                nkill += 1
            else:
                survivors.append(C)
    print(f"  greedy(last-first): {nkill} killed, {len(survivors)} survive",
          flush=True)
    true_frontier = []
    for C in survivors:
        dead = False
        for u in range(nn):
            keep = [i for i in range(nn) if i != u]
            Cd = []
            for i in keep:
                m = 0
                for jp, j in enumerate(keep):
                    if (C[i] >> j) & 1:
                        m |= 1 << jp
                Cd.append(m)
            h = complement(nn - 1, tuple(Cd))
            for hm in all_completions(nn - 1, h):
                c3 = complement(nn - 1, hm)
                if not parents.contains(nn - 1, c3):
                    dead = True
                    break
            if dead:
                break
        if not dead:
            true_frontier.append(C)
    print(f"  full escalation: {len(true_frontier)} true frontier",
          flush=True)
    return dis, true_frontier

def control14(workers):
    corpus14 = parse("aeq_d6_n14.txt")
    corpus15 = parse("aeq_d6_n15.txt")
    with Pool(workers) as pool:
        chunks = pool.map(extensions_of, corpus14)
    allext = [c for ch in chunks for c in ch]
    dis = dedup(allext, 15)
    print(f"control14: {len(allext)} ext, {dis.count} distinct "
          f"(expect 3969)", flush=True)
    missing = 0
    for n, a in corpus15:
        if not dis.contains(15, a):
            missing += 1
    print(f"control14: corpus15 members missing from extension set: "
          f"{missing} (expect 0); counts equal: {dis.count == len(corpus15)}",
          flush=True)
    return dis.count == 3969 and missing == 0 and len(corpus15) == 3969

def control15(workers):
    corpus15 = parse("aeq_d6_n15.txt")
    with Pool(workers) as pool:
        chunks = pool.map(extensions_of, corpus15)
    allext = [c for ch in chunks for c in ch]
    dis = dedup(allext, 16)
    print(f"control15: {len(allext)} ext, {dis.count} distinct "
          f"(expect 18917)", flush=True)
    return dis.count == 18917

def run(workers):
    res = parse("residue644.txt")
    assert len(res) == 644
    parents = dedup([a for n, a in res], 19)
    print(f"run: residue {parents.count} distinct classes", flush=True)
    dis20, tf20 = frontier_step(parents, 19, workers)
    p20 = IsoSet()
    for C in tf20:
        p20.add(20, C)
    dis21, tf21 = frontier_step(p20, 20, workers)
    print(f"RESULT: level20 distinct {dis20.count} truefrontier {len(tf20)}; "
          f"level21 distinct {dis21.count} truefrontier {len(tf21)}",
          flush=True)

if __name__ == "__main__":
    args = sys.argv[1:]
    workers = 24
    todo = args if args else ["--control14", "--control15", "--run"]
    if "--control14" in todo:
        ok = control14(workers)
        print("control14:", "PASS" if ok else "FAIL", flush=True)
    if "--control15" in todo:
        ok = control15(workers)
        print("control15:", "PASS" if ok else "FAIL", flush=True)
    if "--run" in todo:
        run(workers)
