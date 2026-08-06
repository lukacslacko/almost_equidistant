"""Obstruction-coverage pilot: how many deferred level-19 d=6 candidates
contain one of the given 14-vertex patterns as a subgraph (edges into
edges, not induced)? Bounded-backtracking test; timeouts count as
not-covered (conservative).
Usage: python coverage6.py <patterns.txt> <pattern_indices...|-> [sample]
"""
import random, sys, time
from multiprocessing import Pool

def load(fn, n):
    out = []
    for line in open(fn):
        p = line.split()
        if p and p[0] == str(n):
            out.append([int(x) for x in p[1:1 + n]])
    return out

NP = 14
NT = 19
BT_LIMIT = 400_000

def embeds(pat, tgt, deg_p, deg_t):
    order = sorted(range(NP), key=lambda v: -deg_p[v])
    mapping = [-1] * NP
    used = [False] * NT
    nodes = 0
    def bt(i):
        nonlocal nodes
        nodes += 1
        if nodes > BT_LIMIT:
            raise TimeoutError
        if i == NP:
            return True
        v = order[i]
        for w in range(NT):
            if used[w] or deg_t[w] < deg_p[v]:
                continue
            ok = True
            for j in range(i):
                u = order[j]
                if (pat[v] >> u) & 1 and not ((tgt[w] >> mapping[u]) & 1):
                    ok = False; break
            if ok:
                mapping[v] = w; used[w] = True
                if bt(i + 1):
                    return True
                used[w] = False; mapping[v] = -1
        return False
    try:
        return bt(0)
    except TimeoutError:
        return False

_pats = None
def _cover_task(args):
    tgt = args
    deg_t = [bin(a).count("1") for a in tgt]
    for pi, (pat, deg_p) in enumerate(_pats):
        if embeds(pat, tgt, deg_p, deg_t):
            return pi
    return -1

def _init(pats):
    global _pats
    _pats = pats

def main():
    patf = sys.argv[1]
    idxs = sys.argv[2]
    nsample = int(sys.argv[3]) if len(sys.argv) > 3 else 1000
    pats_all = load(patf, NP)
    if idxs == "-":
        sel = list(range(len(pats_all)))
    else:
        sel = [int(x) for x in idxs.split(",")]
    pats = [(pats_all[i], [bin(a).count("1") for a in pats_all[i]])
            for i in sel]
    print(f"{len(pats)} patterns", flush=True)
    # sample deferred targets
    done = set()
    for line in open("killed_d6_n19.log"):
        p = line.split()
        if p: done.add(int(p[0]))
    random.seed(1)
    sample = []
    with open("aeq_d6_n19.txt") as f:
        for idx, line in enumerate(f):
            if idx in done: continue
            p = line.split()
            g = [int(x) for x in p[1:1 + NT]]
            if len(sample) < nsample:
                sample.append(g)
            else:
                j = random.randrange(idx + 1)
                if j < nsample: sample[j] = g
    print(f"{len(sample)} sampled targets", flush=True)
    t0 = time.time()
    hits = 0
    hist = {}
    with Pool(20, initializer=_init, initargs=(pats,)) as pool:
        for r in pool.imap_unordered(_cover_task, sample, chunksize=8):
            if r >= 0:
                hits += 1
                hist[r] = hist.get(r, 0) + 1
    print(f"coverage: {hits}/{len(sample)} = {100*hits/len(sample):.1f}% "
          f"in {time.time()-t0:.0f}s", flush=True)
    top = sorted(hist.items(), key=lambda x: -x[1])[:10]
    print("top patterns:", [(sel[k], v) for k, v in top], flush=True)

if __name__ == "__main__":
    main()
