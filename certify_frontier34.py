"""Certify the 34 level-20 frontier graphs (frontier20_true.txt)
non-realizable in R^6.  Same hotmap strategy as certify_patterns.py: for
each graph try every (seed x circle-vertex) parametrization with a
96-slice tiling; certified when some parametrization's full tiling is
KILLED.  Results appended to frontier34_cert.json.
Usage: python certify_frontier34.py <idx1,idx2,...|all> [workers]
"""
import json, math, os, sys, time
from multiprocessing import Pool
from cdriver6 import _cliques, _layered_order, decide6

TWO_PI = 2 * math.pi
NP = 20
NS = 96
CAP = 3_000_000
GRAPHF = "frontier20_true.txt"
OUTF = "frontier34_cert.json"

def load_pat(i):
    with open(GRAPHF) as f:
        for k, line in enumerate(f):
            if k == i:
                p = line.split()
                assert p[0] == str(NP)
                return [int(x) for x in p[1:1 + NP]]
    raise IndexError(i)

def all_orders(adj, n=NP):
    seeds = _cliques(adj, n, 7) or _cliques(adj, n, 6)
    outs, seen = [], set()
    for seed in seeds:
        placed0 = set(seed)
        c5s = [v for v in range(n) if v not in placed0
               and sum(1 for u in placed0 if (adj[v] >> u) & 1) == 5]
        for cv in (c5s if c5s else [None]):
            if cv is None:
                o = _layered_order(adj, n, seed, lambda c, p: max(c))
            else:
                o = _layered_order(adj, n, seed,
                                   lambda c, p, cv=cv: cv if cv in c else max(c),
                                   force_circle_first=True)
            if o is not None and (seed, tuple(o)) not in seen:
                seen.add((seed, tuple(o)))
                outs.append((seed, tuple(o)))
    return outs

_pat = None
_ords = None
def _init(pi):
    global _pat, _ords
    _pat = load_pat(pi)
    _ords = all_orders(_pat)

def _slice(args):
    di, k = args
    seed, order = _ords[di]
    lo, hi = k * TWO_PI / NS, (k + 1) * TWO_PI / NS
    st, nodes, unres = decide6(_pat, NP, seed=seed, order=list(order),
                               th0=(lo, hi), max_nodes=CAP)
    return di, k, st

def certify(pi, workers):
    adj = load_pat(pi)
    ords = all_orders(adj)
    print(f"[{pi}] {len(ords)} parametrizations", flush=True)
    if not ords:
        print(f"[{pi}] NO parametrization (needs 2-sphere stage)", flush=True)
        return None
    tasks = [(di, k) for di in range(len(ords)) for k in range(NS)]
    hot = {}
    t0 = time.time()
    with Pool(workers, initializer=_init, initargs=(pi,)) as pool:
        for di, k, st in pool.imap_unordered(_slice, tasks, chunksize=4):
            if st != "KILLED":
                hot.setdefault(di, []).append(k)
    clean = [di for di in range(len(ords)) if di not in hot]
    if clean:
        print(f"[{pi}] CERTIFIED non-realizable (dec {clean[0]}, "
              f"{NS}-slice tiling all killed) t={time.time()-t0:.0f}s",
              flush=True)
        return clean[0]
    best = min(range(len(ords)), key=lambda d: len(hot.get(d, [])))
    print(f"[{pi}] NOT certified: best dec {best} has "
          f"{len(hot.get(best, []))} hot slices t={time.time()-t0:.0f}s",
          flush=True)
    return None

def main():
    if sys.argv[1] == "all":
        idxs = list(range(34))
    else:
        idxs = [int(x) for x in sys.argv[1].split(",")]
    workers = int(sys.argv[2]) if len(sys.argv) > 2 else 24
    obs = json.load(open(OUTF)) if os.path.exists(OUTF) else {}
    for pi in idxs:
        if str(pi) in obs:
            print(f"[{pi}] already certified", flush=True)
            continue
        dec = certify(pi, workers)
        if dec is not None:
            obs[str(pi)] = {"dec": dec, "tiling": NS}
            json.dump(obs, open(OUTF, "w"))
    print(f"certified {len(obs)} / 34", flush=True)

if __name__ == "__main__":
    main()
