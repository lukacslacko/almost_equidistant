"""Hot-map: for each remaining graph and EVERY (seed x circle-vertex)
parametrization, coarse-tile the circle (192 slices, small cap) under the
current kernel and count surviving slices. Output: per graph, the list of
(dec-id, hot-count, hot-positions); decs with zero hot slices certify the
graph outright (their tiling is fully KILLED).
Results land in hotmap_results.json; certified graphs are written into
results_n17.json immediately.
"""
import json, math, os, time
from multiprocessing import Pool, cpu_count
from cdriver5 import _cliques, _layered_order, decide5
from reproduce5 import load_level

TWO_PI = 2 * math.pi
NS = 192
CAP = 60_000

_G = None
_ORD = {}

def all_orders(adj, n=17):
    """one circle-early order per (seed, first-circle-vertex) pair"""
    seeds = _cliques(adj, n, 6) or _cliques(adj, n, 5)
    outs, seen = [], set()
    for seed in seeds:
        placed0 = set(seed)
        c4s = [v for v in range(n) if v not in placed0
               and sum(1 for u in placed0 if (adj[v] >> u) & 1) == 4]
        for cv in (c4s if c4s else [None]):
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

def _task(args):
    global _G
    gi, di, k = args
    if _G is None:
        _G = load_level(17)
    adj = _G[gi]
    if (gi) not in _ORD:
        _ORD[gi] = all_orders(adj)
    seed, order = _ORD[gi][di]
    lo, hi = k * TWO_PI / NS, (k + 1) * TWO_PI / NS
    st, nodes, unres = decide5(adj, 17, seed=seed, order=list(order),
                               th0=(lo, hi), max_nodes=CAP)
    return gi, di, k, st

def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=max(1, cpu_count() - 2))
    ap.add_argument("--maxdecs", type=int, default=24)
    args = ap.parse_args()
    graphs = load_level(17)
    results = {int(k): v for k, v in json.load(open("results_n17.json")).items()}
    todo = [gi for gi in range(len(graphs)) if results.get(gi) is None]
    print(f"hot-mapping {len(todo)} graphs", flush=True)
    plans = {}
    tasks = []
    for gi in todo:
        n_orders = min(len(all_orders(graphs[gi])), args.maxdecs)
        plans[gi] = n_orders
        for di in range(n_orders):
            for k in range(NS):
                tasks.append((gi, di, k))
    print(f"{len(tasks)} scan tasks", flush=True)
    hot = {}     # (gi, di) -> [hot ks]
    done_g = set()
    t0 = time.time()
    ndone = 0
    with Pool(args.workers) as pool:
        for gi, di, k, st in pool.imap_unordered(_task, tasks, chunksize=16):
            ndone += 1
            if ndone % 5000 == 0:
                print(f"  scan {ndone}/{len(tasks)} t={time.time()-t0:.0f}s",
                      flush=True)
            if st != "KILLED":
                hot.setdefault((gi, di), []).append(k)
    # summarize
    out = {}
    for gi in todo:
        per = []
        for di in range(plans[gi]):
            ks = sorted(hot.get((gi, di), []))
            per.append((di, len(ks), ks))
        per.sort(key=lambda x: x[1])
        out[gi] = per
        best = per[0]
        if best[1] == 0:
            results[gi] = f"hotmap-dec{best[0]}"
            done_g.add(gi)
            print(f"[{gi}] CERTIFIED outright by dec {best[0]} "
                  f"(zero hot slices at 192-tiling)", flush=True)
        else:
            print(f"[{gi}] best dec {best[0]}: {best[1]} hot slices at "
                  f"{[round(k*360/NS,1) for k in best[2][:6]]} deg", flush=True)
    json.dump({str(k): v for k, v in out.items()},
              open("hotmap_results.json", "w"))
    if done_g:
        json.dump(results, open("results_n17.json", "w"))
    print(f"done in {time.time()-t0:.0f}s; certified outright: {len(done_g)}; "
          f"remaining: {len(todo)-len(done_g)}", flush=True)

if __name__ == "__main__":
    main()
