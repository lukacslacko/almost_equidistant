"""Cheap decomposition-diversity probe for the 34 frontier graphs: for each
parametrization, run a few evenly spaced probe slices at a small node cap
and report which decompositions look globally cheap (all probes KILLED
with small node counts).  Output: per-graph ranking to stdout (append to
probe log), machine-readable summary to frontier34_probe.json.
Usage: python probe_frontier_decs.py <gidx|all> [workers] [cap] [nprobe]
"""
import json, math, os, sys, time
from multiprocessing import Pool
from certify_frontier34 import load_pat, all_orders
from cdriver6 import decide6

TWO_PI = 2 * math.pi
NS = 96

_adj = None
_ords = None
def _init(pi):
    global _adj, _ords
    _adj = load_pat(pi)
    _ords = all_orders(_adj)

def _probe(args):
    di, k, cap = args
    seed, order = _ords[di]
    lo, hi = k * TWO_PI / NS, (k + 1) * TWO_PI / NS
    t0 = time.time()
    st, nodes, unres = decide6(_adj, 20, seed=seed, order=list(order),
                               th0=(lo, hi), max_nodes=cap)
    return di, k, st, nodes, time.time() - t0

def probe(pi, workers, cap, nprobe):
    adj = load_pat(pi)
    ords = all_orders(adj)
    nd = len(ords)
    ks = [round(i * NS / nprobe) % NS for i in range(nprobe)]
    tasks = [(di, k, cap) for di in range(nd) for k in ks]
    res = {}
    t0 = time.time()
    with Pool(workers, initializer=_init, initargs=(pi,)) as pool:
        for di, k, st, nodes, dt in pool.imap_unordered(_probe, tasks,
                                                        chunksize=1):
            res.setdefault(di, []).append((k, st, nodes))
    rank = []
    for di in range(nd):
        r = res.get(di, [])
        nab = sum(1 for _, st, _ in r if st != "KILLED")
        tn = sum(n for _, _, n in r)
        rank.append((nab, tn, di))
    rank.sort()
    print(f"[{pi}] {nd} decs probed ({nprobe} slices @ {cap} cap) "
          f"t={time.time()-t0:.0f}s", flush=True)
    for nab, tn, di in rank[:6]:
        print(f"   dec {di}: {nab}/{nprobe} non-killed, total nodes {tn}",
              flush=True)
    return {"rank": [[nab, tn, di] for nab, tn, di in rank]}

def main():
    arg = sys.argv[1]
    idxs = list(range(34)) if arg == "all" else [int(x) for x in arg.split(",")]
    workers = int(sys.argv[2]) if len(sys.argv) > 2 else 26
    cap = int(sys.argv[3]) if len(sys.argv) > 3 else 30_000
    nprobe = int(sys.argv[4]) if len(sys.argv) > 4 else 6
    fn = "frontier34_probe.json"
    out = json.load(open(fn)) if os.path.exists(fn) else {}
    for pi in idxs:
        out[str(pi)] = probe(pi, workers, cap, nprobe)
        json.dump(out, open(fn, "w"))

if __name__ == "__main__":
    main()
