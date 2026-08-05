"""Scan decompositions across DISTINCT SEED CLIQUES of a stuck graph,
looking for a parametrization with no hot slices (the d=4 'alternate
decomposition' cure). For each 6-clique seed: build one circle-early
layered order, 96-slice coarse scan (small cap); report hot-slice counts.
Usage: python scan_seeds.py <gi> [max_seeds] [workers]"""
import sys, time, math
from multiprocessing import Pool
from cdriver5 import _cliques, _layered_order, decide5
from reproduce5 import load_level

TWO_PI = 2 * math.pi

def scan_one(args):
    gi, seed, order, k = args
    graphs = load_level(17)
    adj = graphs[gi]
    lo, hi = k * TWO_PI / 96, (k + 1) * TWO_PI / 96
    st, nodes, unres = decide5(adj, 17, seed=seed, order=list(order),
                               th0=(lo, hi), max_nodes=250_000)
    return seed, k, st

def main():
    gi = int(sys.argv[1])
    max_seeds = int(sys.argv[2]) if len(sys.argv) > 2 else 48
    workers = int(sys.argv[3]) if len(sys.argv) > 3 else 12
    graphs = load_level(17)
    adj = graphs[gi]
    seeds = _cliques(adj, 17, 6) or _cliques(adj, 17, 5)
    print(f"graph {gi}: {len(seeds)} seed cliques, scanning up to {max_seeds}",
          flush=True)
    # orders per (seed, first-circle-vertex) pair: every choice of the
    # vertex that rides the first circle re-parametrizes the geometry
    decs = []
    for seed in seeds[:max_seeds]:
        placed0 = set(seed)
        c4s = [v for v in range(17) if v not in placed0
               and sum(1 for u in placed0 if (adj[v] >> u) & 1) == 4]
        choices = c4s if c4s else [None]
        for cv in choices:
            if cv is None:
                o = _layered_order(adj, 17, seed, lambda c, p: max(c))
            else:
                o = _layered_order(adj, 17, seed,
                                   lambda c, p, cv=cv: cv if cv in c else max(c),
                                   force_circle_first=True)
            if o is not None and (seed, tuple(o)) not in decs:
                decs.append((seed, tuple(o)))
    seen = set(); uniq = []
    for d in decs:
        if d not in seen:
            seen.add(d); uniq.append(d)
    decs = uniq
    print(f"{len(decs)} orders", flush=True)
    tasks = [(gi, seed, order, k) for seed, order in decs for k in range(96)]
    hot = {}
    t0 = time.time()
    with Pool(workers) as pool:
        for seed, k, st in pool.imap_unordered(scan_one, tasks, chunksize=8):
            if st != "KILLED":
                hot.setdefault(seed, []).append(k)
    print(f"scan done t={time.time()-t0:.0f}s", flush=True)
    clean = [seed for seed, _ in decs if seed not in hot]
    for seed, ks in sorted(hot.items(), key=lambda x: len(x[1])):
        print(f"  seed {seed}: hot {sorted(ks)}", flush=True)
    print(f"CLEAN seeds ({len(clean)}): {clean[:8]}", flush=True)

if __name__ == "__main__":
    main()
