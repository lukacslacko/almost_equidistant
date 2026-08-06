#!/usr/bin/env python3
"""Certified kill campaign for f(6), level 19: decide non-realizability
in R^6 of all 3,971,787 minimal abstract almost-equidistant graphs on 19
vertices (aeq_d6_n19.txt). Killing every one proves f(6) <= 18, hence
with the known 18-point construction f(6) = 18.

Scale-adapted design (vs reproduce5.py): graphs stream from the file
into the worker pool (no per-worker level loads); results are an
append-only log (killed_d6_n19.log: "<idx> <dec>" per line; deferred
graphs go to phase B racing at the end). Resumable.

Usage:
  python reproduce6.py --bulk [--workers N] [--cap M] [--limit K]
"""
import argparse, json, math, os, time
from multiprocessing import Pool, cpu_count

from cdriver6 import decide6, gen_orders, ensure_kernel

TWO_PI = 2 * math.pi
LEVEL = 19
HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)
GRAPHF = "aeq_d6_n19.txt"
KILLLOG = "killed_d6_n19.log"

def ncirc_of(adj, n, seed, order):
    placed = set(seed); c = 0
    for v in order:
        if sum(1 for u in placed if (adj[v] >> u) & 1) == 5: c += 1
        placed.add(v)
    return c

SKIPCIRCLE = os.environ.get("D6_SKIPCIRCLE", "0") == "1"

def _bulk_task(args):
    idx, adj, cap = args
    decs = gen_orders(adj, LEVEL, kmax=1)
    if not decs:
        return idx, "NOORDER", 0
    seed, order = decs[0]
    if SKIPCIRCLE and ncirc_of(adj, LEVEL, seed, order) != 0:
        return idx, "HASCIRCLE", 0
    # a KILLED unrestricted search is a complete certificate (the kernel
    # subdivides circle stages internally); non-kills go to phase B.
    st, nodes, unres = decide6(adj, LEVEL, seed=seed, order=order,
                               max_nodes=cap)
    return idx, st, nodes

def iter_graphs(skip):
    with open(GRAPHF) as f:
        for idx, line in enumerate(f):
            if idx in skip:
                continue
            parts = line.split()
            if not parts:
                continue
            n = int(parts[0]); assert n == LEVEL
            yield idx, [int(x) for x in parts[1:1 + n]], None

_ALL = None
def _load_all():
    global _ALL
    if _ALL is None:
        _ALL = []
        with open(GRAPHF) as f:
            for line in f:
                parts = line.split()
                if parts:
                    _ALL.append([int(x) for x in parts[1:1 + LEVEL]])
    return _ALL

def _slice_task(args):
    gi, seed, order, lo, hi, cap = args
    adj = _load_all()[gi]
    st, nodes, unres = decide6(adj, LEVEL, seed=seed, order=list(order),
                               th0=(lo, hi), max_nodes=cap)
    return lo, hi, st, nodes

def certify_graph(gi, adj, pool, max_orders=12, cap=15_000_000,
                  budget_s=3600, fine=True):
    """Race decompositions; 8-fold split failed slices (per reproduce5)."""
    import queue
    decs = gen_orders(adj, LEVEL, kmax=max_orders)
    if not decs:
        return None
    ncircs = [ncirc_of(adj, LEVEL, *d) for d in decs]
    Q = queue.Queue()
    state = [{"pending": 0, "alive": True} for _ in decs]
    ns1 = 384 if fine else 24
    min_width = TWO_PI / max(ns1, 24) / 8 ** 4

    def submit_range(oi, lo, hi):
        seed, order = decs[oi]
        state[oi]["pending"] += 1
        c = cap * 40 if ncircs[oi] == 0 else cap
        pool.apply_async(_slice_task,
                         ((gi, seed, tuple(order), lo, hi, c),),
                         callback=lambda r, oi=oi: Q.put((oi, r)),
                         error_callback=lambda e, oi=oi: Q.put((oi, e)))

    for oi in range(len(decs)):
        nc = ncircs[oi]
        ns = 768 if nc >= 2 else (ns1 if nc == 1 else 1)
        for k in range(ns):
            submit_range(oi, k * TWO_PI / ns, (k + 1) * TWO_PI / ns)
    deadline = time.time() + budget_s
    while True:
        alive = [oi for oi, s in enumerate(state) if s["alive"]]
        if not alive:
            return None
        for oi in alive:
            if state[oi]["pending"] == 0:
                return oi
        if time.time() > deadline:
            return None
        try:
            oi, r = Q.get(timeout=30)
        except Exception:
            continue
        s = state[oi]
        s["pending"] -= 1
        if isinstance(r, Exception):
            s["alive"] = False
        else:
            lo, hi, st, nodes = r
            if s["alive"] and st != "KILLED":
                if ncircs[oi] == 0 or hi - lo < min_width * 8:
                    s["alive"] = False
                else:
                    w = (hi - lo) / 8
                    for k in range(8):
                        submit_range(oi, lo + k * w, lo + (k + 1) * w)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=max(1, cpu_count() - 4))
    ap.add_argument("--cap", type=int, default=5_000_000)
    ap.add_argument("--limit", type=int, default=None,
                    help="process only the first K uncertified graphs")
    ap.add_argument("--bulk", action="store_true")
    ap.add_argument("--phaseb", action="store_true",
                    help="race the graphs listed in deferred_d6_n19.json")
    args = ap.parse_args()
    ensure_kernel()
    import hashlib
    print(f"kernel: ckernel6 (sha256 "
          f"{hashlib.sha256(open('ckernel6.c','rb').read()).hexdigest()[:12]})",
          flush=True)

    done = set()
    if os.path.exists(KILLLOG):
        with open(KILLLOG) as f:
            for line in f:
                p = line.split()
                if p: done.add(int(p[0]))
    print(f"already certified: {len(done)}", flush=True)

    if args.phaseb:
        targets = [gi for gi in json.load(open("deferred_d6_n19.json"))
                   if gi not in done]
        print(f"phase B: {len(targets)} graphs", flush=True)
        allg = _load_all()
        log = open(KILLLOG, "a", buffering=1)
        t0 = time.time()
        for i, gi in enumerate(targets):
            t = time.time()
            with Pool(args.workers) as pool:
                oi = certify_graph(gi, allg[gi], pool)
            if oi is not None:
                log.write(f"{gi} {oi}\n")
                print(f"[{gi}] CERTIFIED (dec {oi}) ({i+1}/{len(targets)}) "
                      f"t={time.time()-t:.1f}s", flush=True)
            else:
                print(f"[{gi}] UNDECIDED t={time.time()-t:.1f}s", flush=True)
        return

    if not args.bulk:
        print("use --bulk or --phaseb"); return

    t0 = time.time()
    killed = deferred = 0
    deferred_list = []
    log = open(KILLLOG, "a", buffering=1)
    gen = iter_graphs(done)
    if args.limit:
        import itertools
        gen = itertools.islice(gen, args.limit)
    with Pool(args.workers) as pool:
        argsit = ((idx, adj, args.cap) for idx, adj, _ in gen)
        for idx, st, nodes in pool.imap_unordered(_bulk_task, argsit,
                                                  chunksize=4):
            if st == "KILLED":
                killed += 1
                log.write(f"{idx} 0\n")
            else:
                deferred += 1
                deferred_list.append(idx)
            if (killed + deferred) % 20000 == 0:
                el = time.time() - t0
                rate = (killed + deferred) / el
                rem = (3971787 - len(done) - killed - deferred) / rate
                print(f"  bulk {killed+deferred} killed={killed} "
                      f"deferred={deferred} rate={rate:.0f}/s "
                      f"eta={rem/3600:.1f}h", flush=True)
    json.dump(deferred_list, open("deferred_d6_n19.json", "w"))
    print(f"bulk done: killed={killed} deferred={deferred} "
          f"t={(time.time()-t0)/3600:.2f}h", flush=True)
    print(f"deferred graphs saved to deferred_d6_n19.json for phase B")

if __name__ == "__main__":
    main()
