#!/usr/bin/env python3
"""Finisher for the hard tail of level 17.

Strategy (from the diagnosis of graph 11402): each hard graph's circle
parameter dies everywhere except one or two narrow "hot" arcs (the
degenerate-channel phenomenon known from the f(4) campaign). So:

  1. per graph, take the best few decompositions (earliest-circle first);
  2. tile the circle into 96 coarse slices with a small node cap - the
     cold ~95 die immediately and cheaply;
  3. drill ONLY the hot arcs: 16-fold subdivision with an escalating
     node-cap ladder (8M -> 40M -> 200M -> 1G), to a hard width floor;
  4. a graph is certified when one decomposition's tiling of the full
     circle is entirely KILLED (identical certificate semantics to
     reproduce5.certify_graph - a union of closed slices covering
     [0, 2pi], every slice KILLED).

All graphs' slice tasks share one pool: hot arcs from many graphs drill
simultaneously; no graph head-blocks another.

Usage: python finish17.py [--workers N]
Writes results into results_n17.json (merging, resumable).
"""
import argparse, json, math, os, time
from multiprocessing import Pool, cpu_count

from cdriver5 import gen_orders, decide5
from reproduce5 import load_level, ncirc_of

TWO_PI = 2 * math.pi
NS0 = 96                      # coarse scan slices
CAP0 = 400_000                # coarse scan node cap
# splitting, not budget, is what kills channels (the d=4 lesson): keep
# tasks short (<= ~10 min) and keep subdividing
LADDER = (8_000_000, 30_000_000, 80_000_000, 80_000_000)
SPLIT = 16
MIN_W = TWO_PI / NS0 / 16 ** 6        # hard width floor
MAX_DECS = 4

_G = None
def _task(args):
    global _G
    gi, di, lo, hi, cap = args
    if _G is None:
        _G = load_level(17)
    adj = _G[gi]
    decs = gen_orders(adj, 17, kmax=12)
    seed, order = decs[di]
    # deep internal theta floor for drilled (hot) tasks: the engine's
    # adaptive splitting follows only surviving cells, so depth is cheap
    # inside one task, whereas driver-level splitting costs 8^depth.
    tm = TWO_PI / (1 << 18) if cap <= CAP0 else TWO_PI / (1 << 30)
    st, nodes, unres = decide5(adj, 17, seed=seed, order=order,
                               th0=(lo, hi), theta_min=tm, max_nodes=cap)
    return gi, di, lo, hi, cap, st, nodes

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=max(1, cpu_count() - 2))
    args = ap.parse_args()

    graphs = load_level(17)
    resfile = "results_n17.json"
    results = {}
    if os.path.exists(resfile):
        results = {int(k): v for k, v in json.load(open(resfile)).items()}
    todo = [gi for gi in range(len(graphs)) if results.get(gi) is None]
    print(f"remaining graphs: {len(todo)}", flush=True)
    if not todo:
        print("nothing to do"); return

    # per-graph scheduling state
    class GS:
        def __init__(self, gi):
            self.gi = gi
            self.decs = gen_orders(graphs[gi], 17, kmax=12)
            # circle-bearing decs only, capped
            self.dis = [i for i, d in enumerate(self.decs)
                        if ncirc_of(graphs[gi], 17, *d) >= 1][:MAX_DECS]
            if not self.dis:            # circle-free graph: single big task
                self.dis = [0]
            self.cur = 0                # index into self.dis
            self.pending = 0
            self.failed_dec = False
            self.done = False
        def cur_di(self):
            return self.dis[self.cur]

    state = {gi: GS(gi) for gi in todo}
    pool = Pool(args.workers)
    results_q = []

    import queue as _q
    Q = _q.Queue()
    inflight = 0

    def submit(gi, di, lo, hi, cap):
        nonlocal inflight
        state[gi].pending += 1
        inflight += 1
        pool.apply_async(_task, ((gi, di, lo, hi, cap),),
                         callback=lambda r: Q.put(r),
                         error_callback=lambda e, gi=gi, di=di: Q.put(
                             (gi, di, 0.0, 0.0, 0, "ERROR", 0)))

    def start_dec(gi):
        gs = state[gi]
        di = gs.cur_di()
        d = gs.decs[di]
        if ncirc_of(graphs[gi], 17, *d) == 0:
            submit(gi, di, 0.0, TWO_PI, LADDER[-2])
        else:
            for k in range(NS0):
                submit(gi, di, k * TWO_PI / NS0, (k + 1) * TWO_PI / NS0, CAP0)

    for gi in todo:
        start_dec(gi)

    t0 = time.time()
    ncert = 0
    while inflight > 0:
        gi, di, lo, hi, cap, st, nodes = Q.get()
        inflight -= 1
        gs = state[gi]
        gs.pending -= 1
        if gs.done:
            pass
        elif st == "KILLED":
            if gs.pending == 0 and not gs.failed_dec:
                gs.done = True
                results[gi] = gs.cur_di()
                ncert += 1
                json.dump(results, open(resfile, "w"))
                print(f"[{gi}] CERTIFIED (dec {gs.cur_di()}) "
                      f"({ncert}/{len(todo)}, t={time.time()-t0:.0f}s)",
                      flush=True)
        elif st == "ERROR":
            gs.failed_dec = True
        else:
            # hot interval: split and escalate the node cap
            w = hi - lo
            if w > MIN_W:
                nxt = LADDER[-1]
                for c in LADDER:
                    if c > cap: nxt = c; break
                sub = SPLIT if cap <= CAP0 else 8
                for k in range(sub):
                    submit(gi, di, lo + k * w / sub, lo + (k + 1) * w / sub, nxt)
            else:
                gs.failed_dec = True
        # decomposition failed and drained: move to next dec
        if (not gs.done and gs.failed_dec and gs.pending == 0):
            gs.cur += 1
            gs.failed_dec = False
            if gs.cur < len(gs.dis):
                start_dec(gi)
            else:
                print(f"[{gi}] UNDECIDED after {len(gs.dis)} decs", flush=True)
    pool.terminate()
    missing = [gi for gi in todo if results.get(gi) is None]
    print(f"\nfinished: {len(todo)-len(missing)}/{len(todo)} certified; "
          f"missing: {missing}")
    if not missing:
        total = sum(1 for v in results.values() if v is not None)
        print(f"level 17 total: {total}/12654")
        if total == 12654:
            print("==> no 17-point almost-equidistant set in R^5; "
                  "with f(5)>=16: f(5) = 16.")

if __name__ == "__main__":
    main()
