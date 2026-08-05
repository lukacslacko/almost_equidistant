#!/usr/bin/env python3
"""Reproduce the certified decision f(4) = 12.

Steps (see README.md and Appendices A-C of f4_equals_12.pdf):
  1. --enumerate : regenerate the 59 abstract almost-equidistant graphs and
                   verify the counts against BPSSV Table 3 (runs enumerate_aeq).
  2. --controls  : validation controls for the decision engine:
                   G10 (known non-realizable) must be KILLED;
                   the cross-polytope graph (realizable) must NOT be killed.
  3. (default)   : certify all 59 graphs. For each graph, elimination orders
                   are generated (circle-late layered orders); for each order
                   the outer-circle parameter [0, 2pi) is tiled into slices,
                   each decided by the certified C engine; slices that exceed
                   their node budget are split 8-fold and retried. A graph is
                   CERTIFIED as soon as one order's tiling is fully killed.

Requirements: a C compiler (cc), Python >= 3.9, numpy.
The C kernel is compiled automatically on first run.
"""
import argparse, json, math, os, subprocess, sys, time
from multiprocessing import Pool, cpu_count

TWO_PI = 2 * math.pi
HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)

def ensure_kernel():
    dylib, src = "ckernel.dylib", "ckernel.c"
    if (not os.path.exists(dylib)
            or os.path.getmtime(dylib) < os.path.getmtime(src)):
        print("compiling C kernel ...", flush=True)
        subprocess.check_call(["cc", "-O2", "-shared", "-o", dylib, src, "-lm"])

# ---------------------------------------------------------------------------
# elimination orders (layered, circles as late as possible; several tie-breaks)
# ---------------------------------------------------------------------------
def _cliques(adj, n, k):
    res = []
    def ext(clq, cand, start):
        if len(clq) == k: res.append(tuple(clq)); return
        for v in range(start, n):
            if (cand >> v) & 1: ext(clq + [v], cand & adj[v], v + 1)
    ext([], (1 << n) - 1, 0)
    return res

def _layered_order(adj, n, seed, pick):
    placed = set(seed); order = []
    while len(placed) < n:
        while True:
            cands = [(sum(1 for u in placed if (adj[v] >> u) & 1), v)
                     for v in range(n) if v not in placed]
            cands = [(k, v) for k, v in cands if k >= 4]
            if not cands: break
            k, v = max(cands); order.append(v); placed.add(v)
        if len(placed) == n: break
        c3 = [v for v in range(n) if v not in placed
              and sum(1 for u in placed if (adj[v] >> u) & 1) == 3]
        if not c3: return None
        v = pick(c3, placed); order.append(v); placed.add(v)
    return order

def gen_orders(adj, n, kmax=8):
    seeds = _cliques(adj, n, 5) or _cliques(adj, n, 4)
    def closure_size(w, placed):
        P = set(placed); P.add(w); grew = True
        while grew:
            grew = False
            for x in range(n):
                if x not in P and sum(1 for u in P if (adj[x] >> u) & 1) >= 4:
                    P.add(x); grew = True
        return len(P)
    picks = [lambda c3, p: max(c3, key=lambda w: closure_size(w, p)),
             lambda c3, p: min(c3, key=lambda w: closure_size(w, p)),
             lambda c3, p: max(c3), lambda c3, p: min(c3)]
    outs, seen = [], set()
    for seed in seeds:
        for pick in picks:
            o = _layered_order(adj, n, seed, pick)
            if o and (seed, tuple(o)) not in seen:
                seen.add((seed, tuple(o))); outs.append((seed, o))
        if len(outs) >= kmax * 4: break
    def score(so):
        seed, order = so
        placed = set(seed); circ = []
        for i, v in enumerate(order):
            if sum(1 for u in placed if (adj[v] >> u) & 1) == 3: circ.append(i)
            placed.add(v)
        return (len(circ), [-c for c in circ])
    outs.sort(key=score)
    return outs[:kmax]

# ---------------------------------------------------------------------------
# sharded certification of one graph under one decomposition
# ---------------------------------------------------------------------------
_DATA = None
def _slice_task(args):
    global _DATA
    gi, seed, order, lo, hi, cap = args
    from cdriver import decide_c
    if _DATA is None:
        _DATA = json.load(open("aeq_d4_n13.json"))
    st, nodes, unres = decide_c(_DATA["graphs"][gi], 13, seed=seed, order=order,
                                th0=(lo, hi), max_nodes=cap)
    return lo, hi, st, nodes

def certify_graph(gi, adj, pool, max_orders=8, cap=15_000_000,
                  min_width=TWO_PI / (24 * 8 ** 3), budget_s=7200,
                  verbose=True):
    """Race decompositions concurrently: each maintains a frontier of circle
    slices; failed slices split 8-fold; a decomposition whose tiling is fully
    KILLED certifies the graph. Slice budgets/splitting affect only search
    management, never verdict soundness."""
    import queue
    decs = gen_orders(adj, 13, kmax=max_orders)
    Q = queue.Queue()
    state = [{"pending": 0, "alive": True, "submitted": False} for _ in decs]

    def submit_range(oi, lo, hi):
        seed, order = decs[oi]
        state[oi]["pending"] += 1
        pool.apply_async(_slice_task, ((gi, seed, order, lo, hi, cap),),
                         callback=lambda r, oi=oi: Q.put((oi, r)),
                         error_callback=lambda e, oi=oi: Q.put((oi, e)))

    def submit_dec(oi):
        state[oi]["submitted"] = True
        for k in range(24):
            submit_range(oi, k * TWO_PI / 24, (k + 1) * TWO_PI / 24)

    for oi in range(len(decs)):        # all decompositions race from the start
        submit_dec(oi)
    deadline = time.time() + budget_s
    while True:
        alive = [oi for oi, s in enumerate(state) if s["alive"]]
        if not alive:
            return None
        for oi in alive:
            if state[oi]["submitted"] and state[oi]["pending"] == 0:
                return oi                      # fully killed tiling: certified
        if time.time() > deadline:
            return None                        # per-graph wall budget exceeded
        try:
            oi, r = Q.get(timeout=30)
        except queue.Empty:
            continue                           # slices still in flight: wait
        s = state[oi]
        s["pending"] -= 1
        if isinstance(r, Exception):
            s["alive"] = False
        else:
            lo, hi, st, nodes = r
            if s["alive"] and st != "KILLED":
                if hi - lo < min_width * 8:
                    s["alive"] = False         # width floor: give up this order
                else:
                    w = (hi - lo) / 8
                    for k in range(8):
                        submit_range(oi, lo + k * w, lo + (k + 1) * w)

# ---------------------------------------------------------------------------
# validation controls
# ---------------------------------------------------------------------------
def run_controls():
    from cdriver import decide_c
    # G10: BPSSV Lemma 12 graph, proven non-realizable -> must be KILLED
    n = 10; E = set()
    for i in range(8):
        for j in range(i + 1, 8):
            if j != i + 4: E.add((i, j))
    E.discard((2, 4))
    for u in (0, 1, 2, 3): E.add(tuple(sorted((8, u))))
    for u in (2, 3, 4, 5): E.add(tuple(sorted((9, u))))
    E.add((8, 9))
    adj = [0] * n
    for i, j in E: adj[i] |= 1 << j; adj[j] |= 1 << i
    st, nodes, unres = decide_c(adj, n, max_nodes=100_000_000)
    print(f"control G10 (expect KILLED): {st}  nodes={nodes}")
    ok1 = st == "KILLED"
    # cross-polytope graph K_{2,2,2,2}: realizable -> must NOT be killed
    n2 = 8; adj2 = [0] * n2
    for i in range(8):
        for j in range(8):
            if i != j and j != (i + 4) % 8 and i != (j + 4) % 8:
                adj2[i] |= 1 << j
    st2, nodes2, unres2 = decide_c(adj2, n2, theta_min=TWO_PI / (1 << 12),
                                   max_nodes=100_000_000)
    print(f"control cross-polytope (expect SURVIVORS): {st2}  "
          f"nodes={nodes2} unresolved={unres2}")
    ok2 = st2 == "SURVIVORS"
    print("controls:", "PASS" if ok1 and ok2 else "FAIL")
    return ok1 and ok2

# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--enumerate", action="store_true",
                    help="regenerate and verify the 59 graphs")
    ap.add_argument("--controls", action="store_true",
                    help="run the two validation controls")
    ap.add_argument("--graph", type=int, default=None,
                    help="certify a single graph index (0..58)")
    ap.add_argument("--workers", type=int, default=max(1, cpu_count() - 1))
    args = ap.parse_args()
    ensure_kernel()

    if args.enumerate:
        subprocess.check_call([sys.executable, "enumerate_aeq.py"])
        return
    if args.controls:
        sys.exit(0 if run_controls() else 1)

    if not run_controls():
        print("ABORT: validation controls failed"); sys.exit(1)
    data = json.load(open("aeq_d4_n13.json"))
    graphs = data["graphs"]
    targets = [args.graph] if args.graph is not None else range(len(graphs))
    results = {}
    t0 = time.time()
    with Pool(args.workers) as pool:
        for gi in targets:
            t = time.time()
            oi = certify_graph(gi, graphs[gi], pool)
            results[gi] = oi
            status = (f"CERTIFIED non-realizable (decomposition {oi})"
                      if oi is not None else "UNDECIDED - see README")
            print(f"[{gi:2d}] {status}  t={time.time()-t:.1f}s "
                  f"(total {time.time()-t0:.0f}s)", flush=True)
            json.dump(results, open("reproduce_results.json", "w"), indent=1)
    bad = [g for g, o in results.items() if o is None]
    if not bad and args.graph is None:
        print("\nAll 59 graphs certified non-realizable with 13 distinct "
              "points.\n==> no 13-point almost-equidistant set in R^4; "
              "with the known 12-point construction, f(4) = 12.")
    elif bad:
        print(f"\nUndecided graphs: {bad} (raise budgets/orders; see README)")

if __name__ == "__main__":
    main()
