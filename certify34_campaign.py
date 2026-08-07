"""Resumable certification campaign for the 34 level-20 frontier graphs
(frontier20_true.txt).  Strategy per f(5) lessons: parametrization
diversity first, then fine slicing with escalating caps on the best
decomposition only.

Stages (global task streaming across graphs, worker-side caching):
  1. probe: every (graph, dec) x 6 slices @ 30k cap  -> rank decs
  2. tile:  best 2 decs per uncertified graph, full 96 slices @ 150k
            -> certified when a tiling is fully KILLED
  3. ladder: best dec's hot slices, subdivide 8x @ 500k, again 8x @ 2M
            -> certified when the hot set of some dec empties
State in frontier34_state.json (atomic rewrite); log to stdout.
Usage: python certify34_campaign.py [workers]
"""
import json, math, os, sys, time
from multiprocessing import Pool

TWO_PI = 2 * math.pi
NP = 20
NS = 96
GRAPHF = "frontier20_true.txt"
STATEF = "frontier34_state.json"

_cache = {}
def _get(pi):
    if pi not in _cache:
        from certify_frontier34 import load_pat, all_orders
        adj = load_pat(pi)
        _cache[pi] = (adj, all_orders(adj))
    return _cache[pi]

def _task(args):
    pi, di, lo, hi, cap = args
    from cdriver6 import decide6
    adj, ords = _get(pi)
    seed, order = ords[di]
    t0 = time.time()
    st, nodes, unres = decide6(adj, NP, seed=seed, order=list(order),
                               th0=(lo, hi), max_nodes=cap)
    return pi, di, lo, hi, cap, st, nodes, time.time() - t0

def save(state):
    tmp = STATEF + ".tmp"
    json.dump(state, open(tmp, "w"))
    os.replace(tmp, STATEF)

def main():
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 26
    state = json.load(open(STATEF)) if os.path.exists(STATEF) else {
        "probe": {}, "tile": {}, "certified": {}, "ladder": {}}
    ngraphs = 34
    ndecs = {}
    from certify_frontier34 import load_pat, all_orders
    for pi in range(ngraphs):
        ndecs[pi] = len(all_orders(load_pat(pi)))
    print("dec counts:", ndecs, flush=True)

    pool = Pool(workers)
    t00 = time.time()

    # ---- stage 1: probe ----
    ks = [0, 16, 32, 48, 64, 80]
    tasks = []
    for pi in range(ngraphs):
        if str(pi) in state["certified"]:
            continue
        for di in range(ndecs[pi]):
            key = f"{pi}:{di}"
            if key in state["probe"]:
                continue
            for k in ks:
                tasks.append((pi, di, k * TWO_PI / NS, (k + 1) * TWO_PI / NS,
                              30_000))
    print(f"stage1: {len(tasks)} probe tasks", flush=True)
    acc = {}
    done = 0
    for pi, di, lo, hi, cap, st, nodes, dt in pool.imap_unordered(
            _task, tasks, chunksize=1):
        key = f"{pi}:{di}"
        acc.setdefault(key, []).append((st, nodes))
        done += 1
        if len(acc[key]) == len(ks):
            r = acc.pop(key)
            state["probe"][key] = [sum(1 for s, _ in r if s != "KILLED"),
                                   sum(n for _, n in r)]
        if done % 500 == 0:
            save(state)
            print(f"stage1: {done}/{len(tasks)} t={time.time()-t00:.0f}s",
                  flush=True)
    save(state)
    print(f"stage1 done t={time.time()-t00:.0f}s", flush=True)

    # ---- stage 2: full 96-tiles on the best 2 decs ----
    tasks = []
    for pi in range(ngraphs):
        if str(pi) in state["certified"]:
            continue
        rank = sorted((state["probe"][f"{pi}:{di}"], di)
                      for di in range(ndecs[pi])
                      if f"{pi}:{di}" in state["probe"])
        for (_, di) in rank[:2]:
            key = f"{pi}:{di}"
            if key in state["tile"]:
                continue
            for k in range(NS):
                tasks.append((pi, di, k * TWO_PI / NS, (k + 1) * TWO_PI / NS,
                              150_000))
    print(f"stage2: {len(tasks)} tile tasks", flush=True)
    acc = {}
    done = 0
    for pi, di, lo, hi, cap, st, nodes, dt in pool.imap_unordered(
            _task, tasks, chunksize=1):
        key = f"{pi}:{di}"
        acc.setdefault(key, []).append((lo, st))
        done += 1
        if len(acc[key]) == NS:
            r = acc.pop(key)
            hot = sorted(lo for lo, s in r if s != "KILLED")
            state["tile"][key] = hot
            if not hot and str(pi) not in state["certified"]:
                state["certified"][str(pi)] = {"dec": di, "stage": "tile96",
                                               "cap": 150_000}
                print(f"*** [{pi}] CERTIFIED (dec {di}, 96-tile @150k) "
                      f"t={time.time()-t00:.0f}s", flush=True)
            else:
                print(f"[{pi}] dec {di}: {len(hot)} hot slices "
                      f"t={time.time()-t00:.0f}s", flush=True)
            save(state)
        if done % 200 == 0:
            print(f"stage2: {done}/{len(tasks)} t={time.time()-t00:.0f}s",
                  flush=True)
    save(state)
    print(f"stage2 done; certified {len(state['certified'])}/34 "
          f"t={time.time()-t00:.0f}s", flush=True)

    # ---- stage 3: subdivision ladder on the best dec ----
    for rung, (split, cap) in enumerate([(8, 500_000), (8, 2_000_000)]):
        tasks = []
        meta = {}
        for pi in range(ngraphs):
            if str(pi) in state["certified"]:
                continue
            cands = [(len(hot), key) for key, hot in state["tile"].items()
                     if int(key.split(":")[0]) == pi]
            if not cands:
                continue
            cands.sort()
            nhot, key = cands[0]
            di = int(key.split(":")[1])
            prev = state["ladder"].get(key, {"arcs": [
                [lo, lo + TWO_PI / NS] for lo in state["tile"][key]],
                "rung": 0})
            arcs = prev["arcs"]
            meta[key] = {"arcs": [], "rung": rung + 1}
            for (lo, hi) in arcs:
                w = (hi - lo) / split
                for j in range(split):
                    tasks.append((pi, di, lo + j * w, lo + (j + 1) * w, cap))
        print(f"stage3 rung {rung}: {len(tasks)} tasks", flush=True)
        done = 0
        for pi, di, lo, hi, cap, st, nodes, dt in pool.imap_unordered(
                _task, tasks, chunksize=1):
            key = f"{pi}:{di}"
            done += 1
            if st != "KILLED":
                meta[key]["arcs"].append([lo, hi])
            if done % 200 == 0:
                print(f"stage3 rung {rung}: {done}/{len(tasks)} "
                      f"t={time.time()-t00:.0f}s", flush=True)
        for key, m in meta.items():
            pi = int(key.split(":")[0])
            di = int(key.split(":")[1])
            state["ladder"][key] = m
            if not m["arcs"] and str(pi) not in state["certified"]:
                state["certified"][str(pi)] = {"dec": di,
                                               "stage": f"ladder{rung}",
                                               "cap": cap}
                print(f"*** [{pi}] CERTIFIED (dec {di}, ladder rung {rung}) "
                      f"t={time.time()-t00:.0f}s", flush=True)
            elif m["arcs"]:
                print(f"[{pi}] dec {di}: {len(m['arcs'])} hot arcs after "
                      f"rung {rung}", flush=True)
        save(state)
    pool.close()
    pool.join()
    print(f"campaign pass complete: certified {len(state['certified'])}/34; "
          f"remaining need deeper rungs or new ideas "
          f"t={time.time()-t00:.0f}s", flush=True)

if __name__ == "__main__":
    main()
