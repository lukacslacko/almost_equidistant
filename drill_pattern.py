"""Ladder-drill the hot slices of a near-miss obstruction pattern.
Finds hot slices at a small cap, then splits them with escalating caps
and deep theta floors until the full 96-tiling of one parametrization is
KILLED (=> pattern certified) or the width floor is hit.
Usage: python drill_pattern.py <pattern_idx> [dec] [workers]
"""
import json, math, os, sys, time
from multiprocessing import Pool
from certify_patterns import load_pat, all_orders, NP
from cdriver6 import decide6

TWO_PI = 2 * math.pi
NS = 96
LADDER = (300_000, 2_000_000, 8_000_000, 30_000_000, 30_000_000)
MIN_W = TWO_PI / NS / 8 ** 6

_pat = None
_ord = None
def _init(pi, di):
    global _pat, _ord
    _pat = load_pat(pi)
    _ord = all_orders(_pat)[di]

def _run(args):
    lo, hi, cap = args
    seed, order = _ord
    st, nodes, unres = decide6(_pat, NP, seed=seed, order=list(order),
                               th0=(lo, hi),
                               theta_min=TWO_PI / (1 << 30), max_nodes=cap)
    return lo, hi, cap, st

def main():
    pi = int(sys.argv[1])
    di = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    workers = int(sys.argv[3]) if len(sys.argv) > 3 else 12
    t0 = time.time()
    work = [(k * TWO_PI / NS, (k + 1) * TWO_PI / NS, LADDER[0])
            for k in range(NS)]
    pending = len(work)
    failed = False
    with Pool(workers, initializer=_init, initargs=(pi, di)) as pool:
        results = [pool.apply_async(_run, (w,)) for w in work]
        queue = list(results)
        while queue:
            nq = []
            for r in queue:
                if not r.ready():
                    nq.append(r); continue
                lo, hi, cap, st = r.get()
                if st == "KILLED":
                    continue
                w = hi - lo
                if w <= MIN_W:
                    print(f"  floor-park [{math.degrees(lo):.4f},"
                          f"{math.degrees(hi):.4f}] — dec fails", flush=True)
                    failed = True
                    continue
                nxt = LADDER[-1]
                for c in LADDER:
                    if c > cap: nxt = c; break
                for j in range(8):
                    nq.append(pool.apply_async(
                        _run, ((lo + j * w / 8, lo + (j + 1) * w / 8, nxt),)))
            queue = nq
            time.sleep(2)
    if not failed:
        print(f"[{pi}] CERTIFIED via drilled tiling of dec {di} "
              f"t={time.time()-t0:.0f}s", flush=True)
        fn = "obstructions_d6.json"
        obs = json.load(open(fn)) if os.path.exists(fn) else {}
        obs[str(pi)] = {"dec": di, "tiling": f"{NS}+drill"}
        json.dump(obs, open(fn, "w"))
    else:
        print(f"[{pi}] dec {di} fails at floor t={time.time()-t0:.0f}s",
              flush=True)

if __name__ == "__main__":
    main()
