"""8-way parallel probe of one hot slice with a deep internal theta floor.
Usage: python probe_deep.py <gi> <k96> [cap]"""
import sys, time, math
from multiprocessing import Pool
from cdriver5 import gen_orders, decide5
from reproduce5 import load_level

TWO_PI = 2 * math.pi

def sub(args):
    gi, lo, hi, cap = args
    graphs = load_level(17)
    adj = graphs[gi]
    decs = gen_orders(adj, 17, kmax=12)
    seed, order = decs[0]
    t = time.time()
    st, nodes, unres = decide5(adj, 17, seed=seed, order=order, th0=(lo, hi),
                               theta_min=TWO_PI / (1 << 30), max_nodes=cap)
    return lo, hi, st, nodes, unres, time.time() - t

def main():
    gi = int(sys.argv[1]); k = int(sys.argv[2])
    cap = int(sys.argv[3]) if len(sys.argv) > 3 else 40_000_000
    lo0, hi0 = k * TWO_PI / 96, (k + 1) * TWO_PI / 96
    tasks = [(gi, lo0 + j * (hi0 - lo0) / 8, lo0 + (j + 1) * (hi0 - lo0) / 8,
              cap) for j in range(8)]
    t0 = time.time()
    with Pool(8) as pool:
        for lo, hi, st, nodes, unres, dt in pool.imap_unordered(sub, tasks):
            print(f"[{math.degrees(lo):.3f},{math.degrees(hi):.3f}] {st} "
                  f"nodes={nodes} unres={unres} t={dt:.0f}s", flush=True)
    print(f"total t={time.time()-t0:.0f}s", flush=True)

if __name__ == "__main__":
    main()
