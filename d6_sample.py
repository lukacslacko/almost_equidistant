"""Randomized sample of the uncertified d=6 level-19 population:
measure kill-fraction and timing at a small cap to size the campaign."""
import random, time
from multiprocessing import Pool
from reproduce6 import _bulk_task, KILLLOG, iter_graphs

def main():
    done = set()
    for line in open(KILLLOG):
        p = line.split()
        if p: done.add(int(p[0]))
    random.seed(20260806)
    sample = []
    for i, (idx, adj, _) in enumerate(iter_graphs(done)):
        if len(sample) < 600:
            sample.append((idx, adj))
        else:
            j = random.randrange(i + 1)
            if j < 600: sample[j] = (idx, adj)
    print(f"sampled {len(sample)} uncertified graphs", flush=True)
    t0 = time.time()
    stats = {"KILLED": 0, "ABORT": 0, "other": 0}
    with Pool(24) as pool:
        for idx, st, nodes in pool.imap_unordered(
                _bulk_task, [(i, a, 300_000) for i, a in sample],
                chunksize=2):
            stats[st if st in stats else "other"] += 1
            n = sum(stats.values())
            if n % 100 == 0:
                print(f"  {n}/600: {stats} t={time.time()-t0:.0f}s",
                      flush=True)
    print(f"FINAL: {stats} in {time.time()-t0:.0f}s", flush=True)

if __name__ == "__main__":
    main()
