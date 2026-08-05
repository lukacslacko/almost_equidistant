"""Merge shard results and verify campaign completeness.

For each level: every graph index 0..N-1 must have a non-None winning
decomposition recorded (i.e. a fully-killed tiling certificate). Prints
the per-level status and the overall conclusion for f(5).
"""
import glob, json, os, sys

EXPECT = {17: 12654, 18: 8825, 19: 340, 20: 8}

def merge_level(level):
    main = f"results_n{level}.json"
    res = {}
    if os.path.exists(main):
        res = {int(k): v for k, v in json.load(open(main)).items()}
    for f in glob.glob(f"results_n{level}_s*.json"):
        for k, v in json.load(open(f)).items():
            k = int(k)
            if v is not None:
                if res.get(k) is not None and res[k] != v:
                    # both certified via different decompositions: keep main
                    continue
                res[k] = v
    json.dump(res, open(main, "w"))
    return res

def main():
    allok = True
    for level in (20, 19, 18, 17):
        res = merge_level(level)
        n = EXPECT[level]
        missing = [i for i in range(n) if res.get(i) is None]
        ok = not missing
        allok &= ok
        print(f"level {level}: {n - len(missing)}/{n} certified"
              + ("" if ok else f"  MISSING: {missing[:20]}"))
    if allok:
        print("\nAll candidate graphs on 17..20 vertices are certified "
              "non-realizable in R^5 with all points distinct.")
        print("==> f(5) <= 16; with the Larman-Rogers half-cube set, "
              "f(5) = 16.")
    else:
        print("\ncampaign INCOMPLETE")
    return 0 if allok else 1

if __name__ == "__main__":
    sys.exit(main())
