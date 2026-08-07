"""Final formality for the frontier double-check: the C-flow's 34
frontier graphs (frontier20_true.txt) must all be contained, up to
isomorphism, in the INDEPENDENT implementation's extension classes of
residue644.txt.  (Counts already agree at every stage; this pins the
set.)  Proper file script: Windows multiprocessing cannot spawn from
stdin-fed scripts.
"""
from verify_frontier_indep import parse, IsoSet, extensions_of
from multiprocessing import Pool

def main():
    res = parse("residue644.txt")
    with Pool(8) as pool:
        chunks = pool.map(extensions_of, res)
    s = IsoSet()
    for ch in chunks:
        for a in ch:
            s.add(20, a)
    print("indep extension classes:", s.count)
    c34 = parse("frontier20_true.txt")
    missing = [i for i, (n, a) in enumerate(c34) if not s.contains(20, a)]
    print("C-flow frontier graphs missing from indep extension set:", missing)
    print("CONTAINMENT:", "PASS" if not missing and s.count == 117290 else "FAIL")

if __name__ == "__main__":
    main()
