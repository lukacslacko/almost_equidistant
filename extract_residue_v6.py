"""Extract the v6 combined residue (644 graphs) from the codex manifest into
corpus line format ("19 m0 ... m18"), verifying that every manifest adjacency
record is byte-identical to the corresponding line of the level-19 corpus
(aeq_d6_n19.txt) at its claimed index.  Emits residue644.txt (sorted by
corpus index) and residue644_indices.txt.
"""
import json, sys

MAN = "d6_current_residue_manifest_v6.json"
CORPUS = "aeq_d6_n19.txt"

m = json.load(open(MAN))
records = []
for cls in ("K7", "K6_only"):
    for rec in m["classes"][cls]["graphs"]:
        records.append((rec["index"], rec["adjacency"], cls))
records.sort()
idxs = [r[0] for r in records]
assert len(idxs) == len(set(idxs)) == 644, len(idxs)
assert m["combined"]["count"] == 644

want = dict((i, (adj, cls)) for i, adj, cls in records)
found = {}
with open(CORPUS) as f:
    for lineno, line in enumerate(f):
        if lineno in want:
            p = line.split()
            assert p[0] == "19", (lineno, p[0])
            found[lineno] = [int(x) for x in p[1:20]]
assert len(found) == 644, (len(found), "corpus shorter than max index?")

mismatch = 0
for i, (adj, cls) in want.items():
    if found[i] != adj:
        mismatch += 1
        print("MISMATCH at corpus line", i, file=sys.stderr)
assert mismatch == 0, mismatch

with open("residue644.txt", "w") as out:
    for i, adj, cls in records:
        out.write("19 " + " ".join(str(x) for x in adj) + "\n")
with open("residue644_indices.txt", "w") as out:
    for i, adj, cls in records:
        out.write(f"{i} {cls}\n")
print("OK: 644 residue graphs extracted; all adjacencies byte-identical to",
      CORPUS, "at their claimed indices",
      f"({sum(1 for r in records if r[2]=='K7')} K7,",
      f"{sum(1 for r in records if r[2]=='K6_only')} K6-only)")
