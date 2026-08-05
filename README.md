# Reproduction package: f(4) = 12

This package contains exactly what is needed to reproduce the computer-assisted
proof that the maximum size of an almost-equidistant set in R^4 is 12
(Conjecture 1 of Balko–Pór–Scheucher–Swanepoel–Valtr, *Graphs and
Combinatorics* 36 (2020) 729–754). The accompanying note `f4_equals_12.pdf`
states the result; its Appendices A–C document the algorithms implemented here.

## Contents

| file | role |
|---|---|
| `f4_equals_12.pdf` / `.tex` | the result note with method appendices |
| `enumerate_aeq.py` | enumerates the abstract almost-equidistant graphs in R^4 and verifies all counts against BPSSV Table 3 (Appendix A) |
| `aeq_d4_n13.json` | the resulting 59 candidate graphs on 13 vertices (regenerable by the above) |
| `ival.py` | rigorous outward-rounded interval arithmetic (Appendix B.1) |
| `decide.py` | the certified branch-and-prune engine — readable Python reference (Appendix B) |
| `ckernel.c` | the same engine in C (~1000x faster); produces node-for-node identical search trees |
| `cdriver.py` | ctypes driver for the C kernel |
| `reproduce.py` | one-command orchestration: controls, slicing, alternate decompositions (Appendix C) |

## Requirements

POSIX system, a C compiler (`cc`), Python >= 3.9, `numpy`.
The C kernel compiles automatically on first use.

## How to reproduce

```sh
# 1. Regenerate the 59 candidate graphs and check every count against the
#    published tables (~15 s):
python3 reproduce.py --enumerate

# 2. Validation controls (~2 min): the engine must kill BPSSV's known
#    non-realizable graph G10, and must NOT kill the realizable
#    cross-polytope graph:
python3 reproduce.py --controls

# 3. Certify all 59 graphs (hours; uses all cores minus one by default):
python3 reproduce.py
# or a single graph:
python3 reproduce.py --graph 44
```

Most graphs certify in milliseconds to seconds. Three graphs (18, 20, 44 in
this ordering) require the slicing machinery and take the bulk of the time;
graph 18 in addition requires an alternate elimination order (the driver finds
it automatically — historically decomposition 6). A full run took roughly half
a day on a 12-core machine. Per-decomposition budgets and slice caps are
parameters of `certify_decomposition()` in `reproduce.py` and can be raised.

Progress is printed per graph and saved to `reproduce_results.json`. The final
message states the conclusion. The computation is deterministic: identical
hardware/libm produce identical search trees.

## What "certified" means

A graph is certified when, for one elimination order, every slice of a tiling
of the circle parameter is decided KILLED by the interval engine — i.e., every
branch of the placement search is eliminated by an interval that certifiably
excludes a required value, or by one of the two coincidence rules (which only
discard branches that provably cannot contain 13 *distinct* points). Soundness
rests on IEEE-754 correctly rounded arithmetic, libm cos/sin accurate to 8 ulp
(padded), and the engine logic described in Appendix B. A surviving cell would
be reported as SURVIVORS/UNDECIDED — never silently dropped.

## Provenance

Computation carried out 2026-08-04/05. The original run logs are preserved
separately (directory `f4_proof_artifacts` of the same project); this package
is the minimal, cleaned re-derivation path.

## License

MIT — see `LICENSE`.
