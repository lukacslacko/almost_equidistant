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
| `ckernel.c` | the certified branch-and-prune decision engine (Appendix B) |
| `cdriver.py` | compiles and drives the C kernel; constructs the exact seed-simplex enclosures and the elimination orders |
| `ival.py` | minimal outward-rounded interval arithmetic used only for the seed enclosures (Appendix B.1; the full interval engine is in `ckernel.c`) |
| `reproduce.py` | one-command orchestration: controls, slicing, alternate decompositions (Appendix C) |

## Requirements

Python >= 3.9 and a C compiler; the kernel compiles automatically on first
use, and there are no Python package dependencies.

| OS | C toolchain |
|---|---|
| Linux | `cc`/`gcc`/`clang` (e.g. `apt install build-essential`) |
| macOS | Xcode Command Line Tools (`xcode-select --install`) |
| Windows | any one of: LLVM clang (`winget install LLVM.LLVM`), Visual Studio 2022 or Build Tools with the C++ workload (located automatically via `vswhere`), or a `gcc` on `PATH` |

## How to reproduce

The same three commands on every OS — on Linux/macOS write `python3`, on
Windows write `python` (any shell: cmd, PowerShell, or Git Bash):

```sh
# 1. Regenerate the 59 candidate graphs and check every count against the
#    published tables (~10 s):
python3 reproduce.py --enumerate

# 2. Validation controls (~3 min): the engine must kill BPSSV's known
#    non-realizable graph G10, and must NOT kill the realizable
#    cross-polytope graph:
python3 reproduce.py --controls

# 3. Certify all 59 graphs (~10 min on a modern 16-core machine; uses all
#    cores minus one by default):
python3 reproduce.py
# or a single graph:
python3 reproduce.py --graph 44
```

The first command's last line must read `n=17: 0 abstract graphs`, the second
must end `controls: PASS`, and the third must end with

```
All 59 graphs certified non-realizable with 13 distinct points.
==> no 13-point almost-equidistant set in R^4; with the known 12-point construction, f(4) = 12.
```

At startup the driver prints which kernel binary it bound and a hash of
`ckernel.c` — check it against the repository if you have ever built the
kernel elsewhere on the machine.

The driver races all candidate elimination orders concurrently and certifies
with whichever completes first; two-parameter graphs (two circle stages —
graph 44) automatically start at fine initial slicing. Measured on a Ryzen
9950X3D (16C/32T, LLVM clang kernel): 58 of the 59 graphs certify in under a
second of wall time each, graph 44 in ~5 minutes; the whole run, controls
included, takes under 10 minutes. Per-graph budgets and slice caps are
parameters of `certify_graph()` in `reproduce.py` and can be raised.

Progress is printed per graph and saved to `reproduce_results.json`. (The
decomposition index recorded per graph depends on race timing and is not a
stable fingerprint; the verdicts are.) The final message states the
conclusion. The engine is deterministic — and observed to be *bit-identical*
across Apple clang/ARM64, MSVC/x86-64, and LLVM clang/x86-64 on every matched
instance tested, node counts and all.

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

Original computation carried out 2026-08-03/05 on an M2 Pro (the three
hardest graphs by a ~20 h C slicing campaign whose logs are preserved in
directory `f4_proof_artifacts` of the same project); independently
re-certified end-to-end 2026-08-05 on Windows (Ryzen 9950X3D) under both MSVC
and LLVM clang builds of the kernel. This package is the minimal, cleaned
re-derivation path. During development the engine existed in two independent
implementations (a Python reference and the C version); they produced
node-for-node identical search trees on matched instances. The package ships
the C engine.

## License

MIT — see `LICENSE`.
