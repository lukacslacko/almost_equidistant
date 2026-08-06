# Dimension-6 K7 fixed-clique basis-extension pilot

## Scope

This is a bounded exact pilot on the `no near-saturating clique` part of the
258-graph exact K7 residue.  It tests the restricted rational PSD-atom dual
already implemented in `d6_k7_arbitrary_basis_psd_dual.py`; it is not a new
full rejection layer.

The exact replay finds 53 relevant covers in 50 residue graphs.  Every one
has

```text
|N| = 12,       rank(K) <= U = 7,       omega(G[N]) = 5,
zmask = 0.
```

The production pilot takes the first eight covers in the deterministic order
induced by the exact-union residue manifest, K7 seed enumeration, and cover
enumeration.

## Why one fixed maximum clique is enough

Let `C` be any required clique in `G[N]`.  For a feasible normalized Gram
matrix `K`,

```text
K[C,C] = J + diag(e_i),       e_i > 0.
```

Thus `K[C,C]` is positive definite, so the Gram vectors indexed by `C` are
linearly independent.  If `rank(K)=r`, this independent subset extends using
other columns of `K` to an `r`-element basis `B` containing `C`.  Therefore,
after fixing one required maximum clique `C`, every rank in the permitted
range is represented among

```text
B = C union E,       |C| <= |B| <= U,
```

where `E` ranges over subsets of `N-C`.  It is unnecessary to enumerate
principal bases that do not contain `C`.

Here `|C|=5`, `|N-C|=7`, and `U=7`, so each cover has exactly

```text
C(7,0) + C(7,1) + C(7,2) = 1 + 7 + 21 = 29
```

candidate bases, rather than all subsets of sizes five through seven.

For a proposed basis `B`, write `H=K[B,B]^-1` and let `p_y=K[B,y]` for an
outside vertex.  Schur equality gives

```text
p_y^T H p_z = K_yz,           y != z,
p_y^T H p_y = K_yy > 1.
```

The existing rational Farkas certificates contradict this necessary affine
system with a PSD matrix expressed as a nonnegative sum of small integer
rank-one atoms.  The floating-point LP is only a locator; the stored identity
is checked over `Fraction`.

Two direct checks precede the locator:

- if `p_y=0`, Schur equality says `K_yy=0`, contradicting `K_yy>1`;
- if `p_y=p_z`, Schur equality says
  `K_yy=K_yz` and `K_yz` is zero or one, again contradicting `K_yy>1`.

Both are necessary-condition contradictions and do not constrain an omitted
candidate nonedge to be genuinely non-unit.

## Exact pilot result

The eight covers contribute 232 basis candidates.  With the frozen atom
dictionary consisting of primitive `{−1,0,1}` vectors of support at most
three, the exact outcomes are:

| rank | exact PSD-atom certificates | unresolved |
|---:|---:|---:|
| 5 | 8 | 0 |
| 6 | 16 | 40 |
| 7 | 9 | 159 |
| **total** | **33** | **199** |

There were no zero-column or identical-column hits in this sample.  None of
the eight covers was fully eliminated.  `UNRESOLVED` means only that this
restricted certificate dictionary found no proof.

The run used 11 worker processes on the 12-logical-CPU Mac, with BLAS thread
counts fixed at one.  The locator tasks consumed 52.77 aggregate worker
seconds.  End-to-end wall time, including a fresh exact reconstruction of all
53 target covers from the frozen inputs, was 26.87 seconds.  This is cheap
enough for a later full pass, but the 0/8 cover rejection rate says that an
unchanged rollout is not currently worthwhile.

An additional exploratory run on the first cover enlarged the atom
dictionary from support at most three to all `{−1,0,1}` vectors.  It still
certified only 2 of 29 bases while increasing locator time from about 1.50 to
8.13 seconds.  Failure to locate additional certificates is heuristic only;
the exact certificates already found remain valid.

## Reproduction and hashes

```text
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
python3 run_d6_k7_arbitrary_basis_extension_pilot.py \
  --covers 8 --workers 11 --maximum-atom-support 3

python3 verify_d6_k7_arbitrary_basis_extension_pilot.py \
  --report-sha256 \
    2c2dd18ded8a0daf505700001e64eb31fdcbf987f6292f91e58a3c07c8dc8339

python3 -m unittest -v \
  test_d6_k7_arbitrary_basis_extension_pilot.py \
  test_d6_k7_arbitrary_basis_psd_dual.py
```

Artifacts:

```text
d6_k7_arbitrary_basis_extension_pilot_report.json
  2c2dd18ded8a0daf505700001e64eb31fdcbf987f6292f91e58a3c07c8dc8339
d6_k7_arbitrary_basis_extension_pilot_verification.json
  789f53f235599b9b4ce0c4b8b288f7fb1ed767d9038d340dff861d18cfee2d44
```

The independent verifier replays the complete 53-cover population, checks
that the selected prefix and all 232 fixed-clique extension cores are exact,
and verifies all 33 rational certificates.  It returns `PASS`.  Seven focused
and pre-existing tests pass.

## Conclusion and next direction

The fixed-maximum-clique reduction is a useful permanent reduction: it turns
each no-near cover into only 29 bases.  The current arbitrary-basis atom dual,
however, is too weak at ranks six and seven and should not be rolled out
unchanged.

All 53 covers instead expose a common stronger structure.  Fixing the same
required `K5` leaves a Schur complement of rank at most two.  Exact
off-diagonal rank-two Gram identities (the five-index pentad relations) can
use all seven remaining columns at once and avoid enumerating arbitrary
rank-seven bases.  That is the next bounded algebraic pilot.

## Trust scope

The rejection records trust the hash-pinned exact residue and prior campaign,
the exact K7 support/rank replay, integer bit masks, Python arbitrary-precision
integers and `Fraction`, and the documented normalized-Gram reduction.  The
LP locator, IEEE-754 arithmetic, and failure to locate a certificate are
outside the mathematical trust boundary.
