# K6 arbitrary-subset rank Hall: exact production result

This exact layer strengthens support Hall inside a positive-sign,
non-bipartite connected K6 side span.  The frozen hereditary PSD--Z rule can
select only the empty or full set from such a span.  The new rule assigns a
sound rank lower bound to every selected subset and still chooses at most one
subset from each original connected span.

On the 831 graphs surviving the committed hereditary PSD--Z layer, the frozen
production kernel rejects exactly 9 and leaves 822.  A structurally independent
checker recomputed all 831 decisions, replayed all 511 exhaustive certificate
rows, and passed the known realizable 18-point control at all 32 required K6
seeds.  Production source was frozen at commit
`4e1d50492ad76867b67d748792f21ca5c0a6a94f`; the immutable result archive was
committed and pushed at `a21db74bacdf4c9c5c841ac137af326b331d8137`.

## Arbitrary-subset rank theorem

Fix a required unit `K6`, an eligible actual zero-factor set `Z0`, a
bipartite component of the Lorentz graph, and one generic sign orientation.
Let `F` be one connected required-unit graph in a side span.  Its Euclidean
Gram matrix has the frozen form

```text
Q = I + kappa D (I + Adj(F)) D,                         (1)
```

where `D` is nonsingular and `Q` is positive semidefinite.  Candidate
nonedges impose no hidden non-unit constraint: in this derivation they give
the zero off-diagonal entries of the Gram pattern, while every required edge
of `F` gives a nonzero entry.

The previous hereditary rule already handles negative `kappa`, and positive
`kappa` when `F` is bipartite.  Consider the remaining case:

```text
kappa > 0 and F is connected and non-bipartite.
```

For an arbitrary selected subset `S` of `V(F)`, decompose the induced graph

```text
F[S] = H_1 disjoint-union ... disjoint-union H_t
```

into connected components.  The corresponding principal Gram blocks are
mutually orthogonal because there are no `F[S]` edges between components.
For a component `H` of order `h`, production uses

```text
r(H) = max(
    h - Z(H),
    h - n_-(I + Adj(H)),
    h - 1                         if H is bipartite
).
```

Here `Z(H)` is the ordinary zero-forcing number and `n_-` is exact negative
inertia.  The first term is the exact-nonzero-pattern zero-forcing bound.  The
second is the already-proved positive-sign inertia bound applied to the
principal submatrix.  If `H` is bipartite, signature switching turns its
positive edge entries into an irreducible PSD Z-matrix, giving nullity at
most one.  Therefore

```text
rank Q[S,S] >= sum_j r(H_j).                            (2)
```

For `S=V(F)`, formula (2) is asserted to equal the frozen full-span bound;
the implementation checks this equality on every evaluated group.

The discovery probe also tested the optional extension lemma below, which can
raise the last term from `h-1` to `h` in some cases.  It produced no additional
graph rejection, so the lemma and every corresponding code path are omitted
from the production theorem kernel.

## Probe-only one-vertex bipartite extension lemma

This sound but zero-marginal rule is retained only to document the discovery
experiment.  It is not an assumption of the production certificates.  Suppose
a connected bipartite component `H` of `F[S]` has a vertex

```text
v in V(F) - S
```

such that the induced graph `F[V(H) union {v}]` is connected and bipartite.
Restrict (1) to this extension and signature-switch one color class.  After
diagonal congruence it is an irreducible PSD Z-matrix.  The Gram block on `H`
is a proper principal submatrix of that matrix, hence is positive definite by
the Perron--Frobenius proper-principal-submatrix lemma.  Thus

```text
rank Q[V(H),V(H)] = |V(H)|.                            (3)
```

The condition is checked on the induced extension, not merely on the
existence of a boundary edge.  A vertex adjacent to both colors may create an
odd cycle and is not accepted.  A vertex in another component of `F[S]`
cannot connect to `H`, so requiring `v` outside `S` is the correct quantifier.

The synthetic control uses a positive-sign `C5` and a selected induced `P3`.
Without (3), zero forcing, inertia and the bipartite PSD--Z rule give rank 2;
the missing cycle vertex makes an induced `P4`, so (3) raises the exact lower
bound to 3.

## Hall quantifier and no double counting

For every original connected side span, the production kernel chooses exactly
one of:

* the empty set;
* one selected subset `S`, carrying the bound (2), in a positive-sign
  non-bipartite span;
* one hereditary proper subset, carrying rank `|S|`, in a PSD--Z-applicable
  span;
* the full span with its frozen rank lower bound.

It then combines these choices with analogous choices from the opposite
side and with any subset of the singleton `Z0` spans.  Only spans known to be
mutually orthogonal have their ranks added.  Two different subsets of the
same original connected span are alternatives; their ranks are never added,
even when the subsets are disjoint.

For every such choice the necessary inequality is

```text
sum of selected span-rank lower bounds
    <= size of the union of selected allowed coordinate masks.            (H)
```

Allowed masks are upper bounds.  Actual coordinates may be zero, so replacing
actual supports with allowed masks only enlarges the right side of `(H)`.
Failure remains sound without requiring any candidate nonedge to be non-unit.

The optimized checker keeps, for each of the 64 coordinate unions, only the
largest attainable rank after each original span.  A lower rank with the same
union can never produce a later Hall violation that the larger rank misses.
On eight actual input graphs, a raw Cartesian-product checker independently
enumerated every uncompressed one-choice-per-span combination for 4,272
orientation systems and agreed with the optimized decision in every case.

## Pinned input and discovery commands

The input is the ordered 831-graph complement of the committed hereditary
PSD--Z report, with stable index hash

```text
2bcdad095c6bd3038a4bb1d117f2faadfa98113d814633351ec9e79ff553a44b.
```

The extension-enabled command was

```text
/Users/lukacs/claude/opengauss/venv/bin/python3 \
  probe_d6_k6_arbitrary_subset_hall.py \
  --workers 7 --cross-check-graphs 8
```

The exact comparison without the optional extension lemma was

```text
/Users/lukacs/claude/opengauss/venv/bin/python3 \
  probe_d6_k6_arbitrary_subset_hall.py \
  --workers 7 --cross-check-graphs 8 --no-extension \
  --output /private/tmp/d6_k6_arbitrary_subset_hall_no_extension.json
```

Both runs use the same final pilot source, whose SHA-256 is

```text
7570fea6689b2f2ecaedb991415a5a2cc7f9b52c61b81175d4c793ce797eac0e.
```

## Exact pilot result

| quantity | extension enabled | extension disabled |
|---|---:|---:|
| input / rejected / surviving | 831 / 9 / 822 | 831 / 9 / 822 |
| K6 seeds checked / impossible | 25,944 / 9 | 25,944 / 9 |
| `Z0` considered / matchable | 27,626 / 27,626 | 27,626 / 27,626 |
| passing / failing `Z0` | 25,935 / 1,691 | 25,935 / 1,691 |
| bipartite Lorentz components checked / failed | 185,376 / 1,691 | 185,376 / 1,691 |
| frozen-passing components newly failed | 39 | 39 |
| generic orientations / newly failed | 370,752 / 41 | 370,752 / 41 |
| positive-nonbipartite original span groups | 15,047 | 15,047 |
| new proper-subset options | 197,842 | 197,842 |
| one-vertex-extension-upgraded options | 7,183 | 0 |
| raw / coordinate-union-compressed options | 1,457,948 / 1,160,534 | same |
| dominance-DP transitions | 2,176,479 | 2,176,519 |
| measured 7-worker wall seconds | 1.808 | 1.753 |

The extension lemma changes 7,183 option rank bounds, but it adds no graph
rejection and leaves all reported aggregate pass/fail counters unchanged on
this corpus.  All nine graph rejections already follow from arbitrary-subset
zero-forcing, positive-sign inertia and induced-bipartite PSD--Z bounds
without (3).

The rejected indices are

```text
650158, 752132, 1226632, 2554326, 2716230, 2862902, 3178386,
3194938, 3494044.
```

Their stable list hash is

```text
960b3eb759456dd653c2e4d50517f693646fb968b957a787fdd3c2d9addfc448.
```

The nine first impossible seeds contain 511 exhaustive `Z0` rows: four
certificates have 32 rows, four have 64, and one has 127.  As a representative
witness, graph `650158` has impossible seed `[0,2,5,6,14,16]`.  At `Z0=[]`,
the newly failing positive orientation selects rank 3 on the arbitrary subset
`[8,9,17]` of a positive non-bipartite span and rank 2 on hereditary subset
`[7,18]` of the opposite span.  Their allowed-coordinate union has size 4,
so `(H)` fails as `5 > 4`.  The other generic orientation and the lightlike
alternative also fail.

The extension-enabled report is

```text
d6_k6_arbitrary_subset_hall_pilot_report.json
SHA-256 292ca0cac2e19f668bc831e782cb4208694d5264f4b66724fb80b39332d576cc
bytes 3,281,226
```

The immutable no-extension comparison remains in `/private/tmp` with
SHA-256
`d38b87924d7938b961057ab27e710cdf8591e0a3d431cf9e62e1dd52097fd11b`.

## Controls, trust, and nonclaims

The embedded controls establish:

* the `C5/P3` extension rank changes from 2 to 3 exactly;
* frozen full-only and arbitrary-subset-without-extension synthetic Hall both
  pass, while the extension-enhanced synthetic Hall fails;
* optimized and raw Cartesian Hall agree on the synthetic instance and 4,272
  actual orientation systems;
* the known realizable 18-point graph passes all 32 required K6 seeds with
  the extension enabled and disabled;
* every full-subset rank reconstructed by the new kernel equals its frozen
  full-span bound.

All rank, graph, mask and Hall decisions use exact integers.  Floating point
appears only in elapsed-time fields.  The pilot trusts the pinned parent
corpus, the frozen K6 Lorentz-coordinate identity, exact nonzero edge pattern,
ordinary zero forcing, sign-specific inertia, PSD--Z/Perron--Frobenius lemmas,
and Python implementation correctness.

The discovery report by itself is not a theorem certificate.  The production
archive described below supplies the source-bound exact credit for the same
nine rejections.  The 822 survivors are filter non-rejections, not
realizations.  This layer does not settle the K6 branch, the separate K7
branch, or dimension six.

## Source-bound production certification

Production, the independent verifier, and the seven-test control suite were
frozen before the full run in commit

```text
4e1d50492ad76867b67d748792f21ca5c0a6a94f
```

with source hashes

```text
d6_k6_arbitrary_subset_hall.py
  cadc0048fcc7a77299bc2950e2eca57bfed87d54ed6ffd59886d155b2c9a88ff
verify_d6_k6_arbitrary_subset_hall.py
  4ad4f886d67b3d1ebf149277050a639cda243f657ef02b553e6a3fa241f2853e
test_d6_k6_arbitrary_subset_hall.py
  e92896a2d426d9f98a03e210990d76ca314ddfa0f15a56a765f516bfb719cde1
```

`d6_k6_arbitrary_subset_hall.py` contains only the nine-hit theorem kernel:

```text
componentwise max(
    ordinary zero-forcing,
    positive-sign inertia,
    bipartite PSD--Z |H|-1
)
```

on every connected component `H` of an arbitrary induced subset.  It contains
no one-vertex-extension function or rank rule.  It writes a deterministic
report, compact exhaustive first-impossible-seed certificate archive, and an
atomic ordered-prefix checkpoint bound to the committed production-source
hash and the pinned 831-index hash.

`verify_d6_k6_arbitrary_subset_hall.py` imports neither production nor the
discovery probe.  It reconstructs arbitrary induced components and ranks with
the frozen independent SymPy/Sturm inertia and simultaneous zero-forcing
implementations.  Its Hall checker is structurally different from production:
it enumerates every subset alternative inside each original span and all 64
coordinate containers, maximizing one alternative per span.  A violation
exists exactly when one container's summed maxima exceed its size.  The
checker also validates exhaustive `Z0` coverage and every archived witness,
including that no original span is selected twice.  It supports multiple
graph workers for the final replay.

`test_d6_k6_arbitrary_subset_hall.py` has seven source-bound controls:

* production and independent `C5` subset-rank/Hall controls;
* agreement between the independent 64-container criterion and raw Cartesian
  enumeration on a synthetic two-span system;
* the one-subset-per-original-span no-double-count control;
* the pinned 831-index boundary;
* production and independent agreement on fixed pilot rejection `650158`;
* the known realizable 18-point production control, passing all 32 K6 seeds;
* verifier import independence, atomic checkpoint round trip, and an AST
  check that production defines no one-vertex-extension kernel.

The committed-boundary production and verification commands were

```text
/Users/lukacs/claude/opengauss/venv/bin/python3 \
  d6_k6_arbitrary_subset_hall.py \
  --workers 11 --checkpoint-every 50

/Users/lukacs/claude/opengauss/venv/bin/python3 \
  verify_d6_k6_arbitrary_subset_hall.py --workers 11
```

Production completed in 4.675 measured kernel wall seconds and reported these
exact aggregate counts:

| quantity | count |
|---|---:|
| input / rejected / surviving | 831 / 9 / 822 |
| K6 seeds checked / impossible | 25,944 / 9 |
| `Z0` considered / matchable | 27,626 / 27,626 |
| passing / failing `Z0` | 25,935 / 1,691 |
| bipartite Lorentz components checked / failed | 185,376 / 1,691 |
| generic orientations checked | 370,752 |
| positive-nonbipartite original span groups | 15,047 |
| arbitrary proper-subset options | 197,842 |
| raw / coordinate-union-compressed options | 1,457,948 / 1,160,534 |
| dominance-DP transitions | 2,176,519 |
| subset-rank cache entries / hits | 70,380 / 142,509 |

The independent verifier completed in 10.894 wall seconds.  It recomputed all
831 graphs, obtained the same ordered nine-index rejection set, independently
replayed all 511 archived `Z0` rows, and recorded each of the 32 positive-control
K6 seeds individually as passing.  Its inertia calculation uses independent
SymPy/Sturm code and its Hall calculation enumerates subset alternatives via
all 64 coordinate containers rather than importing the production DP.

The post-run source-bound command

```text
/Users/lukacs/claude/opengauss/venv/bin/python3 -m unittest -v \
  test_d6_k6_arbitrary_subset_hall.py
```

passed all seven tests in 2.126 seconds.  Byte-code compilation of production,
verifier and tests also passed.  No production source was altered after the
source-bound run.

The immutable result artifacts were committed and pushed at
`a21db74bacdf4c9c5c841ac137af326b331d8137`:

| artifact | bytes | SHA-256 |
|---|---:|---|
| `d6_k6_arbitrary_subset_hall_report.json` | 614,483 | `bcc4baecf28a95443a5b13e87a796cbba0a97cd42fa76573119bfd81ef3735fc` |
| `d6_k6_arbitrary_subset_hall_certificates.json` | 2,098,623 | `b7fb04f769ddaf5d36630372f5134f2e4a32e572f620b26059667e9cf65f4a23` |
| `d6_k6_arbitrary_subset_hall_checkpoint.json` | 2,871,213 | `43c1b65da336bc0db006313b97630867a660b944b06c814433ee3a8cdc8297e6` |
| `d6_k6_arbitrary_subset_hall_verification.json` | 8,552 | `e6042edfff3683043903f9be544322379e85673207e35ee87349637ff4c73058` |
| `d6_k6_arbitrary_subset_hall_result_manifest.json` | 2,664 | `a245fb4ae68d5c93b75e707799346e360358e78f17cc97c7b8f4ecb4705674bd` |

The exact nine-index list and its stable hash remain

```text
650158, 752132, 1226632, 2554326, 2716230, 2862902, 3178386,
3194938, 3494044

960b3eb759456dd653c2e4d50517f693646fb968b957a787fdd3c2d9addfc448
```

All theorem decisions are exact integer, graph, mask, zero-forcing, inertia,
and Hall calculations; floating point is used only for elapsed-time metadata.
The trust boundary is the pinned parent corpus and K6 Lorentz derivation, the
ordinary zero-forcing and PSD--Z/inertia lemmas, Python/SymPy exact arithmetic,
and correctness of the two independent implementations.  Candidate nonedges
remain unconstrained, allowed defect coordinates may be zero, exactly one
subset is selected per original span, and unresolved cases are never counted
as rejected.  In particular, the optional one-vertex extension lemma is not
used anywhere in this production result.
