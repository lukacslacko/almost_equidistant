# K6 PSD Z-matrix rank refinement

This exact layer rejects 116 of the 977 graphs surviving the frozen K6
side-rank fusion at commit `fc51874`, leaving 861.  Since the parent layer
already rejected 13 of the manifest-pinned 990 K6-only graphs, the cumulative
K6 reduction is now 129 rejected and 861 surviving.

The new ingredient is a finite positive-semidefinite matrix theorem.  On a
connected side-graph block, one Lorentz sign always gives an irreducible PSD
Z-matrix; the other sign does so whenever the block is bipartite.  Such a
matrix has nullity at most one.  This rank bound is fused component by
component with the existing ordinary-zero-forcing and exact-inertia bounds
before the orthogonal block-support inequalities are tested.

An independent checker recomputed all 977 graph decisions and all 4,987
archived `Z0` rows and returned `PASS`.  The known realizable 18-point
configuration passes every one of its 32 required K6 seeds.

## Exact matrix theorem

Fix a required unit `K6`, a possible actual zero-Lorentz-factor set `Z0`, and
a connected bipartite component

```text
C = A union B
```

of the disjoint-defect Lorentz graph `L-Z0`.  In the generic Lorentz case,
all factors on one side are nonzero scalar multiples of one direction.  For
one side, let `F` be its induced required-unit graph.  The coordinate
derivation in `d6_k6_normal_coordinates.md` gives the actual Euclidean Gram
matrix

```text
Q = Gram(u_x : x in F) = I + kappa D (I+Adj(F)) D,       (1)
```

where `D` is a nonsingular real diagonal matrix and `kappa` is nonzero.
Candidate nonedges cause no hidden prescribed non-unit distance here: by
`alpha(G)<=2`, their allowed defect masks are disjoint, hence their defect
inner products are zero.  Required edges give the nonzero off-diagonal terms
in (1), because `kappa` and every diagonal entry of `D` are nonzero.

Let `F_i` be one connected component of `F` and restrict (1) to that block.
Because `Q_i` is a Euclidean Gram matrix, it is positive semidefinite.

If `kappa<0`, positive scaling and congruence by `D_i^{-1}` give

```text
(-kappa)^(-1) D_i^(-1) Q_i D_i^(-1)
    = diag(t) - Adj(F_i).                                (2)
```

The diagonal absorbs the diagonal part of `I+Adj(F_i)`.  Matrix (2) is a PSD
Z-matrix whose off-diagonal graph is connected.  Its nullity is at most one.
For completeness, choose `rho` larger than every diagonal entry and write

```text
B = rho I - (diag(t)-Adj(F_i)).
```

Then `B` is nonnegative and irreducible.  If (2) is singular, `rho` is the
Perron root of `B`; Perron--Frobenius makes that eigenspace one-dimensional.
If (2) is nonsingular, its nullity is zero.  Congruence preserves nullity, so

```text
rank Q_i >= |F_i|-1.                                    (PF-negative)
```

If `kappa>0`, the corresponding congruent PSD matrix has positive edge
entries,

```text
diag(t) + Adj(F_i).
```

When `F_i` is bipartite, switching signs on one color class changes every
edge entry to `-1` and leaves the diagonal fixed.  The switched matrix is
again an irreducible PSD Z-matrix, so

```text
rank Q_i >= |F_i|-1.                                    (PF-positive-bip)
```

No Z-matrix conclusion is imposed on a positive non-bipartite block.

## Componentwise fusion

Write

```text
H_i = I+Adj(F_i),   inertia(H_i)=(p_i,n_i,z_i),
Z(F_i) = ordinary zero-forcing number of F_i.
```

Every actual block Gram matrix also has the exact nonzero graph `F_i`, so
ordinary zero forcing gives rank at least `|F_i|-Z(F_i)`.  The prior exact
inertia argument gives rank at least `|F_i|-n_i` for positive `kappa` and
`|F_i|-p_i` for negative `kappa`.  Consequently the exact component bounds
used here are

```text
r_i(positive) = max(
    |F_i|-Z(F_i),
    |F_i|-n_i,
    |F_i|-1 if F_i is bipartite else 0
),

r_i(negative) = max(
    |F_i|-Z(F_i),
    |F_i|-p_i,
    |F_i|-1
).                                                       (3)
```

Different connected components of `F` are mutually Gram-orthogonal, so the
side-span rank lower bound is the sum of (3), not merely the maximum of
whole-side totals.  This componentwise formulation was checked against the
earlier profiling kernel and gives the same 116 graph rejections.

The two generic Lorentz orientations are still

```text
A positive, B negative;
A negative, B positive.
```

For either orientation, let the resulting side bounds be `r_A,r_B`.  The
actual spans `U_A`, `U_B`, and the one-dimensional spans `U_z` for `z in Z0`
are mutually orthogonal.  Each lies inside the coordinate subspace given by
the union of its allowed defect masks.  Therefore every nonempty subset `J`
of these at most eight blocks must obey

```text
sum_(block in J) rank_lower(block)
    <= |union_(block in J) allowed_coordinates(block)|.  (4)
```

The implementation checks every subset in (4) for both orientations.

### Lightlike alternative

The generic signs and (3) are not applied when the Lorentz direction is
lightlike.  In that case all vectors indexed by `C union Z0` are orthonormal,
so the separate necessary condition is

```text
|C|+|Z0| <= 6
```

together with a matching into the allowed coordinate masks.  A component is
rejected only when both generic orientations and this exact lightlike
alternative fail.

Allowed masks remain one-sided upper bounds: an allowed coordinate may be
zero, and a candidate nonedge may have unit distance.  Replacing actual
supports by allowed masks only weakens (4), so failure with the allowed masks
is a sound obstruction.

## Quantifiers and same-Z0 coupling

For one K6 seed the production evaluator enumerates every eligible `Z0` of
size at most six whose allowed masks admit a matching.  For that same `Z0` it
checks every bipartite Lorentz component with (3)--(4), then applies the
frozen exact non-bipartite light-ray/support system.  A seed is impossible
only if every `Z0` fails.  One impossible required K6 seed rejects its graph.

The production boundary is the ordered 977-index survivor list in
`d6_k6_fused_side_rank_report.json`, with stable hash

```text
27e435506ecaecb936acaef63e5873c3538771be17d1a52f12df82df8dedd947.
```

No graph outside that pinned boundary is silently assumed handled by this
incremental report.

## Exact result

Run:

```text
python3 d6_k6_psd_zmatrix.py --workers 11 \
  --output d6_k6_psd_zmatrix_report.json \
  --certificates d6_k6_psd_zmatrix_certificates.json
```

The source-bound run used 11 workers and completed in 2.48 seconds of
measured production time on macOS 14.5 arm64 with Python 3.11.15.

| quantity | exact value |
|---|---:|
| input parent survivors | 977 |
| newly rejected / surviving | 116 / 861 |
| K6 seeds checked / impossible | 28,597 / 116 |
| `Z0` choices considered / matchable | 35,844 / 35,844 |
| PSD-Z passing / failing `Z0` choices | 28,733 / 7,111 |
| PSD-Z-passing choices failing same-`Z0` non-bipartite support | 252 |
| bipartite Lorentz components checked / failed | 221,715 / 7,111 |
| prior-fused-passing components newly failed | 695 |
| generic orientations checked / newly failed | 443,430 / 1,959 |
| connected side-graph blocks checked | 541,380 |
| blocks where the PSD Z-matrix rule applies | 517,485 |
| strict component-rank improvements from the PSD Z rule | 14,579 |
| prior / new block subsets checked | 1,918,388 / 1,870,804 |
| raw / matchable non-bipartite colorings | 30,317 / 28,481 |
| joint support searches / DFS nodes | 28,481 / 90,731 |

Every one of the 116 rejected graphs has a first impossible seed for which
the PSD Z-matrix bipartite system itself fails for every `Z0`; none requires
the non-bipartite system to complete its rejection.  The archived first-seed
certificates contain 4,987 exhaustive rows: 3,704 are
`nonbipartite_only` and 1,283 are `neither`.  No row passes the PSD Z-matrix
system.

The complete rejection list is stored in the report and archive.  Its stable
SHA-256 is

```text
adb32fa4c00f6209ce7e246fd49337a6c030b20ddd49193765cf3451afb4c36c.
```

## Strict corpus witness

Graph `3278`, seed

```text
[0,1,3,8,12,16],
```

and `Z0` empty have the bipartite Lorentz component

```text
[2,5,6,9,10,11,13,14,17,18]
```

with sides

```text
A=[2,5,9,18],
B=[6,10,11,13,14,17].
```

In the `A_positive` orientation the preceding fused side ranks are `(3,3)`
and pass block support.  The negative `B` block is connected on six vertices,
so the PSD Z-matrix theorem raises its rank to five.  Its allowed coordinate
union has size four, giving the immediate contradiction `5<=4`.

In the `A_negative` orientation the preceding ranks `(2,4)` also pass.  The
negative connected `A` block raises the pair to `(3,4)`, which requires seven
dimensions in a six-coordinate union.  The lightlike alternative would need
ten orthonormal vectors and fails.  Thus both old generic orientations pass
while both strengthened orientations fail: this is a strict witness for the
new theorem, not merely a reordering of an old rejection.

## Independent verification and controls

Run:

```text
python3 verify_d6_k6_psd_zmatrix.py \
  --output d6_k6_psd_zmatrix_verification.json
python3 -m unittest -v test_d6_k6_psd_zmatrix.py
python3 -m py_compile \
  d6_k6_psd_zmatrix.py verify_d6_k6_psd_zmatrix.py \
  test_d6_k6_psd_zmatrix.py
```

The independent verifier imports no new production evaluator or kernel.  It
uses:

- freshly reconstructed side-graph components, bipartiteness, PF
  applicability, sidewise sums, and block subsets;
- exact SymPy characteristic polynomials and Sturm root counts for inertia,
  instead of production rational congruence;
- simultaneous ordinary-zero-forcing closure, instead of production's
  sequential solver;
- direct Hall-subset checks and the frozen independent actual-support CSP,
  instead of production matching/search code.

It recomputed all 977 graph decisions, matched eight core fields and the first
impossible seed for every graph, and matched the full JSON-normalized content
of every one of the 4,987 certificate rows.  The serial verification took
12.40 seconds and returned `PASS`.

The small exact controls are:

| connected graph | positive-side lower | negative-side lower |
|---|---:|---:|
| `K1` | 1 | 0 |
| `K2` | 2 | 1 |
| `P3` | 2 | 2 |
| `C5` | 3 | 4 |

`C5` deliberately checks the sign asymmetry: positive `C5` is non-bipartite,
so it retains only ordinary-ZF/inertia rank three, while negative `C5` gets
the PSD Z-matrix rank four.  The realizable 18-point control passes all 32 K6
seeds in both production and independent implementations.  Seven unit tests
cover these kernels, the strict graph-3278 witness, frozen totals and hashes,
all archived rows, verifier import independence, and the positive control.

## Source-bound hashes

```text
d6_k6_psd_zmatrix.py
  24c2b3abfce5074560183d17f77778f96298bef6e38ece5bc937ccb2a4cba736
d6_k6_psd_zmatrix_report.json
  abb920fbd9eb848502cca37e32c2aa86bf220ab1773e6006ea1fbee4bc2f4e30
d6_k6_psd_zmatrix_certificates.json
  62fc7170d97056941bcaf553c45cfa1cabdf20117ffc69fbd8429b3764ec999c
verify_d6_k6_psd_zmatrix.py
  ca4ca62f75d02846f12e59e16eb9eab4473fe5438fc496e72a25a39db77f9042
d6_k6_psd_zmatrix_verification.json
  8261b58e788f120e3d1eee6dc481a4989e3eb75f9f03babbd552474b4a5302eb
```
