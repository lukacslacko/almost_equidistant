# Exact one-free-neighbour/two-free-center K7 star

## Pre-production theorem statement

Fix a regular unit `K7`, an eligible zero-factor cover `Z`, a labeled family
of actual supports on `Z`, and the propagated masks on its nonzero-factor
complement `N`. Suppose a vertex `x in N` has exactly two coordinates not
pinned by nonbipartite singleton-intersection components, and at least two
required neighbours in `N` have exactly one unpinned coordinate each. For
every correlated assignment of the pinned-component signs, each neighbour
edge gives an affine line in the two free coordinates of `x`. If all those
lines have no common real intersection with the exact diagonal hyperbola of
`x` for every sign assignment, the propagated branch is impossible.

The dry-run corpus replay finds this obstruction in four of the 88 currently
passing families. Under the complete cover/seed quantifier it provisionally
rejects exactly

```text
2592657  3785980  3888410
```

from the frozen v6 19-graph K7 class. These receive theorem-level credit only
after the committed producer runs and the independent committed checker
returns `PASS`.

## Normalization and support semantics

For an outside point define

```text
u_i(x) = ||x-q_i||^2-1,
t_x    = 1+sum_i u_i(x).
```

The zero-factor cover fixes `Z={x:t_x=0}` for the branch. Every point used by
this star lies in the complementary set `N`, so `t_x != 0` and

```text
w_i(x) = sqrt(7) u_i(x)/t_x
```

is legitimate at the center and every neighbour. An allowed defect
coordinate may be zero without changing `t_x != 0`; optional zeros therefore
do not invalidate the normalization.

The propagated masks are supersets of actual support. A coordinate is
deleted only when singleton-overlap orthogonality proves it zero. Every
remaining entry stays optional. A coordinate is treated as a nonzero sign
only after an odd singleton-intersection component forces reciprocal
products around an odd cycle and hence a common value `+1` or `-1`.
Candidate nonedges never enter an equation and may also have unit distance.

## Exact two-free and one-free equations

The normalized diagonal identity is

```text
sum_i w_i^2 - (sum_i w_i)^2
  + 2 sqrt(7) sum_i w_i - 8 = 0.                 (1)
```

If a mask of size `k` has two unpinned coordinates `x,y` and its `k-2`
pinned signs sum to `P`, equation (1) factors as

```text
(x-(sqrt(7)-P))(y-(sqrt(7)-P))
  = (k+4+P^2-2 sqrt(7)P)/2.                     (2)
```

Thus the center lies on `(x-A)(y-A)=R` over `Q(sqrt(7))`. The right side is
nonzero for every admissible integer `k,P`; the tests exhaust these types.

For a one-free neighbour of mask size `l` and pinned-sign sum `Q`, its last
coordinate is uniquely

```text
F(l,Q) = [Q(Q^2-l-5)+(9-l-Q^2)sqrt(7)]/[2(7-Q^2)].   (3)
```

The denominator cannot vanish because `Q` is an integer. Once correlated
component signs are assigned, every coordinate of that neighbour is known
in `Q(sqrt(7))`. Its required edge with the center becomes

```text
L x + M y = D.                                    (4)
```

Pinned shared coordinates contribute to `D`; a center-free coordinate
present at the neighbour contributes to `L` or `M`. No absent coordinate or
candidate nonedge is constrained.

## Complete real-feasibility classification

For each correlated sign assignment the kernel reduces all lines (4) over
the exact quadratic field.

- Rank zero accepts `0=0` and rejects `0=D` only for nonzero `D`.
- Rank one checks coincident and parallel lines exactly. A genuine line is
  substituted into (2). With both coefficients nonzero, feasibility is the
  nonnegativity of a quadratic discriminant in the positive real embedding;
  zero remains feasible as a tangent. With one coefficient zero, the fixed
  coordinate is compared with the hyperbola asymptote, including the `R=0`
  branch defensively.
- Rank two solves the unique point by exact field division, checks every
  remaining line, and substitutes the point into (2) coefficientwise.

The sign of `a+b sqrt(7)` is decided using rational comparisons of `a^2` and
`7b^2`; floating point never enters. Exact zero coordinates are accepted.
The explicit positive control

```text
(x-1)(y-1)=-1,  x=0
```

is retained with witness `(x,y)=(0,2)`.

## Complete dry-run profile

The producer replayed the frozen 19-graph input and the independent checker
separately regenerated the older exact quantifier:

```text
K7 graphs                                           19
K7 seeds                                            55
eligible zero-factor covers                     19,932
inherited-current covers                           108
labeled zero-factor support families             1,589
families passing all older exact layers              88
families rejected by the new star                     4
newly rejected graphs                                 3
surviving K7 graphs                                  16
```

The ordered provisional 16-list is

```text
316173 2581209 2593240 3595554 3648882 3729907 3935560
3936177 3936310 3936435 3945490 3945555 3945557 3945564
3947605 3949382
```

with stable JSON hash

```text
b0235ce7e95918197f0f3b4a57d26c8c83632ff51fd678a989f50197591eee8e
```

The star layer is conditioned only on the already verified v6 boundary and
the already verified one-free exact layer. It neither reads nor depends on
the cap-500,000 interval campaign.

The wider exploratory low-free profile also measured the null cohorts before
promoting the star rule. Across the same 88 families there are 275 required
one-free/one-free edge occurrences and 241 one-free/two-free occurrences.
Testing each edge alone rejects no family, and conjoining all exactly
determined one-free/one-free edges rejects no family. There are 76
two-free/two-free edge occurrences, of which 39 have exactly one shared free
coordinate and no free/pinned cross term. They contain only one alternating
cycle (a triangle); its exact Mobius fixed-point discriminant rejects no sign
assignment. The positive marginal first appears when several one-free edge
lines meet at the same two-free center.

## Frozen inputs

```text
c513bcc40c037af23340dce3cc91d593a101dd1d073f93feab82abd556efadc9
  d6_current_residue_manifest_v6.json
3d11cbd20cb59b42a8c8fd5252fb3ef3dea649576982ea8017622fad3279a7bc
  d6_current_residue_manifest_v6_verification.json
181f63015d373fb69f37f4f390cdd8ad50bc32b5a2d97a36fd0bac876d0877a2
  d6_k7_one_free_conjunction_report.json
ebaacfe663728eff427c9ba2280fc4261d26dfe10fad21e1a00982b9d31ac062
  d6_k7_one_free_conjunction_verification.json
af25b842ca9eaa6e9721f51d0400be2436d930b72cc33c72bc366f0ac5e41290
  ordered v6 K7 input indices
```

## Production commands after source commit

Run the producer only from committed `codex/dimension6` sources with no
modified tracked files:

```text
python3 build_d6_k7_one_two_star_increment.py \
  --workers 9 \
  --output d6_k7_one_two_star_increment_report.json
```

Hash the immutable report and run the independent checker:

```text
python3 verify_d6_k7_one_two_star_increment.py \
  --report d6_k7_one_two_star_increment_report.json \
  --report-sha256 <REPORT_SHA256> \
  --workers 9 \
  --output d6_k7_one_two_star_increment_verification.json
```

The checker imports neither the builder nor `d6_k7_one_two_star.py`. It
independently regenerates all 55 K7 seeds, 19,932 eligible covers, inherited
cover gates, 1,589 labeled families, older exact family failures, all 88 star
inputs, component correlations, sign assignments, affine ranks, and real
feasibility decisions. It also enforces committed source blobs and a clean
tracked launch boundary.

Run the focused controls with

```text
python3 -m unittest -v \
  test_d6_k7_one_two_star.py \
  test_d6_k7_one_two_star_increment.py
```

## Pre-production source hashes

```text
d4cd813b66a0c57dd106bb9987d1417366152f4cd9162cfd9b0486e3fd0d0b5b
  d6_k7_one_two_star.py
3dc6003176903d457e4c2fd697eec9d1a4e036b052afd1359cb8c800b182a41a
  build_d6_k7_one_two_star_increment.py
eafd4f69eb4e3938bbd3ce3c1c2c469ddc38e96d2338ed202398bad5377f1b96
  verify_d6_k7_one_two_star_increment.py
```
