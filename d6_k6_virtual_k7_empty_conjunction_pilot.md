# Virtual-K7 / K6 empty-Hall conjunction: bounded pilot

## Exact finite conjunction

For one required K6 seed `Q`, let `B_Q` be the outside vertices whose
single-apex augmentation survived the exact generic K7 pilot.  The branch
theorem in `d6_k6_virtual_k7_pilot.md` proves that every actual empty K6
defect must lie in `B_Q`.

The current empty-support Hall layer decomposes each bipartite Lorentz
component into two generic sign orientations and a lightlike alternative.  A
generic orientation may make at most one negative-side singleton defect
empty.  The base implementation tests the all-nonempty choice first and then
each possible singleton, retaining the minimum number of empties required by
that component.

This pilot reuses those exact component/Hall APIs but retains the identity of
every passing empty singleton.  For every matchable common `Z0`, it performs a
finite DP over Lorentz components whose state is the sorted tuple of actually
chosen empty vertices.  It keeps exactly the following relaxed states:

```text
(),
(x)       with x in B_Q,
(x,y)     with x,y in B_Q and the double-apex branch still possible.
```

The empty tuple is always offered to every generic orientation and by every
passing lightlike alternative, so the all-nonempty branch is never lost.  A
state of size three is discarded.  Non-bipartite Lorentz components remain
ignored conservatively, exactly as in the current base layer.

This is sound in the necessary direction.  Restricting an empty singleton to
`B_Q` uses an independently exact impossibility statement about that branch;
it does not require a surviving vertex to be empty.  A surviving DP state is
only a necessary relaxation.  A K6 seed is impossible only after every common
`Z0`, component orientation, all-nonempty choice, allowed singleton, and
allowed pair has failed.

## Stronger double-apex branch

The single-apex pilot left 38 candidate-nonedge pairs across six of its 640
K6 seeds.  For each pair `{x,y}`, this conjunction simultaneously adds all
twelve branch-forced pairs

```text
xq and yq,       q in Q,
```

while leaving `xy` a candidate nonedge.  It then reruns the complete generic
graph/rank, labeled-support, propagated-link, and sparse-value quantifier at
both promoted seeds `Q union {x}` and `Q union {y}`.

In an actual two-empty branch the two points are the opposite apices.  Relative
to either promoted K7, the other point is unit to the six old K6 vertices and
has its only allowed defect at the new apex.  Distinctness eliminates the
`-1` seed-collision root, leaving the exact singleton value `4/3`.  Thus the
double augmentation encodes substantially more than intersecting two
single-apex relaxations.  Failure at either promoted K7 rules out the pair;
passing both remains only a necessary condition.

All 38 candidate pairs fail this stronger test.  Equivalently, all 76
promoted-seed checks fail at the generic graph/rank stage:

```text
72 have no eligible zero-factor cover;
 2 have eligible covers only at directly impossible sizes 4--7;
 2 have a small eligible cover but no enhanced-rank passing cover.
```

Therefore this sample's empty-Hall DP has no surviving size-two state.  This
is exact pair-branch information, not a graph rejection.

## Same-20 result

The input is exactly the same deterministic 20-graph selection used by the
single-apex pilot, with ordered selection hash

```text
a45245b9a7fe3b73de131698f1f42d6f68c712e68078a2055c10d2b0d24fd37e.
```

It contains 640 K6 seeds, 136 single-apex survivors, and the 38 pair branches
just described.  The exact conjunction result is

| quantity | count |
|---|---:|
| graphs | 20 |
| K6 seeds | 640 |
| base empty-budget graph rejections | 0 |
| base impossible seeds | 0 |
| conjunction graph rejections | 0 |
| conjunction impossible seeds | 0 |
| marginal graph rejections | **0** |
| marginal impossible seeds | **0** |
| double-apex candidate pairs / survivors | 38 / 0 |

The reason for zero marginal Hall coverage is completely explicit.  For all
640 seeds, the first tested `Z0` is the empty set and every component accepts
the all-nonempty alternative.  Every recorded witness therefore has

```text
chosen_Z0 = [],       chosen_empty_vertices = [].
```

The conjunction evaluated 4,650 bipartite components, ignored 341
non-bipartite components conservatively, tested 9,300 generic orientation
Hall systems, and produced 4,786 component empty options.  Restricting empty
vertices cannot reject a seed whose all-nonempty state already passes, so the
zero marginal count is mathematically expected on this selection.

This negative result changes the recommended sampling strategy.  The
single-apex screen is extremely strong, but its conjunction value can be
measured only on seeds where the base Hall system actually *requires* one or
two empty singletons.  A future bounded target sample should be selected by
positive `minimum_empty_used`, rather than by K6-count strata.  No full
822-graph scan was launched here.

## Commands, reproducibility, and trust

The bounded command was

```text
/Users/lukacs/claude/opengauss/venv/bin/python3 \
  probe_d6_k6_virtual_k7_empty_conjunction.py \
  --workers 11 \
  --output d6_k6_virtual_k7_empty_conjunction_pilot_report.json
```

The report records 0.190 seconds of 11-worker kernel wall time on
`macOS-14.5-arm64-arm-64bit` with Python 3.11.15.  A second 11-worker run
recorded 0.197 seconds.  Removing runtime metadata and per-graph elapsed-time
fields makes the two reports canonically identical, with stable content hash

```text
394cb2daa6f18342ee95e9f074448ede28ce079ecd4a59f4969eae0124df49ab.
```

The report's compact graph decision rows have stable hash

```text
6966e45fcfa262d02bc8a7b775aca426165f801c1ed6d7589e1ea6e22c546caf.
```

The current empty-budget suite passed all nine controls, including its known
18-point positive control, fixed negative witness, checkpoint bindings,
infrastructure-abort semantics, and independent-checker import audit:

```text
/Users/lukacs/claude/opengauss/venv/bin/python3 -m unittest -v \
  test_d6_k6_empty_support_budget.py
```

Additional direct DP controls checked that the empty tuple is retained, a
singleton is admitted without a pair, a second singleton is admitted exactly
when its pair survived the double-apex screen, and no third singleton is ever
admitted.  Byte-code compilation and `git diff --check` also pass.

The read-only pilot is bound to these working sources/artifacts:

```text
d6_k6_empty_support_budget.py
  284b8c3bda4d2a43581665ade3d52d9295b29d460418454e3a760897cabc5fac

probe_d6_k6_virtual_k7.py
  531f4d123f249de03438e411f03fb7c70a591f3398b2a162b92158a31fd655b0

d6_k6_virtual_k7_pilot_report.json
  9f43436fddd5b1aa3511318251cdd992c0621d76807200166ffd49e7fb872148

probe_d6_k6_virtual_k7_empty_conjunction.py
  3a79237751aab9f4710ea85e14994a9ec106cc4f0074a2287b53240ef69d88ab

d6_k6_virtual_k7_empty_conjunction_pilot_report.json
  bytes 234,440
  03db68d823a0de72e8e32c22d857f16775e822747791ea2ce8ec82104ca92320
```

The conjunction is a discovery pilot and receives no theorem credit.  It
does not modify the current empty-budget source or any production artifact.
Candidate nonedges remain unconstrained, allowed defect coordinates may be
zero, all-nonempty remains an explicit branch, and a surviving/unexamined
branch is never turned into a rejection.
