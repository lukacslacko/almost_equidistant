# Odd-cycle/full-support pinning increment

## Exact theorem and measured increment

The triangle double-pin is a special case of a broader exact obstruction.
For a fixed propagated support coordinate, every required edge whose endpoint
masks intersect in exactly that coordinate gives a reciprocal equation.  A
nonbipartite component of this coordinate graph pins every component value to
`+1` or `-1`.  If all coordinates allowed at one vertex are pinned, the
normalized K7 point identity is impossible.

The complete 155-graph measurement rejects 131 graphs: all 130 double-pin
rejections plus graph **3919831**.  To preserve the already frozen and
independently checked 130-graph base, the theorem package records only this
one-graph increment and binds it to the verified base report.

After this increment, the exact K7 residue has 24 graphs:

```text
226183 316173 423661 424226 2581209 2592657 2593240
3595554 3624785 3648882 3729907 3785980 3888410 3935560
3936176 3936177 3936310 3936435 3945490 3945555 3945557
3945564 3947605 3949382
```

Its stable JSON SHA-256 is

```text
71d9ae101445ae90d08699ede81ca405bf5bfe956fda01ee384203fa113d1779
```

Separately verified interval certificates currently reject 3936176, 423661,
and 424226 from this 24-list.  Their union with the exact pinning layers would
leave 21 K7 graphs, but that 21-graph boundary should be published only by a
hash-bound union manifest after the interval artifacts are frozen.

## Proof

Fix a regular unit K7 and a zero-factor cover.  On its nonzero-factor set put

```text
w_i(x)=sqrt(7) u_i(x)/t_x,    t_x=1+sum_i u_i(x) != 0.
```

A required unit edge satisfies `w(x).w(y)=1`.  For coordinate `i`, form the
graph `H_i` whose edges are required edges `xy` for which the propagated
allowed masks intersect in exactly `{i}`.  Every edge of `H_i` gives

```text
w_i(x) w_i(y)=1.                                  (1)
```

Along a path, values alternate between `a` and `1/a`.  An odd cycle forces
`a=1/a`, hence `a^2=1`; every vertex in that connected component is therefore
pinned to a sign.  Equation (1) also proves the coordinate is actually
nonzero, despite candidate nonedges and allowed-mask optionality.

For one point define

```text
T=sum_i w_i,       Q=sum_i w_i^2.
```

The diagonal simplex identity and the definition of `t` give

```text
Q - T^2 + 2 sqrt(7) T - 8 = 0.                    (2)
```

If every bit of the propagated mask is pinned, every allowed coordinate is
actually present and equals a sign.  Hence `T` is an integer and
`Q=|support|<=7`.  Separating the rational and `sqrt(7)` parts of (2) forces
`T=0` and then `Q=8`, a contradiction.

This proof uses only required edges, propagated support supersets, and exact
algebra.  It never requires a candidate nonedge to be genuinely nonunit.

## Graph 3919831 certificate scope

The rejecting seed is

```text
[5,7,10,11,13,16,18].
```

It has three inherited-current covers.  Two have no family surviving the
older propagation/sparse-value conjunction.  Across the remaining branches,
all five older-passing families contain a fully pinned vertex.  The key
propagated mask is `67={0,1,6}`; each of its three coordinate components is
nonbipartite.  Exhausting these three covers rejects the graph.

## Independent verification

`verify_d6_k7_full_pin_increment.py` imports neither the incremental builder
nor the production full-pin locator.  It independently reconstructs graph
3919831 from the v5 corpus, all K7 seeds and eligible covers, the inherited
cover gates, labeled supports, propagation, sparse-value checks, coordinate
graphs, bipartite color conflicts, and the final seed quantifier.  It compares
the complete target record and the 25-to-24 ordered residue update against the
production report.

## Pre-production commands

Freeze and commit the source first.  Then build and independently verify the
130-graph double-pin base as documented in
`d6_k7_two_defect_double_pin.md`.  With the printed base hashes, run:

```sh
python3 build_d6_k7_full_pin_increment.py \
  --base-report-sha256 <double-pin-report-sha256> \
  --base-verification-sha256 <double-pin-verification-sha256> \
  --output d6_k7_full_pin_increment_report.json

python3 verify_d6_k7_full_pin_increment.py \
  --report d6_k7_full_pin_increment_report.json \
  --report-sha256 <full-pin-report-sha256> \
  --base-report-sha256 <double-pin-report-sha256> \
  --base-verification-sha256 <double-pin-verification-sha256> \
  --output d6_k7_full_pin_increment_verification.json

python3 -m unittest -v \
  test_d6_k7_full_pin_odd_cycle.py \
  test_d6_k7_full_pin_increment.py
```

No official increment artifact should be generated before the source commit.
