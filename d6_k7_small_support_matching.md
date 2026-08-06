# Three-disjoint-support obstruction after a `K7` seed

This note proves the exact refinement implemented independently in
`d6_k7_small_support_matching.py`.  It does not alter the frozen sparse-value
layer.

## Normalized unit condition

For the centered regular unit `K7` seed in `R^6`, retain the exact simplex
coordinates

```text
u_i(x) = ||x-q_i||^2-1,       t_x = 1 + sum_i u_i(x),
w_i(x) = sqrt(7) u_i(x)/t_x.
```

For every assigned point whose actual defect support has size one or two,
`t_x` is nonzero:

- With one defect, `t_x=0` gives the root `u_i=-1`, which is the seed point
  `q_i` itself and violates distinctness.  The nonduplicate root is
  `u_i=4/3`, with `t_x=7/3`.
- With two exact nonzero defects, `t_x=0` would give
  `u_i+u_j=-1` and, from the simplex quadratic,
  `u_i^2+u_j^2=1`.  Thus `u_i u_j=0`, contradicting exact two-defect
  support.

The simplex bilinear identity rewrites as

```text
||x-y||^2 - 1 = (t_x t_y/7) (1 - <w(x),w(y)>).
```

Because both factors `t_x,t_y` are nonzero, a unit pair must satisfy
`<w(x),w(y)>=1`.  If the two actual defect supports are disjoint, their
normalized dot product is instead zero, so the pair is necessarily
non-unit.

## Set-packing obstruction

Suppose three assigned one/two-defect points had pairwise-disjoint nonempty
actual supports.  The preceding identity would make all three pairs
non-unit.  That contradicts the defining almost-equidistant condition that
every three points contain a unit pair.  Therefore:

> The assigned one/two-defect supports have set-packing number at most two.

This applies to a partial support assignment: adding more points or enlarging
the still-unassigned part cannot repair the already present non-unit triple.
It is therefore sound to prune as soon as the third pairwise-disjoint actual
small support is assigned.  Candidate nonedges and unassigned supports of
size at least three remain unconstrained.

For exact two-defect points alone, identify support `{i,j}` with edge `ij` in
the simple type graph on the seven seed coordinates.  The theorem says this
graph has matching number at most two.  It rejects, among other patterns,
`3K2`, the six-vertex path `P6`, and the six-cycle `C6`.  In particular it
removes the `C6` case deliberately retained by the order-six Mobius rule.
The stronger formulation also rejects mixed packings such as a singleton
plus two mutually disjoint two-defect types, or two singleton supports plus
one disjoint two-defect type.

## Redundancy in the current `N`-layer CSP

The theorem is useful structurally, but it supplies no new branch rejection
inside the existing frozen `N`-layer support CSP.  The graph `graph_n` there
is an induced subgraph of a candidate unit-edge graph with independence
number at most two.  Consequently any three assigned points contain a
required edge of `graph_n`.  If their supports were pairwise disjoint, the
endpoints of that required edge would have disjoint supports.  The frozen
incremental rule `disjoint_required_small_support` already rejects exactly
that condition as soon as the second endpoint is assigned.

This also explains the apparent six-cycle improvement.  The old unit test's
synthetic `C6` point graph has independence number three: its three
alternating type-points have pairwise-disjoint supports.  Any completion to
independence number at most two must add a required edge within each such
independent triple, but the endpoints' allowed support masks are already
disjoint.  The support-intersection layer then rejects before the Mobius or
set-packing value rules run.

Therefore a zero marginal count is forced for every valid input to this
particular pipeline.  The separate implementation is retained as an exact
audit control and for possible reuse in a context that has not already
encoded `alpha<=2`; no full-corpus post-layer campaign is warranted.
