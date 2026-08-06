# Required-edge K7 Mobius-cycle obstructions

## Result

There are two small, purely required-edge unit-distance obstructions in
dimension six:

- a `K7` plus three outside roles, of order 10 and size 39;
- a `K7` plus four outside roles, of order 11 and size 45.

The outside roles form, respectively, a triangle or quadrilateral of exact
two-defect types.  The order-six Mobius dynamics after a regular `K7`
forbids both cycle lengths.  No absent pattern edge is assigned a distance.

An exhaustive containment screen on the hash-pinned 12,839-graph current K7
selection found zero copies of either obstruction.  Thus these attractive
cores explain heuristic source patterns 905, 936, and 952 but are exactly
redundant after the already completed K7 filters.

## Canonical required-edge cores

Label the seed vertices `q_0,...,q_6`.  For a pair `{i,j}`, a role `x_ij`
has required unit edges to the five seed vertices outside `{i,j}`.  Edges
from `x_ij` to `q_i` or `q_j` are deliberately absent from the pattern and
therefore unconstrained.

The triangle core has roles

```text
x_01, x_12, x_20
```

and required role edges

```text
x_01 x_12,  x_12 x_20,  x_20 x_01.
```

The quadrilateral core has roles

```text
x_01, x_12, x_23, x_30
```

and the four cyclic role edges.  Together with the 21 seed edges and five
seed edges per role this gives 39 and 45 required edges.

The heuristic source motifs canonically relabel as follows.  These labels
are provenance only; the proof and containment scan use the canonical cores.

- pattern 905: seed `[0,1,2,5,10,11,12]`; canonical coordinates
  `(0,1,2)=(2,3,1)` and roles `(x_01,x_12,x_20)=(3,4,9)`;
- pattern 936: seed `[0,1,2,5,11,12,13]`; canonical coordinates
  `(0,1,2)=(2,3,0)` and roles `(3,6,8)`;
- pattern 952: seed `[4,7,9,10,11,12,13]`; canonical coordinates
  `(0,1,2,3)=(1,4,2,3)` and roles `(x_01,x_12,x_23,x_30)=(0,2,6,3)`.

## Why the roles have exact two-defect support

For a point `x` relative to a centered regular unit `K7`, put

```text
u_i(x)=||x-q_i||^2-1,  s_x=sum_i u_i(x),  t_x=s_x+1.
```

The simplex identities are

```text
7 sum_i u_i(x)^2-t_x^2=7,
7 u(x).u(y)=t_x t_y              when xy is a required unit edge.
```

A role `x_ij` is unit from the other five seed vertices, so its actual
defect support is a subset of `{i,j}`.  It is nonempty: support zero would
contradict the first simplex identity.

Every one- or two-defect point distinct from the seed has `t_x != 0`.
For one defect the point equation has roots `u_i=-1,4/3`; the first is the
seed point `q_i` itself and distinctness excludes it.  For two exact defects,
`t_x=0` would give `u_i+u_j=-1` and `u_i^2+u_j^2=1`, hence `u_i u_j=0`, a
contradiction.

Now consider two consecutive roles whose allowed pairs meet only at seed
coordinate `i`.  Their required unit edge and the second identity show that
their actual supports cannot be disjoint: the left side would be zero while
`t_x t_y` is nonzero.  Hence both actual supports contain `i`.  Each role in
the role cycle has two neighbours, one at each endpoint of its allowed pair.
Both endpoint defects are therefore nonzero, making every role an exact
two-defect point.

This step is why the proof is monotone.  Extra unit edges, including an edge
from a role to one of its two omitted seed vertices, do not get interpreted
as non-unit constraints; the cyclic required edges simply rule out the
resulting smaller actual support.

## Exact Mobius contradiction

Put `rho=sqrt(7)` and normalize an exact two-defect point by

```text
w_i=rho u_i/t.
```

If its two components are `A,B`, its point equation gives

```text
(A-rho)(B-rho)=3,
phi(A)=(rho A-4)/(A-rho).
```

Two consecutive roles have exactly one common support coordinate.  Their
required unit edge gives a product of one at that coordinate.  Crossing a
type and then crossing the required role edge therefore applies

```text
T(A)=1/phi(A)=(A-rho)/(rho A-4),
M=[1 -rho; rho -4].
```

Exact arithmetic in `Q(sqrt(7))` gives

```text
M^6=-27 I.
```

The real fixed-point discriminants of `T^k`, for `k=1,...,5`, are

```text
-3, -27, -108, -243, -243.
```

Thus no nonidentity power has a real projective fixed point.  Traversing a
simple type cycle of length `ell` requires a fixed point of `T^ell`, so its
length must be divisible by six.  Lengths three and four are impossible.
The checker retains length six as a negative control for this rule.

This proves that each canonical graph is non-realizable by distinct points
in `R^6` using only its listed unit edges.  Any candidate containing either
graph as a non-induced subgraph is therefore impossible.

## Complete containment enumeration

For each 19-vertex target, the production matcher exhausts every required
`K7`.  Relative to a fixed seed it checks all 35 coordinate triangles and
all 105 unoriented coordinate quadrilaterals.  A role domain consists of
every outside target vertex adjacent to the five required seed images;
adjacency to either omitted seed is allowed.  An injective backtrack then
checks only the cyclic role edges.

An embedding of either canonical core necessarily maps its seven designated
seed vertices to one of the enumerated target `K7`s.  Triangle symmetry
reduces the coordinate choices to the 35 triples.  A four-set has exactly
three Hamilton cycles modulo reversal and rotation, giving `35*3=105`
quadrilaterals.  The enumeration is therefore complete for non-induced
subgraph containment.

On the current selection it checked

```text
14,453 K7 seeds,
505,855 triangle coordinate cycles,
1,517,565 quadrilateral coordinate cycles.
```

Both hit sets are empty.  The independent verifier uses separately written
clique, domain, and direct triangle/quadrilateral assignment routines and
repeats the complete scan.  It also checks the pattern construction and the
Mobius matrix calculation over `Fraction` arithmetic.

## Reproduction and trust scope

```text
python3 -m unittest -v test_d6_k7_mobius_cycle_containment.py
python3 d6_k7_mobius_cycle_containment.py \
  --output d6_k7_mobius_cycle_containment_report.json
sha256sum d6_k7_mobius_cycle_containment_report.json
python3 verify_d6_k7_mobius_cycle_containment.py \
  --report-sha256 <exact report SHA256>
```

The proof trusts Euclidean algebra, exact integer bit masks, Python arbitrary
precision integers, and `Fraction`.  Floating point is used only to report
wall time.  The target TSV and its 12,839-index order are hash-pinned.  A
failed or absent containment search would make no graph-theoretic claim;
here both exhaustive implementations terminate on every target.
