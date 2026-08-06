# Exact graph-only filters for the dimension-6 campaign

This note documents the proof obligation behind `profile_d6.c`.  It is
deliberately independent of the interval realization engine.

Let `G` be the unit-distance graph of an almost-equidistant set in
`R^6`, and let `Q` be a unit clique of size `m`.  The common unit
neighbours of `Q` lie in the intersection of the unit spheres centred
at `Q`.  If `c` is the centroid of `Q`, this intersection is a sphere
of squared radius

```text
(m + 1) / (2m)
```

in the affine space `c + aff(Q-c)^perp`, of dimension `6-m+1`.  The
common-neighbour set is itself almost equidistant.  Consequently:

| seed | link bound | reason |
|---|---:|---|
| `K2` | 16 | it is an almost-equidistant set in an affine `R^5`; `f(5)=16` |
| `K3` | 12 | affine `R^4`; `f(4)=12` |
| `K4` | 10 | affine `R^3`; `f(3)=10` |
| `K5` | 4 | circle argument below |
| `K6` | 2 | a zero-sphere consists of two points |
| `K7` | 0 | the common unit-sphere intersection is empty |

For `K5`, the link circle has squared radius `3/5`.  After normalizing
the circle to radius one, a unit chord has inner product `1/6`, so its
angular step `theta` satisfies `2 cos(theta)=1/3`.  The step is not a
rational multiple of pi: if it were, `2 cos(theta)` would be a rational
algebraic integer and hence an integer, contrary to `1/3`.  Thus the
finite unit-distance graph on the circle has maximum degree two and no
cycle; it is a union of paths.  A path-union on at least five vertices
has an independent set of size at least three, contradicting almost
equidistance.  Hence the link has at most four vertices.

All bounds count graph common neighbours only.  They never assume that
a graph nonedge has non-unit distance.  A violation is therefore a
sound non-realizability certificate consisting of the seed clique and
too many common neighbours.

`profile_d6.c` checks these rules using exact bit operations and reports
individual and cumulative counts.  It also reports whether each graph
has clique number 6 or 7.  The profiler validates symmetry, looplessness,
and `alpha(G) <= 2` before using a graph.

## `K7` facet-reflection rule

Fix a centered regular unit simplex `q_1,...,q_7` in `R^6`.  A point
`x` at unit distance from the six vertices other than `q_i` lies in the
zero-sphere cut out by those six unit spheres.  Its two possibilities
are `q_i` itself and the reflection of `q_i` through the opposite
facet.  Distinctness excludes `q_i`.  Since

```text
||q_i||^2 = 3/7,             <q_i,q_j> = -1/14  (i != j),
reflection(q_i) = -(4/3)q_i,
```

two outside vertices of the same omitted-vertex type are forced to
coincide, while reflections of two different types have squared
distance `16/9`, not one.  Thus a candidate is impossible if, relative
to any `K7` seed,

- two outside vertices each have exactly that same one non-neighbour in
  the seed; or
- two such forced-reflection vertices of different types are joined by
  a required unit edge.

This rule uses only required edges.  A candidate nonedge remains
unconstrained.  The certificate is a `K7` plus the offending one or two
outside vertices.

## Discrete `K7` defect-support CSP

For a point `x` outside a fixed `K7`, write

```text
u_i(x) = ||x-q_i||^2 - 1,
S_x = { i : u_i(x) != 0 }.
```

The actual support `S_x` is a nonempty subset of the seed vertices that
are graph non-neighbours of `x`; allowed non-unit pairs are permitted to
turn out unit.  The profiler applies these additional exact facts:

1. If two actual supports overlap, the corresponding outside vertices
   must be unit: with a shared `q_i`, the other two sides of the triple
   are both non-unit.
2. A one-defect support `{i}` forces the reflection `-(4/3)q_i`.
   Therefore two singleton supports cannot have the same type, and two
   different types can coexist only when their outside vertices are a
   graph nonedge.
3. At most two vertices have one-defect support.  This also follows from
   item 2 and `alpha(G)<=2`.
4. At most two vertices have any fixed two-defect support `{i,j}`.  Such
   vertices are pairwise unit by item 1 and lie on the common-neighbour
   circle of the other five seed vertices.  That circle has squared
   radius `3/5`; three pairwise-unit points would require circumradius
   squared `1/3`, so three are impossible.

For an outside vertex with three or more allowed defects, choosing the
whole allowed set satisfies this discrete layer.  Only vertices with
one or two allowed defects therefore require branching.  The CSP tries
all nonempty subsets of those allowed sets subject to items 2--4.
Overlap across a required graph nonedge cannot occur because the
candidate complement is triangle-free: two adjacent complement
vertices have disjoint complement neighbourhoods in the independent
seed.  CSP infeasibility for any `K7` seed is a rigorous obstruction.
Feasibility is only a screen and makes no realizability claim.
