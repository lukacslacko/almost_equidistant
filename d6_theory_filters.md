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

## Exact simplex matrix identities

For any represented points define

```text
M_xx = -1,
M_xy = ||x-y||^2 - 1  (x != y).
```

Required graph edges give zero off-diagonal entries, while graph
nonedges remain unrestricted and may also give zero.  The matrix has
rank at most eight and at most one positive eigenvalue.

For a centered regular unit `K7`, put

```text
u_i(x) = ||x-q_i||^2 - 1,
s_x = sum_i u_i(x),
c_x = s_x + 1.
```

Direct use of the regular-simplex Gram matrix gives

```text
x = -sum_i u_i(x) q_i,
||u(x)||^2 = 1 + c_x^2/7,
M_xy = c_x c_y/7 - u(x) dot u(y).
```

If `U` has the `u(x)` as columns and `T` is the outside set, the Schur
complement of the seed block `M_QQ=-I_7` is therefore exactly

```text
M_TT + U^T U = c c^T/7.
```

It is positive semidefinite of rank at most one; this is stronger than
a bare rank-eight condition.

## `K7` simplex-clique Hall rule

For each outside graph vertex let

```text
D_x = {i : xq_i is a graph nonedge}.
```

The actual support of `u(x)` is contained in `D_x`.  If `C` is a
required unit clique among outside vertices, then the preceding
identities give

```text
U_C^T U_C = I + c_C c_C^T/7.
```

This matrix is positive definite, so the columns `u(x)`, `x in C`, are
linearly independent.  Sparse-column Hall theory now requires a
matching from every such clique into its allowed seed coordinates.
Equivalently, for every coordinate mask `A subseteq {1,...,7}`,

```text
omega(G[{x : D_x subseteq A}]) <= |A|.
```

The C implementation uses the equivalent direct formulation: it
enumerates every outside clique while accumulating `union D_x`, and
rejects exactly when

```text
|C| > |union_{x in C} D_x|.
```

This checks every Hall subset, not only one maximal clique.  Its witness
packs the failing outside clique into the low 19 bits of `auxiliary`
and its coordinate union into the next seven bits.

On the present corpus this rule is expected to be redundant with the
existing `K8` exclusion.  A failing clique of size at least `|A|+1`,
together with the `7-|A|` seed vertices outside `A`, is a required
`K8`.  Keeping the independent Hall implementation is still useful as
an algebra/control check and for variants not prefiltered by `K8`.

## `K7` disjoint-edge bounded-cover rule

Form a graph `L` on the twelve outside vertices.  Put `xy` in `L` when
`xy` is a required unit edge and `D_x` and `D_y` are disjoint.  The unit
equation then has `u(x) dot u(y)=0`, so

```text
c_x c_y = 0.
```

Consequently `Z={x:c_x=0}` must be a vertex cover of `L`.

Only vertices with at least three allowed defects are eligible for
`Z`.  Indeed, `c_x=0` implies

```text
sum_i u_i(x)=-1,  ||u(x)||^2=1.
```

With at most two nonzero coordinates these equations force
`u(x)=-e_i`, which reconstructs the existing seed point `q_i` and
violates distinctness.

Moreover `|Z|<=7`.  The seed and all points in `Z` lie on the sphere of
squared radius `3/7`.  Sending

```text
p -> (sqrt(2)p, 1/sqrt(7)) in R^7
```

makes them unit vectors and turns unit distance into orthogonality.  If
there are `N` such vectors, the graph of nonorthogonal pairs is
triangle-free.  Bessel's inequality at each vector gives
`tr(B^2)<=2N` for their Gram matrix `B`; since `rank(B)<=7`,

```text
N^2 = tr(B)^2 <= 7 tr(B^2) <= 14N.
```

Thus `N<=14`; the seven seed vectors leave room for at most seven
members of `Z`.  The exact filter rejects precisely when `L` has no
vertex cover of size at most seven contained in
`{x:|D_x|>=3}`.  The implementation uses exhaustive edge-branching on
at most twelve vertices.

The cardinality bound uses almost-equidistance, via the validated
condition `alpha(G)<=2`.  The bounded-cover helper is therefore not
advertised as an obstruction for an arbitrary unit-edge graph lacking
that condition.

### Equality refinement at cover number seven

There is one further exact consequence when the eligible vertex-cover
number of `L` is exactly seven.  The actual set `Z={x:c_x=0}` is an
eligible cover and has size at most seven, so in this case `|Z|=7`.
The lifted seed together with `Z` consists of fourteen unit vectors in
`R^7`.  Equality holds throughout the preceding `N<=2r` proof:

```text
196 = tr(B)^2 <= rank(B) tr(B^2) <= 7*28 = 196.
```

Consequently the frame operator is `2I`.  The seven lifted seed vectors
already form an orthonormal basis and contribute `I`, so the seven
lifted vectors belonging to `Z` form a second orthonormal basis.

In seed-basis coordinates,

```text
<w(x),w(q_i)> = -u_i(x)  for x in Z.
```

Thus the `7 by 7` matrix with columns `u(x)`, `x in Z`, is nonsingular.
Its actual support has a determinant permutation, and hence the allowed
defect masks `D_x` of these seven vertices must admit a perfect matching
to the seven seed coordinates.

The implementation therefore distinguishes three cases for each seed:

1. no eligible cover of size at most seven: reject by the original
   bounded-cover rule;
2. an eligible cover of size at most six: the equality refinement says
   nothing;
3. cover number exactly seven: enumerate **every** eligible size-seven
   cover and reject only if none of their allowed masks has a perfect
   matching.

Checking every minimum cover is essential: the actual `Z` need not be a
canonical cover chosen by the algorithm.  Candidate nonedges remain
unconstrained throughout; the matching is only a necessary support
condition.  A dedicated 15-vertex, `alpha<=2` seed-level control has
relative cover graph `K_{1,7}`, a unique eligible seven-leaf cover, and
leaf masks whose union omits one seed coordinate.  Relative to its
distinguished seed it passes the earlier cover and Hall rules but fails
this refinement.  (It has another `K7` seed that already fails the old
cover rule, so it is not advertised as a graph-level isolation control.)

## `K6` coordinates and virtual-simplex Hall rule

Let `q_1,...,q_6` be a centered regular unit simplex spanning a
five-dimensional space `W`, and let `e` be a unit normal.  Then

```text
q_i dot q_i=5/12,  q_i dot q_j=-1/12,
x=-sum_i u_i(x)q_i+h_x e,
z_x=sqrt(12)h_x,
z_x^2=c_x^2+6-6||u(x)||^2,
M_xy=c_x c_y/6-u(x) dot u(y)-z_x z_y/6.
```

The seed Schur complement is

```text
6(M_TT+U^T U)=c c^T-z z^T.
```

It has rank at most two, at most one positive eigenvalue, and at most
one negative eigenvalue; it is not generally positive semidefinite.

For the first graph-only consequence define the seven-coordinate vector

```text
a_x = (u(x), z_x/sqrt(6)).
```

Its first six coordinates are supported on the graph-allowed defect
set `D_x`, while its last (normal) coordinate is universally allowed.
For a required outside clique `C`, the equations give

```text
a_C^T a_C = I + c_C c_C^T/6,
```

which is positive definite.  Thus a matching into the six defect
coordinates plus the universal normal coordinate is necessary.  Hence,
for every
`A subseteq {1,...,6}`,

```text
omega(G[{x:D_x subseteq A}]) <= |A|+1.
```

Equivalently every required outside clique `C` must satisfy

```text
|C| <= |union_{x in C}D_x|+1.
```

The universal coordinate is used only for this rank/Hall argument.  It
must not be treated as a real vertex in the overlap-forces-unit CSP.

This Hall rule too is redundant after the corpus `K8` filter: a failing
outside clique of size at least `|A|+2` together with the `6-|A|` seed
vertices outside `A` is a required `K8`.
