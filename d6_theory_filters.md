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

### Orthonormal form of the zero-factor set

The bound and equality argument above have a simpler strengthening.  For
every `x in Z`, the simplex identities give `||u(x)||=1`.  For distinct
`x,y in Z` the vectors are always orthogonal.  This follows from the unit
equation when `xy` is a required edge.  When `xy` is a candidate nonedge,
`alpha(G)<=2` makes `D_x` and `D_y` disjoint, since a shared seed
non-neighbour would form an independent triple with `x,y`.  Thus

```text
{u(x):x in Z} is an orthonormal family in R^7.
```

This directly proves `|Z|<=7` and requires the allowed masks of `Z` to
have a matching saturating `Z`.  On the present K8-free corpus this generic
matching condition is new only at size seven.  If a cover of size at most
six failed Hall, a subfamily `Y` would have coordinate union `A` with
`|A|<|Y|<=6`.  Eligibility gives masks of size at least three, so, as
`|A|<=5`, two such masks cannot be disjoint.  The vertices of `Y` are
therefore a required clique.  Together with the `7-|A|` unaffected seed
vertices they make a required K8.  The existing tight-cover check is hence
exactly the nonredundant matching case.

The actual support `S_x=supp(u(x))` carries more information than its
allowed mask.  A necessary support-pattern CSP for any eligible cover is

```text
S_x subseteq D_x,                 |S_x| >= 3,
|S_x intersect S_y| != 1          for distinct x,y in Z,
the family (S_x:x in Z) has a matching saturating Z.
```

The intersection condition is exact: two nonzero vectors supported on sets
meeting in exactly one coordinate cannot be orthogonal.  When `|Z|=7`, the
matrix `U` with these columns is orthogonal.  The equations
`sum_i u_i(x)=-1` say `U^T 1=-1`, and orthogonality then gives `U1=-1`.
Consequently each row support has size one or at least three, never two, and
two row supports cannot meet in exactly one column.  These are support-only
necessary conditions; CSP feasibility does not assert that suitable real
coefficients exist.

### Normalized matrices on the nonzero remainder

For a feasible cover put `N=T\Z`, so every `c_x`, `x in N`, is nonzero, and
normalize

```text
v_x = sqrt(7) u(x)/c_x,       r_x = sqrt(7)/c_x.
```

Their Gram matrix `K` and the shifted matrix `B=K-J` have the exact forms

```text
K_xx = 1+r_x^2,               B_xx = r_x^2 > 0,
K_xy = 1, B_xy = 0            on required edges,
K_xy = 0, B_xy = -1           on candidate nonedges.
```

For a required edge inside `N`, the cover property rules out disjoint
allowed masks.  For a candidate nonedge, `alpha(G)<=2` forces those masks
to be disjoint, which justifies the displayed zero.  Moreover

```text
K is positive semidefinite with rank at most 7;
B has rank at most 8 and at most one negative eigenvalue.
```

There is a crucial dimension drop hidden by these coarse bounds.  Every
`u(z)`, `z in Z`, is orthogonal not only to the other zero-factor vectors but
also to every `u(x)`, `x in N`.  A required pair `zx` gives this from
`M_zx=0` and `c_z=0`; a candidate nonedge has disjoint allowed masks by
`alpha(G)<=2`.  Since the `u(z)` form an orthonormal family,

```text
rank(K) <= 7-|Z|,             rank(B) <= 8-|Z|.
```

If `nu_T` and `nu_N` are the allowed-support term ranks on all outside
vertices and on `N`, respectively, the combined support matrix also gives

```text
rank(K) <= min(7-|Z|, nu_N, nu_T-|Z|).
```

For the fixed `n=19` problem, where a K7 has twelve outside vertices, these
facts force the sharper bound `|Z|<=3`.  Every required clique `C` in `N`
has

```text
K_C = J + diag(r_x^2 : x in C),
```

which is positive definite and has rank `|C|`.  This rules out `|Z|=7`
from the positive diagonal, `|Z|=6` from the required edge forced in the six
vertices of `N`, and `|Z|=5` from the required triangle forced in the seven
vertices of `N` by `R(3,3)=6`.

For completeness, suppose `|Z|=4`.  Then `rank(K)<=3`, `rank(B)<=4`, and
`|N|=8`.  Put `F=overline(G[N])`.  The `F`-neighbours of any vertex form a
required clique; their positive-definite Gram vectors all lie in a
two-dimensional perpendicular space.  Hence `Delta(F)<=2`, so `F` is a
union of paths, cycles, and isolates.  The component/inertia bound below,
together with `nullity(B)>=4`, leaves only `F=4K2` or `F=C4+2K2`.  In the
first case rank at most four forces all four nonzero two-by-two `B` blocks
to be singular positive semidefinite.  Then `B` is positive semidefinite of
rank four, and `K=B+J` has rank at least four, a contradiction.  In the
second case, opposite vertices `b,d` of the `C4` are linearly independent,
while the other opposite pair `a,c` lies in their common one-dimensional
perpendicular.  This makes `a,c` collinear, impossible because `K_ac=1` but
both diagonal entries exceed one.

Thus the actual `Z` must be an eligible vertex cover of `L` of size at most
three.  This exact graph-only rule supersedes the earlier cap seven and its
tight-frame equality refinement, which remain documented above as valid
intermediate milestones.

These statements yield small exact graph tests.  If `Zf(H)` denotes the
ordinary zero-forcing number and
`U_K=min(7-|Z|,nu_N,nu_T-|Z|)`, then

```text
rank(K) >= |N|-Zf(G[N]),       rank(K) <= U_K,
rank(B) >= |N|-Zf(overline(G[N])),
rank(B) <= min(U_K+1,8-|Z|).
```

Two still cheaper consequences use the special diagonal, not merely the
zero--nonzero pattern.  Every required clique `C` in `G[N]` has

```text
K_C = J + diag(r_x^2 : x in C) > 0,
```

so its columns are independent and `omega(G[N])<=U_K`.  Also put
`F=overline(G[N])`.  Because `alpha(G)<=2`, the `F`-neighbours of any
vertex form a required clique.  Their Gram vectors are independent and all
lie in the perpendicular space of that vertex, whence

```text
Delta(F) <= U_K-1.
```

These positive-definite clique and perpendicular-neighbourhood tests are
sound before any zero-forcing search.  They can be strictly stronger than
ordinary zero forcing because zero forcing treats the nonzero entries as
arbitrary, whereas every required off-diagonal entry of `K` is exactly one
and every diagonal excess `r_x^2` is strictly positive.

There is a stronger exact equality-case test whenever `G[N]` contains a
required clique `C` of size exactly `U_K`.  Its vectors are independent, so
they are a basis of the whole normalized-vector span and `rank(K)=U_K`.
Put `R=N\C`, and let

```text
P = K[C,R],                    P_cx in {0,1},
A = adjacency matrix of G[R], K[R,R] = D+A,
D = diag(1+r_x^2 : x in R) > I.
```

Writing the Gram vectors as column matrices gives
`P=V_C^T V_R`.  Since `V_C` is invertible,

```text
ker(P) = ker(V_R) = ker(K[R,R]).
```

Consequently there must be a diagonal `D>I` for which

```text
(D+A) ker(P) = 0.             (basis-kernel compatibility)
```

This has a small exact rational test.  If the columns of `L` form any basis
of `ker(P)`, then independently for each row `x`, the row
`-(A L)_x` must be a common scalar multiple `d_x L_x` with `d_x>1`; when
`L_x=0`, `(A L)_x` must also vanish.  A failure for even one size-`U_K`
clique rules out the cover.  The separate positive-definite clique bound is
still needed when `omega(G[N])>U_K`, because a chosen size-`U_K` clique need
not by itself create a kernel relation in `P`.

The special form of the basis block yields an even cheaper and usually
stronger equality test.  Write

```text
K_C = J+diag(e_i),       x_i=1/e_i>0,       X=sum_i x_i.
```

For `y in N\C`, let `S_y` be its required-neighbour set in `C` and `a_y`
its zero-one indicator.  The zero Schur complement gives

```text
K_yz = a_y^T K_C^{-1} a_z,
K_C^{-1}=diag(x_i)-x x^T/(1+X).
```

It follows immediately that every `S_y` is nonempty and proper, that the
sets are pairwise distinct and pairwise intersecting, and that

```text
yz a candidate nonedge  =>  S_y,S_z are incomparable,
yz a required edge      =>  S_y union S_z is a proper subset of C.
```

Indeed, let `alpha,beta,gamma,delta` be the positive-coordinate sums on the
four Venn regions `S_y intersect S_z`, `S_y\S_z`, `S_z\S_y`, and outside
their union.  Then

```text
K_yz = (alpha(1+delta)-beta gamma)/(1+X).
```

A zero off-diagonal forces `alpha(1+delta)=beta gamma`, hence both set
differences are nonempty; a unit off-diagonal forces
`alpha delta-beta gamma=1+beta+gamma+delta`, hence `delta>0`.  Empty, full,
equal, or disjoint masks contradict the diagonal or the sign of this same
inverse formula.  These saturating-clique Schur mask rules use only finite
bit-set comparisons, but retain information that an arbitrary matrix-pattern
rank bound discards.

The lower bound follows by starting a hypothetical kernel vector at zero on
a zero-forcing set: each color-change row has exactly one still-unknown
coordinate and a nonzero off-diagonal coefficient.  The upper bounds follow
from the support rank of `U` and from adding the rank-one matrix `J`.

Against the old term-rank upper bound, ordinary zero forcing was provably
redundant once `Z` was a cover.  On `N`, a required edge has intersecting
allowed masks (otherwise it
would be an uncovered edge of `L`), while `alpha(G)<=2` makes every candidate
nonedge's masks disjoint.  Thus `G[N]` is exactly the intersection graph of
the `D_x`.  If `W` is their zero-one coordinate-incidence matrix, `W^T W`
has this exact off-diagonal graph and rank at most its term rank `nu_N`.
The universal zero-forcing bound applied to this particular matrix gives

```text
|N|-Zf(G[N]) <= rank(W^T W) <= nu_N.
```

The smaller cross-orthogonal upper bound `U_K` changes this conclusion:
ordinary zero forcing can reject a cover when its lower bound exceeds the
dimension left after removing `span u(Z)`.  The incidence inequality remains
an internal consistency check.

There is also an inertia-sensitive lower bound for `B`.  Let its off-diagonal
graph `F=overline(G[N])` have `p` nontrivial connected components, with
zero-forcing numbers `z_i`.  At most one component block can be indefinite.
Every other nontrivial block, if positive semidefinite, is an irreducible
Z-matrix with nullity at most one: after positive diagonal scaling it is
`I-C` for a connected nonnegative symmetric `C`, and Perron--Frobenius makes
a zero eigenvalue simple.  Isolated positive diagonal blocks are nonsingular.
Therefore

```text
nullity(B) <= max(p, max_i(p-1+z_i)).
```

This component/inertia bound is at least as strong as the ordinary `B`
zero-forcing bound; both are still recorded separately to expose where the
gain comes from.

All these tests quantify existentially over the zero-factor set.  A K7 seed
is rejected only after every eligible cover of `L` of size at most three has
failed the support CSP or one of the normalized-matrix rank tests.

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

## `K6` zero factors and the two-light-ray CSP

Retain the Lorentz factors

```text
ell_x=(c_x,z_x),
<ell_x,ell_y>_L=c_x c_y-z_x z_y,
```

and form `L` on the outside vertices by joining a required edge `xy` when
`D_x` and `D_y` are disjoint.  Such an edge makes `ell_x` and `ell_y`
Lorentz-orthogonal.  The zero pattern alone is not restrictive, because all
nonzero factors may lie on one lightlike projective line.  Coupling it back
to the `u` vectors gives a finite necessary condition.

Let `Z0={x:ell_x=0}`.  For `x in Z0`,

```text
||u(x)||^2=1,             sum_i u_i(x)=-1.
```

As in the K7 case, distinctness requires `|D_x|>=3`, while required edges
give orthogonality and candidate nonedges have disjoint masks.  Hence the
`u(x)`, `x in Z0`, are orthonormal in `R^6`: `|Z0|<=6` and their allowed
masks must match into the six seed coordinates.

For a K6-only graph the matching of `Z0` by itself is automatic, though the
implementation keeps the check as a useful assertion.  A failed Hall
subfamily `Y` has coordinate union `A` with `|A|<|Y|<=6`, hence `|A|<=5`.
All masks in `Y` have size at least three, so each pair intersects and the
corresponding vertices form a required clique.  Adding the `6-|A|`
unaffected seed vertices produces a required K7, contrary to the K6-only
hypothesis.  Matching the larger lightlike color bins below is not redundant.

Remove `Z0` from `L`.  Along an edge, Lorentz orthogonality applies the
orthogonal-line involution on the projective line.  An odd cycle fixes its
direction, so every nonzero factor in a connected non-bipartite component
lies on one of the two lightlike projective lines.  All components choosing
the same line have pairwise Lorentz product zero.  Their `u` vectors are
then pairwise orthonormal: the unit equation handles required pairs and
disjoint allowed masks handle candidate nonedges.  Zero-factor vectors are
orthogonal to both classes.

It follows that, for each K6 seed, the following finite CSP is necessary:

1. choose `Z0` among vertices with `|D_x|>=3`, with `|Z0|<=6`, and require
   its masks to be matchable;
2. find the connected non-bipartite components of `L-Z0`;
3. assign each such component to one of two lightlike colors;
4. for each color, require the masks of its component vertices together
   with `Z0` to match into the six coordinates (and hence total at most six).

Bipartite components are intentionally ignored: their nonzero factors may
alternate between a generic pair of Lorentz-orthogonal projective lines.
The search is existential in `Z0` and in the two-color assignment; a seed is
impossible only when every choice fails.  Vertices in an initially
bipartite component never need to enter `Z0`: deleting them cannot create an
odd component, and omitting them from `Z0` only relaxes both matching bins.
This permits a sound first pruning before the at-most-thirteen-vertex subset
search.

The support-matching helper is also tested outside the K6-only population.
Its synthetic failure is intentionally a kernel-level control rather than a
graph-level isolated rejection, since the preceding argument shows that a
K6-only graph cannot fail on `Z0` alone.  A separate 13-vertex negative
control satisfies `alpha(G)<=2` and has clique number exactly six; relative
to its distinguished seed, all seven outside vertices are ineligible for
`Z0` and form one connected non-bipartite component of `L`, which cannot fit
in either six-dimensional orthonormal bin.

### Bipartite-component zero-forcing rank refinement

The preceding pruning of zero candidates in initially bipartite components
is safe only while those components are ignored.  The following stronger
rule uses them and therefore returns to **every** eligible `Z0` of size at
most six.

Fix a connected bipartite component `C=A union B` of `L-Z0`.  If its Lorentz
direction is generic, connectivity puts the nonzero factors on one
non-lightlike projective line over `A` and its orthogonal line over `B`.
Within `A`, the Gram matrix of the defect vectors has off-diagonal graph
exactly `G[A]`: a required edge has a nonzero Lorentz product, while a
candidate nonedge has disjoint allowed masks by `alpha(G)<=2`.  The analogous
statement holds on `B`.  Every cross-side defect-vector pair is orthogonal,
as are all the unit vectors indexed by `Z0`.  Ordinary zero forcing therefore
gives

```text
|Z0| + (|A|-Zf(G[A])) + (|B|-Zf(G[B])) <= 6.       (K6-bip-rank)
```

If the component direction is lightlike, its orthogonal line is itself.  All
component defect vectors are then orthonormal, giving the stronger bound
`|Z0|+|A|+|B|<=6`, so `K6-bip-rank` remains necessary.  An isolated component
vertex contributes zero because `Zf(K1)=1`.

The combined exact search exhausts every eligible `Z0`, applies this rule to
every bipartite component of `L-Z0`, and applies the earlier two-light-ray
coloring and joint actual-support rules to the non-bipartite components.  On
the complete 1,098-graph residue after the actual-support refinement it
rejects one further graph, corpus index `461363`, leaving 1,097.  The full
proof, all-eligible-`Z0` regression control, source-bound report, and complete
verifier are in `d6_k6_bipartite_rank_audit.md` and the associated
`d6_k6_bipartite_rank_*` files.
