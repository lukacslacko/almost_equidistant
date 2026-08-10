# K6 opposite-ray rank, support, and all-ones conjunction

This note proves the exact necessary condition implemented by
`d6_k6_opposite_ray_rank_conjunction.py` and independently transcribed in
`verify_d6_k6_opposite_ray_rank_conjunction.py`.  The input is the ordered
623-graph residue of the frozen tight-same-`Z0` package.  Every decision uses
exact finite combinatorics and integer bit masks.  Numerical optimization and
floating-point rank tests play no role.

The essential point is quantifier coupling.  A realization has one actual
zero-Lorentz-factor set `Z0`, one assignment of its nonbipartite components to
the two light rays, and one family of actual defect supports.  All the tests
below are imposed on that same choice.  Candidate nonedges remain
unconstrained and every coordinate allowed by a candidate defect mask may
still vanish.

## 1. Exact K6 coordinates

Fix a centered regular unit `K6`, with vertices `q_1,...,q_6`, in a
five-dimensional subspace of `R^6`.  For an outside point `x`, put

```text
u_i(x) = ||x-q_i||^2-1,
c_x    = 1+sum_i u_i(x),
z_x    = sqrt(12) times the remaining normal coordinate.
```

The exact diagonal and pair identities are

```text
z_x^2 = c_x^2+6-6||u_x||^2,

||x-y||^2-1
  = (c_x c_y-z_x z_y)/6-u_x.u_y.                    (1)
```

Write `ell_x=(c_x,z_x)` and use the Lorentz product

```text
<ell_x,ell_y>_L = c_x c_y-z_x z_y.
```

Let `D_x` be the set of seed coordinates where `u_x` is allowed to be
nonzero, and let `S_x=supp(u_x)` be its actual support.  Only

```text
S_x subseteq D_x
```

is known.  In particular, a candidate nonedge is never required to have
non-unit distance.

All candidate graphs have independence number at most two.  Hence a candidate
nonedge `xy` implies `D_x intersect D_y=empty`: otherwise a shared seed
coordinate would give an independent triple consisting of `x`, `y`, and that
seed vertex.  It follows that `u_x.u_y=0` for such a pair, without assigning
its distance.

Define

```text
Z0 = {x : ell_x=0}.
```

For `x in Z0`, (1) gives `||u_x||=1` and `sum_i u_i=-1`.  Distinctness from
the seed implies `|S_x|>=3`.  The `Z0` vectors are orthonormal.

After deleting `Z0` from the disjoint-defect required-edge Lorentz graph, each
connected nonbipartite component lies on one of the two Lorentz light rays

```text
z=+c  or  z=-c.
```

The implementation enumerates every assignment of these components to the
two rays, modulo only simultaneous exchange of the ray names.  Each ray bin,
together with the common `Z0`, is an orthonormal family of defect vectors in
`R^6`.

## 2. Frozen same-`Z0` parent condition

For every possible `Z0`, the producer first invokes the frozen exact
bipartite-component condition from `d6_k6_tight_same_z0.py`.  This checks the
normal inertia, zero-forcing, arbitrary-subset rank, and block-support
conditions on the same `Z0`.  A geometric realization supplies a passing
parent state, so discarding a `Z0` for which none exists is sound.

This layer deliberately does not identify a parent bipartite-state support
witness with the later nonbipartite ray-support witness.  Keeping them
separate is a relaxation, not an extra rejection rule.

## 3. Opposite-ray projection rank

Fix the two nonzero light sets `A` and `B`, with `Z0` included in each
orthonormal bin.  Let `X` be the Gram block from the `A` vectors to the
complete `B union Z0` bin.  Its `Z0` columns vanish.  After naming the rays
oppositely, (1) gives, for `i in A` and `j in B`,

```text
X_ij = c_i c_j/3,  if ij is a required edge,
X_ij = 0,          if ij is a candidate nonedge.     (2)
```

Every displayed `c` is nonzero because the light factor is nonzero.  Project
the orthonormal `A` rows off the opposite complete bin.  The residual Gram
matrix is

```text
R = I-X X^T >= 0,
rank(R) <= 6-|B union Z0|.                            (3)
```

For two distinct `A` rows, diagonal congruence by their nonzero `c` values
makes the off-diagonal entry strictly negative exactly when they share a
required neighbour in `B`.  There is no cancellation: the relevant sum is a
sum of positive squares.  Each nontrivial connected component `H` of this
strict off-diagonal support is an irreducible PSD Z-matrix block.  The
Perron--Frobenius theorem gives

```text
rank(R_H) >= |H|-1.                                  (4)
```

An `A` row having no required neighbour in `B` has a zero row in `X`, hence a
standalone `[1]` residual block and rank contribution one.  Adding these
disjoint contributions gives an exact lower bound on `rank(R)`, which is
compared with (3).  The same check is made in the other direction.

## 4. Both bins cannot be saturated

Suppose both complete bins have six vectors and have a nonempty light part.
Their square cross-Gram matrix is orthogonal.  Equation (2), together with
the no-cancellation observation above, says two source rows cannot share a
nonzero target column.  Applying this in both directions makes the light
cross block a signed permutation.

Let `k=|Z0|`.  Every matched light pair therefore has `u_y=+u_x` or
`u_y=-u_x`.  For a required opposite-ray edge, (1) becomes

```text
u_x.u_y = c_x c_y/3.                                 (5)
```

For a plus match, equal coordinate sums give `c_x=c_y`, so `c_x^2=3`.
The value `c=-sqrt(3)` is impossible because

```text
|sum_i u_i| = sqrt(3)+1 > sqrt(6)||u||.
```

Thus `c=sqrt(3)`, and the squared all-ones coefficient of the common unit
vector is

```text
(sum_i u_i)^2=(sqrt(3)-1)^2=4-2sqrt(3).
```

For a minus match, coordinate sums give `c_x+c_y=2`, while (5) gives
`c_x c_y=-3`.  Hence the pair is `{3,-1}` and its line has all-ones energy
`4`.

The six rows of either bin form an orthonormal basis.  If `p` light matches
have plus sign and `n` have minus sign, Parseval applied to the all-ones
vector gives

```text
6 = k+p(4-2sqrt(3))+4n,
p+n=6-k.
```

Irrationality first forces `p=0`; the remaining integer equations force
`k=6,n=0`, contradicting the assumed light part.  Thus every both-saturated
light coloring is impossible.

## 5. Exhaustive actual supports and tight same-bin closure

For every surviving coloring, the producer exhausts every actual-support
choice

```text
empty != S_x subseteq D_x                 for a nonzero light vertex,
S_x subseteq D_x and |S_x|>=3             for x in Z0.
```

The same `Z0` support is used in both bins.  Within either orthonormal bin,
two supports cannot intersect in exactly one coordinate, and the assigned
supports must pass Hall matching.

There is one further exact closure.  If a subfamily `T` of orthonormal rows
has

```text
|union_{x in T} S_x|=|T|,
```

then those rows form a basis of that coordinate subspace.  Every other row in
the bin is orthogonal to the whole subspace, so its support must be disjoint
from the union.  The depth-first search applies this rule after every support
assignment.

For a required edge `xy` across the two nonzero opposite-ray bins, (5) has a
nonzero right-hand side.  Therefore `u_x.u_y` is nonzero and the actual
supports must overlap.  If both supports are singletons, overlap makes them
the same singleton.  A nonzero-light unit vector on one coordinate cannot be
`-e_i`, because that has `c=0`; hence both vectors are `+e_i` and both have
`c=2`.  Their dot product is then `1`, whereas (5) requires `4/3`.  Thus no
required opposite-ray edge can have two singleton actual supports.  Both
local rules are applied incrementally during the support search.  They have
zero additional graph marginal after the all-ones layer on this corpus, but
they are retained as exact sanity constraints.

## 6. Opposite-ray edge components

Let `E` be the required-edge bipartite graph between the two nonzero light
sets.  By (2), its connected components are exactly the nonzero blocks of the
cross Gram matrix.  For a component `H` with sides `A_H,B_H`, put

```text
t_H = dim(span u(A_H) intersect span u(B_H)).
```

The rows on either side are orthonormal.  After removing the signs of the
nonzero `c` values, the square Gram operator of a connected block is
nonnegative and irreducible.  Perron--Frobenius makes its top singular value
simple.  Therefore

```text
t_H is 0 or 1.                                       (6)
```

The implementation enumerates both hypothetical values whenever they remain
possible.  A realization supplies one of the enumerated choices.

For `t_H=1`, the common singular vectors have every coefficient nonzero.
Let `U_A` and `U_B` be the unions of the actual supports on the two sides.
Every coordinate occurring in exactly one `A_H` support must consequently
occur in `U_B`, and conversely.  The common line has support contained in

```text
I_H = U_A intersect U_B.
```

All selected common lines and the `Z0` lines are mutually orthogonal, so the
mask family consisting of their `I_H` and actual `Z0` supports must pass Hall
matching.

The span contributed by `H` has dimension

```text
|A_H|+|B_H|-t_H.
```

Different component spans, cross-degree-zero light lines, and `Z0` lines are
mutually orthogonal.  Give every such block its exact rank demand and the
union of its actual supports as an upper coordinate mask.  For every block
subfamily `J`, necessity gives the weighted Hall inequality

```text
sum_{H in J} rank_demand(H)
    <= |union_{H in J} support_mask(H)|.              (7)
```

Production checks the equivalent 64 coordinate-container inequalities.  The
independent verifier checks every block subfamily directly.

For an isolated edge component, selecting `t_H=1` makes its two endpoint
vectors parallel.  The private-coordinate rule above forces their actual
supports to be equal, so the exact support of the common line is known.

## 7. Tagged tight-line all-ones obstruction

Consider any collection of selected `t_H=1` isolated-edge lines together
with any collection of `Z0` lines.  These lines are orthonormal.  Tag the
former as edge lines and retain their exact equal endpoint supports.

For an edge line, write `u_y=epsilon u_x`, where `epsilon` is `+1` or `-1`.
Equation (5) gives `epsilon=c_x c_y/3`.  The calculation in Section 4 applies
without requiring either whole bin to be saturated:

```text
epsilon=+1: c_x=c_y=sqrt(3), line energy=4-2sqrt(3),
epsilon=-1: {c_x,c_y}={3,-1}, line energy=4.
```

A `Z0` line has `c=0`, hence `sum_i u_i=-1` and all-ones energy one.

Suppose a tagged subfamily has as many coordinates in its actual-support
union as it has lines.  Its orthonormal lines form a basis of that coordinate
subspace.  Let `z` be the number of `Z0` lines and let `p_+,p_-` count the two
edge-line signs.  Parseval for the all-ones vector gives

```text
z+p_++p_-
  = z+p_+(4-2sqrt(3))+4p_-.
```

After cancelling `z`, irrationality forces `p_+=0`, and then the integer
equation forces `p_-=0`.  This contradicts the presence of a tagged edge
line.  Therefore no coordinate-tight tagged subfamily containing an edge
line is possible.

Production enumerates tagged subfamilies as integer bit masks.  The verifier
uses independently ordered `itertools.combinations` and explicit coordinate
sets.  A singleton endpoint support is the one-line special case.  On this
623-graph corpus that special case already supplies all 87 new graph
rejections; the more general tight-family form has zero additional corpus
marginal, but is retained as the reusable theorem.

## 8. Quantifiers and rejection semantics

For every candidate graph the complete quantifier order is:

1. enumerate every required `K6` seed;
2. enumerate every eligible actual `Z0`;
3. require the frozen bipartite system on that same `Z0`;
4. enumerate every symmetry-reduced two-light-ray component coloring;
5. apply both directed projection-rank inequalities and the saturated-bin
   identity;
6. enumerate every permitted actual-support assignment jointly for both
   bins, including the edgewise overlap and double-singleton exclusions;
7. enumerate every component-intersection bit `t_H` and apply the private,
   ordinary Hall, weighted Hall, tight same-bin, and tagged all-ones rules.

A seed is impossible only after all choices fail.  One impossible required
seed rejects the graph, because any realization of the graph would realize
that seed.  Every aborted or unresolved computation must remain unresolved;
there is no heuristic failure path in this kernel.

## 9. Source-bound discovery census

The frozen parent boundary is

```text
input graphs                                      623
ordered input hash
  b31aeac00b91d0d843ea51909c64f2c4ca3da49e5d2d33792312a45aa9f25e79
```

The pre-all-ones conjunction rejected 285 graphs and left 338.  Its hashes
were

```text
rejected
  db9b04c6a8a599605af88996a8b300d19886844cde020d63517d6c0e95a7d938
residue
  66e695f040675f3726c51ce883547533097f7d3ec313665ac3bb6db7f4c10217
```

The tagged all-ones theorem rejects 87 further graphs.  The source-bound
expected census is therefore

```text
rejected graphs                                  372
surviving graphs                                 251
rejected hash
  968195801d0311b0e43349e6deb8b08f61e5167b76f26edd8905369a6e574c63
residue hash
  a23da7bd9ae21eb7c078a08a7bce92f22d397d5c71f1e4f51d9cf8e59bbfa907
```

The discovery run used 11 workers, checked the complete ordered input, and
passed all 32 required `K6` seeds of the standard realizable 18-point set.
Its report is discovery evidence only.  The theorem-level artifacts must be
generated after all four package sources are committed.

## 10. Independent checking and provenance

The independent checker imports neither this producer nor either discovery
probe.  It uses the frozen independent tight-same-`Z0` verifier, reconstructs
all 623 inputs, enumerates ray colorings itself, constructs projection
components with an edge-list union-find, searches actual supports in a
different order, checks weighted Hall by all block subsets, and checks tagged
tight families by combinations and explicit coordinate sets.  It recomputes
all graph decisions and compares every rejected seed's complete `Z0` and
coloring partition with the compressed certificate archive.

Focused controls cover the PSD Z-matrix component bound, zero-cross-degree
rank, both-saturated arithmetic, tight same-bin support closure, optional
zeros, opposite-edge local support/value rules, tagged tight-line energy, a
fixed old rejection, a fixed new all-ones rejection, and the 32-seed positive
control.

Reproduction after the source-boundary commit is:

```text
python3 -m unittest -v test_d6_k6_opposite_ray_rank_conjunction.py

python3 d6_k6_opposite_ray_rank_conjunction.py --workers 11

python3 verify_d6_k6_opposite_ray_rank_conjunction.py --workers 11
```

The producer refuses to launch unless the producer, verifier, tests, and this
proof are exact blobs at the recorded `codex/dimension6` commit.  Tracked or
staged worktree dirt blocks launch.  Untracked discovery files are permitted,
and their complete porcelain list is hashed into the report.  The verifier
checks the source hashes, re-reads every blob with
`git show <commit>:<path>`, and validates the recorded porcelain binding.

The trust boundary is Python integer and bit-mask logic, exact frozen SymPy
inertia routines in the parent independent package, SHA-256, the Python/Git
runtime recorded in the reports, and the mathematical lemmas proved here.
There is no IEEE-754 or transcendental-library assumption.  The appearances
of `sqrt(3)` are eliminated symbolically by rationality/irrationality and are
never evaluated numerically.

A surviving graph is not asserted realizable.  This package does not settle
the remaining K6 class or prove `f(6)=18`.
