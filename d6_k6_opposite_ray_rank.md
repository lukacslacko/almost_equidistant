# Opposite-light-ray projection rank in the K6 coordinates

This note proves the exact filter in `d6_k6_opposite_ray_rank.py`.  Its
independent checker is `verify_d6_k6_opposite_ray_rank.py`.  The filter is a
standalone necessary condition: it does not rely on numerical optimization,
the current bipartite-component kernel, or the stronger same-coloring
actual-support conjunction.

## 1. K6 identities and optional-zero semantics

Fix a required unit `K6` with vertices `q_1,...,q_6`.  For each of the
thirteen outside points `x`, the exact K6 coordinates give a defect vector
`u_x in R^6` and Lorentz factor

```text
ell_x = (c_x,z_x),
c_x   = 1 + sum_i u_xi,
z_x^2 = c_x^2 + 6 - 6 ||u_x||^2,
M_xy  = ||x-y||^2-1
      = (c_x c_y-z_x z_y)/6 - u_x.u_y.              (1)
```

Let `D_x` be the set of seed coordinates at which `u_x` is *allowed* to be
nonzero.  It is obtained from candidate nonedges between `x` and the seed.
It is only an upper support: a coordinate in `D_x` may still be zero.

Every graph in the corpus has `alpha(G)<=2`.  Therefore, if `xy` is a
candidate graph nonedge, then

```text
D_x intersect D_y = empty.                            (2)
```

Indeed, a shared coordinate `i` would say that `q_i x`, `q_i y`, and `xy`
are all candidate nonedges, an independent triple.  Equation (2) implies
`u_x.u_y=0`; it does **not** prescribe the distance of `x,y`.

Put

```text
Z0 = {x : ell_x=0}.
```

For `x in Z0`, (1) gives `||u_x||=1` and `sum_i u_xi=-1`.  Its actual
support has size at least three.  Support zero is incompatible with the sum;
support one gives `u_x=-e_i`, the existing seed point; and for support two,
`a+b=-1`, `a^2+b^2=1` imply `ab=0`, reducing to support one.  The vectors
`{u_x:x in Z0}` are orthonormal, so `|Z0|<=6`, and their allowed masks must
admit a matching into the six coordinates.  Enumerating every subset of the
vertices with `|D_x|>=3`, up to size six, therefore includes the actual
`Z0` of every realization.

## 2. Nonbipartite Lorentz components force light rays

Let `L` be the graph on the outside vertices in which `xy` is an edge exactly
when

```text
xy is a required unit edge and D_x intersect D_y = empty.
```

For an edge of `L`, (1), (2), and `M_xy=0` give

```text
<ell_x,ell_y>_L = c_x c_y-z_x z_y = 0.               (3)
```

Delete `Z0` from `L`.  In a connected component, Lorentz orthogonality
alternates between one projective line and its orthogonal line.  An odd cycle
forces those lines to coincide.  In the nondegenerate Lorentz plane
`R^(1,1)`, the only self-orthogonal projective lines are

```text
z=+c  and  z=-c.                                      (4)
```

Consequently, every connected nonbipartite component of `L-Z0` lies wholly
on one of the two light rays.  Different components may choose either ray.
The implementation enumerates all choices, fixing the first component to
one ray only to quotient the global exchange of the two ray names.

Within one light ray all Lorentz products vanish.  A required pair is
therefore defect-orthogonal by (1), and a candidate nonedge is
defect-orthogonal by (2).  The diagonal identity and (4) give `||u_x||=1`.
The same statements hold between a light-ray point and `Z0`, and inside
`Z0`.  Thus each ray bin together with the common set `Z0` is an orthonormal
family in `R^6`.  In particular, its allowed masks must pass Hall matching.

## 3. Projection onto the opposite ray

Fix one enumerated ray assignment.  Let `B` be one complete orthonormal bin,
including `Z0`, and put

```text
b = |B|.
```

Let `A` be the nonzero light-ray vertices assigned to the opposite ray; the
shared `Z0` vertices are excluded from `A`.  Write `X` for the `A`-by-`B`
cross Gram matrix of the defect vectors.  The columns indexed by `Z0` are
zero.  After naming the rays `z=+c` on `A` and `z=-c` on `B-Z0`, (1) and (2)
give, for `i in A`, `j in B-Z0`,

```text
X_ij = c_i c_j/3,  if ij is a required graph edge,
X_ij = 0,          if ij is a candidate graph nonedge.       (5)
```

Every nonzero light factor has `c_i!=0`: if `c_i=0`, (4) would also give
`z_i=0`, contrary to `i notin Z0`.

Project the orthonormal vectors indexed by `A` off the span of the
orthonormal bin `B`.  Their residual Gram matrix is

```text
R = I_A-X X^T >= 0,
rank(R) <= 6-b.                                        (6)
```

For distinct `i,k in A`, (5) gives

```text
R_ik = -(c_i c_k/9)
       sum_{j in B-Z0: ij,kj required} c_j^2.          (7)
```

Congruence by the invertible diagonal matrix
`D=diag(1/c_i:i in A)` preserves positive semidefiniteness and rank.  The
matrix

```text
S = D R D
```

has nonpositive off-diagonal entries.  Its strict negative off-diagonal
support is exactly the graph `C_B(A)` joining two vertices of `A` when they
share a required neighbour in `B-Z0`.  There is no cancellation in (7),
because every term is a strictly positive square.

## 4. The irreducible PSD Z-matrix rank lemma

Let `H` be a nontrivial connected component of `C_B(A)`.  The corresponding
principal block `S_H` is a symmetric PSD irreducible Z-matrix.  Choose
`tau >= max_i (S_H)_ii`.  Then

```text
T = tau I-S_H
```

is entrywise nonnegative and irreducible: every edge of `H` becomes a
strictly positive off-diagonal entry.  Perron--Frobenius says the largest
eigenvalue of `T` is simple.  Since the eigenvalues of `T` are
`tau-lambda_i(S_H)`, the smallest eigenvalue of `S_H` is simple.  Hence

```text
nullity(S_H) <= 1,
rank(S_H) >= |H|-1.                                   (8)
```

This also audits the possible zero-diagonal edge case.  In any PSD matrix,
`S_ii=0` forces the whole `i`th row and column to vanish (for example from
every `2`-by-`2` principal minor).  Thus a vertex in a nontrivial strict
off-diagonal component cannot have zero diagonal.  The Perron--Frobenius
argument above remains valid whether or not the whole block is singular.

There is a further exact contribution from some isolated vertices of
`C_B(A)`.  If `i in A` has no required neighbour at all in `B-Z0`, then its
row of `X` is identically zero by (5).  Its row and column in `R` form the
scalar block `[1]`, so it contributes rank one.  An isolated vertex with a
nonempty private target neighbourhood is conservatively assigned rank zero;
its residual diagonal could vanish.

Different components are separate blocks of `S`.  Combining (6)--(8), every
realization must satisfy the directed inequality

```text
#{i in A : i has no required neighbour in B-Z0}
  + sum_{nontrivial components H of C_B(A)} (|H|-1)
  <= 6-b.                                               (9)
```

The filter applies (9) in both directions between the two ray bins.

For example, if `b=6`, the right side is zero: every opposite vertex must
have a target neighbour and no two may share one.  If `b=5`, the common-
neighbour support may contribute rank at most one.  Formula (9) unifies these
special cases without numerical rank tests.

## 5. Exhaustive quantifiers and sound rejection

For every graph, production performs the following exact quantifiers.

1. Enumerate every required `K6` seed.
2. For that seed, enumerate every eligible `Z0` of size at most six.
3. Reject an unmatchable `Z0` by the old orthonormal-support Hall condition.
4. Find every connected nonbipartite component of `L-Z0`.
5. Enumerate all `2^(k-1)` assignments of the `k` components to the two rays
   (one assignment if `k=0`).
6. Require both orthonormal bins to pass allowed-support Hall matching and
   require (9) in both directions.

A seed is impossible only if every `Z0` and every ray assignment fails.  One
impossible required seed rejects the graph.  The certificate archive stores
every `Z0` row and every symmetry-reduced coloring for each rejected seed.

Optional zeros remain optional throughout.  Hall matching is only a
necessary condition on upper supports.  A candidate graph nonedge is used
only through the alpha-two deduction (2), never by requiring its geometric
distance to differ from one.

## 6. Certified K6-623 result and controls

The frozen input is the ordered 623 survivors of the independently verified
strongest same-`Z0` boundary.  The standalone projection-rank filter rejects

```text
171 of 623 graphs,
```

leaving 452.  The stable hashes are

```text
ordered input indices
  b31aeac00b91d0d843ea51909c64f2c4ca3da49e5d2d33792312a45aa9f25e79
rejected indices
  390bca6d3661ba4ab497d2ba23400d217d759702b071ecbf0438a9d03301bdf6
ordered residue indices
  cb62e002449e9be802ae37d66214904168de9528777d582cf6ddf28b350a5656
```

The full producer census records 15,489 required K6 seeds, 25,769 `Z0`
choices, and 36,477 symmetry-reduced ray colorings before graph/seed early
exit.  The rejected-seed archive contains all 9,086 `Z0` rows and all 17,995
ray-coloring rows below them.

The known realizable 18-point half-cube-plus-two-poles construction passes
all 32 of its required K6 seeds.  Focused tests also cover the rank lower
bound on paths and disconnected components, the formerly missed
zero-cross-degree rank-one block, the sharp `b=5` common-neighbour boundary,
a two-component coloring in which only the second assignment passes, a fixed
corpus rejection, deterministic compression, and producer/verifier
agreement.

The independent verifier imports neither the producer nor the producer's
Lorentz/graph helpers.  It reconstructs the graph instances itself and uses
Hall's all-subsets union-cardinality criterion instead of matching DFS.  It
recomputes all 623 decisions and every rejected certificate.  No
floating-point arithmetic or transcendental function occurs in a rejection.

Reproduction commands are

```text
python3 -m unittest -v test_d6_k6_opposite_ray_rank.py

python3 d6_k6_opposite_ray_rank.py --workers 11

python3 verify_d6_k6_opposite_ray_rank.py --workers 11
```

The theorem-level producer refuses to launch unless all four package sources
(producer, verifier, tests, and this proof) are exact blobs at the recorded
`codex/dimension6` commit.  Tracked or staged worktree dirt blocks launch;
untracked discovery files are permitted and their complete porcelain list is
hashed into `execution.git`.  The verifier checks that binding, reads every
recorded source with `git show <commit>:<path>`, and compares its SHA-256 with
the report before replaying a mathematical decision.  Thus output artifacts
made before the source-boundary commit are discovery-only even if their
mathematical contents happen to agree with the later official run.

The official postcommit report, compressed certificate, and verification
JSON bind their own exact artifact and raw-payload hashes.  A survivor is
only a survivor of this necessary condition; it is not asserted realizable.
A separate stronger layer couples this bound to actual-support feasibility on
the identical ray coloring.
