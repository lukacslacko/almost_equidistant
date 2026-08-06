# Saturating-clique Schur rules for the K7 normalized Gram matrix

This note records a cheap exact refinement found after the K7
cover/support/component-rank layer.  It is separate from the current K7
implementation files.

## Setup

Fix a K7 seed, a candidate zero-factor cover `Z`, and its nonzero remainder
`N`.  The normalized vectors have Gram matrix

```text
K_xx = 1+r_x^2,  r_x != 0;
K_xy = 1          on required graph edges;
K_xy = 0          on candidate graph nonedges.
```

Let

```text
U = min(7-|Z|, nu_N, nu_T-|Z|).
```

Then `K` is positive semidefinite and `rank(K)<=U`.

## Positive-definite clique and perpendicular-neighbour rules

For a required clique `C subseteq N`,

```text
K_C = J + diag(r_x^2 : x in C)
```

is positive definite.  Therefore

```text
omega(G[N]) <= U.
```

If `F=overline(G[N])`, the `F`-neighbours of a vertex form a required clique
because `F` is triangle-free.  Their independent Gram vectors all lie in the
perpendicular space of that vertex, so

```text
Delta(F) <= U-1.
```

More generally, if `F` contains a complete bipartite subgraph on `h`
vertices, the corresponding principal submatrix of `K` is the direct sum of
two matrices `J+positive diagonal`, hence has rank `h`; thus `h<=U`.  The
degree rule is the star case.  On the measured sample this biclique extension
did not add a graph rejection beyond the clique rule.

## Exact Schur consequences when a clique saturates the rank

Suppose `C` is a required clique of size exactly `U`.  Its vectors are a basis
of the Gram space.  Put

```text
M = K_C = J + diag(e_i),       e_i=r_i^2>0,
x_i = 1/e_i,                  X=sum_i x_i.
```

For `y in N\C`, let `S_y subseteq C` be its required-neighbour set in `C` and
let `a_y` be its zero-one indicator.  Since the Schur complement of the basis
block is zero, for all `y,z in N\C`, including `y=z`,

```text
K_yz = a_y^T M^{-1} a_z,
M^{-1} = diag(x_i) - x x^T/(1+X).
```

This yields graph-only necessary conditions with no coefficient search:

1. Every `S_y` is nonempty and proper.  For the empty set the displayed
   quadratic form is zero; for all of `C` it is `X/(1+X)<1`; both contradict
   `K_yy=1+r_y^2>1`.
2. The sets `S_y` are pairwise distinct.  Equal indicator vectors determine
   the same vector in the basis, so `K_yz=K_yy>1`, while every off-diagonal
   entry is zero or one.
3. Every two sets `S_y,S_z` intersect.  If they are disjoint and nonempty,
   the inverse formula is strictly negative, again impossible for a zero-one
   off-diagonal entry.
4. If `yz` is a candidate graph nonedge, `S_y` and `S_z` are incomparable.
5. If `yz` is a required edge, `S_y union S_z` is a proper subset of `C`.

For the last two claims, partition `C` into the four Venn regions and write
`alpha,beta,gamma,delta` for the sums of the positive `x_i` over respectively

```text
S_y intersect S_z,  S_y\S_z,  S_z\S_y,  C\(S_y union S_z).
```

Then

```text
a_y^T M^{-1} a_z
  = (alpha(1+delta)-beta gamma)/(1+X).
```

For a nonedge this equals zero, so

```text
alpha(1+delta)=beta gamma.
```

Pairwise intersection gives `alpha>0`, hence both differences are nonempty.
For a required edge the value equals one, which reduces to

```text
alpha delta-beta gamma = 1+beta+gamma+delta;
```

therefore `delta>0`.

These are only necessary conditions.  Passing them does not assert that one
positive vector `(x_i)` satisfies all Schur equations.

## Provisional sample coverage

The current deterministic 512-graph K7 residue sample has 145 survivors after
the existing joint cap-three/support/subspace/component-inertia layer.
Existentially quantifying over every remaining cover gave:

| added rule | rejected among 145 |
|---|---:|
| positive-definite clique cap | 25 |
| perpendicular degree cap alone | 6 |
| clique plus degree | 25 |
| saturating-clique Schur set rules, including the caps | **59** |

Thus the Schur set rules add 34 sample rejections beyond the clique cap and
leave 86 of the 145 sample survivors.  The degree-only graph rejections were
already contained in the clique set, although it did reject extra individual
covers.  The reference measurement examined 549 post-layer cover candidates;
the total 103.8-second wall time was dominated by recomputing the existing
Python zero-forcing layer.  The new work itself is only clique enumeration and
small bit-set comparisons on at most twelve vertices.

The 59 rejected sample indices were:

```text
768969 3396330 1164960 502658 3690509 773478 3436402 3819584
824661 3336291 2565924 862073 762981 31231 3713158 452927
3958222 2644666 3444605 3954034 2573244 861960 2720951 2374276
1796727 3094211 3709819 398872 1167251 837353 1152953 2760238
3704179 1311853 3707945 3834364 3582223 348117 1757396 1796724
364900 636253 181646 278338 3174932 817524 1571972 2228810
2899538 2477214 61094 1049569 917628 2925901 3649340 2238182
156535 3717479 1431273
```

These counts are sample measurements, not a full-population claim.  They
should be rechecked in the production existential-cover pass and against the
realizable positive controls before promotion.
