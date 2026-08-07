# Exact K6 coordinates and the normal-inertia refinement

This note derives the full coordinate system around a centered required unit
`K6` in `R^6`, including the one remaining normal coordinate.  It then proves
the exact bipartite-component inertia screen implemented in
`d6_k6_normal_inertia.py`.  The prototype is layered on the 1,097 exact
survivors of the earlier K6 Lorentz/support/zero-forcing campaign; it does not
alter those frozen references.

## Coordinate identities

Let `q_1,...,q_6` be the centered regular unit simplex in its
five-dimensional span `W`, and choose a unit normal `e`.  Thus

```text
sum_i q_i = 0,
q_i dot q_i = 5/12,
q_i dot q_j = -1/12  (i != j),
sum_i q_i q_i^T = (1/2) I_W.
```

For an outside point `x`, define

```text
u_i(x) = ||x-q_i||^2-1,
s_x    = sum_i u_i(x),
c_x    = s_x+1,
x      = x_W+h_x e,
z_x    = sqrt(12) h_x.
```

Summing the six defect equations and using the tight-frame identity gives

```text
||x||^2 = (2s_x+7)/12,
x_W     = -sum_i u_i(x) q_i,
||x_W||^2 = (1/2)||u(x)||^2-s_x^2/12.
```

Consequently the exact diagonal and pair identities are

```text
z_x^2 = c_x^2+6-6||u(x)||^2,                         (K6-diag)

M_xy  = ||x-y||^2-1
      = (c_x c_y-z_x z_y)/6-u(x) dot u(y).           (K6-pair)
```

Writing

```text
ell_x=(c_x,z_x),
<ell_x,ell_y>_L=c_x c_y-z_x z_y
```

exhibits the remaining coordinate as a two-dimensional Lorentz factor.  The
six-simplex defect coordinates plus `z_x` reconstruct `x` uniquely.

## Actual supports and distinctness

Let

```text
D_x = {i : xq_i is a candidate nonedge},
S_x = {i : u_i(x) != 0}.
```

Only `S_x subseteq D_x` is known.  An allowed coordinate may be zero, and a
candidate nonedge may have unit distance.  In particular, even though a
K6-only candidate has `D_x` nonempty, its actual support may be empty.

The reconstruction gives the precise collision rules:

```text
x=y       iff (u(x),z_x)=(u(y),z_y),
x=q_i     iff u(x)=-e_i and z_x=0.
```

If `S_x` is empty, then `c_x=1` and `z_x=+-sqrt(7)`.  These are exactly the
two common unit neighbours of the K6.  Hence there are at most two distinct
empty-support points, at most one at either sign.  The two signs have
`M_xy=4/3`, or squared distance `7/3`.  This exact fact does not by itself
reject an abstract candidate, because candidate nonedges remain optional.

If `ell_x=0`, then

```text
sum_i u_i(x)=-1,  ||u(x)||^2=1,  z_x=0.
```

Support one forces `u=-e_i`, a seed collision.  For support two, the sum and
norm equations give `2u_i u_j=0`, reducing to the same collision.  Therefore
every zero factor has actual support at least three.  Zero-factor defect
vectors are orthonormal, so their allowed masks must admit a matching into
the six seed coordinates and there can be at most six of them.  These are the
exact eligibility conditions used when enumerating `Z0={x:ell_x=0}`.

For a singleton allowed mask `D_x={i}`, write `u(x)=t_x e_i`; `t_x=0` is
still permitted.  Distinctness from `q_i` gives `c_x=t_x+1 != 0`, and

```text
z_x^2=c_x(12-5c_x).
```

Thus the projective Lorentz ratio `r_x=z_x/c_x` uniquely determines the
point on that coordinate line:

```text
c_x=12/(r_x^2+5),
t_x=(7-r_x^2)/(r_x^2+5),
z_x=12r_x/(r_x^2+5).
```

Two singleton-mask vertices of the same type on the same Lorentz projective
line would therefore collide.  This is sound, but the deterministic sample
has no K6 seed with a repeated singleton allowed-mask type, so it has zero
sample marginal value and is not in the prototype filter.  The basic
cardinality bound for one singleton type is also already subsumed in this
campaign: these vertices and the omitted seed vertex lie in a K5 link; three
outside vertices would be forced to form a required triangle by
`alpha(G)<=2`, producing a previously rejected required K8 with the other
five seed vertices.

For reference, either common K6 neighbour can be used as a virtual seventh
simplex vertex

```text
r_plus  = +sqrt(7/12)e,  w_plus(x)  = (c_x-sqrt(7)z_x)/6,
r_minus = -sqrt(7/12)e,  w_minus(x) = (c_x+sqrt(7)z_x)/6.
```

This converts the formulas to the full K7 simplex identities, but the
virtual apex is not necessarily a point of the almost-equidistant set.
Therefore overlap in its seventh coordinate cannot be used as an
overlap-forces-unit rule.  Treating that coordinate as an actual seed point
would be unsound.

## The new exact inertia bound

For a fixed possible zero-factor set `Z0`, form `L` on the remaining outside
vertices.  A required edge `xy` lies in `L` when `D_x` and `D_y` are
disjoint.  Then `u(x) dot u(y)=0` and `(K6-pair)` gives

```text
<ell_x,ell_y>_L=0.
```

The disjoint-support statement is also valid for every candidate nonedge:
if `xy` is a candidate nonedge and `D_x intersect D_y` contained `i`, then
`q_i,x,y` would be an independent triple, contrary to `alpha(G)<=2`.

Let `C=A union B` be a connected bipartite component of `L-Z0`.  In the
generic case, its Lorentz factors have the form

```text
ell_x=lambda_x ell_A  (x in A),
ell_y=lambda_y ell_B  (y in B),
<ell_A,ell_B>_L=0,
```

with all scalars nonzero.  The two Lorentz norms are nonzero and have
opposite signs.  Put

```text
H_A = I+Adj(G[A]),  H_B=I+Adj(G[B]).
```

If `D_A=diag(lambda_x:x in A)` and
`kappa_A=<ell_A,ell_A>_L/6`, the defect-vector Gram matrix on `A` is exactly

```text
Gram(u(A)) = I+kappa_A D_A H_A D_A.                 (Gram-A)
```

Indeed the diagonal entries follow from `(K6-diag)`.  A required off-diagonal
pair uses `(K6-pair)`, while a candidate nonedge has disjoint allowed masks
and hence zero dot product.  The analogous identity holds on `B`, with the
opposite sign of `kappa`.

Suppose first that `kappa_A>0`.  On the nullspace of `(Gram-A)`, the quadratic
form `D_A H_A D_A` is negative definite.  Its nullity is therefore at most
`n_-(H_A)` by Sylvester inertia.  Hence

```text
rank u(A) >= |A|-n_-(H_A).
```

When `kappa_A<0`, the same argument uses a positive-definite subspace and
gives `rank u(A)>=|A|-n_+(H_A)`.  The sign is opposite on `B`.  Moreover the
`A` and `B` defect spans are orthogonal to one another and to the
`|Z0|`-dimensional orthonormal zero-factor span.  Thus one of the two unknown
Lorentz orientations must satisfy

```text
R_+ = |Z0|+|A|-n_-(H_A)+|B|-n_+(H_B) <= 6,
R_- = |Z0|+|A|-n_+(H_A)+|B|-n_-(H_B) <= 6.
```

Equivalently, the exact necessary condition is

```text
min(R_+,R_-) <= 6.                                  (K6-normal-inertia)
```

If the component direction is lightlike, the orthogonal line is itself and
all defect vectors in `C`, together with `Z0`, are orthonormal.  The stronger
condition `|Z0|+|A|+|B|<=6` then holds, and it implies
`(K6-normal-inertia)`.  The screen is therefore necessary in both cases.

Different connected components cannot be added: their Lorentz projective
directions may differ and their defect spans need not be orthogonal.  The
implementation correctly tests each component separately.  It rejects one
K6 seed only after every eligible, support-matchable `Z0` fails.

All inertias of the integer matrices `I+Adj` are computed by rational
symmetric congruence elimination.  A one-by-one nonzero diagonal pivot or a
two-by-two `[[0,a],[a,0]]` pivot is removed by an exact Schur complement.
There is no eigenvalue tolerance and no floating-point decision.

## Controls and deterministic sample

`d6_k6_normal_inertia_test.py` checks:

- known exact inertias, 160 deterministic random congruences with known
  signatures, and an independent SymPy characteristic-root-isolation check
  of all 1,099 labeled graphs through five vertices;
- a structured pair of side patterns for which either Lorentz orientation
  needs seven dimensions although the older arbitrary-pattern zero-forcing
  lower bounds total only four;
- all 32 K6 seeds in the known realizable 18-point configuration;
- fixed corpus witness `3279`, including exhaustive failure of all 127
  eligible `Z0` choices for seed `[4,11,13,14,15,18]`;
- the frozen sample totals when the sample report is present.

The sample is the 64 records with bytewise-smallest
`SHA256(decimal corpus index)` among the 1,097 exact prior K6 survivors.  The
serial exact run took 0.295 seconds and rejected eight graphs:

```text
95658  420602  502025  1196109
1334813  1860514  2097618  2343525
```

It checked 2,080 K6 seeds, found 13 impossible seeds, enumerated 3,340 `Z0`
subsets, and performed 17,987 bipartite-component checks.  This is exact
marginal sample coverage over the old residue, not a full-residue claim.

Commands:

```text
python3 -m unittest -v d6_k6_normal_inertia_test.py
python3 d6_k6_normal_inertia.py --sample-size 64 \
  --output d6_k6_normal_inertia_sample.json
```

Frozen inputs and sample artifact:

```text
d6_k6_bipartite_rank_input.json
  dfee979e6ca2759e80cf803c3054c8954d98543875d467644b91ee434c22e845
d6_k6_bipartite_rank_report.json
  ede04a71a7f9b6901bbc1e4c198a117797bfd90d4d1366912ca9efd1ca559bdc
d6_k6_normal_inertia_sample.json
  3bc04d4541e11597302631fb30466af19d4c8b18333d4dffabca2f059b06128d
```
