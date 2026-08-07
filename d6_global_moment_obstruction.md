# Global moment and bordered-rank audit in dimension six

This note records exact consequences and a negative boundary for the global
matrix approach.  It does not reject any candidate by itself.

For a 19-point candidate put

```text
M_ii = -1,
M_ij = ||p_i-p_j||^2-1,
A    = M+I.
```

The off-diagonal support of `A` is contained in the triangle-free complement
of the required unit graph.  Allowed entries may be zero.  A realization in
`R^6` has `rank(M)<=8`, at most one positive eigenvalue, and
`-CMC/2` positive semidefinite of rank at most six, where
`C=I-J/19`.

## Trace identity and the positive eigenvalue

Triangle-free support gives `tr(A)=tr(A^3)=0`, regardless of which allowed
entries vanish.  If the nonzero eigenvalues of `M` are `eta_i`, then

```text
sum_i eta_i                         = -19,
sum_i eta_i^2                       = 19+tr(A^2),
sum_i eta_i^3                       = -19-3tr(A^2),
sum_i eta_i^2(eta_i+3)              = 38.             (1)
```

If `M` were negative semidefinite, write `eta_i=-x_i`, `x_i>0`.  Every term
in (1) would be

```text
x_i^2(3-x_i) <= 4.
```

There are at most eight terms, so the right side would be at most 32, not 38.
Thus `M` has a positive eigenvalue; the Euclidean inertia bound makes it
unique.

A required unit `K_k` is the principal block `-I_k`, so it supplies at least
`k` negative eigenvalues.  Together with the positive eigenvalue this proves
`rank(M)>=k+1`.  In particular, every K7 case has exact inertia

```text
(n_+(M),n_-(M),n_0(M)) = (1,7,11).
```

Interlacing with its `-I_7` principal block also puts every negative
eigenvalue of `M` at most `-1`.

These moment facts alone are not contradictory.  For example, let `lambda`
be the real root in `(3,4)` of

```text
16 lambda^3 - 11 lambda^2 - 121 lambda - 264 = 0
```

and put `y=(lambda+11)/7`.  The formal spectrum

```text
A:  lambda, 1 (multiplicity 11), -y (multiplicity 7)
```

has zero first and third moments and gives `M=A-I` the required rank and
inertia.

There is an even stronger exact matrix counterexample.  Let `C_16` be the
Clebsch graph, realized as the Cayley graph on `F_2^4` with connection set

```text
{1000,0100,0010,0001,1111}.
```

Its adjacency spectrum is `5,1^10,-3^5`.  Therefore

```text
A = Adj(C_16) direct-sum Adj(K2) direct-sum [0]
```

has triangle-free support and spectrum

```text
5, 1^11, 0, -1, -3^5.
```

Consequently `M=A-I` has spectrum `4,0^11,-1,-2,-4^5`, hence rank eight and
inertia `(1,7,11)`.  The support graph has independence number seven.  It is
not Euclidean: the vector which is `1` on the Clebsch component, `-8` on both
vertices of the K2, and `0` on the isolate is perpendicular to the all-ones
vector but has quadratic value `x^T M x=64>0`.  This example is a mandatory
negative control for any purported trace/rank/inertia-only obstruction.

## What the centered condition adds

For rank-eight `M` of inertia `(1,7)`, conditional negative semidefiniteness
and centered rank at most six force

```text
1 in range(M),
1^T M^dagger 1 = 0.                                  (2)
```

Conversely, (2) gives the required lightlike hyperplane restriction.  An
equivalent rank statement is

```text
rank [[M,1],[1^T,0]] <= 8.                            (3)
```

With a required K7, the bordered principal block on infinity plus the K7 is
a fixed invertible 8-by-8 core; its Schur complement is exactly the existing
K7 simplex-defect system.  With only a K6, infinity plus the K6 is a fixed
invertible 7-by-7 core; (3) says that the outside Schur complement has rank at
most one, exactly the existing K6 normal/Lorentz system.  Thus the bordered
formulation is clean but does not bypass those seed-coordinate equations.

For a K6 seed define

```text
c_x = 1+sum_i u_i(x),
g_xy = 6 u(x).u(y)-c_x c_y,
h_x = 6-6||u(x)||^2+c_x^2 = z_x^2.
```

On required outside edges, bordered rank one gives

```text
g_xy = -z_x z_y,
g_xy^2 = h_x h_y,
g_ab g_ac g_bc <= 0,
g_ab g_cd-g_ac g_bd = 0.
```

Pure tetrad duals have an exact optional-zero no-go.  Substitute every
allowed coordinate as `u=t-1`, so distinctness gives `t>0`.  At the common
positive point `t=1`, all defects are zero, every required-edge `g=-1`, and
all tetrads vanish.  Hence every polynomial combination of tetrads, including
monomial multiples, vanishes at `t=1`; it cannot be a nonzero
coefficientwise-one-sign polynomial.  A complete 831-graph direct scan had
zero hits, as this argument requires.  The diagonal equations or the global
distinctness multiplicity must be included to escape this no-go.

## Distinctness escape hatch and virtual K7 promotion

An empty K6 defect has

```text
u=0, c=1, z=+sqrt(7) or -sqrt(7).
```

There are exactly two such geometric points and at most one of each sign.
This yields the separate empty-support-budget layer documented in
`d6_k6_empty_support_budget.md`.

It also suggests a stronger finite branch.  If a vertex is declared empty,
it is an actual common unit neighbour of all six seed vertices.  One may
soundly augment the required graph with those six edges, obtaining a virtual
required K7, and apply K7 kernels that depend only on required unit edges and
distinctness (not on corpus minimality).  Two empty vertices must be the
opposite points, so their mutual pair is nonunit and cannot be a required
edge.  A complete K6 disjunction can therefore branch over no empty vertex,
one empty vertex, or an eligible candidate-nonedge pair, promoting every
nonempty branch to the corresponding support constraints and every empty
branch to one or two virtual K7 systems.  This idea is not yet a certified
filter; every reused K7 API must first be audited for assumptions about the
original corpus or maximality.
