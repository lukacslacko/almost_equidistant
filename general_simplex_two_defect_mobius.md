# Two-defect Mobius dynamics after a full simplex seed

This note generalizes the dimension-six two-defect cycle calculation.  It is
independent of the active dimension-six production package.

## Transition in dimension `d`

Fix a regular unit `K_(d+1)` seed in `R^d` and put

```text
n=d+1,   r=sqrt(n),   t_x=s_x+1,   w_i(x)=r u_i(x)/t_x.
```

Here "exact two-defect" refers to the actual support
`{i:u_i(x)!=0}`, not merely the candidate's allowed non-unit incidences.
Such a point has `t_x!=0`: if `t_x=0`, the simplex quadratic and
`u_i+u_j=-1` give `u_i^2+u_j^2=1`, hence `u_i u_j=0`, contrary to exact
two-defect support.  If its two normalized components are `A,B`, the same
quadratic gives

```text
(A-r)(B-r)=d/2.
```

Crossing a two-defect support type therefore applies

```text
phi(A) = (r A-(d+2)/2)/(A-r).
```

At a seed coordinate shared by two distinct support types, the two points
must be unit and their sole normalized product is one.  Crossing a type and
then moving to the next type applies

```text
T(A)=1/phi(A)=(A-r)/(r A-(d+2)/2).
```

Its projective matrix and invariants are

```text
M_d = [ 1   -sqrt(d+1) ]
      [ sqrt(d+1)  -(d+2)/2 ],

trace(M_d) = -d/2,
det(M_d)   =  d/2,
discriminant = trace^2-4 det = d(d-8)/4.
```

Equivalently, form the simple **type graph** on the `d+1` seed coordinates,
putting in edge `{i,j}` when at least one outside point has exact defect type
`{i,j}`.  Choosing one point on each edge of a simple type-graph cycle and
orienting the cycle gives the recurrence `A_(k+1)=T(A_k)`.  Edge
multiplicities are not represented in this simple graph.

Thus `T` is elliptic for `1<=d<8`, parabolic for `d=8`, and hyperbolic for
`d>8`.

## Finite-order classification below dimension eight

For `d<8`, normalize the complex-conjugate eigenvalues to
`exp(+/- i theta)`.  Then

```text
2 cos(theta) = trace/sqrt(det) = -sqrt(d/2),
4 cos(theta)^2 = d/2.
```

If the projective action had finite order, `theta/pi` would be rational.
Then `2 cos(theta)` would be an algebraic integer, so its rational square
`d/2` would be a rational algebraic integer and hence an integer.  For
`1<=d<8`, this leaves only `d=2,4,6`.  Direct exact powers give

```text
d=2: M_d^3 = I,       projective order 3,
d=4: M_d^4 = -4 I,    projective order 4,
d=6: M_d^6 = -27 I,   projective order 6.
```

Every nonidentity elliptic power has no real projective fixed point.  A
simple cycle of `ell` exact two-defect support types requires a real fixed
point of
`T^ell`.  Consequently:

```text
d=2: cycle lengths are divisible by 3,
d=4: cycle lengths are divisible by 4,
d=6: cycle lengths are divisible by 6,
d=1,3,5,7: the exact two-defect type graph is a forest.
```

The forest conclusion for dimension seven is immediate and potentially
useful: after a required K8 seed, any forced cycle among exact two-defect
types is already an exact obstruction, with no coefficient search.

## The parabolic dimension-eight case

For `d=8`, `r=3` and

```text
M_8 = [1 -3; 3 -5] = -2 I + N,   N^2=0.
```

The nontrivial parabolic action has the unique real projective fixed point
`A=1`.  Therefore a simple two-defect type cycle is possible only when every
normalized endpoint encountered around it equals one.  Since `phi(1)=1`,
this value propagates from the cycle through every occupied type in its
connected type-graph component.  Converting back to defects gives

```text
t=3,   u_i=u_j=1
```

for every such point.

Hence dimension eight does not forbid cycles outright, but it makes every
cyclic type-graph component coefficient-rigid.  In particular, two copies
of any occupied support type in such a component would be forced to the same
point and violate distinctness.  The rigid values can also be substituted
directly into constraints involving larger-support points.

As a dimension-independent side condition, three pairwise disjoint
two-defect types have pairwise normalized dot product zero.  The distance
identity (and `t_x t_y!=0`) then makes the three corresponding outside
points pairwise non-unit, impossible in an almost-equidistant set.  Thus the
type graph has matching number at most two.
