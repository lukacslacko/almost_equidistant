# Exact one- and two-defect value constraints after a K7 seed

This note proves the sparse-value refinement implemented in
`d6_k7_small_support_value.py`.  Unlike the preceding labeled-support layer,
these rules use the actual values in the K7 simplex quadratic and bilinear
identities.  They are still finite graph/support tests.

## Normalization

Fix a regular unit K7 seed and write, for an outside point `x`,

```text
u_i(x) = ||x-q_i||^2-1,   s_x = sum_i u_i(x),   t_x=s_x+1.
```

The exact simplex identities are

```text
7 sum_i u_i(x)^2-t_x^2=7,
7 u(x) dot u(y)=t_x t_y                 (xy a required unit edge).
```

For a point with exactly two defects, `t_x` cannot vanish.  Indeed, if its
two nonzero values are `a,b`, then `t_x=0` would give `a+b=-1` and
`a^2+b^2=1`, hence `ab=0`, a contradiction.  We may therefore put

```text
rho=sqrt(7),             w_i(x)=rho u_i(x)/t_x.
```

A required edge becomes `w(x) dot w(y)=1`.

## The two-defect Mobius map

Suppose the exact support is `{i,j}` and abbreviate its two normalized
components by `A,B`.  Dividing the point quadratic by `t_x^2` gives

```text
A^2+B^2 = 1+(rho-A-B)^2,
(A-rho)(B-rho)=3.
```

Thus crossing a two-defect support type sends one endpoint value to the
other by the involution

```text
phi(A) = rho + 3/(A-rho) = (rho A-4)/(A-rho).
```

If two distinct support types meet in exactly coordinate `i`, their points
are both non-unit from `q_i`.  The almost-equidistant condition therefore
forces the two points to be a unit pair.  Their only common nonzero
coordinate gives `A C=1`.  Moving through a support type and then to the next
type at a shared coordinate consequently applies

```text
T = R o phi,       R(A)=1/A,
T(A)=(A-rho)/(rho A-4).
```

The projective matrix of `T` is

```text
M = [ 1    -rho ]
    [ rho   -4  ].
```

Exact multiplication in `Q(rho)`, using `rho^2=7`, gives

```text
M^6 = -27 I.
```

For powers 1 through 5, the discriminants of their real fixed-point
quadratics are respectively

```text
-3, -27, -108, -243, -243.
```

Hence no nonidentity power `T^k`, `k mod 6 != 0`, has a real fixed point.
Poles cause no exception: every normalized component belonging to an actual
point is finite.

Make a graph on the seven seed coordinates, with one edge `{i,j}` for every
outside point whose exact defect support is `{i,j}` (parallel edges are
collapsed for this statement).  Traversing a simple cycle of length `ell`
applies `T^ell` and must return to its starting component value.  We have
proved the exact rule

```text
every simple cycle in the two-defect type graph has length divisible by 6.
```

In particular, triangles, quadrilaterals, pentagons, and seven-cycles are
impossible.  A six-cycle is deliberately retained: the displayed algebra
does not reject it.

## One-defect endpoints

The one-defect point equation has roots `u_i=-1,4/3`.  The first is the seed
point `q_i` itself and is excluded by distinctness.  Thus a one-defect point
of type `i` is unique and has normalized component

```text
w_i=4/rho.
```

There is at most one point of each one-defect type, and at most two
one-defect points altogether: three distinct types would be three pairwise
non-unit points.

If a two-defect point is incident with a one-defect coordinate, its component
there is forced by the unit equation to

```text
h=rho/4,       phi(h)=3/rho.
```

There can be at most one such incident two-defect point.  Two distinct
support types would require `h^2=1` at their sole overlap, but `h^2=7/16`.
Two copies of the same type would have the same two normalized components
and their normalized dot product would exceed one, also contradicting their
forced unit edge.

The exact orbit

```text
h, T(h), ..., T^5(h)
 = rho/4, rho/3, 2rho/5, rho/2, rho, 0
```

then gives two more finite rules.

1. No path of two-defect support types can join two one-defect coordinates.
   A path of length `ell` would require `T^ell(h)=R(h)=4/rho`, which is not
   in the six-element orbit.
2. A component containing a one-defect coordinate cannot contain a vertex
   of degree at least three.  At such a branch, three incident normalized
   components have pairwise product one and hence all equal `1` or all equal
   `-1`; neither value occurs in the propagated orbit from `h`.

There is also an exact constraint between branch coordinates when no
one-defect point is present.  At a coordinate incident with at least three
distinct two-defect types, the endpoint values are a common sign
`sigma in {+1,-1}`.  Along a maximal path of `ell` types from one branch to
another, the endpoint signs must satisfy

```text
T^ell(sigma_left)=sigma_right.
```

The order-six orbit of either sign meets a sign only when `ell=0 mod 3`.
The sign is preserved for `ell=0 mod 6` and reversed for `ell=3 mod 6`.
The checker rejects every other branch distance and checks the resulting
binary sign constraints for consistency.

Finally, a fixed two-defect type has multiplicity at most two.  All such
points are common unit neighbours of the complementary seed K5 and are
pairwise unit.  They lie on its link circle, on which three pairwise unit
chords would require an impossible equilateral central-angle pattern.  This
is the previously used fixed-type multiplicity rule, included here because
the checker exhausts all one/two-support assignments in one place.

If a type occurs twice, neither endpoint can be incident with a different
two-defect type.  Otherwise unit products with a point of that different
type force the two copies to have the same endpoint component.  Applying
`phi` forces their other components to agree as well, making their normalized
vectors equal and contradicting their own required unit edge.  Thus a
parallel two-defect type must be an isolated edge of the simple type graph.

## Checker quantifiers and trust boundary

After a labeled zero-factor support family has been fixed and singleton
orthogonality propagation has produced masks `E_x`, every nonzero-factor
vertex with `|E_x|<=2` has one of the nonempty actual supports contained in
`E_x`.  The checker exhausts all of these choices.  A required edge between
two enumerated points must have intersecting actual supports, because its
normalized dot product is one.  Vertices with masks of size at least three
are left wholly unconstrained, so ignoring them cannot create a rejection.

For each complete choice the checker applies only the exact consequences
proved above.  A propagated support family fails only if every small-support
choice fails; the surrounding cover, seed, and graph quantifiers remain
existential in exactly the same direction as the frozen K7 reference.  All
operations are integer bit masks and finite enumeration.  Floating point is
used only for elapsed-time reporting.

The focused controls independently multiply the matrix `M` in `Q(sqrt(7))`,
check all five negative discriminants, the one-defect orbit, and both full
branch-sign orbits, exercise the forbidden 3-cycle, and retain a 6-cycle as a
negative control.  The sample profiler refuses to run unless the imported
frozen rank reference and labeled-support propagation sources have their
pinned SHA-256 hashes.  Its report records those dependency hashes, every
input/source/proof/test hash, the command and Python/platform data, and the
Git commit, branch, dirty flag, and exact porcelain-status hash.

## Exact sample profile

On the 80 survivors of the committed 512-graph labeled-support sample, a
single-process replay took under nine seconds and rejected 19 graphs:

```text
3946949 2592595 3732572 3333486 2256574 2548658
3330881 3950247 3958767 3336048 2934777 3331537
3328143 3950926 369959 2873261 3330537 2646137 2592168
```

All 239 explored propagation-passing labeled families first passed a control
that used only required-edge support intersection.  The value layer rejected
169 of those families and passed 70.  All 19 graph rejections are additional to
the independent strict-H sample decisions.  These are sample results, not a
full-residue coverage claim.
