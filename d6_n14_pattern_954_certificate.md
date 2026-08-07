# Pattern 954: a 13-vertex required-edge obstruction in dimension 6

## Exact claim

Let `F` be the 13-vertex graph recorded in
`d6_n14_pattern_954_input.json`.  There are no 13 distinct points in
`R^6` for which every edge of `F` has distance one.  No condition is imposed
on a nonedge of `F`: its endpoints may or may not be at distance one.

Consequently, any graph containing `F` as an ordinary, non-induced subgraph
is also impossible as a unit-edge graph on distinct points in `R^6`.

The graph `F` has 54 edges.  It is obtained from graph 954 (zero based) of the
hash-pinned 14-vertex corpus by discarding source vertex 1 and the two
unneeded source edges `(0,3)` and `(5,8)`.  Its reindexed adjacency rows are

```text
2842 6777 3496 8039 4931 7374 8122 5476 7901 7515 7148 5999 4090
```

## Simplex coordinates

Vertices `(3,6,8,9,10,11,12)` of `F` form a unit `K7`; call the corresponding
centered regular-simplex points `q_0,...,q_6`.  For any other point `x`, put

```text
u_i(x) = ||x-q_i||^2 - 1,
s(x)   = sum_i u_i(x),
c_x    = (s(x)+1)/sqrt(7).
```

The Gram matrix of a centered unit `K7` gives the exact identities

```text
||u(x)||^2 - c_x^2 = 1,                                      (1)
||x-y||^2 - 1 = c_x c_y - u(x) dot u(y).                    (2)
```

Thus a required unit edge `xy` gives

```text
c_x c_y = u(x) dot u(y).                                    (3)
```

An edge from `x` to seed vertex `q_i` forces `u_i(x)=0`.  A missing
pattern edge to `q_i` merely allows `u_i(x)` to be nonzero; it never requires
that.  The six outside roles, with their core label, source label, and allowed
support upper bound, are

| role | core | source | `supp(u)` is a subset of |
|---|---:|---:|---|
| A | 0 | 0 | `{1,4,6}` |
| C | 4 | 5 | `{0,4,5}` |
| E | 1 | 2 | `{2,4}` |
| B | 2 | 3 | `{1,3,6}` |
| D | 7 | 8 | `{0,3,5}` |
| F | 5 | 6 | `{2,3}` |

The required outside edges include both triangles `ACE` and `BDF`, and the
bridge `EF`.  Pairwise support intersections in `ACE` are contained in the
singleton `{4}`; those in `BDF` are contained in `{3}`; the `E,F` intersection
is contained in `{2}`.

## Zero factors and distinctness

Suppose `c_x=0`.  Its definition and (1) give

```text
sum_i u_i(x) = -1,       sum_i u_i(x)^2 = 1.                (4)
```

If `u(x)` has support of size at most two, (4) forces `u(x)=-e_i`: for two
putative coordinates `a,b`, `(a+b)^2=a^2+b^2` gives `ab=0`, after which the
remaining coordinate is `-1`.  But `u=-e_i` is exactly the seed point `q_i`.
Distinctness therefore implies

```text
c_x=0  =>  |supp(u(x))| >= 3.                               (5)
```

In particular `c_E` and `c_F` are nonzero because their support upper bounds
have size two.

All four remaining factors are nonzero even if extra target edges have
shrunk the displayed support bounds.  For example, assume `c_A=0`.  By (5),
all three allowed A-coordinates are then actually nonzero, in particular
`u_A4 != 0`.  Equations (3) on edges `AE` and `AC`, whose support overlaps are
contained in `{4}`, force `u_E4=u_C4=0`.  Equation (3) on edge `CE` now gives
`c_C c_E=0`; since `c_E!=0`, this says `c_C=0`.  Applying (5) to C makes all
three of its allowed coordinates nonzero, contradicting `u_C4=0`.  Hence
`c_A!=0`.  Interchanging A and C proves `c_C!=0`.  The identical argument on
triangle `BDF`, with coordinate 3 and the already nonzero `c_F`, proves
`c_B,c_D!=0`.  Thus

```text
c_A,c_C,c_E,c_B,c_D,c_F are all nonzero.                    (6)
```

This argument only used support containment, so it remains valid after any
allowed-support shrinkage.

## Singleton-overlap triangles

Normalize each outside vector by `h_X=u(X)/c_X`.  From (1) and the definition
of `c_X`,

```text
||h_X||^2 - 1 = (sqrt(7) - sum_i h_Xi)^2.                   (7)
```

Every required outside edge now has `h_X dot h_Y=1`.  The three equations on
triangle `ACE` reduce, using the singleton overlap, to

```text
h_A4 h_C4 = h_A4 h_E4 = h_C4 h_E4 = 1.
```

They force

```text
h_A4 = h_C4 = h_E4 = sigma,       sigma in {-1,+1}.         (8)
```

Likewise triangle `BDF` forces

```text
h_B3 = h_D3 = h_F3 = tau,         tau in {-1,+1}.           (9)
```

The bridge `EF` has only coordinate 2 in its possible overlap, so

```text
h_E2 h_F2 = 1.                                               (10)
```

## Two-support factorization and contradiction

Write `e=h_E2`.  Vector `h_E` is supported on coordinates 2 and 4, and its
coordinate 4 is `sigma`.  Substitution in (7) gives

```text
e^2 = (sqrt(7)-e-sigma)^2,
(2e+sigma-sqrt(7)) (sqrt(7)-sigma) = 0.
```

The second factor cannot vanish because `sqrt(7)` is neither `+1` nor `-1`.
Therefore

```text
h_E2 = (sqrt(7)-sigma)/2.                                   (11)
```

The same calculation gives

```text
h_F2 = (sqrt(7)-tau)/2.                                     (12)
```

Equations (10)--(12) would require

```text
(sqrt(7)-sigma)(sqrt(7)-tau) = 4,
3 + sigma*tau = sqrt(7)(sigma+tau).                         (13)
```

If `sigma` and `tau` have opposite signs, (13) reads `2=0`.  If both are
`+1`, it reads `4=2sqrt(7)`; if both are `-1`, it reads
`4=-2sqrt(7)`.  All four cases are impossible.  This proves the claim.

## Why ordinary subgraph containment is sound

Embedding `F` into a target requires only that every edge displayed in `F`
map to a target edge.  Extra target edges to the seed turn some formerly
allowed defect coordinates into forced zeros, so they only shrink the six
support sets used above.  Extra target edges among the outside roles add
equations but remove none of the triangle or bridge equations.  The proof was
written throughout with support *upper bounds* and is therefore monotone
under all such additions.  Pattern and target nonedges are never interpreted
as non-unit distances.

## Reproduction and trust boundary

The extractor pins the source corpus, the current 12,839-target K7 residue,
and all intermediate input hashes.  The independent verifier reconstructs
the core, checks its required-edge hypotheses, and performs the four cases in
the exact quadratic field `Q(sqrt(7))`:

```text
python3 d6_n14_pattern_954_extract.py
python3 d6_n14_pattern_954_verify.py
python3 -m unittest -v d6_n14_pattern_954_test.py
python3 d6_n14_pattern_954_containment.py
```

The final command is preflight-only unless `--run` is supplied.  The C
containment kernel returns a concrete injective mapping for every `HIT`; the
Python orchestrator independently checks every required edge before recording
that hit.  `TIMEOUT` and `INFRA_ERROR` are retained as unresolved and are
never counted as rejections.

The existing interval engine also finds two placement orders, each with two
circle stages, but it is not part of this certificate.  In particular, this
proof does not depend on floating-point rank, numerical optimization,
`nextafter`, transcendental functions, or the engine's conservative 8-ulp
`libm` padding assumption.
