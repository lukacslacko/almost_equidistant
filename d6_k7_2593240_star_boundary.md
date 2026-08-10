# Exact rank-one star obstruction for K7 graph 2593240

## Outcome

The exact exploratory probe
`probe_d6_k7_2593240_star_boundary.py` rejects graph **2593240** at the
independently verified K7 support boundary.  The graph has exactly two
required K7 seeds.  Each seed has 502 raw eligible zero-factor covers.  The
upstream exact layers leave two current covers per seed; the nonzero cover is
already infeasible and `zmask=0` has one surviving propagated-support family.
The two remaining algebraic systems are isomorphic, and the probe gives an
exact contradiction for their common form.

This is not yet included in a production residue manifest.

## Optional-zero audit

Candidate nonedges are unconstrained and may still be unit pairs.  Therefore
the probe does **not** turn a candidate nonedge into a zero normalized-Gram
entry merely because it is a nonedge.  For the sole surviving support family:

- a required graph edge gives normalized entry 1;
- a zero entry is used only when the two independently propagated support
  supersets are disjoint.

Every basis and remainder entry used below is reconstructed by that rule.  No
unknown optional entry enters the certificate.

For the first K7 seed, the local required K6, six remainder masks, and the 15
remainder targets in lexicographic pair order are

```text
C       = (0,1,2,7,9,10)
masks   = (43,30,54,58,46,11)
targets = 011111111110111
```

For the second seed they are

```text
C       = (0,1,3,7,9,10)
masks   = (57,29,43,53,45,11)
targets = 111110111111101
```

These target strings correct one-bit transpositions in the preliminary
numerical handoff.  The probe reconstructs them from the frozen adjacency and
support supersets rather than trusting either string.

The systems are exactly isomorphic.  From the first to the second, the
remainder permutation is `(2,1,3,0,4,5)` and the K6-coordinate permutation is
`(1,0,2,3,4,5)`.

## Rank-one Schur system

Let the required K6 principal block be

```text
M = J + diag(e_0,...,e_5),   e_i>0,
x_i = 1/e_i>0,               T=1+sum_i x_i.
```

It is positive definite.  The normalized Gram matrix has rank at most 7, so
the Schur complement on the six remainder vertices is positive semidefinite
of rank at most one.  For a remainder mask `A_i`, put
`w(A)=sum_{j in A}x_j`.  The known scaled off-diagonal entry is

```text
g_ij = T k_ij - T w(A_i intersect A_j) + w(A_i)w(A_j).
```

Rank one gives `g_ij=h_i h_j`.  In particular, all 30 four-distinct tetrads
vanish:

```text
g_ij g_kl - g_ik g_jl = 0,
g_ij g_kl - g_il g_jk = 0.
```

Write the six positive inverse diagonal weights as `(a,b,c,d,e,f)`.  One
tetrad factors immediately as

```text
-f(a+1)T(c-e)=0,
```

so `e=c`.  Exact rational Groebner reduction of the remaining tetrads after
that substitution gives

```text
(c-f)(2c+1)=0,
(c-d)(c-1)=0,
(c-d)(c+d-2)=0.
```

Because `c>0`, the first equation gives `f=c`.  If `d!=c`, the last two would
give `c=1` and `c+d=2`, hence `d=1=c`, a contradiction.  Thus
`d=e=f=c`.

A second exact elimination gives

```text
(c-1)(c+1)(2c+1)(3c+1)=0,
9a - 6c^3 - 11c^2 + 3c + 8 = 0,
3b - 6c^3 - 11c^2 + 3c - 1 = 0.
```

Strict positivity leaves the unique tetrad solution

```text
(a,b,c,d,e,f) = (2/3,5,1,1,1,1).
```

All ideal membership and factor identities are recomputed over exact
rationals by the probe; numerical optimization is absent.

## The missing three-index obstruction

At the unique positive tetrad point, the nonzero known off-diagonal entries
are

```text
g_01=-8/3,  g_02=8,  g_03=-8/3,  g_04=-8/3,  g_05=-28/3,
```

and every leaf-leaf entry is zero.  Thus the off-diagonal support is a star.
The four-distinct tetrads all vanish on a star, which explains why the prior
tetrad pass could not reject this cover.

But a finite rank-one PSD matrix has `g_ij=h_i h_j`.  Since `g_01` and
`g_02` are nonzero, `h_0,h_1,h_2` are all nonzero, forcing
`g_12=h_1h_2` to be nonzero.  Exactly `g_12=0`, the final contradiction.

This suggests a reusable cheap strengthening: after the tetrad equations,
the graph of proved nonzero rank-one off-diagonal entries must be a clique on
the nonzero coordinates (plus isolated vertices).  Four-index tetrads alone
do not detect a nonzero star.

## Reproduction

```text
python3 probe_d6_k7_2593240_star_boundary.py \
  --output /tmp/d6_k7_2593240_star_boundary_probe.json
python3 -m unittest -v test_probe_d6_k7_2593240_star_boundary.py
```

The probe pins the v6 residue manifest, its independent verification, the
one/two-star support report, and its independent verification by SHA-256.
It trusts exact Python/Sympy rational polynomial arithmetic and those frozen
upstream proof boundaries.  It does not trust floating-point rank, a solver
failure, or a guessed nonedge value.
