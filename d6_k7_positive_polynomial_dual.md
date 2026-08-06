# Exact positive-polynomial dual for the K7 Schur equations

## Why a plain PSD dual does not strengthen strict H

For a saturating required clique `C`, the existing strict-affine-H layer
solves the Schur equalities in the entries of

```text
H = K[C,C]^-1
```

and imposes, simultaneously, negative off-diagonal entries and positive row
sums.  These inequalities already imply that `H` is symmetric, has positive
diagonal, and is strictly diagonally dominant:

```text
H_ii > sum_{j != i} |H_ij|.
```

Thus `H` is positive definite.  A dual certificate separating the affine
slice from the positive-definite cone cannot reject a case that passes these
strict signs.  Positive semidefiniteness alone is not the missing condition
in this clique setup; the missing condition is the special inverse form.

## Sherman--Morrison polynomial system

Write

```text
K[C,C] = J + diag(e_i),       e_i > 0,
x_i = 1/e_i > 0,
T = 1 + sum_i x_i.
```

For an outside basis-neighbour mask `A`, put

```text
w(A) = sum_{i in A} x_i.
```

The inverse formula is

```text
K[C,C]^-1 = diag(x_i) - x x^T/T.
```

Consequently, for every pair of outside masks `A,B`, the zero Schur
complement gives the polynomial equality

```text
f_AB(x)
  = T w(A intersect B) - w(A)w(B) - k_AB T
  = 0,
```

where `k_AB` is the forced zero-or-one entry of the normalized `K` pattern.
These zeros are sound only after fixing the K7 seed and zero-factor cover,
exactly as in `d6_k7_rank_reference.py`.  This does not impose a non-unit
distance on an arbitrary candidate nonedge: the cover decomposition and
`alpha(G)<=2` prove the normalized zero pattern on its nonzero-factor set.

## Farkas certificate on the positive orthant

Let `q_AB(x)` be rational polynomials.  Suppose the exact identity

```text
P(x) = sum_AB q_AB(x) f_AB(x)
```

has a nonzero right-hand side all of whose monomial coefficients are
nonnegative.  Because every `x_i` is strictly positive, every monomial is
strictly positive and hence `P(x)>0`.  But a common zero of the Schur
equations would give `P(x)=0`.  The system is therefore impossible.

The implemented locator restricts each multiplier to degree at most one.
Writing all output coefficients as `Cq`, it solves the linear feasibility
problem

```text
Cq >= 0,       sum(Cq) = 1,
```

with free multiplier vector `q`.  HiGHS supplies only a candidate face of
this rational polyhedron.  The code enforces its zero coefficients and the
normalization over `Fraction`, chooses a nearby rational point, and accepts
only after expanding the complete identity exactly.  Failure to reconstruct
is `UNRESOLVED`, never a rejection.

The serialized certificate stores every nonzero rational equation
multiplier and every positive output coefficient.  The separate verifier
does not invoke HiGHS or any floating-point computation.  It reconstructs
the masks and equations from the graph and checks the polynomial dictionary
identity over arbitrary-precision rational arithmetic.

## Frozen sample result

The input is the 84-graph survivor set in
`d6_k7_strict_h_report.json`, selected from the deterministic 512-graph K7
rank sample.  Degree-one multipliers give:

| quantity | exact count |
|---|---:|
| selected graphs | 84 |
| K7 seeds | 111 |
| eligible covers before prior exact filters | 43,094 |
| enhanced-rank passing covers | 380 |
| strict-H passing covers | 297 |
| strict-H passing saturating cliques | 85 |
| newly failed covers (one retained clique certificate each) | 28 |
| new whole-graph rejections | **2** |
| survivors | 82 |

The newly rejected graph indices are

```text
3649646 3950926
```

All 28 exact identities are stored in
`d6_k7_positive_polynomial_dual_report.json`.  The graph-level count is
smaller because a graph is rejected only when every eligible cover for at
least one required K7 seed fails.  In particular, this distinction avoids
the earlier pitfall of counting failed cliques as failed graphs.

Reproduction:

```text
python3 -m unittest -v test_d6_k7_positive_polynomial_dual.py
python3 d6_k7_positive_polynomial_dual.py \
  d6_k7_rank_sample.json \
  --selection-report d6_k7_strict_h_report.json \
  --multiplier-degree 1 \
  --output d6_k7_positive_polynomial_dual_report.json
python3 verify_d6_k7_positive_polynomial_dual.py
```

## Scope and trust assumptions

- The mathematical certificate uses only rational polynomial identities and
  the strict inequalities `x_i>0`.
- The LP locator is untrusted.  No IEEE-754, `libm`, MPS, or GPU assumption
  enters a rejection.
- The independent checker trusts Python arbitrary-precision integer and
  `Fraction` arithmetic, plus the frozen prior exact seed/cover and strict-H
  implementations whose quantifiers it reruns.
- This is a sample result, not a claim about all 17,161 full strict-H
  survivors.  A full run should be made restartable and committed before it
  is launched.

## Extension to non-clique principal bases

For an arbitrary principal basis `C`, the Sherman variables are unavailable
because `K[C,C]` is no longer `J+diagonal`.  However a complementary exact
PSD dual applies directly to the affine equations in `H=K[C,C]^-1`.
If

```text
<A_i,H> = b_i,                 (outside off-diagonal equations)
<B_y,H> > 1,                  B_y=p_y p_y^T,
H > 0,
```

then rational multipliers `a_i`, nonnegative `lambda_y`, and a rational
positive-semidefinite matrix `Q` satisfying

```text
Q = sum_i a_i A_i - sum_y lambda_y B_y >= 0,
sum_i a_i b_i - sum_y lambda_y <= 0
```

give a contradiction provided

```text
(Q, lambda) != (0, 0).
```

If `Q` is nonzero, positive definiteness gives `<Q,H> > 0`, while the affine
equalities and strict outside-diagonal inequalities force `<Q,H> <= 0`.  If
`Q=0`, the displayed nontriviality condition gives `sum lambda_y>0`, so the
same strict inequalities instead force `0<0`, again a contradiction.  It is
therefore incorrect to require `Q` itself to be nonzero: a strict affine
inequality can already be impossible with the zero PSD matrix.
A restricted but exactly checkable search can represent
`Q=sum_v mu_v vv^T` with rational `mu_v>=0` and a finite dictionary of small
integer vectors `v`.  Its normalization is

```text
sum_y lambda_y + sum_v mu_v = 1,
```

which enforces `(Q,lambda)!=(0,0)` even when all `mu_v` vanish.  To exclude a
rank-at-most-`U` realization without a saturating clique, one must eliminate
every possible principal basis size and subset; eliminating a convenient
subset alone is not sound.

`d6_k7_arbitrary_basis_psd_dual.py` implements this restricted dual and an
exact checker.  On the frozen graph-58186 cover with `zmask=128` (`|N|=11`,
rank upper bound 6, clique lower bound 5), a deterministic pilot over the
first 32 lexicographic cores of each size found exact certificates for all
32 size-five cores and 22 of 32 size-six cores.  These counts include valid
certificates with `Q=0, lambda!=0`; such certificates mean that the proposed
core's strict outside-diagonal conditions are themselves inconsistent.  The
mechanism is nontrivial, but this pilot does **not** reject the cover:
uncertified size-six cores remain.  The corresponding three unit tests
include a real sample-core certificate, a `Q=0,lambda!=0` strict-inequality
certificate, and a feasible negative control.

## Near-saturating clique follow-on

If a required clique has size `U-1`, its Schur complement is positive
semidefinite of rank at most one rather than zero.  Keep the same positive
Sherman variables and define the polynomial numerator

```text
g_yz = T (K_yz - p_y^T H p_z).
```

Writing the Schur complement as `t t^T` shows, for distinct remainder
indices, the exact tetrad identities

```text
g_ij g_kl = g_ik g_jl = g_il g_jk
```

and the sign conditions

```text
g_ij g_ik g_jk >= 0.
```

Moreover the graph of nonzero off-diagonal `g_ij` is one clique together
with isolated vertices.  The tetrads have degree four in the `x_i` and are
natural generators for a next positive-polynomial or exact SOS dual.  This
may be especially useful when `U=5` and a required `K4` is common, but no
coverage claim for this unimplemented follow-on is made here.
