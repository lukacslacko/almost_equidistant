# Degree-two and outside-diagonal preordering dual pilot

## Result in one sentence

A sound exact strengthening of the K7 positive-polynomial dual was built and
checked, but on the complete nine-graph saturating-clique residue of the
frozen 84-graph sample it found **zero** certificates beyond multiplier
degree one.  This is useful negative evidence: the next algebraic effort
should address covers without a saturating clique, especially through the
rank-one Schur complement of a near-saturating clique, instead of launching a
full degree-two run unchanged.

## Extra strict polynomials

For a saturating required clique `C`, retain the notation from
`d6_k7_positive_polynomial_dual.md`:

```text
K[C,C] = J + diag(e_i),  x_i = 1/e_i > 0,
T = 1 + sum_i x_i,
w(A) = sum_{i in A} x_i.
```

For an outside vertex `y` with basis-neighbour mask `A_y`, Schur equality on
the diagonal and distinctness give

```text
K_yy = p_y^T H p_y > 1,
H = diag(x_i) - x x^T/T.
```

Thus the following polynomial is strictly positive:

```text
h_y(x) = T w(A_y) - w(A_y)^2 - T > 0.
```

This is stronger than merely requiring `x_i>0`; it was present only as a
linear inequality in the earlier affine-`H` screen and was absent from the
first positive-polynomial dual.

## Exact truncated-preordering certificate

Let `f_j(x)=0` be the exact off-diagonal Schur equations.  Fix an output
degree `D`.  The locator searches for the exact identity

```text
sum_j q_j(x) f_j(x)
  = P_0(x)
    + sum_y P_y(x) h_y(x)
    + sum_{y<z} P_yz(x) h_y(x) h_z(x),
```

where every coefficient of `P_0`, `P_y`, and `P_yz` is nonnegative and at
least one is nonzero.  Multiplier degrees are truncated so that every term
has degree at most `D`.  The pilot uses `D=4`; hence the equation and
single-generator multipliers have degree at most two, while pair-generator
multipliers are constant.

At a feasible point every `x_i` and every `h_y` is strictly positive.  A
nonzero polynomial with nonnegative coefficients is therefore strictly
positive, as is every retained generator product with a nonzero such
multiplier.  The right side of the identity is positive, while the left side
is zero.  This is a contradiction.

HiGHS is only a locator.  It searches the rational polyhedral cone after
normalizing the total positive coefficient mass to one.  The code freezes
the numerical zero face, reconstructs a rational point in that face, and
accepts only after expanding the entire identity with `Fraction`.
`d6_k7_dual_degree2_verify.py` independently rebuilds every `f_j`, every
`h_y`, and every product before checking a serialized identity.  It never
imports the new numerical locator.

Original candidate nonedges are not constrained.  As in the frozen K7 rank
layer, the zero-factor cover first proves the exact normalized zero/one
off-diagonal pattern on `G[N]`; the dual is applied only after that step.

## Deterministic pilot

The pilot selected all nine degree-one survivor graphs in the frozen
84-graph sample that still had at least one strict-H-passing saturating
clique:

```text
369959 2592168 2681964 3328143 3331537
3939503 3950247 3958767 3959628
```

Single-process reproduction:

```text
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  /Users/lukacs/claude/opengauss/venv/bin/python3 \
  d6_k7_dual_degree2_pilot.py \
  --indices 369959 2592168 2681964 3328143 3331537 \
            3939503 3950247 3958767 3959628 \
  --output-degree 4 --maximum-generator-order 2 \
  --output d6_k7_dual_degree2_pilot_report.json

/Users/lukacs/claude/opengauss/venv/bin/python3 \
  d6_k7_dual_degree2_verify.py
```

Measured on the local Apple-silicon Mac while the full degree-one campaign
owned the main CPU budget:

| quantity | count |
|---|---:|
| graphs | 9 |
| strict-H-passing covers | 62 |
| covers with no saturating clique | 26 |
| degree-one-surviving covers with a saturating clique | 15 |
| covers killed by degree one | 21 |
| new pure degree-two cover certificates | **0** |
| new degree-four preordering cover certificates | **0** |
| new graph rejections | **0** |
| one-process wall time | 17.11 s |

The report is `d6_k7_dual_degree2_pilot_report.json`.  No negative claim is
made for a failed LP search: all 41 covers not killed at degree one remain
unresolved.  In particular, “zero new certificates” is a measured locator
result, not a proof that no degree-four identity exists.

## Structural blocker and recommended next extension

The sample separates two limitations:

1. Twenty-six of the 62 covers have no saturating clique at all.  No increase
   in multiplier degree can make the current rank-zero Schur equations apply
   to them.
2. On the remaining degree-one-passing saturating-clique covers, neither pure
   degree-two multipliers nor the new strict-diagonal preordering produced a
   certificate.  A full run over thousands of graphs is therefore not yet
   justified by this pilot.

For a required clique of size `U-1`, the Schur complement is positive
semidefinite of rank at most one rather than zero.  If

```text
g_yz = T (K_yz - p_y^T H p_z),
```

then four distinct remainder vertices satisfy the exact tetrads

```text
g_ij g_kl - g_ik g_jl = 0,
g_ij g_kl - g_il g_jk = 0,
```

and triples satisfy

```text
g_ij g_ik g_jk >= 0.
```

These equations apply precisely where the saturating-clique dual often does
not.  A sensible next pilot is a positive-polynomial/SOS dual generated by
the tetrads and rank-one sign polynomials on a stratified sample of the 26
no-saturating-clique covers.  Any numerical solution or dual must again be
converted to an exact polynomial identity before it rejects a cover.

## Trust and scope

- Exact rejections would trust Python integer/`Fraction` arithmetic and the
  frozen seed/cover quantifiers.  IEEE-754 and HiGHS are discovery aids only.
- The independent checker validates emitted rational identities, not the
  claim that the LP locator exhausts all identities at the chosen degree.
- The pilot emitted no new real certificate, so it changes no current
  dimension-six residue count.
- The synthetic infeasible and feasible controls, plus a real previously
  certified clique and a tamper rejection, are in
  `d6_k7_dual_degree2_test.py`.
