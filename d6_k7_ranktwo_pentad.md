# Exact rank-two Schur pentads over the K7 no-near covers

## Reduction to rank two

Every exact `no near-saturating clique` cover in the current K7 residue has

```text
|N| = 12,       rank(K) <= 7,       omega(G[N]) = 5,       zmask = 0.
```

Fix one deterministic required `K5`, denoted `C`.  Its normalized Gram block
is

```text
K[C,C] = J + diag(e_i),       e_i > 0,
```

so it is positive definite.  The Schur complement after `C` is positive
semidefinite and has rank at most `7-5=2`.  This uses all seven remaining
vertices at once and is stronger than enumerating arbitrary principal bases.

Set

```text
x_i = 1/e_i > 0,       T = 1 + sum_i x_i,
H = K[C,C]^-1 = diag(x_i) - x x^T/T.
```

For remaining vertices `y,z`, their zero-one columns to `C` are `p_y,p_z`.
The known off-diagonal Schur entry and its polynomial numerator are

```text
S_yz = K_yz - p_y^T H p_z,
g_yz = T S_yz.
```

Each `g_yz` is an exact polynomial of degree at most two in the five strictly
positive `x_i`.

## Derivation of the pentad

Let `h_ij=<v_i,v_j>` be the off-diagonal Gram entries of five vectors in a
two-dimensional real space.  If `v_1` is nonzero, rotate coordinates and put
`t=||v_1||^2>0`.  For indices `j,k` in `{2,3,4,5}`, define

```text
r_jk = t h_jk - h_1j h_1k.
```

These are, up to the common positive scale `t`, products of the second
coordinates.  They therefore obey the residual rank-one tetrads

```text
r_23 r_45 - r_24 r_35 = 0,
r_23 r_45 - r_25 r_34 = 0.
```

After the constant terms cancel and the common factor `t` is removed, these
two equations are

```text
t A + B = 0,       t C + D = 0,
```

where

```text
A = h_23 h_45 - h_24 h_35,
B = -h_12 h_13 h_45 + h_12 h_14 h_35
    +h_13 h_15 h_24 - h_14 h_15 h_23,

C = h_23 h_45 - h_25 h_34,
D = -h_12 h_13 h_45 + h_12 h_15 h_34
    +h_13 h_14 h_25 - h_14 h_15 h_23.
```

Eliminating `t` gives the degree-five identity `A D - C B=0`.  Its exact
12-term expansion is

```text
 h12 h13 h24 h35 h45 - h12 h13 h25 h34 h45
-h12 h14 h23 h35 h45 + h12 h14 h25 h34 h35
+h12 h15 h23 h34 h45 - h12 h15 h24 h34 h35
+h13 h14 h23 h25 h45 - h13 h14 h24 h25 h35
-h13 h15 h23 h24 h45 + h13 h15 h24 h25 h34
+h14 h15 h23 h24 h35 - h14 h15 h23 h25 h34 = 0.
```

If `v_1=0`, every term contains two entries incident with vertex 1 and the
identity is immediate.  Thus the pentad holds without a nonzero-vector
assumption.  The independent checker also substitutes five generic 2D
vectors and expands all 384 raw products to zero.

The pentad is homogeneous of degree five.  Consequently

```text
P(g_yz) = T^5 P(S_yz) = 0.
```

If either `P(g)` or `-P(g)` has at least one term and every expanded
coefficient is positive, it is strictly positive for all `x_i>0`, an exact
contradiction.  There is no numerical locator: this is just integer/rational
polynomial expansion.

For one unordered five-subset, all 120 orderings give the same pentad up to
sign (60 with each sign).  The independent symbolic checker exhausts that
orbit, so the full scan needs one increasing ordering for each of the
`C(7,5)=21` subsets, not 120 redundant copies.

## Frozen target manifest and execution discipline

`build_d6_k7_ranktwo_pentad_targets.py` independently replays the hash-pinned
K7 seed, cover, rank, and clique quantifiers and writes the 53-cover target
manifest.  The full runner requires the exact manifest hash.

Sparse degree-ten multiplication creates many short-lived dictionaries.
Although the final polynomial has at most `C(15,5)=3003` monomials, retaining
the corpus replay and long-lived worker processes caused macOS infrastructure
terminations during discovery.  The production runner therefore:

- loads only the small frozen target manifest;
- keeps at most four short-lived cover children live;
- exits each child after one atomic JSON record;
- closes the scheduler every 16 covers;
- atomically checkpoints each completed batch;
- distinguishes infrastructure exit from `UNRESOLVED`.

The independent full verifier separately replays the 53 targets from the
original frozen corpus in a short-lived gate, derives the pentad without
importing the production kernel, checks the generic rank-two identity and all
120 labelings, reconstructs all 1,113 cover pentads in fresh children, checks
every stored positive polynomial, and runs the realizable 18-point control.
That control has no required `K7`, so the layer must report it as passing and
not applicable.

## Pre-freeze discovery observation (not a final production claim)

A discovery run before the source-only Git boundary found one-sign pentads on
36 of 53 covers.  Exact cover-level conjunction with the committed
joint-support checkpoints would reject a complete K7 seed in 34 of the 189
joint-support survivors.  Those 34 are disjoint from the 69 joint-support
graph rejections, suggesting a combined residue of 155 graphs.

These counts are deliberately labeled discovery only.  They must be
reproduced from the committed/pushed source boundary and independently
verified before entering the theorem-level residue manifest.

## Trust scope

The mathematical layer trusts the normalized-Gram/Schur reduction, the
hash-pinned exact target corpus, Python arbitrary-precision integers and
`Fraction`, and the independently expanded pentad identity.  It does not use
floating-point rank, an optimizer, a failed search, `libm`, an epsilon, or a
requirement that an omitted candidate nonedge be genuinely non-unit.
