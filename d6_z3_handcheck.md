# Hand-verification of the |Z| <= 3 bounded-cover refinement (2026-08-07)

Scope: the claim in `d6_theory_filters.md` that for a K7 seed of a
level-19 candidate (twelve outside vertices), the zero-factor set
`Z = {x : c_x = 0}` — an eligible vertex cover of the disjoint-edge graph
`L` — has size at most three.  This note records an independent
line-by-line verification of that derivation, performed while the
in-house frontier computation ran.  Verified by the Claude track; the
derivation is the codex track's.

## Verified chain

1. **Setup identities** (from the centered regular unit K7,
   `||q_i||^2 = 3/7`, `<q_i,q_j> = -1/14`): `x = -sum u_i(x) q_i`,
   `||u(x)||^2 = 1 + c_x^2/7`, `M_xy = c_x c_y / 7 - u(x).u(y)`.
   These were previously machine-checked over random configurations
   (verify_codex_identities.py).  For `z in Z`: `||u(z)|| = 1`, and
   pairwise orthogonality of `{u(z)}` holds: required pairs via
   `M = 0` with `c_z = 0`; candidate nonedges have disjoint supports
   (a shared seed non-neighbour of a non-adjacent pair would be an
   independent triple, contradicting `alpha(G) <= 2`).  Likewise every
   `u(x)`, `x in N = T \ Z`, is orthogonal to every `u(z)`.  Hence
   `rank(K) <= 7 - |Z|` and `rank(B) = rank(K - J) <= 8 - |Z|`. VERIFIED.
2. **Required cliques have PD Gram** `K_C = J + diag(r_x^2)`, `r_x^2 > 0`:
   PD as PSD + positive definite diagonal shift. VERIFIED (elementary).
3. **|Z| = 7 impossible**: `rank(K) <= 0` forces `K = 0`, but diagonals
   are `1 + r_x^2 > 0` and `N` (5 vertices) is nonempty. VERIFIED.
4. **|Z| = 6 impossible**: `rank(K) <= 1`; `|N| = 6`; `alpha(G) <= 2`
   forces a required edge inside `N` (six vertices cannot be pairwise
   candidate-nonedges), giving a PD 2x2 block, rank >= 2. VERIFIED.
5. **|Z| = 5 impossible**: `rank(K) <= 2`; `|N| = 7`; colour the pairs of
   any six vertices of `N` by required/candidate; `R(3,3) = 6` gives a
   monochromatic triangle; a candidate triangle is an independent triple
   (forbidden), so a required triangle exists: PD 3x3 block, rank >= 3.
   VERIFIED.
6. **|Z| = 4 impossible**: `rank(K) <= 3`, `rank(B) <= 4`, `|N| = 8`,
   `F := complement of G[N]` (candidate-nonedge graph).
   - `F`-neighbourhoods are required cliques (independent-triple
     argument); their Gram vectors are independent (2.) and lie in the
     perpendicular of the vertex inside a span of dimension <= 3, so
     `Delta(F) <= 2`; `F` is triangle-free (a candidate triangle is an
     independent triple), so `F` is a union of paths, cycles of length
     >= 4, and isolated vertices. VERIFIED.
   - **Component/inertia bound**: `B` is block-diagonal over
     `F`-components; `K` PSD implies `B = K - J` has at most one
     negative eigenvalue, so at most one non-PSD block; a nontrivial PSD
     irreducible block is, after positive diagonal scaling, `I - C` with
     `C >= 0` irreducible symmetric, and Perron-Frobenius makes a zero
     eigenvalue simple (nullity <= 1); the possibly-indefinite block has
     nullity bounded by its zero-forcing number (standard); isolated
     blocks are nonsingular.  Hence
     `nullity(B) <= max(p, max_i(p - 1 + z_i))` with `Zf(path) = 1`,
     `Zf(cycle) = 2`. VERIFIED.
   - `nullity(B) >= 8 - rank(B) >= 4` forces, on eight vertices with the
     above component menu, exactly `F = 4K2` (p = 4) or `F = C4 + 2K2`
     (p = 3 with the cycle indefinite-or-not, 2 + 1 + 1). All other
     shapes (paths only, C4+K2+isolateds, C5/C6/C8 variants) give
     nullity bound <= 3. VERIFIED by exhaustion.
   - `F = 4K2`: each 2x2 block `[[r_x^2, -1], [-1, r_y^2]]` must be
     singular (four blocks, nullity 1 each); a singular such block has
     positive trace and zero determinant, hence is PSD; so `B` is PSD of
     rank 4, and for PSD summands `ker(B + J) = ker B intersect ker J`
     gives `rank(K) = rank(B + J) >= 4 > 3`. VERIFIED.
   - `F = C4 + 2K2`: on the cycle `a-b-c-d`, `K`-zeros make `v_b, v_d`
     both orthogonal to `v_a, v_c`; `v_b, v_d` are linearly independent
     (collinearity contradicts `K_bd = 1` against diagonals exceeding
     one: if `v_d = t v_b` then `K_dd = 1/(1 + r_b^2) < 1`); inside a
     rank-<=3 span their common perpendicular is one-dimensional, so
     `v_a, v_c` are collinear, and `K_ac = 1` gives the same diagonal
     contradiction. VERIFIED.

## Consequence

The `--zcap3` option of `indep_profile6.c` (reject a K7 seed whose
disjoint-edge graph has no eligible vertex cover of size at most three;
n = 19, twelve outside vertices) is now backed by an independently
hand-verified proof, upgrading its 5,144 additional exact kills
(K7-class union 3,685,523; first-generation residue 116,693 -> 111,551)
from reproduction tier to independently-verified tier.

## Not yet verified

The remaining second-generation machinery on feasible covers of size
<= 3: the actual-support CSP (`|S_x intersect S_y| != 1`, saturation
matchings), the normalized-matrix zero-forcing/inertia rank tests as
implemented, the saturating-clique Schur mask rules, the K7 tetrad
identities, and the K6-side chains (PSD Z-matrix, Hall, empty-support,
virtual-K7, tight-Hall, singleton-fan, repeated-arm).  These remain the
conditional tier.
