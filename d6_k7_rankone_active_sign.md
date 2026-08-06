# Exact active/sign layer for a near-saturating clique

## Scope and outcome

This note isolates a graph-only consequence of the rank-one Schur system used
by `d6_k7_rankone_tetrad.py`.  It is cheaper than a polynomial dual and uses
no floating point.  For each off-diagonal Schur numerator it decides whether
the numerator is always positive, always negative, identically zero, or
sign-variable (and zero-capable) as the positive Sherman variables vary.

The resulting discrete abstraction has exactly two rejection mechanisms:

1. an identically-zero pair whose two endpoints are forced active, either by
   forced-nonzero pairs or by the diagonal extreme-mask rule below; or
2. an unbalanced signed graph of forced-nonzero pairs.

This is a necessary-condition layer only.  Passing it does not assert that a
single positive vector satisfies all pair equations, tetrads, diagonal Schur
conditions, or the geometric realization problem.

## Rank-one setup

Let `K` be the normalized Gram matrix on a fixed post-cover nonzero-factor set
`N`, with

```text
K_xx > 1,
K_xy in {0,1} for x != y,
rank(K) <= U.
```

Take a required clique `C` of size `U-1`.  Its principal block and inverse are

```text
K[C,C] = J + diag(e_i),       e_i > 0,
x_i = 1/e_i > 0,             T = 1 + sum_i x_i,
H = K[C,C]^-1 = diag(x_i) - x x^T/T.
```

For a remainder vertex `y`, let `A_y` be its required-neighbour mask on `C`,
let `p_y` be the indicator of that mask, and write

```text
w(A) = sum_{i in A} x_i.
```

For distinct remainder vertices define the off-diagonal numerator

```text
g_yz = T(K_yz - p_y^T H p_z)
     = K_yz T - T w(A_y intersect A_z) + w(A_y)w(A_z).
```

The Schur complement is positive semidefinite of rank at most one.  Hence it
is `t t^T`; after putting `a_y=sqrt(T)t_y`,

```text
g_yz = a_y a_z.
```

Define the actual active set

```text
A = {y : a_y != 0}.
```

Every pair inside `A` has nonzero numerator with sign
`sign(g_yz)=sign(a_y)sign(a_z)`.  Every pair touching the complement of `A`
has zero numerator.  Thus the nonzero graph of the off-diagonal numerators is
one clique, possibly empty or a singleton, together with isolated vertices.

### The diagonal implication, used only in the inactive direction

If `y` is inactive, then its diagonal Schur entry is zero.  Hence

```text
K_yy = p_y^T H p_y > 1
```

because distinctness gives `K_yy>1`.  Multiplication by `T` says that an
inactive vertex must satisfy

```text
h_y = T w(A_y) - w(A_y)^2 - T > 0.
```

For an empty mask, `h_y=-T`; for a full mask, with
`X=sum_i x_i` and `T=1+X`,

```text
h_y = TX-X^2-T = -1.
```

Thus empty-mask and full-mask vertices are forced active.  For a nonempty
proper mask, writing `a=w(A_y)` and `b=w(C\A_y)` gives

```text
h_y = b(a-1)-1,
```

whose sign is variable over positive `a,b`; this discrete layer makes no
activity inference in that case.  Conversely, **no** inequality `h_y>0` may
be imposed on an active vertex: its diagonal Schur entry is positive, not
zero.  Rank zero remains permitted when no diagonal or off-diagonal rule
forces activity, and non-extreme singleton active sets remain permitted.

## Exact Venn-region classification

For two masks `A_y,A_z`, partition the clique coordinates into four regions
and use the same letters for the sums of their positive variables:

```text
alpha = w(A_y intersect A_z),
beta  = w(A_y \ A_z),
gamma = w(A_z \ A_y),
delta = w(C \ (A_y union A_z)).
```

A region sum is zero exactly when its region is empty; otherwise it is
strictly positive.

For a zero target `K_yz=0`, expansion gives

```text
g_yz = beta gamma - alpha(1+delta).
```

Therefore:

| Venn condition | exact classification |
|---|---|
| `alpha=0`, `beta>0`, `gamma>0` | forced positive |
| `alpha=0` and (`beta=0` or `gamma=0`) | identically zero |
| `alpha>0` and (`beta=0` or `gamma=0`) | forced negative |
| `alpha>0`, `beta>0`, `gamma>0` | flexible |

In the last row, fixing positive `alpha,beta,delta` and varying `gamma`
crosses the unique zero

```text
gamma = alpha(1+delta)/beta.
```

Thus “flexible” means both signs and zero are possible at the level of this
one pair.  It does not claim that arbitrary choices for different pairs are
simultaneously attainable.

For a one target `K_yz=1`, expansion gives

```text
g_yz = 1 + beta + gamma + delta + beta gamma - alpha delta.
```

Therefore:

| Venn condition | exact classification |
|---|---|
| `alpha=0` or `delta=0` | forced positive |
| `alpha>0` and `delta>0` | flexible |

There is no forced-negative or identically-zero required-edge case.  In the
second row, varying `alpha` crosses the positive zero

```text
alpha = (1+beta+gamma+delta+beta gamma)/delta.
```

The special empty-mask case is worth making explicit.  Two disjoint masks do
not give a forced-positive zero-target numerator unless both masks are
nonempty.  If either mask is empty, the numerator is identically zero.

## Complete criterion for the discrete abstraction

Form a signed graph `F` on the remainder vertices.  Include precisely the
forced-positive and forced-negative pairs, labelled respectively `+1` and
`-1`.  Let `W` contain the endpoints of `F` **and** every remainder vertex
whose mask is empty or full.

Every vertex of `W` must be active, by either an off-diagonal forced-nonzero
numerator or the diagonal implication.  Consequently an identically-zero
pair contained in `W` is a contradiction.

If no such zero pair exists, the forced signs require vertex signs
`sigma_y in {+1,-1}` satisfying

```text
sigma_y sigma_z = sign(g_yz)       for every edge yz of F.
```

These equations are feasible exactly when every cycle of `F` has positive
edge-sign product, equivalently when the signed graph is balanced.  Breadth-
first sign propagation either constructs the signs or exposes a conflicting
edge.

This criterion is complete for the deliberately relaxed discrete layer.  If
`F` is balanced and no identically-zero pair lies in `W`, choose the relaxed
active set to be exactly `W` and use the propagated signs.  Then:

- every forced pair lies inside `W` and has the required sign;
- every pair touching the inactive complement is either identically zero or
  flexible, hence is allowed to be zero in the abstraction;
- every pair inside `W` is nonzero-capable, and a flexible pair permits the
  sign dictated by the chosen vertex signs; and
- every inactive vertex has a nonempty proper mask, for which `h_y>0` is
  individually attainable.

If `W` is empty, the empty active set is the canonical relaxed survivor.
Singleton active sets are also harmless at this layer because there is no
off-diagonal active-active constraint; an extreme-mask singleton can be
forced by the diagonal rule and is still consistent with this abstraction.
Rank zero is permitted whenever `W` is empty.  None of these survivor
statements is a realization claim; they only prevent an unsound rejection.

## A no-go theorem for the no-saturating-clique stratum

There is a useful campaign-specific simplification that eliminates the need
for a broad run on one whole stratum.

Assume the post-cover required-edge graph still has independence number at
most two, let `C` have size `U-1`, and suppose the graph has no required
`K_U`.  Then the active/sign layer can never reject at `C`.

Indeed, a full mask would extend `C` to a required `K_U`, so no full mask
occurs.  If remainder vertices `y,z` are nonadjacent, every `q in C` must be
adjacent to at least one of them; otherwise `{y,z,q}` is an independent
triple.  Thus

```text
A_y union A_z = C.
```

An identically-zero nonedge would then have one empty mask and the other full,
which is impossible.  A forced-negative nonedge has comparable overlapping
masks; because their union is `C`, the larger one is full, again impossible.
Required edges are only forced positive or flexible.  Consequently there are
no identically-zero pairs, every forced edge is positive, and its signed graph
is automatically balanced.  Empty masks may be diagonally forced active, but
without an identically-zero pair or a negative sign they create no conflict.

This proves zero coverage—not merely measured zero coverage—for the
no-saturating-clique covers processed by the degree-four tetrad layer.  It
does **not** cover degree-one survivors that contain a `K_U`; near cliques
inside that stratum can have full masks and forced-negative pairs.

## Candidate-nonedge and cover semantics

An original candidate nonedge is not prescribed to have non-unit Euclidean
distance.  The zero target used here is valid only inside the exact K7
zero-factor-cover construction documented in `d6_k7_rank_reference.py`:
after fixing the zero-factor set, the normalized Gram support on the
nonzero-factor vertices is exactly known.  The active/sign classifier accepts
an already constructed post-cover graph and never adds a zero merely because
two vertices were nonadjacent in the original candidate graph.

One contradictory near clique eliminates that cover, because every required
`K_(U-1)` in a hypothetical realization must obey the rank-one Schur
condition.  It does **not** by itself eliminate the original graph.  A graph
rejection requires reconstructing the full quantifiers: for some required K7
seed, every eligible zero-factor cover surviving the earlier exact layers must
be rejected.  Pilot output must therefore report marginal near-clique and
cover coverage separately and call graph coverage unproved unless those
quantifiers have actually been replayed.

## Implementation and independent checking

`d6_k7_rankone_active_sign.py` is the locator.  It serializes the complete
abstract system—variable count, remainder labels, basis masks, and every
zero/one pair target—plus the classification profile and a contradiction
witness.  It emits certificates only for rejections.

`verify_d6_k7_rankone_active_sign.py` does not import the locator.  It
independently:

1. validates the masks and the complete pair-target table;
2. reconstructs every Venn-region classification;
3. reconstructs the diagonal extreme-mask set, `W`, and the forced signed
   graph;
4. checks the recorded zero or parity contradiction; and
5. checks the recorded profile and witness data against its reconstruction.

For a theorem-level graph rejection, a wrapper checker must additionally pin
the corpus and upstream decision hashes, independently rebuild each K7 seed
and every eligible cover, rerun the earlier exact filters in their documented
order, confirm that the certificate's near clique is present and has size
`U-1`, and verify that no passing cover was omitted.  That outer quantifier
checker is intentionally not inferred from a small pilot.

The focused synthetic controls are in
`test_d6_k7_rankone_active_sign.py`.  They cover every row of both Venn tables,
empty and singleton active-set safety, a zero contradiction, a negative-sign
cycle, a balanced signed survivor, flexible-pair non-forcing, empty/full masks
forced by the diagonal, an extreme-mask zero contradiction, and rejection of
a tampered certificate by the independent checker.

## Stratified 24-cover pilot

`run_d6_k7_rankone_active_sign_pilot.py` reconstructed upstream state for a
deterministic one-process sample of 12 saturating and 12 no-saturating covers
among the completed degree-four campaign's 937 graph survivors.  It stopped
at 24 covers and therefore did not reconstruct complete graph quantifiers.

| stratum | covers | near cliques | rejected covers | rejected near cliques |
|---|---:|---:|---:|---:|
| contains a saturating `K_U` | 12 | 268 | 0 | 0 |
| no saturating `K_U`, tetrad survivor | 12 | 49 | 0 | 0 |
| total | 24 | 317 | **0** | **0** |

Across the 317 near cliques the classifier saw 1,420 forced-positive pairs,
306 forced-negative pairs, 3,029 flexible pairs, and no identically-zero
pairs.  The diagonal rule forced one vertex active in 156 near cliques, all
in the saturating stratum.  All forced signed graphs were balanced.  The 306
negative pairs all came from the saturating stratum, as predicted by the
no-go theorem.

The machine-readable report is
`d6_k7_rankone_active_sign_pilot_report.json`.  Its graph-coverage field is
explicitly `NOT_ASSESSED` with zero certified graph rejections.  The pilot's
zero hits do not justify a broad active/sign-only run; a stronger coupled
algebraic layer would need a separate measured pilot.
