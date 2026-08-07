# Remote steering after the K7 cover milestone

Review boundary: branch `codex/dimension6` at
`7b29c58691b690caf49d1575db4cec16943f1626`, with implementation commit
`0380903e1b848fbd3607e8bb91b612d30de75fa3`.

Read this after `AGENTS.md`, `STATUS.md`, and `REMOTE_STEERING.md`.  It
supersedes the immediate-action section of the earlier steering note, but not
its collaboration protocol or proof-engineering requirements.

## 1. Review of the completed milestone

I reviewed the new handoff, `d6_theory_filters.md`, `d6_profile.json`, the
reference implementation, the controls, and the obstruction-pilot summary.
I found no mathematical soundness defect in the implemented bounded-cover or
tight-cover rules.

The measured progress is substantial:

```text
old deferred population             916,313
exactly rejected by the new union    627,358
remaining K7 graphs                  113,136
remaining K6-only graphs             175,819
total current residue                288,955
```

The zero coverage of both clique-Hall rules is also explained correctly: on
this already K8-free corpus, each such Hall failure literally supplies a
required K8.

One conceptual simplification below is worth documenting, but it does **not**
by itself add coverage beyond the current size-seven matching rule.

## 2. The K7 zero-factor set is an orthonormal family

Fix a required K7 seed `Q={q_1,...,q_7}` and use the existing notation

```text
u_i(x) = ||x-q_i||^2 - 1,
c_x    = 1 + sum_i u_i(x),
D_x    = {seed coordinates allowed to be nonzero}.
```

Let

```text
Z = {outside x : c_x = 0}.
```

For `x in Z`, the simplex identity gives

```text
||u(x)||^2 = 1.
```

For any two distinct `x,y in Z`, in fact

```text
u(x) dot u(y) = 0.
```

There are only two cases:

1. If `xy` is a required graph edge, then `M_xy=0`, and the K7 identity
   `M_xy=c_x c_y/7-u(x) dot u(y)` gives orthogonality.
2. If `xy` is a candidate graph nonedge, `alpha(G)<=2` implies
   `D_x intersect D_y` is empty: otherwise a shared seed vertex would form an
   independent triple with `x,y`.  The actual supports are subsets of the
   allowed supports, so their dot product is again zero.

Thus the columns `u(x)`, `x in Z`, are orthonormal in `R^7`.  This directly
implies

```text
|Z| <= 7
```

and the allowed masks of `Z` must admit a matching saturating `Z`.

### Why a generic all-cover matching pass would be redundant

Do not spend a milestone merely replacing the current tight size-seven test by
"matching for every cover size".  On this K8-free corpus, every eligible cover
of size at most six is automatically matchable.

Indeed, a failed matching gives a Hall subset `Y` whose coordinate union `A`
has `|A|<|Y|<=6`.  Every eligible mask has size at least three.  Since
`|A|<=5`, no two masks in `Y` can be disjoint, so `Y` is a required clique.
Together with the `7-|A|` seed vertices outside `A`, it forms a required clique
of size at least eight, impossible in the corpus.

For a size-seven cover the only new Hall obstruction can occur at the full
seven-vertex set, exactly the equality case already implemented.  This also
answers the earlier near-equality question: the six-vector case is full-rank,
but its pure matching consequence is already forced by K8-freeness.

## 3. A genuinely stronger K7 support-pattern filter

The orthonormality statement does yield a new combinatorial layer beyond
matching.  Let the **actual** support of a zero-factor vertex be

```text
S_x = {i : u_i(x) != 0} subseteq D_x.
```

Necessary conditions are:

```text
|S_x| >= 3,
|S_x intersect S_y| != 1  for every distinct x,y in Z.
```

The first condition is the existing distinctness argument: with `c_x=0`, a
support of size at most two collapses to `u=-e_i`, hence to an existing seed
point.  The second condition follows because two orthogonal vectors cannot
have exactly one common nonzero coordinate.

This gives a small exact CSP for every eligible cover `Z`:

1. each variable chooses a subset `S_x subseteq D_x` of size at least three;
2. every pair of chosen supports has intersection size zero or at least two;
3. the chosen support family must still have a matching saturating `Z`.

For `|Z|=7`, the `7 x 7` matrix `U` is orthogonal and satisfies

```text
U^T 1 = -1,
U 1   = -1.
```

Therefore the row supports give additional necessary conditions:

- every row support has size one or at least three, never exactly two;
- two row supports cannot intersect in exactly one column.

Passing this support CSP is only a screen, not a realization claim.  Start with
a transparent Python implementation on the committed reference sample and on
a deterministic sample of the 113,136 K7 residue.  Promote it to the full C
profile only if it has measurable incremental coverage.

## 4. A second K7 reduction on the nonzero remainder

For a feasible cover `Z`, write `N=T\Z`; hence `c_x != 0` for `x in N`.  Define

```text
v_x = sqrt(7) u(x) / c_x,
r_x = sqrt(7) / c_x.
```

Let `K` be the Gram matrix of the `v_x`.  Exact algebra gives

```text
K_xx = 1 + r_x^2,
K_xy = 1  if xy is a required graph edge,
K_xy = 0  if xy is a candidate graph nonedge.
```

The last line uses `alpha(G)<=2`, which makes the allowed defect masks of a
candidate nonedge disjoint.  The cover property excludes a required edge with
disjoint masks from remaining wholly inside `N`.

Consequently

```text
K is positive semidefinite,
rank(K) <= 7.
```

Now put

```text
B = K - J.
```

Then

```text
B_xx = r_x^2 > 0,
B_xy = 0   on required graph edges,
B_xy = -1  on candidate graph nonedges,
rank(B) <= 8,
B has at most one negative eigenvalue.
```

Equivalently, the off-diagonal graph of `B` is the candidate nonedge graph
`F=overline(G[N])`, which is triangle-free.  This is a sparse fixed-weight
matrix problem and is a better target than the original dense coordinate
system.

### 4.1 Cheap exact rank tests

For each feasible cover, calculate the following exact bounds.

**Ordinary zero forcing.**  If `Zf(H)` is the exact ordinary zero-forcing
number of a graph `H`, every symmetric matrix whose off-diagonal graph is
exactly `H` has nullity at most `Zf(H)`.  The proof needed here is elementary:
if a kernel vector vanishes on a zero-forcing set, each color-change row, whose
unique uncolored neighbor has a nonzero matrix entry, forces the next
coordinate to vanish.

Apply this to both matrices:

```text
rank(K) >= |N| - Zf(G[N]),
rank(B) >= |N| - Zf(F).
```

Since `|N|<=12`, exact zero forcing by subset enumeration is tiny.

**Support upper bound.**  Let `nu_N` be the maximum bipartite matching from
vertices of `N` into their seven allowed defect coordinates.  Since the
columns `v_x` retain the zero pattern of `u(x)`,

```text
rank(K) <= nu_N,
rank(B) <= nu_N + 1.
```

Thus reject a cover if either lower rank bound exceeds respectively
`min(7,nu_N)` or `min(8,nu_N+1)`.

**Component/inertia refinement for B.**  The matrix `B` is block diagonal over
connected components of `F`.  At most one block can be indefinite.  Every
other connected nontrivial block is a positive-semidefinite irreducible
Z-matrix with positive diagonal and therefore has nullity at most one.  A
self-contained proof is obtained by scaling it to `I-C`, where `C` is a
connected nonnegative symmetric matrix; if singular and positive
semidefinite, eigenvalue one is its Perron root and is simple.  An isolated
vertex block has nullity zero.

If `p` is the number of nontrivial components of `F`, and `z_i` is the exact
ordinary zero-forcing number of component `F_i`, then

```text
nullity(B) <= max(
    p,
    max_i (p - 1 + z_i) over nontrivial components i
).
```

The first term covers the all-PSD case; the second allows component `i` to be
the sole indefinite block.  Convert this to a rank lower bound and compare it
with `min(8,nu_N+1)`.

### 4.2 Existential search over covers

These are necessary conditions for the **actual** zero-factor cover, not for a
canonical minimum cover.  For one K7 seed:

1. enumerate eligible vertex covers `Z` of the disjoint-edge graph `L` with
   `|Z|<=7`;
2. discard covers failing the exact support-pattern CSP above;
3. on `N=T\Z`, apply the `K` and `B` rank tests;
4. the seed is contradictory only if every possible cover is discarded.

Before a full run, report the histogram of feasible cover sizes and the
individual/marginal coverage of each rank rule on a deterministic residue
sample.

## 5. The useful finite K6 Lorentz rule

For a K6 seed retain the existing factors

```text
ell_x = (c_x,z_x) in R^(1,1),
<ell_x,ell_y>_L = c_x c_y - z_x z_y.
```

A required outside edge with disjoint allowed defect masks imposes Lorentz
orthogonality.  Let `L` be the graph of those edges.

A warning first: the Lorentz zero pattern of `L` **alone** is vacuous.  All
nonzero projective factors can be placed on one lightlike line, so do not build
an elaborate ratio-propagation system that ignores the coupled `u` vectors.
The following coupled consequence is nontrivial.

### 5.1 Zero factors and odd components

Let

```text
Z0 = {x : ell_x = 0}.
```

For `x in Z0`, the K6 quadratic identity gives

```text
||u(x)||^2 = 1,
sum_i u_i(x) = -1.
```

Distinctness again implies `|D_x|>=3`.  Exactly as in the K7 argument, the
vectors `u(x)`, `x in Z0`, are pairwise orthonormal in `R^6`.  Hence

```text
|Z0| <= 6
```

and their allowed masks must be matchable into the six seed coordinates.

Remove `Z0` from `L`.  In a connected non-bipartite component, every remaining
Lorentz factor is nonzero.  Orthogonality along an odd cycle forces its
projective direction to be fixed by the orthogonal-line involution in
`R^(1,1)`, hence to be lightlike.  Connectedness then places the entire
component on one of the two lightlike projective lines.

For every vertex on such a lightlike line,

```text
||u(x)||^2 = 1.
```

Moreover, any two vertices on the same lightlike line have Lorentz product
zero.  If their candidate pair is an edge, the unit equation forces their
`u`-dot-product to zero; if it is a nonedge, their allowed defect masks are
disjoint and the dot product is again zero.  Thus all `u` vectors in one
lightlike class are pairwise orthonormal.

Different non-bipartite components may choose either of the two lightlike
lines.  The zero-factor vectors are orthogonal to both classes.  Therefore a
necessary finite CSP is:

1. choose an eligible, support-matchable set `Z0`, with `|Z0|<=6`;
2. compute the non-bipartite connected components of `L-Z0`;
3. assign each such component one of two colors, corresponding to the two
   lightlike projective lines;
4. for each color, the union of that color's component vertices with `Z0`
   must have allowed masks admitting a matching into the six seed
   coordinates.  In particular its total size is at most six.

If no choice of `Z0` and no two-color assignment satisfies these matching
conditions, the K6 seed is impossible.  Bipartite components are deliberately
left unconstrained: they may use a generic non-lightlike pair of orthogonal
projective lines.  This makes the rule conservative and sound.

This is the recommended first exact attack on the 175,819 K6-only graphs.  It
uses only bitsets, connected components, bipartiteness, small matchings, and an
exhaustive search over at most thirteen outside vertices.

The same support-intersection condition can strengthen each orthonormal bin:
actual supports of two orthogonal nonzero vectors cannot intersect in exactly
one coordinate.  Implement that only after measuring the pure matching CSP.

### 5.2 Required controls

At minimum include:

- the realizable 18-point lower-bound configuration, tested over every K6 seed;
- a synthetic case with a robust odd component too large for either lightlike
  bin and with no eligible zero-factor deletion;
- a case with two odd components that fit only when assigned opposite colors;
- a case in which a candidate zero-factor set has enough vertices but fails
  support matching;
- C/Python agreement on a committed sample of K6-only residue graphs.

Do not use a synthetic negative graph as a theorem-level control unless all
preconditions used by the proof, especially `alpha(G)<=2`, are checked.

## 6. Obstruction patterns 954 and 952

Do not certify these patterns before remeasuring them against the stronger
direct filters above.  Their present sampled value is real but based on only
21 post-filter K7 targets, and a direct algebraic filter may absorb much of it.

After the new exact K7 profile:

1. recompute the deterministic incremental containment sample;
2. if pattern 954 still has substantial marginal coverage and remains usable
   by the existing certified engine, attempt it first;
3. try 952 only if it still adds material marginal coverage;
4. keep every containment claim heuristic until the pattern itself has a
   rigorous unit-edge-only non-realizability certificate.

The current 89-pattern family can never touch the K6-only population because
every pattern contains a K7.  In parallel with the exact K6 CSP, construct a
separate K6-only obstruction pilot: sample K6-only graphs at orders 14--16,
rank them numerically only for triage, measure containment on the 175,819
residue, and rigorously certify only the highest-marginal candidates.

## 7. Immediate local milestone

Before any new full interval campaign:

1. Add the orthonormal-zero-factor derivations and the normalized `K/B`
   matrices to `d6_theory_filters.md`.
2. Implement independent Python references and controls.
3. Profile, in this order:
   - K6 odd-component/two-light-ray CSP on a deterministic sample and then the
     175,819 K6-only graphs if the sample shows coverage;
   - K7 support-pattern CSP;
   - K7 `K/B` zero-forcing and component/inertia rank filters.
4. Report each rule separately and cumulatively, plus feasible-cover-size
   histograms and measured runtime.
5. Re-evaluate patterns 954 and 952 only after the direct exact union is
   updated.
6. Push a focused implementation commit and a separate handoff commit with the
   exact hashes, controls, counts, unresolved split, and specific questions.

Do not restart the 288,955-case interval grind unchanged.  Any long run still
needs atomic checkpoints and a manifest as specified in `AGENTS.md`.
