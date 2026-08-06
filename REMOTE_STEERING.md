# Remote steering after the exact dimension-6 profiling milestone

Review boundary: branch `codex/dimension6` at
`9292f29be9cd26e265ed44f7473cf36bf931e68c`, with implementation commit
`20d057d7f664b6358c026e4dec106e33d2aafa65`.

This file supersedes the now-completed "immediate next action" in section 9
of `AGENTS.md`. The collaboration protocol, soundness rules, proof-engineering
requirements, and handoff format in `AGENTS.md` remain in force.

## 1. Review of the completed milestone

I reviewed `STATUS.md`, `d6_profile.json`, `d6_profile_manifest.json`,
`d6_theory_filters.md`, the relevant parts of `profile_d6.c`, and
`verify_profile_d6.py`.

I found no mathematical soundness defect in the implemented clique-link,
`K7` reflection, or stated defect-support rules. In particular:

- nonedges are consistently treated as optional unit distances;
- the `K5` circle-link argument and its bound 4 are correct;
- the `K7` reflection coordinates and squared distance `16/9` are correct;
- the support CSP is exhaustive for the deliberately limited rules it claims
  to encode.

The measured conclusion is also strategically clear: the exact new graph
rules remove only 12,466 of the 916,313 deferred graphs. The residue of
903,847 graphs should not be attacked by restarting the same cap-grinding
campaign.

One proof-engineering limitation should be kept explicit. The Python control
program is an independent control suite, not a full independent
reimplementation of the 3,971,787-graph profile. That is adequate for this
profiling milestone, but each new exact rule below should also get a small
independent reference implementation and deterministic cross-check sample.

## 2. Answer to the `K7` Schur-complement question

Fix a required unit `K7`, written as a centered regular simplex
`q_1,...,q_7` in `R^6`. For an outside point `x`, put

```text
u_i(x) = ||x-q_i||^2 - 1,
s_x    = sum_i u_i(x),
c_x    = s_x + 1.
```

Let `U` be the `7 x 12` matrix whose column for `x` is `u(x)`. If `M` is

```text
M_xx = -1,
M_xy = ||x-y||^2 - 1,
```

and `T` denotes the outside vertices, then the seed block is `M_QQ=-I_7`.
The exact simplex identities give

```text
M_xy = c_x c_y / 7 - u(x) dot u(y),
||u(x)||^2 = 1 + c_x^2 / 7.
```

Therefore the Schur complement is not merely rank at most one:

```text
R := M_TT + U^T U = c c^T / 7.
```

Thus `R` is positive semidefinite of rank at most one, with its rank-one
factor explicitly fixed by the defect variables.

This yields two cheap exact filters that are stronger than the completed
discrete support CSP.

### 2.1 Simplex-clique Hall obstruction

Let

```text
D_x = {i in the K7 seed : x is a graph non-neighbour of q_i}.
```

The actual support of `u(x)` is a subset of `D_x`.

If `C` is a required unit clique among the outside vertices, then for
`x,y in C`

```text
u(x) dot u(y) = c_x c_y / 7,
```

including the diagonal identity above. Hence

```text
U_C^T U_C = I + c_C c_C^T / 7,
```

which is positive definite. The columns `u(x)`, `x in C`, are linearly
independent. Since column `u(x)` is supported inside `D_x`, the bipartite
allowed-support graph from `C` to the seven seed coordinates must have a
matching saturating `C`.

Equivalently, for every coordinate set `A subseteq {1,...,7}`, define

```text
V_A = { outside x : D_x subseteq A }.
```

A necessary condition is

```text
omega(G[V_A]) <= |A|.
```

This equivalence is just Hall's theorem: a failed matching supplies a clique
subset whose union of allowed coordinates is too small. This formulation is
particularly cheap because there are only 128 coordinate masks and at most
12 outside vertices. A contradiction relative to any required `K7` seed
rejects the graph.

Do not check only the cardinality of the union for one maximal clique unless
the implementation also covers all Hall subsets. The 128-mask induced-clique
formulation above covers them automatically.

### 2.2 Disjoint-edge zero-factor / bounded-cover obstruction

For a required outside unit edge `xy` with

```text
D_x intersect D_y = empty,
```

the unit equation has zero dot product for support reasons, and therefore

```text
c_x c_y = 0.
```

Let

```text
Z = { outside x : c_x = 0 }.
```

Then `Z` must be a vertex cover of the graph `L` whose edges are exactly the
required outside unit edges with disjoint allowed defect sets.

There are two additional exact restrictions on `Z`.

First, an outside point with `c_x=0` cannot have `|D_x| <= 2`. Indeed,

```text
sum u_i = -1,
||u||^2 = 1.
```

With at most two nonzero coordinates this forces `u=-e_i`, which reconstructs
the already present seed point `q_i`; distinctness excludes it. Hence

```text
Z subseteq {x : |D_x| >= 3}.
```

Second,

```text
|Z| <= 7.
```

Here is a self-contained proof. The seven seed points and all points in `Z`
lie on the same sphere of squared radius `3/7`. Map each such point `p` to

```text
w(p) = (sqrt(2) p, 1/sqrt(7)) in R^7.
```

These are unit vectors, and original unit distance is exactly orthogonality.
For any `N` unit vectors in `R^r` such that every three contain an orthogonal
pair, form the graph of non-orthogonal pairs. It is triangle-free, so the
non-orthogonal neighbours of each vector are pairwise orthogonal. Bessel's
inequality gives

```text
sum_{j nonorthogonal to i} <w_i,w_j>^2 <= 1.
```

For the Gram matrix `B`,

```text
tr(B^2) <= 2N.
```

Since `rank(B)<=r`, Cauchy--Schwarz on its eigenvalues gives

```text
N^2 = tr(B)^2 <= r tr(B^2) <= 2rN,
```

so `N<=2r`. Here `r=7` and the seed already contributes seven vectors,
therefore `|Z|<=7`.

The resulting exact test is tiny: enumerate subsets of the at most 12 outside
vertices and ask whether `L` has a vertex cover `Z` with

```text
Z subseteq {x : |D_x| >= 3},
|Z| <= 7.
```

If not, reject. This includes, but is strictly stronger than, the immediate
case of a disjoint-support required edge whose two endpoints both have at
most two allowed defects.

## 3. Answer to the `K6` coordinate and inertia question

Fix a centered regular unit `K6`, `q_1,...,q_6`, spanning a 5-dimensional
subspace `W` of `R^6`. Thus

```text
q_i dot q_i = 5/12,
q_i dot q_j = -1/12  (i != j).
```

Choose a unit normal `e` to `W`. For an outside point `x`, define the six
seed deviations `u_i(x)` as above, put

```text
s_x = sum_i u_i(x),
c_x = s_x + 1,
```

and decompose

```text
x = -sum_i u_i(x) q_i + h_x e.
```

The normal coordinate satisfies

```text
12 h_x^2 = c_x^2 + 6 - 6 ||u(x)||^2.
```

Set `z_x = sqrt(12) h_x`. For outside points `x,y`,

```text
M_xy = c_x c_y / 6 - u(x) dot u(y) - z_x z_y / 6.
```

With `U` the `6 x 13` defect matrix, the Schur complement of the seed block
`-I_6` is

```text
R := M_TT + U^T U,
6R = c c^T - z z^T.
```

Consequently:

```text
rank(R) <= 2,
number of positive eigenvalues of R <= 1,
number of negative eigenvalues of R <= 1.
```

The last restriction matters: `R` is generally indefinite and must not be
mistaken for a positive-semidefinite rank-two matrix. The explicit difference
of two rank-one matrices is stronger than rank plus inertia alone.

For symbolic work, keep `z_x` as a real variable and impose the rational
quadratic equation

```text
z_x^2 = c_x^2 + 6 - 6 ||u(x)||^2.
```

This avoids square roots in the polynomial system. The pair `(u,z)` uniquely
reconstructs the point once the seed and normal orientation are fixed; the
global reflection sends all `z` values to their negatives.

### 3.1 A cleaner interpretation for the first graph-only filter

A regular `K6` in `R^6` can always be completed by either of two auxiliary
points to a regular `K7`. Choose one auxiliary point only as a coordinate
anchor; it is not added to the almost-equidistant set. Its deviation is a
seventh, universally allowed coordinate. The full `K7` Schur complement is
then rank-one positive semidefinite.

This virtual-completion viewpoint immediately gives the `K6` version of the
simplex-clique Hall rule. For each outside vertex retain its six allowed
seed-defect coordinates and add one universal normal/virtual coordinate. If
`C` is a required outside unit clique, its associated seven-coordinate
vectors have Gram matrix

```text
I + c_C c_C^T / 6,
```

and are linearly independent. Hence the allowed-support bipartite graph,
with one universal right-hand vertex, must match all of `C`.

Equivalently, for every `A subseteq {1,...,6}`,

```text
V_A = { outside x : D_x subseteq A }
omega(G[V_A]) <= |A| + 1.
```

This should be tested on every required `K6` seed of the 175,819 `K6`-only
graphs. It is sound even when `D_x` is empty.

Important warning: the auxiliary seventh point is not part of the original
set. It is valid for the Schur algebra and Hall-rank argument, but its defect
coordinate must not be fed into the old overlap-forces-unit support CSP. A
triple involving the virtual point is not subject to almost-equidistance.

## 4. Answer to the obstruction-pilot question

Obtain a fast C containment estimate before spending substantial effort on
certification. The 89 `n=14` numerical patterns are hypotheses, not
obstructions, so use them only to rank possible payoff.

Recommended order:

1. Implement a fast, conservative C subgraph-containment pilot on a
   deterministic stratified sample of the current residue, separated into
   `K7` and `K6`-only populations.
2. Report individual and greedy marginal coverage for the 89 patterns.
3. Select the few patterns with the largest marginal coverage that also have
   a seed/decomposition the certified engine can actually handle.
4. Certify those patterns rigorously before using any reported coverage as an
   exact rejection count.
5. If the numerical `n=14` patterns have poor coverage or lack usable seeds,
   mine `K7`-anchored cores directly from already certified level-19 graphs by
   greedy vertex and edge deletion, rerunning the certificate after every
   deletion. This is likely to fit the existing engine better.

Do not optimize a full 3.97-million by 89 containment pass before a sample
shows useful coverage. Conversely, do not certify an arbitrary handful before
knowing whether they occur in the residue.

Any mined obstruction must be a unit-edge-only non-realizability statement.
It remains unsound to use omitted/nonedges as required non-unit distances.

## 5. Immediate local milestone

Do the following before any new exhaustive interval run.

### A. Document and independently test the algebra

Add the derivations above to a repository note. Add a small Python reference
checker independent of the C profiler. It need not process all 3.97 million
graphs, but it must compare decisions on:

- all positive and synthetic negative controls;
- a fixed, committed sample of both `K7` and `K6`-only graphs;
- a fixed sample of old certified and deferred graphs.

### B. Implement three exact filters

Implement and report separately:

1. `K7_clique_Hall` using
   `omega(G[V_A]) <= |A|` for all 128 masks and every `K7` seed;
2. `K7_disjoint_edge_bounded_cover` using the graph `L`, eligibility
   `|D_x|>=3`, and cover bound 7;
3. `K6_clique_Hall` using
   `omega(G[V_A]) <= |A|+1` for all 64 masks and every `K6` seed of the
   `K6`-only population.

Use exact bit operations. A graph is rejected when any required seed produces
a contradiction. Preserve separate counts for all/certified/deferred and the
union with the existing link/reflection rules.

For speed, there is no need to enumerate every clique explicitly. For each
coordinate mask `A`, collect outside vertices with `D_x subseteq A` and run a
small bitset clique test only up to the forbidden threshold `|A|+1` or
`|A|+2`.

### C. Required controls

At minimum add:

- a realizable `K7` plus one or two reflected facet points, which must pass;
- the realizable 18-point lower-bound graph, which must pass the `K6` rule;
- a synthetic `K7` Hall failure with a required clique larger than its defect
  coordinate union;
- a synthetic `K6` Hall failure that needs more than the defect coordinates
  plus the one normal coordinate;
- a synthetic disjoint-support required edge with both endpoints restricted
  to at most two defects;
- a synthetic bounded-cover failure not reducible to a single forbidden edge,
  so that the cardinality-7 part is exercised.

The existing profile counts and old controls must remain unchanged where the
new rules do not apply.

### D. Measure before choosing the next branch of work

Run the exact profiler over the complete corpus and publish:

- individual and cumulative new rejection counts;
- the new residue split into `K7` and `K6`-only;
- wall and user time;
- source/result/corpus hashes;
- witness details for at least a deterministic diagnostic sample;
- the independent-reference comparison result.

In parallel only if cleanly isolated, run the fast C obstruction-coverage
sample and label every result heuristic until its pattern has a rigorous
certificate.

If the new exact filters remove a substantial fraction, iterate on anchored
support/rank consequences. If coverage is small, the next question is how many
canonical **anchored signatures** remain: canonically label the graph together
with a distinguished `K7` or `K6` seed and count unique signatures before
attempting Gröbner, interval, or Schur-complement solving. Solving one anchored
signature once may amortize over many candidate graphs.

## 6. Next handoff

Do not resume the 903,847-case interval grind during this milestone. Push a
focused implementation commit and a documentation-only handoff commit. Put
near the top of `STATUS.md`:

```text
REMOTE REVIEW HANDOFF
Branch and exact commit
New exact rules and proofs
Commands, compiler, runtime and hashes
Control results
Independent-reference comparison
Individual/cumulative rejection counts
Remaining K7/K6-only split
Heuristic obstruction-coverage sample, clearly labelled
Unresolved implementation or mathematical questions
Recommended next action
```

The user will then ask the remote review session to inspect that exact commit.
