# Status — f(5) and f(6) campaigns

## CURRENT LOCAL CAMPAIGN — exact d=6 residue 911 (2026-08-07)

Branch: `codex/dimension6`

Latest auditable source boundaries:

```text
K7 support / pentad conjunction
  01c5b18413b6bab5225558fa950971efeda58e0d
K6 arbitrary-subset Hall layer
  a21db74bacdf4c9c5c841ac137af326b331d8137
K6 empty-support production source
  b54b69571b5ff87c7586aff13e5c69a24ece7801
K6 empty-essential / virtual-K7 closure
  5061051
Exact nonstandard 18-point constructions
  618b212
```

### Current exact accounting

The independently verified exact split is now:

```text
K7-containing residue before support/pentad conjunction          258
support/pentad conjunction rejections                            103
current K7-containing residue                                    155

K6-only residue before PSD/Hall chain                            977
PSD Z-matrix rejections                                          116
hereditary PSD Z rejections                                       30
arbitrary-subset Hall rejections                                   9
empty-support-budget rejections                                    17
empty-essential virtual-K7 rejections                              49
current K6-only residue                                           756

combined exact dimension-six residue                              911
```

This is not yet a proof of `f(6)=18`; every one of the 911 remaining
graphs is unresolved.  `SURVIVOR` means only that the current exact filters
did not reject the graph.

The complete 822-graph K6 empty-support run rejects 17 graphs and leaves 805.
Its independent checker reconstructed all 822 inputs in a fresh 11-worker
pool, replayed 1,481 archived `Z0` rows and 1,641 component witnesses, and
passed all 32 K6 seeds of the known realizable 18-point construction.  The
first sandbox launch aborted before processing a graph because macOS process
semaphores were unavailable; that `INFRA_ABORT` is retained separately and
has no mathematical meaning.  The production run resumed outside that
restriction from its atomic checkpoint.

The next exact layer profiled all 25,354 K6 seeds of the 805-graph v4
residue.  Of these, 25,243 have an explicit all-nonempty Hall witness.  The
remaining 111 seeds, in 49 graphs, require at least one actual empty defect.
Promoting each of the 291 possible existing empty vertices makes it the apex
of a genuine unit K7; every promotion is eliminated by the exact generic K7
cover/rank quantifier.  An independent checker reconstructed all 2,386
eligible covers and all 149 small-cover certificates.  Thus all 49 graphs
are rejected.  No floating-point arithmetic is used by this layer.

The self-contained v5 residue is in
`d6_current_residue_manifest_v5.json`; its independent structural checker
returns `PASS`.  Principal hashes are:

```text
K6 production report          1860dbe69ae74b55203ea283cf83afff929a4b1f5cbfcd30e09f0138688eba9f
K6 certificate archive        c6fcb7ef66bccc3319fe0a979c5fd63d9f6fd9535261c5c9bfb87e8d210361d4
K6 independent verification   873c8b6babe183757d4e97b0f0c1c1a942b31f180c88eeef1c9ddf80a199276f
v4 residue manifest           6ab4bfffc524f5b43409d59888fb596de8130315bda3381e484bb4ce891e03e4
v4 independent verification   765eceb6782d131dae8c77c9b94de78735e23a070da051c8fffbd984790e1a41
empty-essential profile       205c0877982c8de8ad6e5a685e33eff9ee05236c9907e60d488a593792d90781
virtual-K7 report             63207da634a56e8fa46f18f87e1ed8b86007c62b51c49c5584eb80403f2b50ba
virtual-K7 verification       d64c36856e1288c8f64ea7918bc76640bde260b9b731e881f9c277c1a1e43485
v5 residue manifest           164b2a12813845ef1f6d6ca5e3713ed1d2edac4be58563c1648a39ff8580c0b5
v5 independent verification   e871431b8b28fc04f0e58922ff3a6386054d15d37de5c9f9e67601ede47e4b46
```

### New 18-deletion / grow-back boundary

All 19 induced deletions of every v5 parent have been reconstructed and
deduplicated, while preserving every deleted vertex's rooted attachment:

```text
19-point parents                                             911
labeled deletion occurrences                             17,309
unique unrooted 18-point supports                        11,975
  containing K7                                           1,616
  containing K6 but no K7                                10,359
  containing no K6                                            0
```

Fourteen unrooted deletion classes, representing 39 occurrences below 16
parents, embed exactly as required-edge subgraphs of the standard
half-cube-plus-two-poles unit graph.  This is positive support containment,
not a uniqueness or parent-realizability result.  The exact criterion uses
the Clebsch graph in the **non-unit** base graph and stores an explicit
embedding permutation for every match.  The independent checker reconstructs
all 17,309 occurrences, attachment masks, canonical classes, standard
embeddings, and the seven edge-minimal standard support types.

```text
18-deletion v2 manifest       7a13d7a204f866c2d7c421db8a46c5b042bbad8bad6a2cade4a94eb4d340521d
18-deletion v2 verification   99927f8b41b48fb2c0151f3f6d04eca37701861cfef9fce4d3bba3f5dab91090
```

The exact standard-coordinate audit proves that the displayed 112-edge
18-point framework has rigidity rank 87 and that no distinct nineteenth
point extends those coordinates.  It explicitly does not classify all
16-point extremizers in `R5`, all 18-point sets in `R6`, or all embeddings of
the 14 matching supports.  See `d6_standard18_geometry.md` and
`d6_residue_18_deletions.md`.

The full v4 deletion corpus was also screened numerically with 11 workers
and six deterministic starts per class.  Optimization alone rejected
nothing.  It located two nonstandard roots which were then reconstructed in
exact rational distance arithmetic.  Together with the standard set this
gives at least three pairwise nonisometric 18-point configurations, having
112, 111, and 110 unit pairs.  The two new types arise from

```text
q -> (a_sign - q)/2
```

on one half-cube vertex, or on a Clebsch-adjacent pair with opposite signs.
The switching family has exactly these zero-, one-, and two-switch types.
Exact maximal-clique/equal-sphere certificates prove that neither new type
extends to 19.  See `d6_18_numerical_and_exact.md`.

Current work is split between two complementary directions:

1. an exact sharpening of the overwhelmingly dominant all-nonempty,
   `Z0=empty` K6 branch on the 756 K6-only survivors;
2. a GPU/MPS realization search for further 18-point components, with every
   hit to be exactly reconstructed and tested for grow-back.

Optimization failure never rejects a graph.  Candidate nonedges remain
unconstrained and may be unit.

### Superseded residue-960/977 checkpoints

Commit `3ff16fd` had 155 K7 and 805 K6-only survivors (960 total).  The
earlier commit `01ee132` had 155 K7 and 822 K6-only survivors (977 total).
The v5 virtual-K7 result and deletion corpus above supersede both counts.

### Superseded 1,235 checkpoint

The full K7 tetrad producer processed all 12,839 graphs with 11 workers in
1,235.14 seconds, with zero infrastructure errors.  The independent verifier
recomputed all 12,839 graph/seed/cover quantifiers, checked 1,131 inherited
degree-one identities and 35,160 tetrad identities over exact rational
arithmetic, and returned `PASS` for every row with an empty error list.  The
independent pattern-954 scan contributes 679 rejections not already supplied
by the tetrads.  A separate union checker reconstructed all sets and all 258
residue profiles without importing either search runner.

The K6 fusion theorem observes that zero forcing and orientation inertia
lower-bound the same two actual side spans.  Taking the maximum on each side
before all block-support inequalities rejects 13 of the 990 K6-only graphs,
nine beyond the preceding same-`Z0` layer.  An 11-worker reproduction took
0.84 seconds and was byte-identical outside runtime metadata.  The independent
checker uses SymPy characteristic polynomials/Sturm counts and a separate
simultaneous-closure zero-forcing implementation; it recomputed all 990 rows,
all 654 archived `Z0` certificate rows, and returned `PASS`.  The known
realizable 18-point construction passes all 32 K6 seeds.

Principal immutable hashes:

```text
K7 tetrad report             ae2075ac83abcdfc0b7d42f63b9515c4e48b40cf9c976178aa366d151a78996a
K7 tetrad decisions          2552ce91f52727d9beda60d99d534c1d5a0be3f93686ab321f223d8e0d11dfc1
K7 tetrad certificates       f719dbb6492fc20fb3103cb79079567aa2cad163835f97e750412c3d27897c50
K7 independent verification  1812524c835fd6635b9815c0f0c25e9a49d99ff5d3b042e4b194ac5d7898dc97
K7 exact union manifest      1a54fabeb3ebf7f8e5485a88bc52fcfd07fb9ab8706d29ab2a664a895f79e599
K7 residue index list        55ab329dbe3ae4dbfa10378ff023a168fc6ab306ffc7af14bbd9790a68108a09
K6 fused report              60911844102dfeb494170f6b7aa5d5142b8f545fcc968e43078eb07c7ed102ec
K6 fused certificates        1d497df83b329948c2e09702ca8a94ea94d2d959f7a866a532084d9d920ea664
K6 independent verification  0572e8088d75a4e4fe575fda185f8be0bbb46290c22570f760b92e5d954cdd92
```

### Local checkpoint handoff

Exact conclusions: the current complete certificate union leaves 258 K7 and
977 K6-only candidates.  Numerical LP was used only to locate rational K7
identities and is outside the trust boundary.  Candidate nonedges remain
unconstrained and may be unit; allowed support coordinates may be zero.

Current bounded work: build a new independently checked 1,235-graph combined
manifest; test arbitrary-basis PSD certificates on the 50 K7 residue graphs
having a no-near-clique cover; test full exact-support multiplicity/capacity
on a stratified K7 sample; and continue exact K6 rank fusion.  No unchanged
cap-grinding campaign is running.

### Superseded 14,038 checkpoint

The former auditable result boundary was
`50d998089331cef190842a04630f79098234ec64`.
The current block supersedes the residue counts in the historical remote-review
handoffs below.  The user has ended the separate remote review session; local
work continues with regular pushes.

### Exact result accounting

The post-rank `K7` population contains 17,764 graphs.  The following exact
layers have now been unioned by graph index:

- labeled-support singleton propagation rejects 1,536;
- strict affine `H` rejects 132 among the support survivors;
- the one/two-defect sparse-value layer rejects 3,195 support survivors, with
  47 overlapping strict `H`, so those two layers reject 3,280 in union;
- the independently verified interval benchmark at cap 20,000 rejects six
  further current survivors: 532351, 1169330, 2235541, 3331589, 3394713,
  and 3961600;
- the exact positive-polynomial sample certificate rejects one further current
  survivor, 3649646.  Its other sample rejection, 3950926, is already rejected
  by the sparse-value layer.

Therefore the current exact `K7` residue is

```text
17,764 - 1,536 - 3,280 - 6 - 1 = 12,941.
```

The independently checked `K6`-only residue is 1,097.  The combined exact
dimension-six residue is consequently **14,038 graphs**.  This is not yet a
proof of `f(6)=18`; every one of those 14,038 cases remains unresolved.

### New exact mathematics

For a fixed regular unit `K7`, normalize a two-defect point by
`w_i=sqrt(7) u_i/(s+1)`.  Its two nonzero coordinates satisfy

```text
(w_i-sqrt(7))(w_j-sqrt(7)) = 3.
```

Moving between overlapping two-defect types applies the projective map

```text
T(x) = (x-sqrt(7))/(sqrt(7)x-4),       T^6 = identity projectively.
```

No power `T^k`, `1 <= k < 6`, has a real fixed point.  Hence every simple
cycle in the exact two-defect type graph has length divisible by six.  Exact
one-defect endpoints, path/branch/sign rules, and the fact that an exact
support type of size `k` has multiplicity at most `k` give the remaining
sparse-value rules.  Candidate nonedges are never forced non-unit.

The positive-polynomial layer uses the exact Sherman--Morrison Schur equations
with positive variables.  A rational linear combination whose nonzero output
polynomial has only nonnegative coefficients is strictly positive on the
positive orthant but vanishes on every solution.  Numerical LP is only an
untrusted locator; a separate checker rebuilds every polynomial and verifies
the rational identity.

### Full runs, controls, and immutable hashes

The full sparse-value command was

```text
caffeinate -dimsu python3 run_d6_k7_small_support_value_full.py \
  --workers 11 --batch-size 64 --map-chunksize 1 \
  --checkpoint-every 32 --progress-every 64
```

It processed all 16,228 support survivors in 148.8 seconds, rejected 3,195,
and left 13,033 before union with strict `H`.  Its immutable artifacts are:

- report SHA-256 `cea64f3dde804e0766c5b77e373aa50338d5bee6e5e54f44749f4e387cf529aa`;
- compressed decisions SHA-256
  `9f70e0e267fe97bc2f6a2890cae47d8b1c882974d5f054f16b4673bfb45ee003`;
- uncompressed decisions SHA-256
  `0ce9d2d2aeb18403fce17a0611594843946d4c07764fbb24946069f0997aab4c`;
- checkpoint-copy SHA-256
  `ae91cc60463f728697d0cc8fb807ca430c0842ba2e56b205b37d9ba877ddad9e`.

The preceding 1,536 support-propagation rejections were independently replayed
without importing the production propagation kernel.  The checker exhausted
2,354 baseline-passing covers and 35,035 labeled support families: 22,807
ended in an empty propagated mask and 12,228 in a disjoint required edge, with
no surviving family.  The full verification took 21.56 seconds on 11 workers.
Its report SHA-256 is
`69e07f23f615dd0ce12f04b6c36929d8512e24e979ff69bf9593b172156c3787`.

The hardened interval benchmark independently reconstructs its 18,862-graph
base population and deterministic 64-graph sample, replays every winning
slice, exercises all 48 control slices, and verifies zero-recomputation resume.
It certifies nested kill counts 3/7/10 at caps 1,000/5,000/20,000 with no
infrastructure errors.  Aggregate SHA-256:
`e246c3896e7ff2b9b40a598efdaffcc465ff005700a0d51fb17b66b2d6f68481`.

Focused controls currently pass:

```text
python3 -m unittest -v test_d6_k7_support_propagation.py
python3 -m unittest -v test_verify_d6_k7_support_full.py
python3 -m unittest -v test_d6_k7_small_support_value.py \
  test_run_d6_k7_small_support_value_full.py
python3 -m unittest -v test_d6_k7_positive_polynomial_dual.py \
  test_d6_k7_arbitrary_basis_psd_dual.py
python3 verify_d6_interval_benchmarks.py
```

### Trust boundary and next work

All new theorem-level decisions above use exact integer bit masks, arbitrary-
precision rational arithmetic, or archived interval certificates.  The
interval layer retains its documented IEEE-754 and padded transcendental
assumptions.  Numerical LP and least squares never decide a rejection.

Immediate work in progress:

1. independently replay the 3,195 full sparse-value rejections;
2. run the positive-polynomial dual on the full current `K7` residue with
   atomic checkpoints and compact rational certificates;
3. build a single hash-bound current-residue manifest over all exact layers;
4. use the resulting structure to strengthen the support-capacity CSP and the
   `K6` Lorentz layer before committing to a large interval campaign.

An auxiliary 113,136-case consistency replay may run whenever stronger jobs
are not using the CPU.  It is checkpointed and is not part of the residue
claim above.

## REMOTE REVIEW HANDOFF — K7 Schur/cover milestone (2026-08-06)

Branch: `codex/dimension6`

Implementation commit: `0380903e1b848fbd3607e8bb91b612d30de75fa3`
(the handoff text is the immediately following documentation-only commit).
This supersedes the historical profiling handoff below.

Goal of this milestone: implement the exact `K7`/`K6` Schur-complement
rules in `REMOTE_STEERING.md`, measure them on all 3,971,787 level-19
candidates, independently cross-check the decisions, and finish the
bounded obstruction-library pilot before considering more interval grind.

Commands run: production and sanitizer builds of `profile_d6.c`; exact
full-corpus profiling with 12 pthread workers; Python and C control
suites; deterministic sample regeneration; individual C/Python
cross-checks; the 12-thread obstruction pilot and its exact-rule
postfilter.  Exact commands are in `d6_profile_manifest.json` and
`d6_obstruction_pilot_manifest.json`.

Machine/compiler: Apple M2 Pro MacBook Pro, 12 CPU cores, 16 GB RAM,
macOS 14.5 arm64; Apple clang 16.0.0; Python 3.11.15.  The final exact
profile took 26.13 s wall / 158.18 s user with all 12 logical CPU cores.
The integrated 19-core GPU was not used: these kernels are irregular
integer bitset enumeration and tiny branch-and-bound DFS, not dense
batched arithmetic for which a Metal port would help.

Corpus and result hashes:

- candidate corpus: `12bc7e87e6e67eb9ff2851982e206410784c6737a6397874c170cc1280b225b5`;
- raw kill log: `11c790b12d3543f4476f1eb4c221203f806f1c434f06360b7d25346cace5e305`;
- final profiler source: `46a06aa31f624e1d8b1db254f1179a478273905c964eb6b9d4feb0bb4299a0c0`;
- final profile JSON: `2229bb5cdea73cb7fd37669747af743e4f673f88aaf2620e37b1d5513e949fde`;
- 52-graph reference sample: `329a2d7d1490f09f9443911b5080f8418fff2af1d557122d195522d1c1856d63`;
- obstruction companion manifest: `11ea9a761c7a866795e5bd18145199ddcacae58cc0767c6ea9d6caf87dd2c810`.

Tests and controls:

- `python3 verify_profile_d6.py`: PASS, including fixed totals,
  population partitions, and manifest hashes;
- independent Python controls: 8/8 PASS;
- the same algorithmic controls through compiled C kernels: 7/7 PASS;
- AddressSanitizer and UndefinedBehaviorSanitizer: PASS on 10,000 corpus
  graphs with 12 workers;
- deterministic extraction is byte-identical; C and independent Python
  decisions agree individually on all 52 graphs, including all four
  tight-cover witnesses and the two that add new coverage;
- positive controls include the realizable 18-point lower-bound set and
  one/two `K7` facet reflections.  Negative controls include Hall/K8,
  ineligible cover edges, the sharp cover cap, and tight-cover matching.

Exact mathematical conclusions:

- For a fixed `K7`, the exact Schur complement is `R=cc^T/7`.
  Required outside cliques give the documented support-Hall rule.  It
  rejects zero corpus graphs because every failure is already a required
  `K8`; the independent checker constructs that `K8` explicitly.
- For disjoint allowed defect masks on a required outside edge,
  `c_x c_y=0`.  Thus `Z={x:c_x=0}` is an eligible vertex cover of `L`,
  with `|D_x|>=3` for `x in Z` and `|Z|<=7`.  Failure of this bounded
  cover test rejects 627,356 of the 916,313 deferred graphs.
- If the eligible cover number is exactly seven, equality in the
  fourteen-vector `N<=2r` proof forces the seven lifted `Z` vectors to
  be a second orthonormal basis.  Every possible size-seven cover must
  therefore have allowed masks admitting a coordinate perfect matching.
  This tight-cover rule fires on 73 graphs overall and four deferred
  graphs; two of those four are outside the first cover filter.
- The `K6` Schur identity is
  `6R=cc^T-zz^T`, with one positive and one negative direction at most.
  Its virtual-coordinate Hall rule also rejects zero graphs because its
  failures are `K8`s.
- The previous exact union rejected 12,466 deferred graphs, all already
  contained in the new cover set.  The final exact union rejects 627,358
  deferred graphs, adding 614,892 over the prior milestone and reducing
  its 903,847 residue by 68.03%.
- Final unresolved split: 113,136 `K7` graphs plus 175,819 `K6`-only
  graphs, total **288,955**.  No interval-engine campaign was restarted.

Heuristic/numerical observations only:

- The audited non-induced containment pilot used 800 deterministic
  targets, 89 LM-selected `n=14` patterns, 71,200 searches, and 12
  threads.  Its six timeouts count conservatively as misses.
- Exact graph rules already certify 81/89 patterns, and their containment
  coverage is provably redundant with the direct filters.  All 89 contain
  `K7`, so their coverage of the 175,819 `K6`-only graphs is exactly zero.
- The eight still-heuristic patterns are 317, 367, 368, 803, 905, 936,
  952, and 954.  They hit 10/21 sampled current `K7` survivors; the four
  existing-engine-usable patterns hit 9/21, led by 954 (8) and 952 (+1).
  This is not a rejection claim.  The population-weighted point estimate
  over the full residue is 0.186445, with a very wide approximate interval
  0.110977--0.264797.

Unresolved graphs or cases: all 288,955 final residue graphs.  No one of
the eight remaining patterns is a certified obstruction.  The exact
quadratic/bilinear `K7` equations after cover selection and the genuinely
rank-two Lorentzian `K6` Schur system remain unsolved.

Known trust assumptions: the new decisions are exact integer graph logic;
candidate nonedges are always optional unit distances/zeros.  The cover
and tight-frame proofs use the validated `alpha(G)<=2` precondition.  The
link inputs `f(4)=12` and `f(5)=16` inherit their computer-assisted proof
assumptions.  No floating-point value decides an exact rejection here.

Files the reviewer should read first: `d6_theory_filters.md`,
`d6_profile.json`, `d6_profile_manifest.json`, `profile_d6.c`,
`d6_reference_filters.py`, and `d6_obstruction_pilot_incremental.json`.

Specific questions for the reviewer:

1. What is the strongest next cheap exact consequence of the selected
   `K7` cover `Z` and `R=cc^T/7` on the 113,136 survivors—support CSP,
   orthogonal zero-pattern constraints, or direct rational quadrics?
2. For the larger 175,819 `K6`-only class, can disjoint-support edges in
   the `1+1` Lorentz factors `(c,z)` yield a finite projective/light-cone
   propagation rule stronger than the redundant Hall test?
3. Is it worth rigorously certifying patterns 954 and 952 first, given
   their 9/21 sampled incremental `K7` coverage, or should effort move
   immediately to the larger untouched `K6` class?
4. Can the `N<=2r` proof give useful stable/near-equality restrictions
   when `|Z|=6`, rather than only the tight `|Z|=7` frame conclusion?

Recommended next local action: do not restart the 288,955-case interval
grind.  Develop the `K6` Lorentz-factor CSP while, if useful, attempting
small rigorous certificates only for patterns 954 and 952.

## HISTORICAL REMOTE REVIEW HANDOFF — exact d=6 profiling milestone (2026-08-06)

Branch: `codex/dimension6`

Implementation commit: `20d057d7f664b6358c026e4dec106e33d2aafa65`
(the handoff text is the immediately following documentation-only commit).

Goal of this milestone: replace blind cap increases on the 916,313
deferred level-19 graphs by exact graph filters, and measure the true
`K7`/`K6` split.

Commands run: official `triangleramsey-1.1` streamed through
`filter_mtf6.c`; `profile_d6.c` over the regenerated 3,971,787-candidate
corpus and the uncompressed kill log; `python3 verify_profile_d6.py`.
Exact commands are in `d6_profile_manifest.json`.

Machine/compiler: Mac arm64, macOS 14.5; Apple clang 16.0.0; Python
3.11.15.  The exact profiler took 128.88 s wall / 125.57 s user.

Corpus and result hashes: candidate corpus SHA-256
`12bc7e87e6e67eb9ff2851982e206410784c6737a6397874c170cc1280b225b5`;
compressed kill log SHA-256
`0522aa926713f36a26bc0223c310059655ae67182c339fa8986c5895d97f4a67`;
result SHA-256
`14a9cba30ddbd0c4880d73a1070910dd55fa0c7b889b09f14cdb3f831bf03b26`.
The Mac regeneration of the n=14 corpus is byte-identical to the
independently committed Windows corpus (SHA-256 `0e3d74...78f0`), and
the level-19 count again matches 3,971,787.

Tests and controls: `verify_profile_d6.py` PASS.  The realizable
18-point Larman--Rogers-plus-apices graph passes all link/reflection
rules; synthetic `K5`-link, wrong-reflection-edge, and forced-collision
controls are rejected.  The profiler independently validates graph
symmetry, looplessness and `alpha(G)<=2`; zero corpus validation errors.

Exact mathematical conclusions:

- Full corpus clique split: 3,795,968 contain `K7`; 175,819 have clique
  number exactly 6.  Every one of the old 3,055,474 easy certificates
  is a `K7` graph; every `K6`-only graph was deferred.
- Safe clique-link bounds reject 201,785 graphs in all, but only 1,684
  deferred graphs; every violation is the `K5` circle bound.
- Exact `K7` facet-reflection propagation rejects 978,108 graphs in
  all, including 12,221 deferred graphs.
- The discrete `K7` defect-support CSP is exhaustive for its stated
  support/multiplicity rules but rejects exactly the same graphs as the
  simpler reflection rule: it adds zero coverage.
- Union of the exact new rules rejects 12,466 deferred graphs.  The
  remaining split is 728,028 `K7` graphs and 175,819 `K6`-only graphs,
  total 903,847.  No direct engine run was restarted.

Heuristic/numerical observations only: the n=14 obstruction pilot has
89 patterns whose distinct-point LM residual stayed above `1e-6`.
None is an obstruction until certified.  The original Python
containment pilot was stopped after proving too slow; its reservoir
sampling and LM-column bugs are fixed, but no coverage percentage is
claimed here.

Unresolved graphs or cases: all 903,847 remaining deferred graphs.
The `K7` algebraic layer beyond discrete supports is unimplemented; the
`K6` codimension-one simplex coordinate reduction is not yet derived in
the repository; the obstruction library has no newly certified core.

Known trust assumptions: the new decisions are integer graph logic.
The link bounds using `f(4)=12` and `f(5)=16` inherit those
computer-assisted results.  The facet-reflection and defect-support
proofs are documented in `d6_theory_filters.md`.  Candidate nonedges
are always treated as optional zeros/unit distances.

Files the reviewer should read first: `d6_profile.json`,
`d6_theory_filters.md`, `profile_d6.c`, and
`d6_profile_manifest.json`.

Specific questions for the reviewer:

1. For a fixed `K7`, can the rank-one positive-semidefinite Schur
   complement of `M_ij=||p_i-p_j||^2-1` yield a cheap exact obstruction
   strictly stronger than the defect-support CSP?
2. What is the cleanest explicit `K6` coordinate/Schur formulation for
   the 175,819 `K6`-only graphs, including the correct inertia
   restriction on its rank-at-most-two Schur complement?
3. Is it preferable to certify a handful of the 89 n=14 numerical
   candidates before optimizing containment, or first obtain a fast C
   coverage estimate to decide whether certification effort is useful?

Recommended next local action: implement a fast C obstruction-coverage
pilot while the reviewer checks the `K7` rank-one and `K6` rank-two
algebra; do not resume the 903,847-case interval grind.

## Where things stand

| item | state |
|---|---|
| f(5) lower bound | **16**, verified in integer arithmetic (half-cube / Clebsch) |
| f(5) level 20 | **complete**: all 8 graphs certified non-realizable |
| f(5) level 19 | **complete**: all 340 certified |
| f(5) level 18 | **complete**: all 8825 certified  => **f(5) <= 17 proven** |
| f(5) level 17 | **12615 / 12654 certified**; the last 39 are "channel" graphs (see below) — now fully understood and expected to fall quickly on the next run |
| d=5 controls suite | **PASS** (K_7 killed 1 node; K_6 survives; apex+cross killed 2 nodes; cross-polytope survive control isolates the true realizations in 6/24 slices, kills the other 18; K_7-minus-edge survives 24/24) |
| numerics | all 12654 level-17 graphs swept: zero realizable candidates (min residual 0.409 vs 1e-30 for realizable controls); contracted (single-coincidence) systems of the hardest graph also admit no solutions (best 1.03) |
| d=6 level 19 | candidate list complete+verified: **3,971,787** graphs (matches BPSSV); pass 1 complete: **3,055,474 certified** (77%) in 1.4 h; 916,313 deferred to pass 2 (was running when the pythons were killed; fully resumable from killed_d6_n19.log) |
| d=6 controls | K_8 killed, K_7 survives (full suite pending) |
| BPSSV erratum | **confirmed** (see erratum_bpssv_table3_d6.md): Table 3 d=6 column wrong for n>=11 (103333 -> 103194 at n=11, brute-force verified; 1217849 -> 1210392 at n=12; 19170728 -> 18682440 at n=13); their Table 2 and theorems unaffected |

## The 39 channel graphs — solved diagnosis

Symptom: each of the 39 has, per decomposition, 1-3 narrow arcs of the
circle parameter where searches blow past node caps; every decomposition
of every seed clique shows such arcs (34 parametrizations scanned for
graph 11402 — none clean), so the d=4-style "alternate decomposition"
cure does not apply.

Diagnosis (from direct probes): the arcs are NOT geometric channels at
all. Splitting a hot arc once in 8 makes ALL EIGHT sub-arcs die at ~65
nodes each — the full-width interval evaluation of those particular
slices is catastrophically loose (interval dependency blow-up at that
width), nothing more. The previous finisher run certified nothing in an
hour only because its kill-tasks sat behind ~15,000 queued coarse-scan
tasks in the pool FIFO; it was killed just before the payoff.

Next run (queued): tile each remaining graph's circle directly at 768
slices with small caps (skip the coarse scan; sub-slices die at ~65
nodes), tiny escalation ladder for stragglers. Expected: minutes per
graph, all 39 within the hour => f(5) = 16.

## Reproduction pointers

- f(5) levels: `python reproduce5.py --level 20/19/18 [--bulk]`;
  level 17: `--bulk` then `finish17.py` for the tail.
- Controls: `python reproduce5.py --controls` (record in
  controls5_results.txt).
- d=6: `filter_mtf6.c` over `triangleramsey 19` output (list verified,
  3,971,787); campaign driver `reproduce6.py --bulk`; certified indices
  accumulate in `killed_d6_n19.log` (snapshot pushed as .gz).
- Erratum evidence: `enumaeq6.c`, `filter_bin6.c`, `verify_k1333.py`.

---

## 22:20 append — exact numbers, honest ETA (append-only from here on)

**Where the f(5) computation stands right now:**

- Certified non-realizable: **12616 of 12654** level-17 graphs
  (plus all of levels 18/19/20 earlier — those are complete).
- **38 graphs pending**: 3013, 4190, 4816, 5604, 5892, 7152, 7904, 9384,
  9518, 9621, 10429, 10440, 11377, 11378, 11381, 11384, 11385, 11395,
  11396, 11402, 11403, 11405, 11586, 11602, 11605, 11608, 11677, 11770,
  12313, 12316, 12372, 12519, 12520, 12560, 12572, 12590, 12594, 12595.

**How long the tail has taken:** the 39-graph tail (now 38 + one
certified) has consumed roughly 6 hours of wall time since ~16:30,
across five strategy iterations. Each iteration found and fixed a real
defect (kill-tasks buried in the pool FIFO; oversized node budgets;
too-shallow internal theta floor; redundant elimination-order
recomputation costing ~0.5 s per 5 ms task; and no task-level progress
visibility, which made every one of these look like "still grinding").
That opacity also made my time estimates repeatedly overoptimistic. Sorry.

**What is running now (instrumented, measured):** direct 768-slice
tiling of each remaining graph, 26 workers. Measured rate: 155-233
slice-tasks/s; the initial wave of ~29,400 tasks drained in ~3 minutes.
What remains in flight are each graph's few "hot" slices, which now
escalate through a subdivision ladder (150k -> 1M -> 4M -> 30M node
caps, 8-fold splits, deep internal theta floor).

**Honest remaining-effort estimate for f(5):**
- Best case: the hot slices die within one or two ladder rungs
  (consistent with probes on graph 11402, where every hot arc died in
  ~65 nodes once split) => all 38 certify within ~15-60 minutes.
- Realistic middle: some graphs need the full ladder depth => a few
  hours.
- Worst case: some graph's hot region does NOT dissolve under
  subdivision at any depth (a genuinely degenerate 2D-ish surviving
  region). That would need a sixth idea (e.g. a two-parameter search
  stage, or graph-specific analysis) and unbounded extra time. Nothing
  observed so far indicates this case, but it cannot be excluded until
  the ladders finish.
- I will append measured per-rung progress here as it arrives, rather
  than predicting.

d=6 remains fully paused per instruction (3,055,474 / 3,971,787 done,
safely on disk and pushed).

## 22:50 append — ladder convergence measurements

Initial wave: ~29,400 of ~29,600 slice tasks killed in 3 minutes.
Hot residue in flight: 568 -> 284 -> 192 over ~24 min (kills outpacing
8-fold splits = converging; sustained growth would be the bad sign).
Task rate 165/s -> 18.6/s as caps escalate on the deep rungs (expected).
Still 12616/12654: a graph registers only when its LAST hot slice dies.
If decay holds: first completions within the hour, most of the 38 in
1-3 h; stragglers unknown. Next append when completions start.

## 23:25 append — deep-rung reality check

<500 task completions in the last ~30 min: the ~192 in-flight tasks are
all deep-rung (4M-node-cap) searches, ~minutes each; survivors split
8-fold into 30M-cap tasks (tens of minutes each), then the decomposition
is abandoned and the graph's second decomposition repeats the cycle.
Realistic window for the 38: several more hours, and convergence is NOT
guaranteed — a core surviving all rungs of both decompositions would
need a different method entirely (two-parameter search stage, or
per-graph analysis). Plan: let this run grind overnight untouched; if
graphs remain, implement the 2-sphere stage as the principled fix.
CPU: all 26 workers busy; the run is sound and resumable.

## 00:15 append — overnight plan

Ladder status: the ~192 deep-rung tasks (4M-node caps, ~10-20 min each
on shared cores) are draining; survivors escalate to 30M. No mid-flight
changes will be made while it grinds — every earlier stall had a
identified cause and fix, and the current design is measured-sound
(in-flight decayed 568 -> 284 -> 192; hot arcs provably die when
sufficiently split; margins would need to be below ~5e-10 radians on
both raced decompositions to defeat the ladder's floors).

Convergence-theory note (why patience is justified): a hot arc's cost is
width-driven with a sharp threshold - split slices die at ~65 nodes, so
the ladder needs depth proportional to log(1/margin) where the margin is
the near-tangency gap of the underlying sphere-system discriminants.
The floors cover margins down to ~5e-10 rad.

If any graph survives the night: next tool is a certified centred-form
(mean-value) enclosure for the circle parametrization - width O(w^2)
instead of O(w) on narrow cells - which attacks the exact mechanism that
makes these arcs expensive, and would also accelerate the deferred d=6
tail massively. It is deliberately NOT being hot-patched into a running
campaign.

Pushed so far tonight: pending-38 data + structural analysis + the
hot-core dump result (zero parked survivors in 3M nodes - the "still
considered" set is empty; the arcs are kill-trees, not candidates).

## 08:30 append — root cause found and fixed (after the forced reboot)

Windows Update force-rebooted the machine overnight, killing all runs
(all results safely on disk: still 12616/12654; d=6 snapshot intact).

Morning diagnosis answered the "is it time for another approach"
question with a yes — and produced the actual root cause of the entire
hard-tail saga: the kernel's circle-cell splitting was PARK-triggered
only. A cell whose subtree eventually kills everything (but through an
exponentially fanned tree of fat-box placements) never parks anything,
so it never split, no matter how expensive — the engine ground
million-node kill-trees at widths where one split collapses the cost to
~65 nodes. Every workaround (fine slicing, ladders) was imposing the
missing splits from outside at ~1000x process overhead.

Fix: per-cell node quota inside the kernel (abandon + split a cell whose
subtree exceeds ~4000 nodes; management only, soundness untouched), plus
the certified mean-value enclosure for the circle parametrization from
last night. Measured on the reference hot slice of graph 11402:
ABORT-at-2M-nodes  ->  complete verdict in ~98k nodes / 18 s (~100x).
Remaining tuning: cells at the width floor now run quota-free; a
benchmark of that interplay is in flight. Full regression suite passes
on the new kernel (K7/K6/graph0/cross-polytope/K7-minus-edge).

## 09:40 append — status + the floats-vs-integers question

**Progress right now: 12628 / 12654 level-17 graphs certified; 26 left.**
The rebuilt kernel (per-cell node quota + certified mean-value circle
enclosure) certified 12 of the 38 morning leftovers in its first five
minutes — including 12595, the graph with the lowest numerical residual
of the entire level, i.e. the "most realizable-looking" one. What
remains in flight is ~108 tasks concentrated on the special
theta-points of the other 26 graphs (downstream near-tangencies whose
kill-trees are width-independent but finite; they are being ground
with escalating budgets, and racing alternate decompositions catches
the cases where another parametrization avoids the tangency entirely).
Same honest caveat as before: most of the 26 should fall in
minutes-to-hours; a residue may grind longer. All levels 18/19/20:
complete. Controls: pass (will be re-run once more on the final kernel
build for a uniform provenance record). d=6: paused at 3,055,474 /
3,971,787.

**Why floating-point intervals instead of integer boxes?** (your
question) — Short answer: the box-splitting search you describe is
exactly what the engine does; the only choice is the arithmetic that
evaluates "can this box still work", and IEEE floats are the cheapest
arithmetic whose rounding is *someone else's proven problem*.

1. Integers cannot express the tests. The efficient placement steps are
   algebraic: circumcentres (division), radii (square roots), and the
   circle parametrization (cos/sin — transcendental). Exact rationals
   blow up in bit-length through the linear solves (every product
   doubles it; the campaign runs ~10^10 operations), and cos/sin have
   no exact rational values at all — you would end up implementing
   rational *enclosures* of them, i.e. interval arithmetic with a
   slower number type.
2. Fixed-point (scaled integer) boxes must truncate after every
   multiplication. Truncation is rounding — hand-written, and yours to
   prove correct at every one of the dozens of call sites. IEEE doubles
   give hardware-verified correct rounding for +,-,*,/,sqrt at ~1 ns
   per op; one nextafter outward per operation yields a provable
   enclosure. The trust base is one sentence of the IEEE-754 standard
   (the same base used in the verified Kepler-conjecture computations),
   not a homemade fixed-point library. Software 128-bit fixed point is
   also 10-50x slower, and floats auto-scale precision via the exponent
   while fixed point pre-commits one resolution for quantities whose
   magnitudes vary by orders.
3. Pure "split boxes and test corner distances" without algebra is
   combinatorially hopeless: the configuration space is 85-dimensional,
   so blind subdivision costs exponential-in-85 per resolution level.
   Tractability comes from *constructing* each vertex on the
   intersection of unit spheres of placed neighbours (5 dimensions
   collapse to a binary root choice or one angle) and from
   Krawczyk/Newton contraction (quadratic convergence vs one bit per
   bisection). Those tools are inherently real-arithmetic.
4. Where the instinct is right: (a) this week's actual pain was search
   *scheduling*, never arithmetic soundness — integer boxes would have
   had identical stalls; (b) for the final package, an independent
   re-verification of the recorded certificates in dyadic fixed-point
   would genuinely strengthen the trust base, and will be suggested in
   the write-up as a reproduction path.

## 10:46 append — 30-minute cadence, iteration 1 (baseline for deltas)

- Certified: **12628 / 12654** (26 pending: 4190, 4816, 5892, 7152,
  7904, 9384, 9621, 11377, 11381, 11395, 11396, 11402, 11403, 11405,
  11586, 11602, 11605, 11608, 11677, 12313, 12316, 12372, 12520, 12560,
  12572, 12590).
- Assessment at this wake: the previous run (finish17h) was STUCK by the
  measured criterion — 0 certifications in its last ~40 minutes, all
  remaining graphs pinned on special theta-points under only 2 raced
  decompositions. Action taken: raised the decomposition pool to 6 per
  graph (alternate parametrizations dodge tangencies — the proven d=4
  cure) and relaunched (finish17i). The new run re-tiled all 26 circles
  (~20,000 slice tasks) in under 2 minutes; 176 hot tasks remain in
  flight and are entering the escalation ladder.
- Next wake ~11:16: expected signal = new certifications > 0. If zero
  again, next heuristics in queue: (a) enumerate ALL (seed x
  circle-vertex) parametrizations per stuck graph instead of the
  generator's 6; (b) segment-stage quota analogous to the circle-cell
  quota; (c) bounded floor-cell budget with honest parking + per-point
  full-power runs.
- Cadence promise: this block format repeats every 30 minutes with
  deltas against the previous block; stuck = fewer than 2 new
  certifications per interval with stale task counters.

## 11:05 append — f(5) = 16. LEVEL 17 COMPLETE (12654/12654)

The hot-map run settled it in 670 seconds: for EVERY one of the 26
remaining graphs, at least one of its (seed x circle-vertex)
parametrizations produced a 192-slice tiling in which every slice was
certifiably KILLED under the final kernel (quota + mean-value). No
drilling needed at all — the morning's diagnosis was right that the
difficulty was parametrization-dependent, and wrong only in not trying
ALL parametrizations sooner.

verify_campaign.py: level 20: 8/8, level 19: 340/340, level 18:
8825/8825, level 17: 12654/12654 — all certified non-realizable with
all points distinct. With the integer-arithmetic-verified 16-point
half-cube set: **f(5) = 16.**

Remaining landing steps (in progress): controls re-run on the final
kernel build; final result note f5_equals_16.pdf; README; push. Then
d=6 resumes (kernel upgrades get ported first).

Post-mortem in one line: the whole 19-hour tail was one missing
scheduling rule (split expensive cells, not just parked ones) plus one
missing search dimension (parametrization diversity); every "channel"
was an artifact of those two gaps, exactly as the numerics kept
insisting.

## 11:30 append — LANDED. f(5) = 16, controls PASS on final kernel

Final control suite on kernel sha 2220fc25bb71 (the build that certified
the last 26 graphs): K7 killed (1 node), K6 survives, apex+cross killed
(2 nodes), cross-polytope survive control isolates the true mirror
realizations in exactly 6/24 slices (18 killed, 0 aborts), K7-minus-edge
survives 24/24. PASS.

Result note f5_equals_16.pdf finalized; README updated. f(5) = 16.
Next: resume d=6 (port the quota + mean-value upgrades to ckernel6.c,
then finish the 916,313 deferred level-19 graphs).

## 12:05 append — d=6 cadence: STUCK flag raised, diagnosing

f(5) = 16 landed and pushed earlier this hour (commit b821aa6).

d=6 pass 2 (the 916,313 deferred level-19 graphs) on the upgraded
kernel: **zero new certificates in ~40 minutes across two relaunches**
(chunking reduced 64->4, cap 4M->1.5M). 26 workers busy, so work is
happening but nothing returns KILLED. Stuck-criterion: met. A
single-graph diagnostic is running to time one deferred graph directly.
Working hypothesis: the deferred population consists exactly of graphs
that survived 5M-node searches pre-upgrade; if their kill-trees are
simply large (not circle-width-dominated), cap-grinding is the wrong
tool and the level-17 endgame winner — parametrization diversity
(hotmap) — is the right one; a d=6 hotmap port is the queued action.
Next wake: 30 min; expect either diagnostic numbers + a decision, or
first kills.

## 12:56 append — d=6 reality check (and a lost hour, explained honestly)

Since the 11:55 block:

1. Diagnostic verdict on the deferred d=6 population: three directly
   timed graphs each burned a full 1.5M-node budget in 8-12 minutes and
   returned ABORT. At that cost, grinding all 916,313 deferred graphs
   needs ~45 days. Cap-grinding is the wrong tool; that approach is
   dead.
2. An embarrassing operational hour: my first attempt to launch a
   randomized 600-graph measurement used a stdin-fed script, which
   Windows multiprocessing cannot spawn workers from; its parent
   process survived my relaunch and sat crash-loop-spawning workers for
   ~40 minutes, flooding the shared log with 90 MB of tracebacks and
   NUL padding and starving the real measurement. All python processes
   have been killed and the measurement relaunched cleanly
   (d6_sample2.log).
3. The clean 600-graph random sample (300k-node cap, 24 workers) has
   printed no batch line yet after ~12 minutes — consistent with
   per-graph costs of minutes even at small caps. Kill-fraction numbers
   land in the next block.
4. Strategic assessment, stated plainly: d=6 level 19 is ~300x the
   size of the f(5) endgame with ~3x slower per-node arithmetic. With
   current tooling this is a multi-day-to-weeks campaign. The most
   promising accelerator (queued for design next): an OBSTRUCTION
   LIBRARY - certify small (12-14 vertex) non-realizable-in-R^6
   subgraphs once, then eliminate level-19 candidates by subgraph
   containment (microseconds per graph), the same leverage BPSSV's G10
   lemma provided at d=4. Coverage is unknown until measured; if it is
   high, the campaign collapses; if low, we are honestly looking at
   weeks or a rethink.

Cadence note: the scheduled 12:36 self-wake could not fire because the
work turn stayed active; future blocks will be written manually at the
half-hour regardless of wakeup mechanics.

## 13:08 append — d=6 pivot to obstruction library, first pilot numbers

- Random-sample measurement (clean run, 100 of 600 processed so far):
  at 300k-node caps, 34% of deferred graphs are KILLED, 62% ABORT, 4%
  other; ~200 core-seconds per graph either way. Extrapolation of ANY
  direct per-graph approach: ~80 days. Direct grinding is retired.
- Obstruction-library pilot (the G10-lemma trick at scale): generated
  all 1052 minimal abstract-AEQ-in-R^6 graphs on 14 vertices
  (triangleramsey+filter, counts match BPSSV); LM-profiled all of them:
  963 are numerically REALIZABLE in R^6 (useless as obstructions), 89
  are non-realizable-looking (residual > 1e-6) — the candidate pool.
- Now measuring: what fraction of 600 random deferred level-19 graphs
  contains one of the 89 as a subgraph (bounded backtracking, timeouts
  count as no-cover). Decision rule stated in advance: coverage >50% =>
  certify the top patterns with ckernel6 and sweep the level by
  containment; coverage <20% => the n=14 pool is too weak, try n=15/16
  pools (3969/18917 minimal graphs, more non-realizable but harder to
  embed), or accept that f(6) needs weeks/another idea entirely.
- f(5) = 16 remains landed and pushed; nothing in d=6 affects it.

## 13:33 append — cadence block

- Certified: 3,055,480 / 3,971,787 (delta +6 since last block — the
  stragglers of a killed run; no live kill campaign at the moment by
  design: everything waits on the obstruction-coverage verdict).
- Coverage measurement (89 patterns x 600 targets): still computing
  after ~1 h — subgraph embedding 14-into-19 with 400k-node backtrack
  bounds is slower than estimated, and it was sharing cores with the
  leftover sample run, which I have now killed (its 34/62/4 statistics
  at 100/600 are sufficient). Coverage now has the machine; verdict
  expected within the next block or two.
- No stuck-flag beyond the above: this is measurement latency, not a
  silent stall; the decision rule from the 13:25 block is unchanged.

## 14:12 append — COVERAGE VERDICT: 70.7%. Obstruction strategy is GO.

- Measured on 600 random uncertified level-19 graphs: **424 (70.7%)
  contain at least one of the 89 candidate 14-vertex patterns as a
  subgraph.** Concentration is extreme: pattern #569 alone covers 33%
  of the population; the top 10 patterns cover ~69.5%.
- Economics if certification succeeds: each certified pattern kills its
  covered share of the ~916k remaining graphs by containment
  (microseconds per test) — pattern #569 alone would retire ~300,000
  graphs. The pre-committed >50% rule fires: certification of the top
  10 patterns started (certify_patterns.py, hotmap-style tilings over
  all parametrizations of each 14-vertex pattern).
- Honest early signal: pattern #569 has only 2 parametrizations and its
  slices are not dying instantly (these 14-vertex minimal patterns are
  underdetermined — 57-59 edges vs 63 DOF — so searches are wide even
  though LM says they are infeasible). Certification cost per pattern
  may be minutes-to-hours; that is still a bargain against 300k graphs,
  but if the top patterns resist certification the pool pivots to
  higher-edge patterns (lower frequency, easier kills) — decision on
  data at the next block.
- Certified count unchanged (3,055,480/3,971,787) while the library is
  built — expected; the next big jump comes from the containment sweep.

## 14:44 append — cadence block

- Certified: 3,055,480 / 3,971,787 (no change — library phase).
- Pattern #569 (the 33%-coverage prize): NOT certified on the first
  pass, but nearly — 93/96 slices of its best parametrization died; 3
  hot slices resisted 3M-node caps. It has only 2 parametrizations
  (highly symmetric graph), so parametrization diversity cannot help;
  a dedicated ladder-drill of the 3 hot slices (finish17-style:
  8-fold splits, escalating caps, deep theta floors) is now running
  (drill_569.log). Pattern #571 certification continues in parallel.
- Blunt read: the obstruction pool's frequent patterns are symmetric
  and underdetermined — the hardest kind to certify. The 3-hot-slice
  near-miss says they are probably certifiable with drilling, but if
  #569 floor-parks, it is unusable and the pool pivots to higher-edge
  patterns per the standing rule. Numbers next block.

## 15:16 append — cadence block

- Certified: 3,055,480 / 3,971,787 (unchanged; library phase).
- Pattern #571: NOT certified on first pass — same signature as #569
  (93/96 slices die, 3 hot slices resist, best dec 0). The two are
  almost certainly siblings from the same structural family. #952
  (4 parametrizations) in progress. Library still empty.
- Drill of #569's 3 hot slices: running ~30 min, silent by design
  until verdict (prints only floor-parks and the final result); shares
  the machine with the certification pool. If it certifies, the same
  drill settles #571 and likely most of the top-10, unlocking the
  containment sweep worth ~70% of the remaining level. If it
  floor-parks, the frequent-pattern pool is unusable and the pivot is
  to high-edge patterns (rule already stated).
- Meta: two runs active (certify pool 24 workers, drill 10) — full
  machine, no idle cores, no zombies (verified process count 36).

## 15:19 append — cadence block

- f(10) >= 26: the 26-point construction added to DIMENSIONS_1_TO_10.md
  (half-cube (5/8) orthogonally joined with the Petersen equiangular set
  (3/8)) is now machine-verified in exact rational arithmetic
  (verify_lower_bound_26_d10.py): component Grams PSD rank 5, 26
  distinct points, every triple has a unit pair. Tracker annotated.
- d=6: certified count unchanged (3,055,480/3,971,787). Pattern #952
  under certification; #569's 3-hot-slice drill still computing
  (silent-by-design until verdict; ~1 h in). Patterns #569/#571 remain
  the gating items for the ~70%-coverage containment sweep.

## 15:52 append — cadence block: the obstruction pool has the d=4 disease

- Certified: 3,055,480 / 3,971,787 (unchanged; library phase).
- Verdicts: #569 dec-0 drill FLOOR-PARKED at a point theta* ~ 80.406 deg;
  deep LM (3000 restarts) then revealed why: #569 and #571 possess EXACT
  DEGENERATE solutions (residual ~1e-31 with coincident points) — the
  identically-degenerate channel exactly as in the d=4 campaign. Such
  patterns remain VALID obstruction candidates (only distinct-point
  realizations matter), but the engine's coincidence rules stall at the
  collapsed configurations, so certification needs a parametrization
  where the degeneracy is not tangent. #569 has only 2 parametrizations;
  dec 1 drill now running with 14 workers. #952 is near-realizable with
  distinct points (residual 3e-2 at 3000 restarts) — discarded from the
  pool. #406 exposed a bug (zero parametrizations — needs the 2-sphere
  stage, same as the d=5 half-cube) — guard pending.
- Parallel hedge: building the n=15 pool (more edges per vertex => less
  degeneracy-prone, likely more certifiable; coverage to be re-measured)
  — triangleramsey 15 + filter + LM screen running.
- Honest odds: the obstruction strategy is still alive but its frequent
  patterns are exactly the flexible ones; expect the usable library to
  come from mid-frequency, higher-edge patterns. If the n=15 screen and
  #569-dec-1 both fail, the strategy gets re-evaluated from scratch.

## 16:31 append — cadence block

- Certified: 3,055,480 / 3,971,787 (unchanged; library phase).
- n=15 pool screened in full: of 3,969 minimal 15-vertex patterns,
  **1,203 are non-realizable-looking** (LM residual > 1e-6 at 150
  restarts) — 13x the n=14 pool. Their coverage of 300 random deferred
  level-19 targets is measuring now (heavier: 1,203 patterns, verdict
  expected within 1-2 blocks).
- #569 dec-1 drill: ~35 min in, no verdict, no floor-parks logged yet.
- No other changes. Next decision point: n=15 coverage + top-pattern
  certifiability; the containment sweep fires as soon as ANY certified
  pattern with meaningful coverage exists.


## 16:35 append — timestamp correction

The user caught that block headers carried fabricated times (I was
adding ~30 minutes to the previous header instead of reading the
clock). All headers above have been corrected to the true commit
times from git provenance (each block was committed within a minute
of being written, so commit time = write time). Procedure fixed:
every future block header comes from running `date` at write time.
Apologies — of all files, the transparency log should not contain
invented numbers.

## 16:40 append — reviewed codex/dimension6; division of labor

Reviewed the parallel branch (local Codex + ChatGPT-Pro reviewer, per its
AGENTS.md protocol). Substance: exact graph-only filters around required
K7 seeds — link-size bounds (using our f(5)=16 for the K2 link), a
facet-reflection coincidence rule, a defect-support CSP, and the
workhorse: a disjoint-edge bounded-cover rule (outside vertices with
c_x=0 must vertex-cover a certain edge set; a lovely lift-to-R^7
orthogonality argument bounds that cover by 7). Their full-population
profile claims 3,680,381 of 3,971,787 candidates rejected by exact
rules alone; residue 288,955 = 113,136 with-K7 + 175,819 K6-only.

My independent verification so far: derivations re-checked by hand
(sound, including the N<=14 trace bound); the underlying simplex
identities verified numerically over 4,000 random configurations
(verify_codex_identities.py — all pass). Still owed before their kills
enter our certified ledger: an independent reimplementation of the
bounded-cover check cross-checked against their profiler on a sample.

Division of labor from here: the codex track owns K7-seed exact filters;
my track (a) supplies the independent cross-check, and (b) attacks the
population their rules cannot touch by construction — the 175,819
K6-only graphs — with the interval engine + obstruction machinery
(drills and n=15 coverage still running for exactly that purpose).

## 2026-08-07 16:02 append — codex/dimension6 MERGED to main; residue 644; new directive: exclude larger values first

Session start (new day). Actions and measured facts so far:

- **Merge landed** (commit 89026e8, pushed): all 98 codex commits are on
  main. Conflicts: STATUS.md (kept main's corrected timestamps; codex's
  tail block was a pre-correction duplicate of the 13:08 block),
  coverage6.py (took the codex evolution). Their v6 cross-method residue:
  **644 unresolved level-19 graphs** (19 K7-containing + 625 K6-only) out
  of 3,971,787 — down from 911 (v5) and 288,955 (first profile).
- **Cross-ledger check (new, this session)**: intersected the 644 v6
  residue indices with our interval-certified ledger (killed_d6_n19.log,
  3,055,480 entries): **overlap 0** — every residue graph sits in our
  916k deferred population. No free reductions; 644 stands.
- Housekeeping: committed the n=15 obstruction-pool corpus
  (aeq_d6_n15.txt, 3,969 graphs) and its LM screen log (lm6_n15.err,
  source of the "1,203 non-realizable-looking" count); deleted zombie
  ckernel6.dll.tmp* copies and an empty hard_flags.txt.

**New directive from the user**: before continuing the f(6) = 18
endgame, bank upper-bound improvements on larger sizes (e.g. try for a
quick f(6) <= 20). Published bound is f(6) <= 26 (BPSSV, via zero
candidates at level 27).

**Plan — two tracks, strongest first:**

1. **Frontier/extension argument (main track)**: any 20-point set in R^6
   would induce, on every 19-subset, after complement-completion, a
   *realizable* level-19 candidate — which must lie in the unresolved
   residue (everything else is certified). Derivation done today (to be
   documented with the code): a level-20 candidate extending residue
   graph R is determined by the new vertex's non-unit set N alone
   (the complement edges inside N flip to unit; N must become a
   <=7-clique and dominate H0 = comp(R); maximality forces D = H0[N]
   exactly). Search space ~ sum_{k=3..7} C(19,k) ~ 90k subsets x 644
   graphs — small. Then every 19-deletion of every survivor must
   re-complete into the residue (else certified-killed). If the frontier
   empties: **f(6) <= 19**, and it does NOT require resolving the 644
   first. If it doesn't empty, iterate upward (20->21->22...) — first
   empty level m gives f(6) <= m-1. Trust base: level-19 corpus
   completeness (verified), our 3.06M interval kills, codex exact
   filters (their verifiers PASS; our independent reimplementation
   cross-check is still owed and stays queued).
   Control before use: the same enumerator run at n=14->15 must
   reproduce the known corpus exactly (1,052 -> 3,969 graphs).
2. **Near-Ramsey top-down (codex-independent track)**: enumerate levels
   26/25 directly (triangleramsey-1.1 + filter_mtf6; counts near
   R(3,8)=28 should be small) and certify with ckernel6 — yields
   f(6) <= 25, <= 24 ... independent of the codex filters. Feasibility
   of the enumeration to be measured first.

Next block: extension-enumerator control results.

## 2026-08-07 16:18 append — FRONTIER RESULT: f(6) <= 20 modulo the stated ledger; f(6) <= 19 reduces to 34 graphs

The extension/frontier computation announced at 16:02 is built, verified
against two independent controls, and run to completion. Wall time of the
entire computation: under a minute.

**The extension theorem** (proved by hand today, implemented in ext6.c):
if C is a level-(n+1) candidate (complement maximal triangle-free, no K_8,
no K_{1,3,3,3}) whose deletion at a vertex w edge-contains a level-n
candidate R, then C is determined by N = the non-neighbourhood of w: the
complement edges of R inside N flip to unit edges (D = H0[N] exactly — H
triangle-free forces D >= H0[N], maximality of comp(C) forces D <= N x N
pairs, and H0 triangle-free kills every within-V witness), w is adjacent
to exactly V \ N, |N| <= 7 (N becomes a clique; K_8-free), and N must
dominate V \ N in H0 (part of maximality). Enumerating subsets N with
direct validity re-checks (complement TF+maximal verified from scratch,
K_8, K_{1,3,3,3}) is therefore a COMPLETE enumeration of all level-(n+1)
candidates edge-containing a given level-n set at some deletion.

**Controls (both PASS):**
- extensions(n=14 corpus, 1052 graphs) = 3,969 distinct-up-to-iso graphs,
  and canonical SET EQUALITY with the independently generated
  triangleramsey+filter n=15 corpus (ext15.log).
- extensions(n=15 corpus) = 18,917 distinct — exactly the documented
  n=16 pool count (ext16.log).
- canonical-form self-test: invariance under random relabelings over both
  corpora (canontest; mini-nauty with exact signature refinement, trace
  pruning, and verified-automorphism orbit pruning; worst case 1,561
  search nodes).

**The run** (all artifacts committed):
- residue644.txt: the codex v6 combined residue (19 K7 + 625 K6-only),
  extracted from d6_current_residue_manifest_v6.json with every adjacency
  verified byte-identical to aeq_d6_n19.txt at its claimed index
  (extract_residue_v6.py).
- 19 -> 20: extensions of the 644 = 223,847 with multiplicity, **117,290
  distinct** level-20 candidates. Greedy re-completion of every one of
  their 20 single-vertex deletions against the residue kills **117,256**
  (a deletion re-completes to a certified level-19 graph). 34 survive.
- Escalation (escalate_frontier.py): for each of the 34, enumerate ALL
  maximal TF completions of every deletion (branching: first addable pair
  in / blocked-by-cherry; deduped) — every completion of every deletion
  of all 34 lands inside the 644 residue. **The level-20 true frontier is
  exactly these 34 graphs** (frontier20_true.txt). Their deletions have
  essentially unique completions (20-24 per graph).
- 20 -> 21: extensions of the 34 = **3,602 distinct** level-21
  candidates; greedy re-completion against the 34-set kills **all
  3,602**. Kill soundness at this level: a re-completed 20-vertex
  candidate outside the true frontier violates a NECESSARY condition for
  realizability (all its own deletion-completions inside the certified
  ledger's complement — and extensions(residue) provably contains every
  candidate satisfying that condition, by the same extension theorem).
  **The level-21 frontier is empty.**

**Consequences, with the trust base stated exactly:**
1. **f(6) <= 20**: no 21-point almost-equidistant set exists in R^6 —
   REST ON: (a) completeness of the 3,971,787-graph level-19 corpus
   (verified against BPSSV and enumaeq6 cross-checks), (b) the certified
   status of every level-19 candidate outside the 644 residue = our
   3,055,480 interval certificates PLUS the codex exact-filter kills for
   the deferred rest (their independent verifiers PASS; OUR independent
   reimplementation cross-check is still owed and remains queued), (c)
   the extension theorem + ext6.c implementation (two corpus controls),
   (d) hereditarity (any m>21-point set contains a 21-point set).
2. **f(6) <= 19 reduces to certifying 34 explicit 20-vertex graphs**
   non-realizable (frontier20_true.txt) — the interval engine's next
   targets. If all 34 die: f(6) <= 19 on the same trust base.
3. The published bound was 26. This session's provisional improvement is
   6 (or 7 pending the 34).

Next actions (in order): LM numerical screen of the 34 (any residual ~ 0
would instead signal f(6) >= 20); independent Python reimplementation of
the extension+frontier pipeline as a cross-check; ckernel6 certification
of the 34; then the owed codex-filter cross-check to make the whole chain
independently verified. DIMENSIONS_1_TO_10.md gets a clearly-marked
provisional note only — the main table keeps 26 until the trust base is
fully in-house.

## 2026-08-07 16:52 append — cadence: 34-graph certification campaign launched; independent checks green so far

- **LM screen of the 34 frontier graphs (500 restarts each): every best
  residual is far from zero** — min 0.2199, typical 0.4-1.5, vs ~1e-30
  for realizable controls (lm_frontier34.err). No 20-point set is hiding
  in the frontier; f(6) <= 19 via certification of the 34 remains the
  expected outcome.
- **Independent re-implementation of the whole frontier pipeline**
  (verify_frontier_indep.py: own subgraph tests with different loop
  structures, own iso decision via profile buckets + backtracking — no
  canonical labeling shared with ext6.c, own completion enumeration,
  greedy completion deliberately last-pair-first where ext6 is
  first-pair-first): **control14 PASS and control15 PASS with exact
  multiplicity agreement** (250,339 and 962,217 raw extensions — equal
  to ext6.c's counts). The 19->20->21 re-run is in progress (currently
  in its single-threaded dedup phase); its verdict on 117,290 / 34 /
  3,602 / 0 lands in a later block.
- **The 34 are hard for the interval engine** (measured): graph 0 dec 0
  slice 0 ABORTs at 3M nodes / 310 s; 8-fold subdivision does NOT
  collapse the cost (all 8 sub-slices ABORT at 150k) — width-independent
  kill-trees, the same disease as the deferred level-19 class, most
  likely near-realizable 18-point cores inside the graphs. Probe of all
  46 decs of graph 0: best decs have 2/6 probe slices hot at 30k cap; no
  globally cheap parametrization.
- **Campaign launched** (certify34_campaign.py, 24 workers, resumable
  state in frontier34_state.json): stage 1 probes every (graph, dec) at
  6 slices @ 30k; stage 2 full-96-tiles the best 2 decs @ 150k; stage 3
  runs an 8-fold subdivision ladder @ 500k then 2M on the best dec.
  Certifications are logged as they land; a monitor reports milestones.
- Next work item while it grinds: the owed independent cross-check of
  the codex exact filters (starting with the workhorse disjoint-edge
  bounded-cover rule), which is what upgrades f(6) <= 20 from
  "modulo their ledger" toward fully in-house.
