# Status — f(5) and f(6) campaigns (2026-08-05, ~21:00 push)

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

## 11:55 append — d=6 cadence: STUCK flag raised, diagnosing

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

## 13:00 append — d=6 reality check (and a lost hour, explained honestly)

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

## 13:25 append — d=6 pivot to obstruction library, first pilot numbers

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
