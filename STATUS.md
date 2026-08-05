# Status — f(5) and f(6) campaigns (2026-08-05, ~21:00 push)

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
