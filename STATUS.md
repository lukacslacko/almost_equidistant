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
