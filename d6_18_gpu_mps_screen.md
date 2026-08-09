# Dimension-6 18-deletion MPS realization screen

## Status and scope

`run_d6_18_gpu_mps_screen.py` is a **heuristic discovery screen**, not a
non-realizability checker.  It searches all 12,712 isomorphism classes obtained
by deleting one vertex from the surviving level-19 graphs, with particular
interest in the 12,698 classes that do not embed in the known 18-point
halfcube-plus-two-poles configuration.

Failure in every restart proves nothing.  A low residual only nominates a
graph and coordinates for exact or interval reconstruction.  Candidate
nonedges are never constrained; in particular, a nonedge is allowed to end up
at unit distance.  The report and every checkpoint use the status
`HEURISTIC_ONLY_NO_REJECTIONS` and record zero mathematical rejections and zero
realizability conclusions.

## Gauge and objective

Every deletion graph contains a unit `K6`.  The corrected producer chooses the
lexicographically first `K7` when one exists, otherwise the lexicographically
first `K6`.  There are 1,616 `K7`-gauge classes and 11,096 `K6`-only classes.
The 14 standard-compatible controls are all in the `K6`-only class, as must be
the case because the standard-18 unit graph has clique number six.

One fixed Helmert realization of a regular `K7` in `R^6` supplies a common
coordinate frame.  A `K6` gauge uses its first six vertices; a `K7` gauge also
fixes the seventh.  Thus a `K6` run optimizes 72 coordinates and a `K7` run
optimizes 66, while padded tensors allow both types in one MPS batch.

```text
q_1,...,q_7 = rows of a Helmert basis for 1^perp in R^7, divided by sqrt(2).
```

This removes translation and rotation freedom.  For each required unit edge, the loss
contains `(squared_distance - 1)^2`.  A small smooth short-distance repulsion
is added as a search bias, and the output separately records:

- the best edge RMS even if vertices collide;
- the best edge RMS among endpoints whose minimum pair distance is at least
  `0.01`.

Thus a collapsed solution is never reported as a distinct candidate.  The
repulsion is not asserted as a geometric condition and may cause genuine
solutions to be missed.

## Exact seed-sphere starts, LM refinement, and controls

Random directions use NumPy PCG64.  The stream for a graph and restart wave is
determined independently of batch order by

```text
(seedbase*1000003 + corpus_index*7919 + wave*104729) mod 2^64.
```

For each of the 14 standard-compatible classes, start zero in wave zero is
replaced by an exact standard-18 realization aligned to the fixed `K6`.
These are positive controls for input semantics, gauge alignment, MPS tensor
layout, and residual calculation.  The random-start total in the final report
subtracts these 14 seeded starts explicitly.

The corrected initializer does not waste iterations repairing edges to the
fixed simplex.  If an outside vertex has required seed-neighbour set `A` of
size `k`, it starts at

```text
m_A + sqrt((k+1)/(2k)) v,
m_A = (1/k) sum_{i in A} q_i,
```

where `v` is a random unit vector orthogonal to the affine span of the seed
face.  Every required seed edge is therefore unit at step zero.  This is the
exact common-sphere identity; it is not a numerical approximation.

For every graph, the producer retains the best distinct MPS endpoint even if
it is far from zero.  A configurable top `k` of the ranked distinct endpoints
is refined by SciPy float64 Levenberg--Marquardt with an explicit analytic
Jacobian.  The LM objective includes a fully disclosed short-distance hinge,
zero outside its collision radius, because edge-only LM frequently finds
exact but coincident graph homomorphisms.  This hinge is a discovery bias, not
a permissible proof constraint.  The raw checkpoint stores every LM attempt's
input rank, residual, separation, status, and evaluation counts.

## First sweep: stopped weak calibration

The first launched configuration was:

```text
caffeinate -dimsu \
  /Users/lukacs/claude/opengauss/venv/bin/python3 \
  run_d6_18_gpu_mps_screen.py \
  --full --device mps --chunk-size 64 \
  --total-restarts 128 --restarts-per-wave 32 \
  --steps 1000 --cpu-refine-max-nfev 2000
```

Its configuration hash is
`6f2d86c247bb39e5badda1aa7e3df71e8bab4c415958f247d21ccc22a1521173`.
Atomic gzip checkpoints live in
`d6_18_gpu_mps_checkpoints/6f2d86c247bb39e5/`.  A restart with the exact
same command verifies and skips complete chunks.  Temporary files are fsynced
and renamed into place; partial temporary files are not treated as chunks.

This run was intentionally stopped after 33 complete chunks, covering graphs
`[0,2112)`.  An unseeded calibration on the 14 known standard-compatible
classes showed that 32 random starts and 1000 Adam steps recovered none of
them: the best edge RMS values were only `0.1046--0.1352`.  Because the first
producer retained coordinates only below RMS `3e-4`, completing that run would
have yielded little information about promising LM basins.  The preserved
partial output is therefore labeled `HEURISTIC_WEAK_CALIBRATION`, not a
completed screen.  `termination_state.json` records the exact boundary and
reason without modifying any checkpoint.

The producer source hash frozen into this run is
`ee7a1f30abf345f06008a0682262c7204872e552f6031b543370ca8629bc9f6d`.
The pinned corpus hash is
`9d08f9ec579434c9601ad3908bd2ae51f10e09e9d06a7f38fbc25794a634b732`.

## Benchmark

On the Apple MPS device, a realistic 64-graph, 32-start, 100-step benchmark
took 1.575 seconds end to end; its optimizer wave took 0.417 seconds and
processed about 491,000 graph-start-steps per second.  All 14 seeded controls
passed before optimization.  Full 64-graph chunks settled at roughly
12.5--13.2 seconds each, projecting about 42 minutes for 199 chunks before
final assembly and independent verification.

That throughput benchmark describes only the stopped first producer.  The
corrected producer deliberately adds CPU LM work, which dominates wall time
when only two workers are available.

## Corrected-producer calibration

The controls were always seeded in ordinary benchmarks to validate tensor
semantics, but the following calibrations disabled that exact seed to measure
genuine basin recovery:

- Gaussian MPS starts followed by LM recovered only 1 of 14 controls in 32
  starts.
- Exact seed-sphere starts followed by LM recovered 3 of 14 in 32 starts.
- With 128 seed-sphere starts, collision-aware LM recovered 9 of 14.  The
  recovery-by-ranked-endpoint curve was `4,5,7,8,9` at top
  `1,2,4,8,16`; refining ranks 17 through 128 found no additional control.
  The four MPS waves took 3.28 seconds, whereas 1,792 LM attempts took 70.53
  seconds on two CPU workers.
- A complementary all-`K6`-gauge pilot used 8--12 gauges per control, 16
  starts per gauge, and top-8 LM.  It recovered 7 of 14, including one control
  missed by the single-gauge run, but was not uniformly better.

The frozen machine-readable artifacts are
`d6_18_gpu_mps_final_unseeded_control_calibration.json` (the 128-start/top-16
result, 9 of 14 controls in 15.57 seconds) and
`d6_18_gpu_mps_all_k6_gauge_calibration.json`.  Both match the final producer
hash exactly.  Earlier local Gaussian and 32/128-rank exploratory JSON files
are chronological diagnostics, not part of the frozen package.  None of these
heuristic calibrations is evidence that an unrecovered graph is unrealizable.

No corrected full run should start until the finalized producer and checker
are committed.  The control curve supports top-16 LM as the measured
cost/coverage knee, but its 9/14 recovery rate is an explicit sensitivity
warning, not a completeness claim.

The final 64-graph audit-boundary benchmark used 128 starts, top-16 LM, two
CPU workers, and included both gauge types (58 `K6`, 6 `K7`).  It took 56.08
seconds: four MPS waves took 13.27 seconds total and 1,024 LM attempts took
40.39 seconds.  This projects 11,140 seconds (3.09 hours) over the corpus with
only two LM workers; near-linear CPU scaling after the separate interval job
finishes would reduce the LM component substantially.  None of the 50 sampled
nonstandard classes became an LM candidate.  Its artifact is
`d6_18_gpu_mps_corrected_benchmark.json`, configuration hash
`3249d4d3a4afd60ddeb921de9a4bf1402e04cc7ee2a673d10dc28f22b44bd311`,
and producer hash
`cd57a6193fba1a9b05b88be566598413bd7ff79a13cdc7b5530df90a14293627`.

## Independent checks

`test_d6_18_gpu_mps_screen.py` checks all 12,712 deterministic mixed `K7`/`K6`
gauges, the 14 exact aligned standard controls, the exact common-sphere
initializer, batch-independent random streams, atomic gzip payload hashing,
and a seeded CPU tensor smoke test.

After assembly, `verify_d6_18_gpu_mps_screen.py` independently checks:

- corpus, producer, configuration, result, and checkpoint hashes;
- complete ordered row coverage and graph metadata;
- each checkpoint's internal payload hash and equality with the compact TSV;
- every deterministic gauge seed and every retained endpoint's recomputed
  residual, maximum error, and separation;
- all LM-attempt aggregation, input ranks, collision-aware candidate flags,
  and termination counts;
- all threshold and candidate summary arithmetic;
- the exact standard-18 control embeddings; and
- edge RMS, maximum error, and minimum pair distance from every retained
  candidate coordinate witness.

The verifier also emits zero mathematical rejections and zero realizability
conclusions.  Any promising nonstandard witness must be reconstructed and
certified by a separate exact or rigorous interval method.

The complete first replay exposed one checker-only precision issue: retained
MPS endpoints store binary32 coordinates, while their archived metrics were
reduced by Torch/MPS in binary32 and the checker recomputes them in Python
binary64.  The largest observed minimum-distance discrepancy was below
`8e-7`.  The corrected checker uses `rel_tol=1e-5, abs_tol=1e-6` only for
those MPS endpoint metrics; CPU-LM witnesses and summaries retain the tighter
binary64 tolerance.  This is an audit tolerance for a heuristic artifact and
does not enter a geometric rejection.
