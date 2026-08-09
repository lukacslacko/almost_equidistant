# All-K6-gauge GPU triage of the 181-class n=18 cover

`run_d6_18_gpu_mps_cover_v7.py` is a focused **heuristic discovery tool** for
the production cover selected by `run_d6_interval_18_cover_v7.py`.  It is not
a certificate generator and cannot reject a graph.  Exhausting every gauge
and start proves nothing; a small residual only identifies a candidate for an
exact reconstruction or a rigorous interval check.

## Frozen input and gauge coverage

The wrapper rebuilds the pinned v7 residue and its dense deletion cover rather
than accepting an ad hoc graph list.  The resulting campaign has 181 canonical
18-vertex classes: the 179-class greedy cover, the known standard-realizable
positive control, and two non-positive actionable backups.  It enumerates
every required-edge K6 in every class:

- 4,099 K6 gauges in total;
- 500 gauges among the 12 K7-population classes;
- 3,599 gauges among the 169 K6-population classes;
- 9 gauges for the standard-realizable class at deletion-manifest ordinal 441.

The complete gauge manifest and variant-id list are source-pinned.  A variant
id is `128 * class_index + gauge_ordinal`, where gauges are lexicographically
enumerated.  This makes NumPy random streams independent of MPS batch and
checkpoint boundaries.

## Semantics

Only entries present in the graph adjacency matrix enter the unit-edge
residual.  A candidate nonedge is unconstrained and may also end up at unit
distance.  The optimizer has an all-pair short-distance repulsion, but this is
only a soft discovery bias against coincident numerical roots; it does not
ask a nonedge to be non-unit.  The detailed output and summary report both
record zero mathematical rejections and zero realizability conclusions.

Every numerical witness is restored to original vertex order and retained
with its class ordinal, class id, K6 seed, gauge ordinal, residuals, and
distinctness metric.  MPS endpoints and LM candidates remain floating-point
triage data.

## Current resource-balanced command

First audit the exact selection without touching the GPU:

```sh
python3 run_d6_18_gpu_mps_cover_v7.py --selection-only
```

While eleven CPU threads are occupied by the interval run, use:

```sh
caffeinate -dimsu python3 run_d6_18_gpu_mps_cover_v7.py \
  --device mps \
  --total-restarts 8 \
  --restarts-per-wave 8 \
  --steps 800 \
  --lm-top-k 2 \
  --cpu-refine-max-nfev 1000 \
  --cpu-workers 1 \
  --class-chunk-size 16
```

This is 32,792 MPS optimizer starts and at most 8,198 sequential LM attempts.
An out-of-sandbox MPS pilot on the largest class (56 gauges, 448 starts and 112
LM attempts at exactly these settings) took 15.3 seconds and preserved all 56
distinct endpoints.  Linear scaling predicts about 19 minutes; allow roughly
18--30 minutes for the full run because MPS batch scaling, graph difficulty,
and the sequential LM tail vary by chunk.  A one-chunk pilot can be stopped
safely with `--max-new-chunks 1`; the same command without that option resumes
from the matching configuration directory.

The default seeds every gauge of the known standard realization with aligned
standard coordinates and requires the completed report to audit that control.
Use `--no-seed-standard-controls` only for a separately labelled unseeded
recovery calibration.

Atomic checkpoints live below
`.runs/d6_18_gpu_mps_cover_v7/<config-prefix>/`.  A complete run writes a
compact JSON report plus a gzip-compressed detail payload containing all
per-gauge rows, retained endpoint coordinates, LM attempts, and candidate
witness coordinates.

## Fast controls

```sh
python3 -m unittest -v test_run_d6_18_gpu_mps_cover_v7.py
```

These controls verify the frozen 181/4,099 counts and hashes, every K6 seed,
all nine exact positive-control alignments, chunk-independent random streams,
checkpoint integrity, one-worker defaults, and the optional-nonedge objective
semantics without requiring MPS.
