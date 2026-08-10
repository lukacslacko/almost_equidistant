# Dimension-6 level-19 intrinsic-manifold GPU screen

`run_d6_19_gpu_mps_residue_v1.py` is a restartable numerical discovery
screen. It is **not** a rejection certificate and it cannot certify a
realization. A failure in all gauges has no mathematical meaning; every small
residual still requires exact reconstruction or rigorous interval checking.

## Frozen boundary and gauges

The runner pins four committed inputs by SHA-256: the verified v8 residue and
verification, followed by the exact graph-3936435 report and its independent
verification. All 261 v8 records remain physically present. Graph `3936435`
is tagged as an exact-negative numerical calibration, so it is never counted
among the 260 active unresolved graphs. Prospective source-only reductions are
deliberately not removed.

The full deterministic census is:

```text
active K7 graphs / gauges                 11 / 34
exact-negative K7 calibration graphs       1 / 3
active K6-only graphs / gauges           249 / 7,684
physical total graphs / gauges           261 / 7,721
```

Gauge-manifest SHA-256:
`a545686fc30b230473bcfcc091a826bd3d781fe0079163d4f1a8a31d04c10ca3`.
The all-K7 policy is used for the K7 class and the all-K6 policy for K6-only.

## Intrinsic MPS parameterization

For an outside vertex adjacent to `k >= 1` fixed simplex vertices, Adam does
not optimize a free Cartesian point. It optimizes a padded coefficient vector
on the common-unit-sphere manifold

```text
x = center + radius * B * (z / ||z||),
radius^2 = (k + 1) / (2k).
```

Thus every required edge to the K6/K7 gauge remains unit by construction.
After each Adam step, every active coefficient block is projected back to unit
norm and padded coordinates are zeroed. Vertices with `k=0` remain Cartesian.
The active intrinsic degrees of freedom range from 25–31 for K7 gauges and
27–36 for K6-only gauges, rather than 72 and 78 free Cartesian coordinates.

For a K7 gauge and `k=6`, the zero-sphere consists of the missing simplex
vertex and its reflection through the K6 face. Distinctness excludes the
missing-vertex collision, so the unique reflected root is fixed. Six physical
K7 gauges contain such a root. Duplicate forced roots fail closed, as does a
K8. The subsequent SciPy float64 LM stage intentionally returns to free
Cartesian coordinates as an independent numerical cross-check.

Only adjacency-one entries enter the unit-edge residual. Candidate nonedges
are unconstrained and may also have distance one. The all-pair term is solely
a soft short-distance collision bias.

## Controls and resources

Every configuration first runs three seeded controls through the same final
kernel:

- the standard 18-point configuration plus a distinct lightly tethered point,
  in its first and last K6 gauges;
- the standard 18-point configuration plus a distinct reflected apex over a
  standard K6, giving a K7 gauge.

The exact-negative graph `3936435` is part of the physical sweep and is tagged
separately in the final report. A numerical candidate there is an alert, not a
reversal of its exact certificate.

Defaults use MPS with 16 simultaneous restarts, 1,000 Adam steps, and the best
two endpoints sent to 1,500-evaluation float64 LM. A production configuration
is blocked unless all three controls both enter with the planted low residual
and recover as numerical candidates after Adam plus LM. Runs with LM disabled
are labeled `PASS_SEEDED_INPUT_ONLY`, not full final-kernel control passes.
Ten LM processes are spawned only after all BLAS/OpenMP thread-limit variables
have been set to one; PyTorch itself is limited to one CPU thread while MPS is
active.

The default checkpoint chunk is eight physical graphs. Across this frozen
corpus that means 22--291 gauges per chunk, avoiding the 116--560 range of the
older 16-graph draft while retaining thousands of simultaneous MPS starts.

## Commands and restart behavior

Accelerator-free boundary tests and a tiny final-kernel control smoke test:

```bash
python3 -m unittest -v test_run_d6_19_gpu_mps_residue_v1.py
python3 run_d6_19_gpu_mps_residue_v1.py --selection-only
python3 run_d6_19_gpu_mps_residue_v1.py --control-only --device cpu \
  --total-restarts 1 --restarts-per-wave 1 --steps 2 \
  --lm-top-k 0 --cpu-refine-max-nfev 0 --cpu-workers 1 \
  --checkpoint-root /tmp/d6_n19_gpu_control_smoke
```

Full production command (not launched while preparing this source boundary):

```bash
caffeinate -dimsu python3 run_d6_19_gpu_mps_residue_v1.py
```

Atomic checkpoints live below
`.runs/d6_19_gpu_mps_residue_v1/<config-hash-prefix>/`. Resume with the exact
same command. The source hash is part of the configuration hash, so edited
code cannot silently consume old checkpoints. Completed report or detail
paths also fail closed rather than overwrite an artifact carrying another
configuration hash. `--max-new-chunks 1` provides a checkpoint-safe pilot;
`--index` may be repeated for a pinned subset.
