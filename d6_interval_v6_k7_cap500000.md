# Pre-production cap-500,000 interval campaign on the v6 K7 residue

## Scope

This package freezes a theorem-level interval campaign on exactly the 19 K7
graphs in the independently verified v6 current-residue manifest.  It is a
launch plan, not a result: no graph receives interval credit until the frozen
runner has completed and the independent checker has replayed every winning
slice.

The pinned inputs are:

```text
c513bcc40c037af23340dce3cc91d593a101dd1d073f93feab82abd556efadc9
  d6_current_residue_manifest_v6.json
3d11cbd20cb59b42a8c8fd5252fb3ef3dea649576982ea8017622fad3279a7bc
  d6_current_residue_manifest_v6_verification.json
af25b842ca9eaa6e9721f51d0400be2436d930b72cc33c72bc366f0ac5e41290
  ordered list of the 19 selected K7 indices
```

The search uses nine native worker threads, four non-bulk candidate orders,
24 certified angular slices, and a 500,000-node cap for each kernel call.
Every graph is checkpointed atomically.  Only `KILLED` is a mathematical
conclusion; `ABORT`, `UNRESOLVED`, and `INFRA_ERROR` remain unresolved.

## Frozen launch and resume commands

Run only after the source package has been committed on
`codex/dimension6` and `git status --porcelain=v1 --untracked-files=all`
shows no modified tracked files.  Unrelated untracked research files are
recorded in the launch provenance and are permitted.

```text
caffeinate -dimsu python3 run_d6_interval_v6.py \
  --workers 9 --orders 4 --slices 24 --cap 500000 \
  --checkpoint-dir .runs/d6_interval_v6_k7_cap500000 \
  --output .runs/d6_interval_v6_k7_cap500000.json \
  --decisions .runs/d6_interval_v6_k7_cap500000.tsv \
  --checkpoint-every 1 --progress-every 1 \
  --outer-launch-command \
  'caffeinate -dimsu python3 run_d6_interval_v6.py --workers 9 --orders 4 --slices 24 --cap 500000 --checkpoint-dir .runs/d6_interval_v6_k7_cap500000 --output .runs/d6_interval_v6_k7_cap500000.json --decisions .runs/d6_interval_v6_k7_cap500000.tsv --checkpoint-every 1 --progress-every 1'
```

The same command resumes a clean interruption from the atomic checkpoints.
If and only if the completed report records an `INFRA_ERROR`, rerun the same
command with `--retry-infra-errors`; that flag changes checkpoint scheduling,
not the hashed mathematical configuration.  Never edit a checkpoint or log
in place.

The exploratory observation that index `3595554` may be killed at this cap is
not a certificate and is deliberately excluded from the v6 manifest and this
pre-production claim.  It must be rediscovered by the committed campaign and
replayed by the checker.

## Independent verification

After completion, first freeze the report hash, then run:

```text
python3 verify_d6_interval_v6.py \
  --report .runs/d6_interval_v6_k7_cap500000.json \
  --expected-report-sha256 <SHA256_OF_COMPLETED_REPORT> \
  --expected-cap 500000 --expected-workers 9 \
  --output /tmp/d6_interval_v6_k7_cap500000_verification.json
```

The checker imports neither `run_d6_interval_v6.py` nor its selection logic.
It reconstructs the 19 graphs directly from the pinned v6 manifest, checks
the committed source blobs and clean tracked launch boundary, checks the
report, decision TSV, campaign state, sessions, and all 19 atomic result
files, reruns the positive and negative controls, and exactly replays every
slice in every claimed `KILLED` witness.  A report is publishable only after
this command returns `PASS` with replay enabled.

## Trust boundary

The interval conclusions assume correctly rounded IEEE-754 binary64 basic
operations and square root, with every interval endpoint expanded using
`nextafter`.  The macOS `libm` cosine assumption is an error of at most eight
ulps; the kernel pads endpoints by eight `nextafter` steps in both directions
and includes all interior extrema.  Candidate nonedges are unconstrained and
may also be at unit distance.

The wrapper and checker additionally bind the exact 264 exact-algebra/graph
rejections and three incremental interval rejections used to construct v6.
No non-`KILLED` interval status is used as a rejection at either layer.

## Frozen source hashes

The final pre-production source hashes are recorded below after the wrapper
hash is embedded in the independent checker:

```text
509ec8055542bb60863756d1ce64a24cccfe57ccc9f56f13925edf3290ac38be
  run_d6_interval_v6.py
bbfe851282934980a1f0343f80296b8c8884d4ca305b6babcff09ede214501d4
  verify_d6_interval_v6.py
58fbbbc548ec3cc6b9d5b231a6790249b02a75d7bfd82f1afcd7fc77fce59404
  run_d6_interval_residue.py
8f65948c8bc95242d0fde50a8cbcb798cb0e225b78a8ce626382d00135b95173
  cdriver6.py
383ef7a17c328869c338f64a16f2b90875ff968fe2fc20e9a5c3094c830e0874
  ckernel6.c
524e41e0d0637a5352c59ec998009b968f2c9f5e7f31b7f7927e8bc36c2e6fa0
  ival.py
84ec7ae7e682de660329fd1fb02bd88d774bde542032ece95a7eb7d3c54ece4a
  verify_d6_interval_benchmarks.py
```
