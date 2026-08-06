# Full K7 positive-polynomial dual campaign

This package productionizes the frozen exact certificate locator documented in
`d6_k7_positive_polynomial_dual.md`; it does not change that mathematical
kernel or the frozen 84-graph sample package.

## Exact scope and claim

`build_d6_k7_positive_dual_selection.py` independently reconstructs the
current K7 residue from the 17,764 rank survivors.  Its hash-pinned manifest,
`d6_k7_positive_dual_selection.json`, removes the union of all existing exact
layers in this order:

| existing exact layer | reported rejections | new in this union | residue |
| --- | ---: | ---: | ---: |
| support propagation | 1,536 | 1,536 | 16,228 |
| strict-H | 603 | 132 | 16,096 |
| sparse value | 3,195 | 3,148 | 12,948 |
| interval benchmark, cap 20,000 | 10 | 6 | 12,942 |
| frozen sample positive dual | 2 | 1 | **12,941** |

The six interval-only indices are `532351, 1169330, 2235541, 3331589,
3394713, 3961600`; the sample-dual-only index is `3649646`.  The selected-index
hash is
`a07bafe688c45ff132921519320a19687d2917a79a6c9009bda51a19a6b7f290`,
and the complete selection manifest hash is
`814c80651746ce05048f72f0cfa7e49a921554abaf167bd5c8cf4063a34d35c0`.

Only a `REJECTED` row is a new mathematical conclusion.  For such a row, one
K7 seed has no surviving cover: each eligible cover either failed an earlier
exact layer or carries an exact rational positive-polynomial identity.  A
`SURVIVOR` row means only that the degree-one locator did not establish this
obstruction.  `INFRA_ERROR` makes no claim and forces exit status 2.

## Production command and restart

Commit the exact source and start from a stable worktree before running.  The
runner pins 11 worker processes (one below the 12 logical CPUs), one numerical
thread per worker, and user-initiated macOS scheduling QoS.  `caffeinate`
prevents an intentional run from sleeping.

```sh
caffeinate -dimsu python3 run_d6_k7_positive_polynomial_dual_full.py \
  --input .runs/d6_k7_rank_survivors.json \
  --selection-report d6_k7_positive_dual_selection.json \
  --selection-report-sha256 814c80651746ce05048f72f0cfa7e49a921554abaf167bd5c8cf4063a34d35c0 \
  --workers 11 \
  --max-inflight 22 \
  --checkpoint-dir .runs/d6_k7_positive_dual_full \
  --decisions-archive d6_k7_positive_dual_full_decisions.tsv.gz \
  --certificate-archive d6_k7_positive_dual_full_certificates.jsonl.gz \
  --checkpoint-copy d6_k7_positive_dual_full_checkpoint.json \
  --report d6_k7_positive_dual_full_report.json \
  --pid-file .runs/d6_k7_positive_dual_full.pid.json \
  --outer-command 'caffeinate -dimsu python3 run_d6_k7_positive_polynomial_dual_full.py (full explicit argv is recorded separately)'
```

After interruption or reboot, rerun exactly the same command from the same Git
state.  Each graph has an atomic JSON file below
`.runs/d6_k7_positive_dual_full/results/`; the compatible campaign hash is
checked before loading any result.  Complete checkpoints are independently
checked on load, including every rational identity, and are not recomputed.
Use `--retry-infra-errors` only to retry rows explicitly marked
`INFRA_ERROR`.  The directory is protected by a nonblocking process lock, and
a new campaign refuses to overwrite existing final artifacts.

The final decision TSV and certificate JSONL are deterministic gzip streams.
Certificates are retained for every failed cover, including cover failures on
graphs that are not rejected; this preserves enough evidence for independent
checking without embedding graph records or duplicate first certificates.

## Independent verification

The full archive checker does not import or call the numerical LP locator or
the production runner.  It rebuilds the selection, checks report/source/archive
hashes and complete graph coverage, expands every rational polynomial identity
over `Fraction`, and independently repeats the exact seed/cover quantifiers for
every graph claimed rejected:

```sh
python3 -m unittest -v test_d6_k7_positive_polynomial_dual.py \
  test_run_d6_k7_positive_polynomial_dual_full.py
python3 verify_d6_k7_positive_polynomial_dual.py
python3 verify_d6_k7_positive_polynomial_dual_full.py \
  --selection d6_k7_positive_dual_selection.json \
  --report d6_k7_positive_dual_full_report.json
```

The focused production tests include byte-deterministic gzip, explicit
selection-hash rejection, no-work restart, overwrite refusal, and a known
rejection (`3649646`) round-tripped through the compact archive and rebuilt by
the independent graph quantifier.

## Bounded benchmark

Three 64-graph strata at offsets 0, 6,438, and 12,877 ran with 11 workers.
All 192 completed without infrastructure errors; three cover certificates were
emitted and independently checked, while no graph in these bounded strata was
rejected.  The measured wall times were 12.18, 18.68, and 11.89 seconds.
Together they sustained 4.49 graphs/second and project about 48 minutes for the
12,941-graph residue.  The per-stratum throughputs give 40--63 minutes; use
0.7--1.5 hours as a conservative planning range.

Observed checkpoint payload extrapolates to 14.4 MB.  Because the targeted
frozen sample had much higher certificate incidence, reserve 15--150 MB for
logical checkpoint/certificate data, plus filesystem directory-block overhead.
The measurements and report hashes are in
`d6_k7_positive_dual_full_benchmark.json`.  No full campaign was launched in
this milestone.

This layer is CPU-bound in many small HiGHS linear programs followed by exact
rational reconstruction.  A GPU would add transfer/setup complexity without a
suitable audited solver path and is not useful for this campaign.
