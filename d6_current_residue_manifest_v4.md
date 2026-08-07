# Exact dimension-six residue v4

This boundary incorporates the complete K6 empty-defect support-budget run
into the previously verified v3 union.  It is a rejection boundary, not a
realizability classification.

## Exact accounting

```text
v3 K7-containing residue                              155
v3 K6-only residue                                    822
new K6 empty-support-budget rejections                 17
v4 K6-only residue                                    805
v4 combined residue                                   960
```

The 17 new rejected corpus indices are

```text
202556 2175302 2275063 2593296 2672875 2673439
3015139 3015146 3334740 3340669 3555775 3655287
3724295 3968267 3968279 3968280 3968436
```

For each rejected graph the production run archived the K6 seed and every
zero-factor/component row needed by the decision.  The independent checker
reconstructed all 822 inputs, replayed its own 64-state Hall dynamic program,
validated 1,481 archived `Z0` rows and 1,641 component witnesses, and retained
all 32 K6 seeds of the known realizable 18-point control.

The failed sandbox launch in
`d6_k6_empty_support_budget_abort.json` processed zero graphs and has status
`INFRA_ABORT`; it contributes no mathematical conclusion.  The production
run was then restarted outside that macOS semaphore restriction from the
atomic checkpoint.

## Reproduction boundary

Production source commit:
`b54b69571b5ff87c7586aff13e5c69a24ece7801`.

```text
d6_k6_empty_support_budget_report.json
  1860dbe69ae74b55203ea283cf83afff929a4b1f5cbfcd30e09f0138688eba9f
d6_k6_empty_support_budget_certificates.json
  c6fcb7ef66bccc3319fe0a979c5fd63d9f6fd9535261c5c9bfb87e8d210361d4
d6_k6_empty_support_budget_checkpoint.json
  da13eecacc47559f617924a593692f42ed6d05b49d4e5deeed7b0987ebbd0bd0
d6_k6_empty_support_budget_verification.json
  873c8b6babe183757d4e97b0f0c1c1a942b31f180c88eeef1c9ddf80a199276f
d6_current_residue_manifest_v4.json
  6ab4bfffc524f5b43409d59888fb596de8130315bda3381e484bb4ce891e03e4
d6_current_residue_manifest_v4_verification.json
  765eceb6782d131dae8c77c9b94de78735e23a070da051c8fffbd984790e1a41
```

The committed verification artifact is immutable runtime provenance.  A
fresh replay should use a separate output path rather than overwrite it:

```text
python3 verify_d6_k6_empty_support_budget.py \
  --expected-report-sha256 \
  1860dbe69ae74b55203ea283cf83afff929a4b1f5cbfcd30e09f0138688eba9f \
  --workers 11 --output /tmp/d6_k6_empty_support_replay.json

python3 verify_d6_current_residue_manifest_v4.py \
  --expected-manifest-sha256 \
  6ab4bfffc524f5b43409d59888fb596de8130315bda3381e484bb4ce891e03e4 \
  --output /tmp/d6_v4_replay.json
```

Candidate edges remain required unit distances.  Candidate nonedges remain
unconstrained and may also be unit.  Every v4 survivor is unresolved.
