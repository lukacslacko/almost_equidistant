# Certified interval pass on the exact d=6 residue v5

## Scope and result

The frozen launcher at commit
`2b97bcec90a032070369e37b78a2ca94e7639b26` reconstructed all 911 graphs in
`d6_current_residue_manifest_v5.json` and ran the outward-rounded dimension-6
placement kernel with 11 native threads, four candidate orders, 24 certified
circle slices, and a 20,000-node cap per slice.

The complete status partition is

```text
KILLED          22
ABORT          888
UNRESOLVED       1
INFRA_ERROR      0
```

All 22 certified graphs are in the K7 stratum; the K6 stratum has zero
interval rejections at this cap.  The exact rejected indices are

```text
379078   2280137  3330235  3331585  3342165  3368321
3394720  3643248  3723506  3726282  3727225  3793488
3910558  3910562  3910564  3936176  3936987  3939486
3939503  3945374  3950483  3960292
```

`KILLED` means that every one of the 24 slices for one archived placement
order was eliminated by a rigorous interval enclosure.  `ABORT` and
`UNRESOLVED` reject nothing.  In particular, this pass does not claim that
any of its 889 non-KILLED graphs is realizable.

The production command was

```text
caffeinate -dimsu python3 run_d6_interval_v5.py \
  --workers 11 --orders 4 --slices 24 --cap 20000 \
  --include-circle-orders \
  --checkpoint-dir .runs/d6_interval_v5_full_cap20000 \
  --output .runs/d6_interval_v5_full_cap20000.json \
  --decisions .runs/d6_interval_v5_full_cap20000.tsv \
  --progress-every 50 \
  --outer-launch-command \
  'caffeinate -dimsu python3 run_d6_interval_v5.py --workers 11 \
   --orders 4 --slices 24 --cap 20000 --include-circle-orders \
   --checkpoint-dir .runs/d6_interval_v5_full_cap20000 \
   --output .runs/d6_interval_v5_full_cap20000.json \
   --decisions .runs/d6_interval_v5_full_cap20000.tsv --progress-every 50'
```

The run took 537.43 wall seconds.  Its maximum single-graph time was 184.44
seconds.  All 911 graph results were written atomically, and no infrastructure
failure occurred.

## Independent checking

`verify_d6_interval_v5.py` imports neither the v5 launch wrapper nor its
selection code.  It independently reconstructs the 155 K7 and 756 K6-only
graphs from the pinned manifest, validates the launch commit and source
blobs, checks every report/checkpoint/TSV counter, and reconstructs each
winning order.  It then reruns the frozen interval kernel on every winning
slice.

The independent replay returned `PASS`:

```text
winning slices replayed       528
winning nodes replayed     223854
status/node/cell matches       528 / 528
positive/negative controls      48 / 48
```

The verification command was

```text
python3 verify_d6_interval_v5.py \
  --report .runs/d6_interval_v5_full_cap20000.json \
  --expected-report-sha256 \
    b8f38e765a7212a6bb8a30f138800c170a03992d76f08eb920e67effa0c00c73 \
  --output /tmp/d6_interval_v5_full_verification.json
```

To reproduce from the committed archive, extract
`d6_interval_v5_full_cap20000_checkpoints.tar.gz` at the repository root and
run the same verifier against the restored `.runs/...json` path.

## Trust assumptions and hashes

The theorem credit assumes correctly rounded IEEE-754 binary64 basic
operations and square root.  Every endpoint is expanded with `nextafter`.
For macOS `libm` cosine endpoints, the kernel assumes an error below 8 ulps
and pads both directions by 8 `nextafter` steps while including all interior
extrema.  Candidate nonedges are unconstrained and may have unit distance.

```text
b8f38e765a7212a6bb8a30f138800c170a03992d76f08eb920e67effa0c00c73
  d6_interval_v5_full_cap20000_report.json
9b0f19181ca4f7f50f48d09bd3c9d6e0326c3abe397ae5a3ed67123c89e0be66
  d6_interval_v5_full_cap20000_decisions.tsv
d080140f45942022899a9f59a2849babd028a48c91f7a16924fbfee682848ab8
  d6_interval_v5_full_cap20000_verification.json
a961b7de8b0f2faeb2503d25ee02583d3f3149453be4c7f13f77d02f293d2af7
  d6_interval_v5_full_cap20000_checkpoints.tar.gz
```
