# Certified cap-100,000 interval pass on the K7 residue

## Scope and exact result

The frozen residue-v5 interval launcher was run on all 155 K7-containing
graphs with 11 workers, four candidate orders, 24 certified circle slices,
and a 100,000-node cap per kernel call.  The launch report records commit
`6cf6c53a83127ba68cf188f99f508a1930ec524b`, with every tracked source clean;
the working tree contained unrelated untracked research artifacts.

The complete status partition is

```text
KILLED          24
ABORT          131
UNRESOLVED       0
INFRA_ERROR      0
```

Only `KILLED` is a non-realizability conclusion.  Every `ABORT` remains
unresolved.  The certified indices are

```text
379078   423661   424226   2280137  3330235  3331585
3342165  3368321  3394720  3643248  3723506  3726282
3727225  3793488  3910558  3910562  3910564  3936176
3936987  3939486  3939503  3945374  3950483  3960292
```

This contains all 22 cap-20,000 K7 kills and adds exactly `423661` and
`424226`.  The complete run used 1,523.92 wall seconds, with 12,357.60 summed
graph-seconds and an 883.56-second maximum graph.  All 155 per-graph results
were written atomically; there were no infrastructure failures.

The production command was

```text
caffeinate -dimsu python3 run_d6_interval_v5.py \
  --workers 11 --population K7 --orders 4 --slices 24 --cap 100000 \
  --include-circle-orders \
  --checkpoint-dir .runs/d6_interval_v5_k7_cap100000 \
  --output .runs/d6_interval_v5_k7_cap100000.json \
  --decisions .runs/d6_interval_v5_k7_cap100000.tsv \
  --progress-every 20 \
  --outer-launch-command \
  'caffeinate -dimsu python3 run_d6_interval_v5.py --workers 11 \
   --population K7 --orders 4 --slices 24 --cap 100000 \
   --include-circle-orders \
   --checkpoint-dir .runs/d6_interval_v5_k7_cap100000 \
   --output .runs/d6_interval_v5_k7_cap100000.json \
   --decisions .runs/d6_interval_v5_k7_cap100000.tsv --progress-every 20'
```

## Independent replay

The independent verifier reconstructs the exact 155-graph selection from
the verified v5 manifest, checks all checkpoint/report/TSV counters and
hashes, reconstructs every winning seed and placement order, and reruns the
frozen interval kernel on every winning slice.  It imports neither the v5
launcher nor its graph-selection implementation.

The replay returned `PASS`:

```text
winning witnesses                 24
winning slices replayed          576
winning nodes replayed        323608
positive/negative controls        48 / 48
```

The verification command was

```text
python3 verify_d6_interval_v5.py \
  --report .runs/d6_interval_v5_k7_cap100000.json \
  --expected-report-sha256 \
    85b175f3b3348204c3ba9065b30432458528bdd2e22b471e64fc56341480bac3 \
  --expected-cap 100000 \
  --output /tmp/d6_interval_v5_k7_cap100000_verification.json
```

To reproduce from Git, extract
`d6_interval_v5_k7_cap100000_checkpoints.tar.gz` at the repository root and
run the same command against the restored `.runs` report.

## Relation to the exact K7 pinning layers

The exact double-pin/full-pin algebra rejects 131 of the 155 K7 graphs.
Twenty-one interval kills overlap that exact set.  The three interval-only
indices are

```text
423661 424226 3936176
```

Thus the union of the independently checked exact-pinning and interval sets
rejects 134 graphs and leaves 21 K7 graphs.  This accounting must be bound by
a separate combined-residue manifest before it becomes the input to another
theorem-level campaign.

## Trust assumptions and hashes

The interval credit assumes correctly rounded IEEE-754 binary64 basic
operations and square root.  Every endpoint is expanded with `nextafter`.
For macOS `libm` cosine endpoints, the kernel assumes an error below eight
ulps and pads both directions by eight `nextafter` steps while including all
interior extrema.  Candidate nonedges are unconstrained and may have unit
distance.

```text
85b175f3b3348204c3ba9065b30432458528bdd2e22b471e64fc56341480bac3
  d6_interval_v5_k7_cap100000_report.json
3eddc672b0facc6cf1ee170394c1de798592597430fbab0281b0aafb0341c3fb
  d6_interval_v5_k7_cap100000_decisions.tsv
a56fb1561b778c55bffb433bf0eea90e5cca823573d0717e74f2a73528bc1b10
  d6_interval_v5_k7_cap100000_verification.json
07632a7ad9bd903270612d68af6f156c01283ff5d2fcf40ee8e47d199127d023
  d6_interval_v5_k7_cap100000_checkpoints.tar.gz
```

This pass does not claim that any non-killed graph is realizable, and it does
not by itself settle the K7 class or dimension six.
