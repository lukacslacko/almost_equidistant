# Partial certified interval increment from deletion class 2100

## Result and scope

The atomic result for canonical 18-vertex deletion class `2100` in the live
v7 dense-deletion campaign is self-contained: its second generated placement
order is `KILLED` on every one of the 24 closed theta slices covering
`[0, 2*pi]`.  Under the same interval trust boundary as the campaign, this
proves that the required-edge graph has no realization by distinct points in
`R6`.

Class 2100 has 106 required unit edges and canonical graph6 string

```text
QTm}BxB{nTzVnwX}i}y^^^^]}~o
```

Its two induced occurrences are obtained from K7 parent `3945564` by deleting
vertex 4 or vertex 15.  Therefore a realization of that 19-point parent would
restrict to a realization of class 2100, which is impossible.  Candidate
nonedges are never constrained and may also have distance one.

This is a valid partial increment before the other 180 deletion classes have
finished.  The campaign file is immutable after initialization and each
result file is written atomically.  Neither an unfinished campaign nor any
`ABORT` status is used as a rejection.

## Exact frozen evidence boundary

```text
raw campaign.json SHA-256
  70ba451094eb4eb5d039fc059ab8e3c351ebeed157e2eae39cc35fe4dcb3e348
raw result_000032_0002100.json SHA-256
  0e246dc1e2a84d5fafb3746366e5f1b06fb3b2368c0ad1b61ac12446d835815a
configuration SHA-256
  82c8d332af4ac4e3984949b98cde969f5b239eadb640e9633f704f8ec88fbcba
campaign executable-source commit
  22be775b1dafde30122ab579f66a7a539581ba73
kernel binary SHA-256
  ec93d78302d6fcf0f3a21ebcdb65333004eb7a406ab4ada8445231536fed07bd
```

The producer writes a deterministic gzip JSON evidence object containing the
exact UTF-8 bytes of both raw files, not merely a transcription of the
winning witness.  The report binds that evidence archive by SHA-256.  It also
pins the five committed campaign/kernel sources, the v7 parent boundary, the
independently verified deletion corpus, all positive-control inputs, and the
current v8 residue.

## Search witness

The checker independently regenerates exactly two placement orders.  Both use
the K7 seed

```text
1 7 8 13 14 15 16
```

The first order certifies six slices and then returns `ABORT` at its node cap
on slice 6.  That order makes no rejection.  The second order is

```text
12 6 17 11 10 9 5 4 3 2 0
```

and has one circle placement, at its first vertex.  All 24 slices return
`KILLED`.  The complete two-order replay has 31 kernel calls: 30 `KILLED`, one
`ABORT`, 108,355 nodes, zero survivors, and zero unresolved cells.  The
checker requires exact agreement of every status, node count, and cell count.

## Independent quantifier checks

`verify_d6_interval_18_class2100_increment.py` imports only the low-level
kernel ABI wrapper.  It does not import the producer, launch wrapper, campaign
manager, or existing interval verifier.  It independently performs all of
the following:

1. Decodes the graph6 string and checks its SHA-derived class ID and adjacency.
2. Rebuilds all 263 v7 parents, all dense deletion choices, the 179-step greedy
   cover, the two actionable backups, and the sorted 181-class production run.
3. Scans the v7 parent deletions with exact graph isomorphism and reconstructs
   precisely the two rooted occurrences above, including attachment masks.
4. Reimplements the clique/closure placement-order generator and replays the
   losing prefix and all 24 winning slices through the pinned binary.
5. Replays all 48 flexible-circle/K8 slice controls and all four independently
   exact known-realizable 18-point controls.
6. Scans every deletion of every graph in the 261-case v8 residue for class
   2100 as a non-induced required-edge subgraph.  There are 4,959 deletions,
   of which 1,085 have at least 106 edges.  The checker enumerates all 17,344
   ways to delete zero through three edges from the 47-edge complement,
   filters by independently computed 1-WL invariants, and finishes with exact
   isomorphism backtracking.  Only the same two deletions of parent 3945564
   occur; no additional current parent is rejected by this obstruction.

The producer's containment scan uses a different exact algorithm: direct
spanning-subgraph domain propagation on complements.  Agreement between the
two methods is required.

## Trust assumptions

The conclusion is interval-dependent rather than exact algebraic.  It assumes
correctly rounded IEEE-754 binary64 basic operations and `sqrt`; every interval
endpoint is expanded with `nextafter`.  On the recorded macOS platform it also
assumes `libm` cosine endpoint error is at most eight ulps; the kernel pads
both directions by eight `nextafter` steps and includes all interior extrema.
The source and Mach-O binary hashes are fixed above.  A different platform or
binary needs its own validated trust boundary.

## Source-bound artifact commands

First commit these four source files.  Then package the immutable result
without touching the live checkpoint directory:

```text
python3 build_d6_interval_18_class2100_increment.py \
  --campaign .runs/d6_interval_18_cover_v7_cap100000_w3/campaign.json \
  --result .runs/d6_interval_18_cover_v7_cap100000_w3/results/result_000032_0002100.json \
  --evidence d6_interval_18_class2100_increment_evidence.json.gz \
  --output d6_interval_18_class2100_increment_report.json
```

Run the independent default full replay:

```text
python3 verify_d6_interval_18_class2100_increment.py \
  --report d6_interval_18_class2100_increment_report.json \
  --evidence d6_interval_18_class2100_increment_evidence.json.gz \
  --output d6_interval_18_class2100_increment_verification.json
```

Focused source tests omit the approximately 30-second kernel replay by
default.  Set the environment variable to exercise it:

```text
python3 -m unittest -v test_d6_interval_18_class2100_increment.py
D6_INTERVAL_REPLAY=1 python3 -m unittest -v \
  test_d6_interval_18_class2100_increment.LiveEvidenceIntegrationTests.test_full_kernel_and_control_replay
```

No report, evidence archive, or verification JSON is part of the source-only
boundary.  Those artifacts must be generated only after the package sources
are committed, because the producer fails closed unless all four files equal
their blobs in the launch commit.
