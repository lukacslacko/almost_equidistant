# Completed 181-class dimension-six 18-point interval campaign

`verify_d6_interval_18_cover_v7_completion.py` is a kernel-free independent
completion audit for the raw campaign in
`.runs/d6_interval_18_cover_v7_cap100000_w3`.

It independently rebuilds the 179-class maximum-edge greedy deletion cover
of the frozen 263-graph v7 residue, adds the two frozen non-positive backups,
and obtains the same ordered 181-class production selection.  It pins and
checks the campaign sources, compiled kernel, all eight input/control
artifacts, campaign configuration, 181 atomic checkpoints, raw report, and
decisions TSV.  The TSV is reconstructed byte for byte from the checkpoints.

The completed raw accounting is:

```text
KILLED          1   class 2100, population K7
ABORT         179   no rejection claim
UNRESOLVED      1   class 7259, no placement order, no rejection claim
INFRA_ERROR     0
```

Only `KILLED` has certificate semantics.  The sole such deletion class is
2100; its required-edge subgraph consequence for parent 3945564 is handled
by the separate class-2100 producer/checker package.  This completion audit
does not rerun the interval kernel and does not turn any numerical failure,
`ABORT`, or `UNRESOLVED` result into evidence of non-realizability.
Candidate nonedges remain unconstrained and may also be unit distance.

Raw SHA-256 boundaries:

```text
campaign   70ba451094eb4eb5d039fc059ab8e3c351ebeed157e2eae39cc35fe4dcb3e348
report     0a8dc0d875f3f3f547c2feec3992152d2b841c4a4ba4d6b1e2dad97397ab0a53
decisions  f7416adec9ad442d4749eef3f5afd3e1a4ec25972267d432302dc8cf68a5aab8
selection  8c624c969b517bb5663dc32bdf577ab3addbe8a108458b2803443c17f294d8e6
results    cf5807d8408a76dda996ca5bbfaccbb76260d45c2311cce9f0ac3c9989060b5c
```

For durable review, byte-identical copies of the completed raw report and
decision table are committed as
`d6_interval_18_cover_v7_completed_report.json` and
`d6_interval_18_cover_v7_completed_decisions.tsv`.  Their hashes are the
`report` and `decisions` values above.  The local `.runs` directory retains
the 181 atomic checkpoint envelopes used by the independent audit.

Run the focused tests without calling the interval kernel:

```sh
python3 -m unittest -v test_d6_interval_18_cover_v7_completion.py
```

After committing the three package source files, create a source-bound
machine-readable verification artifact with:

```sh
python3 verify_d6_interval_18_cover_v7_completion.py \
  --source-commit "$(git rev-parse HEAD)" \
  --output d6_interval_18_cover_v7_completion_verification.json
```

The verifier fails closed if any raw hash, source/input hash, selection root,
checkpoint identity, status counter, decisions byte, or package source blob
differs from the pinned boundary.
