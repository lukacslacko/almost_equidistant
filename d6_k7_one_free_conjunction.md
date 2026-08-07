# K7 one-free-edge corpus conjunction

## Frozen input boundary

This production layer starts from the independently verified 24-graph exact
K7 residue produced by the double-pin and odd-cycle/full-support-pin layers.
Its four required upstream artifact hashes are:

```text
d6_k7_double_pin_conjunction_report.json
  7c90c9a518243de4095f4ec4394bb7b1d22c5889844e90cd7fe4d895f22ac873
d6_k7_double_pin_conjunction_verification.json
  890568a80b8363de84997db86ab5271fee13f707544f91145d15f4976512a33b
d6_k7_full_pin_increment_report.json
  1d2c3a315aaaa8e7536c532d86935f528a1c80ca05159977377ca795a61409ab
d6_k7_full_pin_increment_verification.json
  cacab17a687a7666824bf8548fa052dd76be477e6a37b3da976999972b9a68cc
```

The input 24-list has stable JSON SHA-256

```text
71d9ae101445ae90d08699ede81ca405bf5bfe956fda01ee384203fa113d1779
```

## Quantifier replay and coverage accounting

For every input graph, every unit `K7` seed, every inherited-current
zero-factor cover, and every labeled support family, the builder independently
replays support propagation and the sparse-value screen.  A surviving family
then passes through four parallel pipelines:

```text
old         = double-pin + full-support odd-cycle pinning
singleton   = old + singleton-overlap one-free edge
correlated  = old + pinned-product correlated one-free edge
combined    = old + both new one-free layers
```

A cover passes if at least one family passes, a seed is infeasible if all its
current covers are infeasible, and a graph is rejected if any of its required
`K7` seeds is infeasible.  The separate singleton and correlated pipelines
measure individual marginal coverage; the combined pipeline also detects any
seed-level synergy between them.  The report records rejected and surviving
ordered lists and stable hashes for all four pipelines, individual graph
overlap, combined-only synergy, exact family counts, all new certificates,
and the first passing witness for each cover/pipeline.

The `old` replay must leave exactly the frozen 24-list.  Any old-layer
rejection is a hard consistency failure and aborts production.

## Independent verification

`verify_d6_k7_one_free_conjunction.py` does not import the production builder
or either production one-free locator.  It separately reconstructs:

- the 24 input graphs from the v5 corpus and all inherited-current covers;
- labeled supports, propagated masks, and the sparse-value screen;
- nonbipartite coordinate components and full-support pinning;
- the singleton one-free formula and exact reciprocal classification;
- correlated component-sign variables and the complete local
  `Q(sqrt(7))` equation enumeration;
- every cover, seed, graph, pipeline, count, certificate, and summary hash.

It requires an explicit SHA-256 for the production report and for all four
upstream artifacts, and writes a separate atomic PASS artifact.

## Production commands after source commit

The Mac currently reserves CPU for two MPS least-squares workers, so this
campaign uses nine process workers once no earlier exact run is active:

```sh
python3 build_d6_k7_one_free_conjunction.py \
  --base-report-sha256 7c90c9a518243de4095f4ec4394bb7b1d22c5889844e90cd7fe4d895f22ac873 \
  --base-verification-sha256 890568a80b8363de84997db86ab5271fee13f707544f91145d15f4976512a33b \
  --full-pin-report-sha256 1d2c3a315aaaa8e7536c532d86935f528a1c80ca05159977377ca795a61409ab \
  --full-pin-verification-sha256 cacab17a687a7666824bf8548fa052dd76be477e6a37b3da976999972b9a68cc \
  --workers 9 \
  --output d6_k7_one_free_conjunction_report.json
```

Hash the report, then run:

```sh
python3 verify_d6_k7_one_free_conjunction.py \
  --report d6_k7_one_free_conjunction_report.json \
  --report-sha256 <printed-report-sha256> \
  --base-report-sha256 7c90c9a518243de4095f4ec4394bb7b1d22c5889844e90cd7fe4d895f22ac873 \
  --base-verification-sha256 890568a80b8363de84997db86ab5271fee13f707544f91145d15f4976512a33b \
  --full-pin-report-sha256 1d2c3a315aaaa8e7536c532d86935f528a1c80ca05159977377ca795a61409ab \
  --full-pin-verification-sha256 cacab17a687a7666824bf8548fa052dd76be477e6a37b3da976999972b9a68cc \
  --workers 9 \
  --output d6_k7_one_free_conjunction_verification.json
```

Finally run all kernel and artifact tests:

```sh
python3 -m unittest -v \
  test_d6_k7_one_free_edge.py \
  test_d6_k7_correlated_one_free_edge.py \
  test_d6_k7_one_free_conjunction.py
```

No conjunction artifact should be generated before the builder, independent
verifier, tests, and this document are committed.
