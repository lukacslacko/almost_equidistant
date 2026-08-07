# Dimension-six current residue v6

## Result boundary

The v6 manifest is a cross-method certificate union over the independently
verified v5 residue.  It reconstructs and embeds every surviving 19-vertex
adjacency record.

```text
v5 K7-containing graphs                              155
exact K7 algebra/graph-logic rejections              133
exact K7 residue                                      22
verified interval KILLED records                      24
interval KILLED already covered exactly               21
new interval-dependent K7 rejections                   3
v6 K7 residue                                          19

v5 K6-only graphs                                     756
tight-Hall exact rejections                           107
singleton-fan exact rejections                         15
repeated-two-support-arm exact rejections               9
v6 K6-only residue                                    625

v6 combined residue                                   644
```

Survival means only that the current certificate union did not reject the
graph.  It is not a realizability claim and does not yet prove `f(6)=18`.

The ordered K7 residue is

```text
316173 2581209 2592657 2593240 3595554 3648882 3729907
3785980 3888410 3935560 3936177 3936310 3936435 3945490
3945555 3945557 3945564 3947605 3949382
```

with stable JSON SHA-256

```text
af25b842ca9eaa6e9721f51d0400be2436d930b72cc33c72bc366f0ac5e41290
```

The ordered 625-graph K6-only residue has stable JSON SHA-256

```text
04d9475217e294efed0c6eac34a7fad88616c90ad7d37a2525776a713e6dd0e1
```

The combined class-ordered 644-index list has stable JSON SHA-256

```text
907b6fe460846ab613e7d35a94dbfc6eb1eb94322a36f10b4500a8ce4fcc6ccf
```

## Trust tiers

The manifest deliberately does not blur exact and interval conclusions.

### Exact algebra and graph logic

The exact K7 chain rejects 133 graphs and leaves 22.  The exact K6 chain
rejects 131 graphs and leaves 625.  These 264 rejections use exact integer,
rational, algebraic, polynomial, and graph-logic certificates with
independently replayed quantifiers.  Floating point does not enter these
rejections.

### Interval certificates

The cap-100000 interval campaign has 24 independently replayed `KILLED`
records.  Twenty-one were already rejected by the exact K7 chain.  The three
new rejections are

```text
423661 424226 3936176
```

`ABORT`, `UNRESOLVED`, and `INFRA_ERROR` reject nothing.  The three
incremental conclusions retain the interval verifier's explicit assumptions:

- IEEE-754 binary64 basic operations and `sqrt` are correctly rounded, with
  every endpoint expanded using `nextafter`;
- macOS `libm` cosine endpoints are within 8 ulps, and the kernel pads both
  ways by eight `nextafter` steps while including interior extrema;
- candidate nonedges remain unconstrained and may also be unit distances.

Thus the exact-only residue has 647 graphs; the mixed certificate union has
644.

## Audited source graph

The builder pins these eight primary artifacts:

```text
d6_current_residue_manifest_v5.json
  164b2a12813845ef1f6d6ca5e3713ed1d2edac4be58563c1648a39ff8580c0b5
d6_current_residue_manifest_v5_verification.json
  e871431b8b28fc04f0e58922ff3a6386054d15d37de5c9f9e67601ede47e4b46
d6_k7_one_free_conjunction_report.json
  181f63015d373fb69f37f4f390cdd8ad50bc32b5a2d97a36fd0bac876d0877a2
d6_k7_one_free_conjunction_verification.json
  ebaacfe663728eff427c9ba2280fc4261d26dfe10fad21e1a00982b9d31ac062
d6_interval_v5_k7_cap100000_report.json
  85b175f3b3348204c3ba9065b30432458528bdd2e22b471e64fc56341480bac3
d6_interval_v5_k7_cap100000_verification.json
  a56fb1561b778c55bffb433bf0eea90e5cca823573d0717e74f2a73528bc1b10
d6_k6_repeated_two_support_arm_report.json
  ea8d9438d062b92857a1057b950e76b8ae08bbe1bf97b4327130fd4058171ebe
d6_k6_repeated_two_support_arm_verification.json
  3abf48bab092671a9dff03c2dc47c0063a75515cd219e1d5a4bfa9a5fa797b5b
```

It also follows and checks every named transitive hash in the K7 algebraic
report and the three-stage K6 chain, including production sources,
independent verifications, and raw certificate archives.  It checks that the
K6 ordered partitions concatenate exactly:

```text
v5 756 -> tight Hall 649 -> singleton fan 634 -> repeated arm 625.
```

The source-independent v6 checker imports neither the builder nor any K6/K7
production kernel.  Starting from the pinned v5 adjacency records and result
lists, it independently reconstructs both set differences and their union,
checks every embedded 19-row adjacency matrix, verifies independence number
at most two, verifies the K7 versus K6-only clique classification, recomputes
all stable hashes, audits positive controls, and rejects trust-tier tampering.

## Production after source commit

Do not generate official v6 artifacts before the builder, verifier, tests,
and this document are committed.  After that source boundary is frozen:

```sh
python3 build_d6_current_residue_manifest_v6.py \
  --output d6_current_residue_manifest_v6.json

shasum -a 256 d6_current_residue_manifest_v6.json

python3 verify_d6_current_residue_manifest_v6.py \
  --manifest d6_current_residue_manifest_v6.json \
  --expected-manifest-sha256 <printed-manifest-sha256> \
  --output d6_current_residue_manifest_v6_verification.json

python3 -m unittest -v test_d6_current_residue_manifest_v6.py
```

The official JSON files are immutable proof artifacts.  Later replays should
write to separate paths.
