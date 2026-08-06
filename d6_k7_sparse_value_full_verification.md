# Independent verification of the full K7 sparse-value layer

## Exact result

The restartable checker `verify_d6_k7_sparse_value_full.py` independently
re-certified all 3,195 sparse-value rejections in the archived 16,228-graph
K7 support-layer residue.  It found no independently surviving assignment.

For the recorded failing K7 seed of every rejected graph, the checker
exhausted:

- 209,519 eligible covers;
- 14,292 covers passing the frozen baseline rank layer;
- 1,506,596 labeled zero-factor support families; and
- 15,893 families reaching the one/two-support value stage.

Singleton propagation rejected 1,018,325 families by an empty mask and
472,378 by a disjoint required edge.  Every one of the remaining 15,893
families was rejected by an exact one/two-support rule.  The independent
search visited 153,159 small-support search nodes; every branch was pruned
before a complete admissible assignment, so the complete-assignment count is
zero.

The per-rule independent prune totals were:

```text
disjoint required small supports                 246670
duplicate one-defect type                         25967
forbidden two-defect cycle                        15228
incompatible branch distance                        663
one-defect component branches                      3749
two one-defects joined by two-defect types          7986
parallel two-defect type touches another               2
three one-defect points                              274
two-defect types repeat at a one-defect endpoint    21644
```

These totals exactly reproduce every nonzero sparse-value branch total in
the production archive, despite the independent implementation and search
driver.  The two production branch classes with total zero (three copies of
one two-defect type and inconsistent branch signs) also remained zero.

## Independence boundary

The verifier never imports `d6_k7_small_support_value.py`,
`profile_d6_k7_small_support_value.py`, or
`run_d6_k7_small_support_value_full.py`.  It reads them only to check the
hashes frozen in the production report.  It imports:

- `d6_k7_rank_reference.py` for classification by the already verified
  baseline cover filters; and
- `verify_d6_k7_support_full.py`, the separately committed independent
  labeled-support/singleton-propagation kernel.

The new verifier separately implements the actual one/two-support domain
enumeration, required-edge intersection search, simple-cycle enumeration,
one-defect uniqueness and cardinality rules, two-defect multiplicities,
component and branch rules, and branch-sign consistency through a parity
union-find.

As a mandatory exact algebra control it works directly in `Q(sqrt(7))` and
recomputes

```text
M^6 = -27 I
```

as well as the five fixed-point discriminants

```text
-3, -27, -108, -243, -243
```

and the one-defect orbit coefficients

```text
1/4, 1/3, 2/5, 1/2, 1, 0.
```

All decision operations use exact integers, rational numbers, and bit masks.
There are no floating-point decisions, numerical rank tests, random choices,
or transcendental functions.

## Quantifiers and soundness

For each recorded failing required K7, the checker reconstructs every
eligible zero-factor cover.  Covers eliminated by the frozen exact rank layer
remain eliminated there.  For every other cover, it enumerates every labeled
actual zero-factor support family and computes singleton-overlap closure.  If
that closure passes, every nonzero-factor vertex whose propagated mask has
size at most two receives every nonempty actual sub-support allowed by its
mask.  Vertices with masks of size at least three are deliberately ignored,
which can only make the screen weaker.

Required edges between selected vertices must have intersecting actual
supports.  Every remaining assignment is subjected only to the exact
one/two-support consequences proved in `d6_k7_small_support_value.md`.
Because every enumerated possibility fails, the recorded K7 seed—and hence
the whole required-edge graph—has no realization.  Candidate nonedges are
never required to be non-unit.

## Population and archive checks

Before doing mathematical work, the checker validates:

- the 17,764-graph rank-residue input and complete support archive;
- the exact ordered population of 16,228 support survivors;
- all 16,228 production sparse-value rows and their order;
- the 3,195 rejection indices and recorded failing seeds;
- every production report scalar and decision-column total;
- production input, archive, report, checkpoint, source, proof, profiler,
  runner, frozen-reference, and independent-support-verifier hashes; and
- the production checkpoint/configuration equality and all compressed and
  uncompressed artifact bindings.

The independent output is a deterministic 3,195-line JSONL archive.  Each
line records its graph index, failing seed, cover/family populations, exact
failure totals, deletion count, and search-node totals.

## Reproduction

The nine exact controls were run with:

```sh
D6_ALLOW_BACKGROUND_TEST_ONLY=1 \
  /Users/lukacs/claude/opengauss/venv/bin/python3 \
  -m unittest -v test_verify_d6_k7_sparse_value_full.py
```

The full replay used:

```sh
caffeinate -dimsu \
  /Users/lukacs/claude/opengauss/venv/bin/python3 \
  verify_d6_k7_sparse_value_full.py \
  --workers 11 \
  --checkpoint-every 32 \
  --progress-every 64
```

It completed with `PASS` in 95.316057333 seconds of internal wall time on the
machine recorded in the report.  The checkpoint was fsynced after every 32
results.  `--resume` validates all compatibility hashes, trims any uncommitted
suffix, and continues the exact ordered population.  The resume/finalization
path was separately exercised on an eight-case bounded control.

## Artifact hashes

```text
verifier source
6c0012e66d0b89e280f3cebed7a61f24cc1c6e128e726663a7b246e7da0edbbc

verification report
5f9ec90305e51c78a1ca1e827e354e7e51b5617f383a011eb554e7f68ae5de7e

compressed independent decision archive
7ff59f5cbcf4471195793e99d1d33369d50e26cf9949e78fb4d8e71a2f356219

uncompressed independent decision JSONL
cce8779f0353b81ddf6e9715c31693d904be952cae8407223957478d4df9e794

tracked completion checkpoint
3fe3878932c7a4fc69a6fea8bd390adc80df5071b4f889b10ac7074c89c3e669

production sparse-value report
cea64f3dde804e0766c5b77e373aa50338d5bee6e5e54f44749f4e387cf529aa

production compressed decision archive
9f70e0e267fe97bc2f6a2890cae47d8b1c882974d5f054f16b4673bfb45ee003

production uncompressed decisions
0ce9d2d2aeb18403fce17a0611594843946d4c07764fbb24946069f0997aab4c
```

The remaining trust assumptions are Python's exact integer/rational,
container, JSON/CSV/gzip, and SHA-256 implementations; the proved support and
sparse-value lemmas; and the frozen rank reference for already-existing
baseline cover classifications.  No IEEE-754 or `libm` assumption enters the
decisions.
