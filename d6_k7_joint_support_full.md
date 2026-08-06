# Full exact K7 cover/support conjunction

## Certified conclusion

The frozen rank-one-tetrad/pattern-954 complement contains 258 dimension-six
level-19 candidate graphs.  The full conjunction campaign proves that 69 of
these graphs are impossible and leaves 189 graphs unresolved:

| status | graphs |
|---|---:|
| `JOINT_PRE_CAPACITY_REJECTED` | **69** |
| `SURVIVOR` | **189** |
| capacity-only rejection | 0 |
| unresolved | 0 |
| infrastructure error | 0 |

The ordered 258-graph input has index-list SHA-256

```text
55ab329dbe3ae4dbfa10378ff023a168fc6ab306ffc7af14bbd9790a68108a09
```

and is the exact complement recorded in
`d6_k7_rankone_pattern_union.json` (SHA-256
`1a54fabeb3ebf7f8e5485a88bc52fcfd07fb9ab8706d29ab2a664a895f79e599`).
The full campaign source was committed and pushed before execution at

```text
2265af8aab199e806a98f75f678a0bc9a251549a
```

The 69 rejected graph indices are explicit in the production report and
decision archive.  The certificate archive contains 90 failing K7 seeds and
the first labeled-support failure for every one of their 112 jointly relevant
covers.

`SURVIVOR` is not a realizability claim.  It says only that this exact
conjunction did not reject the graph.

## Exact quantifiers and theorem

Fix a required K7 seed in a hypothetical realization.  The frozen K7
coordinate reduction associates to every outside point a seven-coordinate
defect vector.  Its zero-factor vertices form an eligible cover `Z` of the
disjoint-required-edge graph.  Consequently every realization determines:

1. one eligible cover `Z`;
2. one labeled family of exact zero-factor supports allowed by the candidate
   defect masks; and
3. propagated nonzero-factor support masks satisfying all required-edge and
   rank consequences.

For each K7 seed, the full campaign exhausts every eligible cover.  A cover is
discarded only when it fails one of the already exact, independently checked
cover layers:

- the frozen baseline rank/support rules;
- strict-H at a saturating clique;
- a recorded exact degree-one Schur certificate; or
- a recorded exact degree-four rank-one tetrad certificate.

For every cover surviving those gates, the campaign exhausts every labeled
zero-factor support family.  It then applies the frozen singleton-overlap
propagation and one/two-defect sparse-value rules.  A graph is labeled
`JOINT_PRE_CAPACITY_REJECTED` only if at least one required K7 seed has no
cover/support family surviving this complete conjunction.

This is sound because a realization would have to induce one of the exhausted
cover/support choices at every required K7 seed.  If every choice at one seed
fails a necessary exact condition, no realization exists.

The graph semantics are unchanged throughout:

- every candidate edge is required to have distance one;
- candidate nonedges are unconstrained and may also have distance one;
- all represented points must be distinct;
- no epsilon separation or numerical failure is used;
- an unresolved computation is never converted into a rejection.

## Complete production accounting

The 258 graphs contain 390 required K7 seeds and 145,209 eligible covers.  The
cover partition is exact:

```text
baseline-rejected covers                         143259
strict-H-rejected covers                            376
degree-one-certificate-rejected covers              159
tetrad-certificate-rejected covers                  677
jointly relevant covers                             738
                                                   ------
total                                             145209
```

The 738 relevant covers generated 23,350 labeled zero-factor support
families:

```text
empty propagated mask                            10489
disjoint required edge after propagation         12147
families passing propagation                       714
  sparse-value failures                             36
  pre-capacity survivors                            678
```

At the cover level, 342 of the 738 relevant covers lose every labeled family;
396 retain at least one.  At the graph level this exact intersection rejects
69 graphs and leaves 189.

## Independent verification

`verify_d6_k7_joint_support_full.py` imports neither the production full
runner, the production pilot evaluator, nor the support-capacity kernel.  It
validates the complete ordered decision archive, compressed and uncompressed
hashes, the certificate population, checkpoint/configuration binding, source
hashes, and inherited artifacts.  It then uses the previously independent
support and sparse-value implementations to replay all 90 recorded failing
seeds.

The independent replay obtained:

```text
joint graph rejections verified                      69
failing K7 seeds replayed                             90
jointly relevant covers reconstructed                112
labeled support families exhausted                  2355
empty propagated masks                              1413
disjoint required edges                              938
propagation survivors                                  4
sparse-value failures among those four                 4
independent survivors                                  0
errors or unresolved cases                             0
```

The verification result is `PASS`.  The production run used 11 workers and
completed in 33.014 seconds; the independent 11-worker replay completed in
3.878 seconds.  Elapsed times are measurements, not mathematical decisions.

## Capacity layer: measured but zero credit

The support-type multiplicity theorem says that at most `|S|` outside points
can have the same exact nonempty defect support `S`.  The production evaluator
also ran the finite capacity-aware CSP on every one of the 678 families that
reached it.  All 678 were feasible, covering all 396 pre-capacity-surviving
covers.  There were:

```text
capacity-only graph or cover rejections                 0
capacity-infeasible families                            0
capacity-unresolved families or covers                  0
capacity search nodes                                7742
capacity Hall-flow checks                            7742
capacity Hall-flow prunes                               0
```

Therefore none of the 69 graph rejections is credited to the capacity layer.
They are entirely cover/support cross-layer conjunctions.

There is one known bookkeeping defect in the frozen production report.  Its
totals field named `capacity_maximum_depth_observed` is `2436`, because the
full runner aggregates the 258 per-graph maxima with `Counter.update` and thus
sums them.  A read-only scan of all 258 hash-bound atomic checkpoint records
shows that the true maximum capacity-search depth is **12**; the per-graph
maxima sum to 2436.  This mislabeled diagnostic is not read by any decision,
certificate, or independent verification path.  The frozen source and output
artifacts were not edited after discovery.

## Controls and trust boundary

The focused suite passed all 14 controls before production.  These include
120 deterministic tiny random capacity instances cross-checked against full
Cartesian enumeration, exact positive and negative capacity cases, explicit
candidate-nonedge optionality, checkpoint/configuration binding, and separate
joint-versus-capacity decision columns.

The known realizable 18-point lower-bound graph has no required K7, so this
K7-conditioned layer correctly records it as
`NOT_APPLICABLE_NO_K7` (`18` vertices, `0` K7 seeds).  It is not claimed as a
direct positive control for a K7-conditioned rule; the cover-level synthetic
controls exercise that kernel.

The mathematical conclusion trusts:

- the documented K7 coordinate, rank, strict-H, Schur/tetrad, labeled-support,
  singleton-propagation, and sparse-value lemmas;
- the independently verified inherited certificates and their SHA-256
  bindings;
- Python arbitrary-precision integer, bit-mask, container, JSON/CSV/gzip, and
  SHA-256 implementations; and
- the independently replayed finite quantifiers described above.

No floating-point rank test, failed numerical optimization, IEEE-754/libm
bound, or transcendental computation enters a rejection.

## Exact commands and restart procedure

Production command:

```sh
caffeinate -dimsu \
  /Users/lukacs/claude/opengauss/venv/bin/python3 \
  run_d6_k7_joint_support_full.py \
  --workers 11 \
  --max-inflight 22 \
  --node-limit 200000 \
  --checkpoint-dir .runs/d6_k7_joint_support_full \
  --decisions d6_k7_joint_support_full_decisions.tsv.gz \
  --certificates d6_k7_joint_support_full_certificates.jsonl.gz \
  --checkpoint-copy d6_k7_joint_support_full_checkpoint.json \
  --report d6_k7_joint_support_full_report.json \
  --progress-every 16
```

Independent verification command:

```sh
caffeinate -dimsu \
  /Users/lukacs/claude/opengauss/venv/bin/python3 \
  verify_d6_k7_joint_support_full.py \
  --report d6_k7_joint_support_full_report.json \
  --report-sha256 \
    b3fd2752a371f8909b79b50b8e2b6135f93423d0adf7286aee26d4a8c0643326 \
  --workers 11 \
  --output d6_k7_joint_support_full_verification.json
```

To resume an interrupted production run, repeat the production command with
`--resume`.  Each graph result is an atomically replaced, configuration-bound
file in `.runs/d6_k7_joint_support_full/`.  Resume rejects source, selection,
node-limit, or schema mismatches; validates all completed graph identities;
and processes only the missing ordinals.  Final archives are deterministic
gzip files.  Infrastructure exceptions are durable `INFRA_ERROR` records and
make finalization exit nonzero.

## Artifact and source hashes

```text
source boundary commit
  2265af8aab199e806a98f75f678a0bc9a251549a

run_d6_k7_joint_support_full.py
  db74a8ac4635f2acbf36672d813340a9b94e78e875be5cf370063574f9dab493
verify_d6_k7_joint_support_full.py
  28f05687ea8a5f1f4be51723e727716f4e8b4e7bb910c4772a1dbcb2ab76d74a
d6_k7_support_capacity.py
  df8c010faac37f9de5481316cc8ccf58e8df46230bd9d1e8542f7a4322b64996
run_d6_k7_support_capacity_pilot.py
  9af542217d178bec2a71cb4c30faabc3dea8c8279053d729d924d9ca35eb128a
independent bounded support verifier used by the full verifier
  d16e7e9037d14aea797ff90df51df2ca1bb583b412a3f65992d995964263da6b

d6_k7_joint_support_full_report.json
  b3fd2752a371f8909b79b50b8e2b6135f93423d0adf7286aee26d4a8c0643326
d6_k7_joint_support_full_decisions.tsv.gz
  b6040e7796e3d7df511a4e460ad71082e63c074ac2cf2713c19d8f3b772bcc90
uncompressed decisions TSV
  d884eec3e09c73ad8cb970f6148e8a962d262a46e420616cff997868da0f791a
d6_k7_joint_support_full_certificates.jsonl.gz
  774d3d2119d5ce5eb965a689853d85704c44e081b12ff0947f360b4571d00522
uncompressed certificate JSONL
  583cc33eb0b0e898c800ccf462f7cb8a52e136b20225c532e0c11b4bcfc41ca6
d6_k7_joint_support_full_checkpoint.json
  91638dbcde9e18feb598521c801578a5e95a6a3c6fe48a1568c079d1c1d8ef3a
d6_k7_joint_support_full_verification.json
  41132fb0cb6d7fcb9c02d4e5171de9c0b8afd9c7eab4da2a2776f2518a1facd1

production configuration SHA-256
  73b95fc921805be327e827d0cfc65f5632c0d79cd7ac38ef3630e8a2dc3ec52a
```
