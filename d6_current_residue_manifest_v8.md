# Dimension-six current residue v8

## Result boundary

The v8 manifest is an exact two-graph delta on the frozen, independently
verified v7 residue.  It consumes only the immutable v7 manifest boundary and
the source-bound saturated-singleton-basis report, certificate archive, and
independent verification.

```text
v7 K7-containing residue                            12
new K7 rejections                                    0
v8 K7-containing residue                            12

v7 K6-only residue                                 251
exact saturated-singleton-basis rejections           2
v8 K6-only residue                                 249

v8 combined residue                                261
```

The two removed K6-only indices are

```text
3138618 3673988
```

Their stable JSON SHA-256 is

```text
45ce66ccf4db5fe1a565a576cd16fc719e51b20c6ea9fee858dc3682463e9e84
```

The K7 class is copied value-for-value from v7, including all 12 embedded
graphs and its exact/interval accounting.  Its ordered-index hash remains

```text
e545ea5a8fc9691ed98e885919e149f84dc7e8b4796d0fe9e6cb3a519df6fd8c
```

The ordered 249-index K6-only residue hash is

```text
f274aab3f2d408fce34c7cd7eb4518621fa9d634a8943e3ee1a38d76df9f3d31
```

and its embedded-graph hash is

```text
92894350a8a3f76e791cd419049e34a52903c0ff6a4a45ed7a262ab801ca630e
```

With class order `K7, K6_only`, the combined 261-index hash is

```text
2391a93a3629517363106603bdad00be9b6f960966d9089ef37d93be2988c213
```

and the hash of the numerically sorted combined indices is

```text
2ac404f7260a130549f572ba5c34812b95b9df96bab7ab546b2f3f4c6656ce65
```

Survival means only that the present certificate union has not rejected the
graph.  It is not a realization claim and does not prove `f(6)=18`.

## Exact increment

The new K6 layer starts from the exact ordered 251-index K6-only list embedded
in v7.  The saturated-singleton report declares that same input list and
partitions it into the two rejected indices above and the ordered 249-index
residue.  Its independent verifier replayed all 251 graph decisions and both
complete rejected seed-certificate trees.  Producer and verifier both pass
all 32 K6 seeds of the known realizable 18-point configuration.

The manifest additionally checks:

- every report, verification, and certificate-archive hash;
- the complete ordered input, rejection, and residue partitions;
- all 251 per-graph result rows;
- the certificate archive's two-index coverage;
- all producer and verifier source hashes against the recorded launch commit;
- the exact semantics that candidate nonedges remain unconstrained, allowed
  defect coordinates may vanish, and the new rule prunes only one complete
  actual-support leaf at a time.

The cumulative post-v6 K6 rejection list is stored in sequential evidence
order: the 374 rejections inherited from v7 followed by the two v8 increment
indices.  Its stable hash is

```text
621e3d8f9859cbab831f307ef20d57a0a1ffeba988b5c3178a0352f26ea7fa0b
```

## Trust-tier accounting

Starting from the 911-graph v5 residue, exact certificates now reject 645
graphs.  The five inherited interval-only rejections remain separate and
disjoint.  Their union rejects 650 graphs and leaves 261.

No floating-point computation enters the two new rejections.  The inherited
interval tier retains the v7 assumptions verbatim:

- IEEE-754 binary64 basic operations and `sqrt` are correctly rounded, with
  interval endpoints expanded by `nextafter`;
- macOS `libm` cosine endpoints are within eight ulps, with the recorded
  outward padding and interior-extremum checks;
- candidate nonedges are unconstrained and may also be unit distances.

Only independently replayed `KILLED` interval records reject anything.
Inherited `ABORT`, `UNRESOLVED`, and `INFRA_ERROR` records reject nothing.

## Independent verification

`verify_d6_current_residue_manifest_v8.py` imports neither the v8 builder nor
any K6 production or verification module.  It independently:

- checks the five immutable primary artifact hashes;
- validates the frozen v7 manifest, its `PASS` verification, source package,
  counts, hashes, and embedded graphs;
- validates the saturated-singleton partition, controls, source boundary,
  and certificate coverage;
- proves the two new indices are members of v7 `K6_only`, absent from `K7`,
  and exactly the records removed from the embedded K6 graph list;
- proves the entire K7 object is unchanged;
- validates every surviving graph as symmetric and loopless, with
  independence number at most two and the declared K7/K6-only clique class;
- recomputes class, combined, sorted, tagged-index, and tagged-graph hashes;
- checks exact/interval accounting and rejects graph, evidence, source-root,
  control, or accounting tampering.

## Production after source commit

Do not create the official JSON artifacts until the builder, independent
verifier, tests, and this note have been committed.  Then run:

```text
python3 -m unittest -v test_d6_current_residue_manifest_v8.py
python3 build_d6_current_residue_manifest_v8.py
python3 verify_d6_current_residue_manifest_v8.py \
  --expected-manifest-sha256 <printed-manifest-sha256>
```

The builder and verifier require all four v8 package sources to equal blobs
in the recorded source-boundary commit.  The generated manifest and
verification JSON may then be committed without changing that boundary.
