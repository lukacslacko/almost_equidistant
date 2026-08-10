# Dimension-six current residue v7

## Result boundary

The v7 manifest is the complete certificate union currently available after
the independently checked v6 boundary.  It embeds every surviving adjacency
record and keeps exact and interval-dependent conclusions in separate trust
tiers.

```text
v6 K7-containing residue                              19
exact one-two-star rejections                          3
exact 2593240 star-boundary rejection                  1
exact 3949382 Schur rejection                          1
interval cap-500000 rejection                          1
focused interval cap-2000000 rejection                 1
v7 K7 residue                                         12

v6 K6-only residue                                   625
exact tight-same-Z0 rejections                          2
exact opposite-ray conjunction rejections            372
v7 K6-only residue                                   251

v7 combined residue                                  263
```

The ordered K7 residue is

```text
316173 2581209 3648882 3729907 3935560 3936310
3936435 3945490 3945555 3945557 3945564 3947605
```

with stable JSON SHA-256

```text
e545ea5a8fc9691ed98e885919e149f84dc7e8b4796d0fe9e6cb3a519df6fd8c
```

The ordered K6-only residue hash is

```text
a23da7bd9ae21eb7c078a08a7bce92f22d397d5c71f1e4f51d9cf8e59bbfa907
```

and the class-ordered combined 263-index hash is

```text
a50baf54bfff642e7c1711fe7041b5184df12f90de2cf52557cd086f10875127
```

Survival is only non-rejection by the present certificate union.  It is not a
realizability claim and does not by itself prove `f(6)=18`.

## Trust-tier accounting

Starting from the 911-graph v5 residue, exact algebra and graph logic now
reject 643 graphs.  Five additional, disjoint rejections depend on rigorous
interval arithmetic.  Thus the union rejects 648 and leaves 263.

The exact tier uses integer, rational, algebraic, polynomial, and graph-logic
certificates.  Floating point does not enter those 643 rejections.  The two
new K7 exact methods preserve optional-nonedge semantics: a candidate nonedge
is never required to have non-unit distance.

The interval tier uses only independently replayed `KILLED` records.
`ABORT`, `UNRESOLVED`, and `INFRA_ERROR` reject nothing.  Its five incremental
conclusions retain these assumptions:

- IEEE-754 binary64 basic operations and `sqrt` are correctly rounded, with
  endpoints expanded by `nextafter`;
- macOS `libm` cosine endpoint values are within eight ulps, with eight-step
  outward padding and interior extrema included;
- candidate nonedges are unconstrained and may also be unit distances.

## Post-v6 evidence union

The K7 evidence consists of four independently checked roots:

- the full 19-graph cap-500000 campaign, whose only `KILLED` graph is
  `3595554`;
- the exact one-two-star report, rejecting `2592657`, `3785980`, and
  `3888410`;
- the focused cap-2000000 replay rejecting `3936177`;
- exact source-bound certificates rejecting `3949382` and `2593240`.

The five exactly rejected graphs and two interval-rejected graphs are pairwise
distinct.  The manifest independently subtracts their union from the ordered
v6 K7 list.

The K6 evidence is sequential:

```text
v6 625 -> tight same-Z0 623 -> opposite-ray conjunction 251.
```

Both reports and their independent verifiers use exact arithmetic, retain
optional nonedges, pass the known realizable 18-point control, and expose
their complete ordered input/rejection/residue partitions.

The manifest pins every report, verifier output, decision table, certificate
archive, and checkpoint archive.  It also preserves every source hash map in
the constituent reports.  Where an interval wrapper later changed, the old
source is checked directly against its recorded Git commit blob rather than
against the current working file.  The official manifest and its checker also
archive their exact commands, UTC timestamps, platform, architecture, Python
version, logical CPU count, Git commit, branch, and working-tree digest.

## Independent verification

`verify_d6_current_residue_manifest_v7.py` imports neither the builder nor any
K6/K7 producer, kernel, or verifier.  It independently:

- checks all primary and transitive artifact hashes;
- validates every exact partition and every interval selection/KILLED list;
- verifies source maps against current files or the recorded historical Git
  blobs;
- checks control outcomes and interval trust assumptions;
- proves the rejection layers are disjoint where required;
- reconstructs all 263 embedded graphs from v6;
- checks every adjacency matrix is symmetric and loopless, has independence
  number at most two, contains a K6, and has the declared K7/K6-only class;
- recomputes every class and combined hash;
- rejects graph, source-root, or trust-tier tampering.

## Production after source commit

Do not create the official JSON artifacts until the builder, verifier, tests,
and this note have been committed.  Then run:

```text
python3 -m unittest -v test_d6_current_residue_manifest_v7.py
python3 build_d6_current_residue_manifest_v7.py
python3 verify_d6_current_residue_manifest_v7.py \
  --expected-manifest-sha256 <printed-manifest-sha256>
```

The builder and checker both require their four-file production package to
equal blobs in the source-boundary commit.  Commit the generated manifest and
verification JSON without changing that boundary.
