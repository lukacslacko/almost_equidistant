# Exact K7 support-or-pentad conjunction

## Result

This layer combines two already certified statements on the frozen ordered
258-graph K7 residue:

1. The full joint support campaign rejects 69 graphs.
2. The rank-two Schur pentad campaign rejects 36 labeled K7 covers on 34
   further graphs.

The conjunction rejects all 34 further graphs, disjointly from the earlier 69.
Thus this layer rejects exactly 103 of the 258 graphs and leaves the following
ordered residue size:

```text
258 - 69 - 34 = 155.
```

This is an exact reduction of the K7 class.  It is not by itself a proof of
the dimension-six theorem, because the 155 K7 graphs and the separate K6-only
class still require exact treatment.

## Seed-local quantifier

Fix a required unit `K7` seed.  The inherited exact cover filters enumerate a
finite set of still-possible zero-factor covers.  For each such cover, the
joint propagation/sparse-value layer either proves that no labeled support
family exists or leaves the cover alive.

The new graph-level rule is:

```text
There exists a K7 seed such that its inherited-passing cover set is nonempty,
and every cover in that set is either
  (a) joint-pre-capacity support-infeasible, or
  (b) rejected by an exact rank-two pentad certificate.
```

Any realization would induce one of the enumerated covers at every K7 seed, so
exhausting the covers at one seed is enough to reject the graph.  Candidate
nonedges remain optional throughout.

The builder replays the full cover partition, not only the older aggregate
`has_no_near_clique_cover` graph profile.  It reconstructs:

```text
K7 seeds                                      390
eligible covers                           145,209
inherited-passing covers                      738
joint-pre-capacity support-infeasible          342
pre-capacity support-passing                   396
```

For this particular 36-certificate pentad library, the full OR rule and the
more conservative rule requiring pentad certificates on every surviving cover
happen to reject the same 34 graphs.  There are 36 qualifying seeds: graphs
561011 and 3337921 each have two, and every other new graph has one.  All 36
certified covers are support-passing, so the support disjunct does not enlarge
the graph count in this frozen pilot.  Retaining the full quantifier is still
important for correctness and for future obstruction libraries.

Graph 3968822 is the key aggregate-profile control.  It has two K7 seeds and
two current covers in total, but one seed has exactly one current cover, which
the pentad rejects.  A graph-level marginal summary misses this; the seed-local
quantifier catches it.

## Artifact boundary

The builder and checker read only committed root artifacts.  In particular,
they do not read the ignored `.runs/d6_k7_rank_survivors.json` or per-graph
joint checkpoints.  The 258 adjacencies come from
`d6_current_residue_manifest_v2.json`, and the following upstream roots are
hash-pinned in both programs:

- current v2 residue manifest;
- full joint checkpoint copy, decision archive, certificate archive, report,
  and independent verification;
- full pentad report and independent verification;
- inherited degree-one and tetrad certificate/decision archives.

The production builder uses the frozen production propagation and sparse-value
kernels.  The independent checker does not import the builder or its production
kernel.  It has its own size-at-most-seven cover enumerator, replays all 738
current covers using the independently transcribed support/sparse-value
implementations, reconstructs every seed quantifier, and re-expands all 36
pentad certificates in fresh bounded subprocesses.

## Commands

Freeze the source before producing theorem artifacts.  Then run:

```bash
python3 build_d6_k7_pentad_joint_conjunction.py --workers 11
shasum -a 256 d6_k7_pentad_joint_conjunction_report.json
python3 verify_d6_k7_pentad_joint_conjunction.py \
  --report d6_k7_pentad_joint_conjunction_report.json \
  --report-sha256 <printed-report-sha256> \
  --workers 11 --pentad-workers 4
python3 -m unittest -v test_d6_k7_pentad_joint_conjunction.py
```

The 11 graph workers use the 12-logical-core Mac while leaving one logical core
for the OS.  Pentad re-expansion is limited to four fresh processes because its
exact polynomial dictionaries have materially higher peak memory.

The report stores the exact ordered 103-graph rejection list and its SHA-256,
the exact ordered 155-graph residue and its SHA-256, all 34 graph witnesses,
and all 36 qualifying seed/cover references.  The verifier output binds the
report by an explicit command-line SHA-256.
