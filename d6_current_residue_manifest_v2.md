# Dimension-six current exact residue manifest v2

This package records the current exact level-19 residue at commit boundary
`fc51874458196c0391530d0e98968f530ab738c9` without changing the frozen
`d6_current_residue_manifest.json` used by earlier tools.

## Exact boundary

The independently verified K7 tetrad/pattern-954 union rejects 12,581 of its
12,839 ordered inputs. Its ordered complement has 258 indices and hash
`55ab329dbe3ae4dbfa10378ff023a168fc6ab306ffc7af14bbd9790a68108a09`.

The independently verified fused K6-only layer rejects 13 of the frozen 990
inputs. Its ordered complement has 977 indices and hash
`27e435506ecaecb936acaef63e5873c3538771be17d1a52f12df82df8dedd947`.
The known realizable 18-point positive control passes the exact fused K6
checker for all 32 K6 seeds.

The classes are disjoint. In the declared class order `K7`, `K6_only`, the
combined boundary therefore has exactly 1,235 graphs. Its principal hashes are:

- ordered concatenated indices:
  `62d73f5604e6ce06de72335ea77f8e612e6b79027ee328518f91dc1203328c4d`;
- class-tagged ordered indices:
  `dff478080c33a71142c77b44dc43033bd755f407c3a9fa08801c73ea52a4a49c`;
- class-tagged ordered adjacency records:
  `00ba3bbbb74d11314ae0c50e8c1cadbf75905682552659145b3f279f3d5f6e94`.

The machine-readable manifest embeds every 19-row adjacency bitmask in residue
order. Downstream tools can iterate
`manifest["classes"][class_name]["graphs"]` in `manifest["class_order"]`
without access to historical transient files. The original adjacency sources
and their SHA-256 hashes remain recorded for provenance and independent replay.

## Rebuild and verification

Use the pinned Python 3.11 environment (ordinary Python 3 also suffices because
the package uses only the standard library):

```text
/Users/lukacs/claude/opengauss/venv/bin/python3 build_d6_current_residue_manifest_v2.py
/Users/lukacs/claude/opengauss/venv/bin/python3 verify_d6_current_residue_manifest_v2.py
/Users/lukacs/claude/opengauss/venv/bin/python3 -m unittest -v test_d6_current_residue_manifest_v2.py
```

The builder is deterministic: a fresh rebuild is byte-identical to
`d6_current_residue_manifest_v2.json`, whose SHA-256 is
`961dac1f9b44bb541e2c5bd1626027826eeea7bc628bf5ff0a85ab26afa949f4`.
The independent checker imports neither the builder nor any production engine.
It separately reconstructs both ordered complements, checks the K7 set
partition, checks the K6 independent report and positive control, validates all
1,235 embedded symmetric loopless adjacencies, and verifies the source hashes.

Rebuilding or independently replaying adjacency provenance requires the
historical `.runs/d6_k7_rank_survivors.json` file with SHA-256
`7a0a350142930e4e13830ea44c4fd217f25aeb07bbd6286e51e316f72712c4f9`.
The committed v2 manifest itself embeds the 258 K7 records, so downstream
residue analysis does not require that ignored historical source.

## Semantics and nonclaims

Each candidate edge is a required unit distance. Candidate nonedges remain
unconstrained and may also be unit. Points must be distinct, and allowed defect
coordinates may be zero. Manifest membership uses only exact set filtering and
integer adjacency copying; it makes no floating-point decision.

Survival means only that the graph was not rejected by the exact layers listed
here. It is not a realization claim. In particular, the combined count 1,235
is not yet a proof that `f(6) = 18`, and this boundary artifact does not replace
the upstream rejection certificates.
