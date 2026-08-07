# Dimension-six current exact residue manifest v3

This package advances the self-contained level-19 boundary from the committed
v2 manifest using two independently verified result commits:

- K7 support-or-pentad conjunction:
  `01c5b18413b6bab5225558fa950971efeda58e0d`;
- K6-only arbitrary-subset Hall chain:
  `a21db74bacdf4c9c5c841ac137af326b331d8137`.

It creates new v3 files and does not modify the frozen v1/v2 manifests used by
earlier tools.

## Exact boundary

The K7 class changes as follows:

```text
258 v2 graphs - 69 joint-support - 34 pentad increment = 155.
```

The K6-only chain is kept explicit rather than jumping directly to its final
report input:

```text
977 v2 graphs - 116 PSD Z-matrix = 861
861           -  30 hereditary PSD Z = 831
831           -   9 arbitrary-subset Hall = 822.
```

All stage rejection sets are exact, ordered, and disjoint. The v3 boundary is
therefore `155 + 822 = 977` graphs. Its ordered class-concatenated index hash
is:

```text
71af8031e783b8a87d710156cae9cc4a19c575190a7c164a65935061f8ac1d7d
```

The manifest reconstructs both ordered subsets from the adjacencies embedded
in `d6_current_residue_manifest_v2.json`; it embeds all 977 final 19-row
adjacency records itself. Every record is checked for symmetry, looplessness,
and independence number at most two. Every K7-class record contains a `K7`.
Every K6-only record contains a `K6` and contains no `K7`.

## Positive control

The known realizable 18-point construction is stored in a separate top-level
control record, not among the 977 level-19 candidates. Its adjacency hash is

```text
7697e049810093251c017328fa1043647cc3bf12df7c61ecd5433225154497b5
```

and it has exactly 32 `K6` seeds and no `K7` seed. The upstream K6 exact
checkers also pass all 32 seeds.

## Rebuild and independent verification

After committing the v3 source boundary, run:

```bash
python3 build_d6_current_residue_manifest_v3.py
shasum -a 256 d6_current_residue_manifest_v3.json
python3 verify_d6_current_residue_manifest_v3.py
python3 -m unittest -v test_d6_current_residue_manifest_v3.py
```

The builder is deterministic and uses standard-library exact set/bitmask
operations. The checker imports neither the builder nor any K7/K6 production
kernel. It independently reconstructs the 155 and 822 ordered subsets from the
hash-pinned upstream reports/manifests, rechecks the K6 three-stage chain, uses
a separate clique-search implementation, and compares every embedded graph to
v2.

All directly used upstream sources, reports, verifications, result manifests,
certificate archives, and checkpoints are SHA-256 pinned. The manifest also
records the exact source commits underlying the two final result packages:
`4982e564eecf552f83160644a61fb78662b13e96` for K7 and
`4e1d50492ad76867b67d748792f21ca5c0a6a94f` for K6.

## Semantics and nonclaims

Edges are required unit distances. Candidate nonedges are unconstrained and
may also be unit. Points must be distinct, and allowed defect coordinates may
be zero. The K7 reduction is the full seed-local support-or-pentad disjunction;
the K6 arbitrary-subset layer chooses one subset per original span and does not
use the one-vertex extension.

Survival means only that no listed exact filter rejected the graph. It is not a
realizability claim, and the 977-graph boundary is not yet a proof that
`f(6)=18`.
