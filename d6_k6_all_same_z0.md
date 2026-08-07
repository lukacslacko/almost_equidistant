# All K6 systems on one common Z0

This exact layer applies every current K6 necessary system to the same actual
zero-Lorentz-factor set.  On the manifest-pinned 990-graph K6-only residue it
rejects four graphs:

```text
652900, 2301548, 2842523, 3289061.
```

The first and third were already rejected by the two-system same-`Z0` layer.
Coupling the older ordinary zero-forcing rank bound exposes the two new
indices `2301548` and `3289061`.  An independent checker recomputed all 990
decisions and every archived per-`Z0` certificate and returned `PASS`.

## Exact quantifier argument

Fix a required unit `K6` seed in a putative realization and let

```text
Z0 = {x : ell_x = 0}
```

be its actual zero-Lorentz-factor set.  This is one geometrically determined
set.  It is eligible, has size at most six, and its allowed defect masks admit
a matching.  For this same `Z0`, every realization satisfies three exact
necessary systems.

### Ordinary zero-forcing rank

For each connected bipartite component `C=A union B` of `L-Z0`, the defect
Gram matrices on the two sides have the exact off-diagonal graphs `G[A]` and
`G[B]`.  Ordinary zero forcing therefore gives

```text
rank u(A) >= |A|-Zf(G[A]),
rank u(B) >= |B|-Zf(G[B]).
```

The two spans and the `|Z0|` orthonormal zero-factor vectors are mutually
orthogonal, so

```text
|Z0| + |A|-Zf(G[A]) + |B|-Zf(G[B]) <= 6.          (ZF)
```

The lightlike component case is stronger and still implies `(ZF)`.  The exact
derivation, including why the Gram pattern has required nonzero entries where
ordinary zero forcing needs them, is frozen in
`d6_k6_bipartite_rank_audit.md`.

### Normal-inertia block support

The exact inertias of `I+Adj(G[A])` and `I+Adj(G[B])` give rank lower bounds
for both generic Lorentz orientations.  The `A` span, `B` span, and every
one-dimensional `Z0` block are mutually orthogonal and lie in their allowed
coordinate unions.  Every subset of blocks must therefore obey its exact
support-Hall dimension inequality.  A bipartite component must pass one of
the two generic orientations or its separate lightlike dimension/matching
alternative.  This is system `(B)`, proved in
`d6_k6_normal_block_support.md`.

### Non-bipartite light-ray support

Every non-bipartite component of `L-Z0` lies on one of two Lorentz light
rays.  Some symmetry-reduced coloring of those components must have
allowed-mask matchings in both ray bins and a joint actual-support assignment.
The support of every `Z0` vertex is one shared choice in both bins.  Call this
system `(N)`.

The earlier exact filters established the three weaker statements

```text
exists Z0: ZF(Z0),
exists Z0: B(Z0),
exists Z0: N(Z0)
```

separately.  A realization requires

```text
exists Z0: ZF(Z0) and B(Z0) and N(Z0).             (all-same-Z0)
```

The production evaluator exhausts every eligible subset of size at most six
and rejects a K6 seed only if `(all-same-Z0)` fails.  One impossible required
K6 seed rejects the graph.  This is a finite exact quantifier argument, not a
heuristic search failure.

Candidate nonedges remain unconstrained and may have unit distance.  Allowed
defect masks are upper bounds on actual support; allowed coordinates may be
zero.  No condition requires a candidate nonedge to be genuinely non-unit.

## Exhaustive certificates

`d6_k6_all_same_z0_certificates.json` stores all permitted `Z0` choices for
the first impossible seed of each rejected graph.  Every row includes the
three system Booleans, its class, exact failed zero-forcing components, exact
failed block components, and the non-bipartite coloring/support decision.
The archive contains 255 rows in total.

| index | first impossible K6 seed | `Z0` choices |
|---:|---|---:|
| 652900 | `[1,3,10,11,14,17]` | 64 |
| 2301548 | `[2,4,10,12,14,18]` | 32 |
| 2842523 | `[0,3,4,14,15,17]` | 127 |
| 3289061 | `[4,11,12,14,15,18]` | 32 |

The exact pass-set histograms are:

| index | systems passed by one `Z0` | count |
|---:|---|---:|
| 652900 | block + zero forcing | 32 |
|  | non-bipartite only | 22 |
|  | zero forcing + non-bipartite | 10 |
| 2301548 | block + non-bipartite | 6 |
|  | non-bipartite only | 25 |
|  | zero forcing + non-bipartite | 1 |
| 2842523 | block + zero forcing | 96 |
|  | non-bipartite only | 15 |
|  | zero forcing + non-bipartite | 16 |
| 3289061 | block + non-bipartite | 6 |
|  | non-bipartite only | 25 |
|  | zero forcing + non-bipartite | 1 |

No row passes all three systems.

The two new graphs have especially transparent certificates.  Each has five
eligible zero-factor vertices, hence all 32 subsets are checked, and the
non-bipartite system passes for every subset.  The block system passes six
subsets, the zero-forcing system passes one, and those witness families are
disjoint.  For `2301548`, the empty set passes block but fails zero forcing;
`Z0=[6,11]` passes zero forcing but fails block.  For `3289061`, the analogous
representatives are the empty set and `Z0=[7,16]`.

At the empty set, the failed component in each new graph has zero-forcing
side lower ranks `3` and `4`, requiring seven dimensions.  At the unique
zero-forcing-passing representative, both exact normal-inertia orientations
also require seven dimensions and the lightlike component has ten
orthonormal vectors.  The archive records the full labeled components and
all numerical fields needed to check these statements.

## Production result

The source-bound serial run used

```text
python3 d6_k6_all_same_z0.py --workers 1 \
  --output d6_k6_all_same_z0_report.json \
  --certificates d6_k6_all_same_z0_certificates.json
```

It completed in 7.07 seconds on macOS 14.5 arm64 with Python 3.11.15.  One
worker was used intentionally while the K7 production campaign owned the
other cores.

| quantity | exact value |
|---|---:|
| input graphs | 990 |
| graphs rejected / surviving | 4 / 986 |
| K6 seeds checked | 31,590 |
| impossible seeds | 4 |
| `Z0` choices considered and matchable | 32,143 |
| zero-forcing-passing / failing `Z0` | 32,044 / 99 |
| subsequently block-passing / failing `Z0` | 31,997 / 47 |
| subsequently non-bipartite-failing `Z0` | 411 |
| zero-forcing components checked | 219,896 |
| block components checked / failed | 219,663 / 47 |
| generic orientations checked | 439,326 |
| proper-subset-only orientation failures | 2 |
| block subsets checked | 1,348,101 |
| raw / allowed-mask-matchable ray colorings | 33,635 / 31,586 |
| joint actual-support searches / DFS nodes | 31,586 / 95,593 |

The known realizable 18-point construction passes all 32 required K6 seeds.

## Independent verification

`verify_d6_k6_all_same_z0.py` imports no production evaluator or K6
production engine.  It reuses the frozen independent reconstruction of the
block and non-bipartite systems from `verify_d6_k6_same_z0.py` and implements
ordinary zero forcing afresh:

- it enumerates every initial black subset in increasing cardinality;
- it computes closure by simultaneous valid color changes, rather than the
  production solver's sequential deterministic forces;
- it reconstructs every induced side graph and the full inequality `(ZF)`;
- it compares every failed labeled zero-forcing component in the certificate
  archive, not only graph-level answers.

The serial command

```text
python3 verify_d6_k6_all_same_z0.py \
  --output d6_k6_all_same_z0_verification.json
```

recomputed all 990 graphs in 11.61 seconds and returned `PASS`.  It matched
ten core decision fields for every graph, all four first impossible seeds,
every archived three-system Boolean, every exact failed zero-forcing
component record, and the complete four-index rejection set.  Its independent
18-point positive control also passes.

## Controls and hashes

Run the controls with

```text
python3 -m unittest -v test_d6_k6_all_same_z0.py
python3 -m py_compile \
  d6_k6_all_same_z0.py verify_d6_k6_all_same_z0.py \
  test_d6_k6_all_same_z0.py
```

The seven tests cover frozen totals and hashes, every certificate histogram,
the two new disjoint-witness certificates, independent zero-forcing controls
on empty, isolated, path, cycle, and edgeless graphs, verifier import
independence, the frozen verification artifact, and the realizable positive
control.  All seven pass.

```text
d6_k6_all_same_z0.py
  d15fe376993393a5106e5d01f8c0b2316f43945c5fd9ee01b1ebdabfb4fe2b94
d6_k6_all_same_z0_report.json
  8fa616455d930a1ac3c2255ecb123d8095f7731780f6b6d0bec1ba71c173edeb
d6_k6_all_same_z0_certificates.json
  bce608e936cc6826962ecae8ea26baf3b9162072a861f568fc15570350e3f401
verify_d6_k6_all_same_z0.py
  9ad84985cffbc6f53562730f1aa536a6d7046709c21af0dfd1bfbab907a32162
d6_k6_all_same_z0_verification.json
  a503ad16025bdc5c252fc9fc4a871a0319293dd73ef850f535aad5abe77c9203
test_d6_k6_all_same_z0.py
  7e2916b0d5c2e9e3635c83bac0e8e8d196fa16cda40aecbc9a1eed83a40e2a4a
```
