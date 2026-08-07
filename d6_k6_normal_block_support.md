# K6 normal-coordinate block-support refinement

This note records an exact support-Hall strengthening of the K6
normal-inertia rule.  It was profiled on the exact 990-index K6-only class in
`d6_current_residue_manifest.json`.  It gives genuine extra pruning for two
individual Lorentz-orientation cases, but no whole component, zero-factor
choice, K6 seed, or graph is eliminated.  Thus its standalone marginal graph
coverage is zero.

## Exact block-subset inequality

Fix a required unit `K6` seed and write `u_x in R^6` for the six normal defect
coordinates of an outside vertex.  Let `D_x` be the seed positions at which
`u_x` is *allowed* to be nonzero.  These are upper bounds on actual support:
a candidate nonedge may still have unit distance, and an allowed coordinate
may still vanish.

Fix a possible actual zero-Lorentz-factor set `Z0` and a connected bipartite
component `C=A union B` of the disjoint-defect Lorentz graph after deleting
`Z0`.  The derivation in `d6_k6_normal_coordinates.md` gives two generic
Lorentz orientations.  Put

```text
H_A = I + Adj(G[A]),        H_B = I + Adj(G[B]).
```

If the direction on `A` has positive Lorentz norm, exact inertia gives

```text
r_A = |A| - n_-(H_A),       r_B = |B| - n_+(H_B).
```

If the direction on `A` has negative norm, the signs swap:

```text
r_A = |A| - n_+(H_A),       r_B = |B| - n_-(H_B).
```

In either orientation, define the Euclidean subspace blocks

```text
U_A = span{u_x : x in A},
U_B = span{u_x : x in B},
U_z = span{u_z}             for each z in Z0.
```

The K6 pair identity and the Lorentz-component structure make these blocks
mutually orthogonal.  For `z in Z0`, the diagonal identity gives
`||u_z||=1`, so each `U_z` has dimension one.  The inertia argument gives
`dim U_A >= r_A` and `dim U_B >= r_B`.

Every block also lies in a known coordinate subspace:

```text
U_A subset span{e_i : i in union_(x in A) D_x},
U_B subset span{e_i : i in union_(x in B) D_x},
U_z subset span{e_i : i in D_z}.
```

Consequently, for every subset `J` of these mutually orthogonal blocks,

```text
sum_(j in J) rank_lower(j)
    <= | union_(j in J) allowed_coordinates(j) |.       (block-support)
```

Indeed, the selected blocks have an orthogonal direct sum, whose dimension is
at least the left side, and that direct sum lies in the coordinate subspace
on the right.  This proves the inequality without a generic-position or
nonzero-coordinate assumption.  There are at most `|Z0|+2 <= 8` blocks, so
the implementation checks all nonempty subsets directly for both Lorentz
orientations.

The old normal-inertia rule is the full-block inequality after replacing its
coordinate-union capacity by the ambient bound six.  Therefore
`(block-support)` can be strictly stronger on a proper block subset.

## Lightlike alternative and rejection semantics

If the component direction is lightlike, all vectors indexed by `C union Z0`
are individually orthonormal.  The implementation retains this as a separate
alternative and requires both

```text
|C| + |Z0| <= 6
```

and a full matching from those vectors into their allowed coordinate masks.
It does not reject a component unless both generic orientations and this
lightlike alternative fail.

For candidate nonedges, `alpha(G)<=2` implies that their two allowed defect
masks are disjoint: a shared seed nonneighbor would form an independent
triple.  Thus all orthogonality statements used above remain valid whether a
pair is a required edge or a candidate nonedge.  The public graph evaluator
checks `alpha(G)<=2` explicitly.

The mask semantics are deliberately one-sided.  Replacing an actual support
by its allowed mask only enlarges every right side in `(block-support)`, so a
failure with allowed masks is still a sound impossibility certificate.
Enlarging allowed masks can never create a rejection.  No candidate nonedge
is required to be genuinely non-unit.

Non-bipartite Lorentz components are ignored by this prototype.  This loses
possible coupled coverage but cannot create an invalid rejection.  For one
K6 seed, every eligible support-matchable `Z0` is tried.  A graph would be
rejected only if at least one of its required K6 seeds had no surviving
choice.  A reported survivor is only a filter non-rejection, not a
realizability claim.

All matrix inertias are computed by exact rational symmetric congruence in
`d6_k6_normal_inertia.py`; all remaining decisions are integer bitset and
matching calculations.  Floating point is used only for elapsed-time
reporting.

## Full 990-graph result

The exact run completed as follows:

| quantity | exact value |
|---|---:|
| input K6-only graphs | 990 |
| graphs rejected | 0 |
| graphs surviving | 990 |
| K6 seeds checked | 31,654 |
| `Z0` choices considered | 31,673 |
| support-matchable `Z0` choices | 31,673 |
| failing `Z0` choices | 19 |
| bipartite components checked | 218,377 |
| generic orientation cases checked | 436,754 |
| old ambient-passing generic cases | 436,666 |
| all failing generic cases | 90 |
| ambient-passing cases newly failed by support | 2 |
| block subsets checked | 1,310,900 |
| passing lightlike cases | 211,499 |
| old-passing components newly failed | 0 |
| old-passing `Z0` choices newly failed | 0 |

The two newly failing visited orientation cases occur once each in corpus
graphs `2552186` and `2870953`.  This count concerns the deterministic choices
visited before each seed's first passing `Z0`; it is not a count of all latent
orientation failures after that short circuit.

A concrete strict-strengthening witness is graph `2552186`, seed
`(2,5,7,10,13,17)`, `Z0` empty, with component

```text
A = [0,9,12,14,16,18],      B = [4,8,15].
```

In the `A_positive` orientation the rank lower bounds are `(5,1)`, so the old
ambient total is exactly six and passes.  But `U_A` is allowed in only seed
coordinates `{2,3,4,5}`.  The one-block subset therefore violates `5 <= 4`.
The `A_negative` orientation has ranks `(3,2)` and passes all block subsets,
so the component and graph are not rejected.

## Controls and reproduction

`test_d6_k6_normal_block_support.py` checks:

- a synthetic proper-subset failure whose full ambient inequality passes;
- monotonicity under enlargement of allowed masks and an optional-zero case;
- the eight-block limit;
- the fixed graph-`2552186` strict-strengthening witness;
- the `alpha(G)<=2` precondition;
- all 32 K6 seeds in the known realizable 18-point construction;
- every frozen aggregate and the two graph-level orientation counters in the
  full report.

The known realizable 18-point positive control passes: no seed is rejected
and no generic orientation is newly failed by this support rule.

Reproduce the report and controls with

```text
python3 d6_k6_normal_block_support.py \
  --output d6_k6_normal_block_support_report.json
python3 -m unittest -v test_d6_k6_normal_block_support.py
python3 -m py_compile \
  d6_k6_normal_block_support.py test_d6_k6_normal_block_support.py
```

The frozen report records the command, Python version, platform, wall time,
input-index hash, dependency hashes, and complete per-graph decisions.

## Artifact hashes

The final SHA-256 values are recorded after the source-bound report is
regenerated:

```text
d6_k6_normal_block_support.py
  14faafc6b2f9c3c34c954bb9cbaa21e9a7f646cfb76fcdfb30dfe7af742cf9cf
d6_k6_normal_block_support_report.json
  e8fb18ba12dc5a211fface7b504e1f1119f4bfa5c1db0bb8360953787504c92d
test_d6_k6_normal_block_support.py
  ca68f6126550c18559f1ea6e27811921748ac24714e28956dc5cb2bcc2dbc976
```
