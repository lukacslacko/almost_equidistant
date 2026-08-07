# Exact K6 same-Z0 conjunction

This layer removes two more graphs from the exact 990-index dimension-six
K6-only residue.  The new fact is a quantifier coupling: a geometric
realization has one actual zero-Lorentz-factor set `Z0`, so its bipartite
normal-coordinate constraints and its non-bipartite light-ray constraints
must hold for the same choice.

The complete production rejection set is

```text
652900, 2842523.
```

An independent checker rebuilt both constraint systems and all 990 graph
decisions separately and returned `PASS`.

## Why the conjunction is necessary

Fix a required unit `K6` seed.  Every outside vertex has a Lorentz factor
`ell_x` and a normal defect vector whose allowed coordinate mask is `D_x`.
In a putative realization define the actual set

```text
Z0 = {x : ell_x = 0}.
```

This is one geometrically determined set, not a fresh existential choice for
each lemma.  It has all of the following necessary properties.

1. Every `x in Z0` has `|D_x|>=3`, the orthonormal `Z0` vectors fit in six
   dimensions, and their allowed masks admit a matching.  Thus `Z0` is among
   the exhaustively enumerated eligible subsets of size at most six.
2. Every bipartite component of `L-Z0` satisfies one of its two exact generic
   normal-coordinate block-support systems or the separate lightlike system.
   The proof of those inequalities is in `d6_k6_normal_block_support.md`.
3. Every non-bipartite component of `L-Z0` lies on one of the two Lorentz
   light rays.  Some symmetry-reduced two-ray coloring must pass the two
   allowed-mask Hall tests and the joint actual-support CSP.  Each `Z0`
   support is one shared choice in both light-ray bins.

Write `B(Z0)` for property 2 and `N(Z0)` for property 3.  Checking the two
filters separately establishes only

```text
(exists Z0: B(Z0)) and (exists Z0: N(Z0)).
```

A realization requires the strictly stronger statement

```text
exists Z0: B(Z0) and N(Z0).                         (same-Z0)
```

The production search enumerates every eligible `Z0` and rejects a K6 seed
only when `(same-Z0)` fails.  One impossible required K6 seed rejects the
whole graph.  This preserves the existential quantifier exactly; no heuristic
or optimization failure is used.

Candidate nonedges remain unconstrained and may also have unit distance.
Every `D_x` is only an upper bound on actual support, and any allowed
coordinate may be zero.  The block inequalities enlarge actual supports to
these upper masks, while the joint-support search explicitly ranges over
subsets of them.  Neither subsystem requires a candidate nonedge to be
genuinely non-unit.

## Exact certificates for the two graphs

For corpus graph `652900`, the first impossible seed is

```text
[1,3,10,11,14,17].
```

It has six eligible zero-factor vertices, hence all `2^6=64` subsets are
enumerated.  The full cross-classification is

| status of one `Z0` | count |
|---|---:|
| block system passes, non-bipartite system fails | 32 |
| non-bipartite system passes, block system fails | 32 |
| both pass | 0 |
| neither passes | 0 |

The empty set is a representative of the first class; `Z0=[4]` represents
the second.

For corpus graph `2842523`, the first impossible seed is

```text
[0,3,4,14,15,17].
```

It has seven eligible vertices.  The full seven-element set cannot be `Z0`
in six dimensions, so all `sum_(k=0)^6 binom(7,k)=127` permitted subsets are
enumerated.  Their cross-classification is

| status of one `Z0` | count |
|---|---:|
| block system passes, non-bipartite system fails | 96 |
| non-bipartite system passes, block system fails | 31 |
| both pass | 0 |
| neither passes | 0 |

The empty set represents the first class; `Z0=[6,18]` represents the second.

Thus each constituent filter has witnesses when considered separately, but
no common witness exists.  The report stores every `Z0`, both Boolean system
decisions, its class, and the detailed failure used by the main search.

For both rejected seeds, every non-bipartite failure already has zero
allowed-mask-matchable light-ray colorings.  Consequently no negative claim
for these certificates depends on exhaustion of the deeper actual-support
DFS.  Every block failure records all failed components, both generic exact
inertia orientations, their first Hall-subset violation, and the lightlike
dimension/matching failure.

The newer proper-subset block inequality is present in this layer, but it is
not what eliminates these two seeds: their certificate counters contain zero
generic cases that pass the old ambient inequality and fail only the new
proper-subset inequality.  The two rejections are genuinely due to preserving
the common `Z0` quantifier across previously separate exact systems.

## Production result

The serial source-bound run completed on the exact manifest-pinned K6 residue
with the following totals:

| quantity | exact value |
|---|---:|
| input graphs | 990 |
| graphs rejected | 2 |
| graphs surviving | 988 |
| K6 seeds checked | 31,602 |
| impossible seeds | 2 |
| `Z0` choices considered and matchable | 32,093 |
| `Z0` choices passing the block system | 32,011 |
| `Z0` choices failing the block system | 82 |
| block-passing `Z0` choices failing the non-bipartite system | 411 |
| bipartite components checked | 219,831 |
| bipartite components failed | 82 |
| generic orientations checked | 439,662 |
| proper-subset-only generic failures | 2 |
| block subsets checked | 1,358,802 |
| raw light-ray colorings | 33,649 |
| allowed-mask-matchable colorings | 31,600 |
| joint actual-support searches | 31,600 |
| actual-support DFS nodes | 95,607 |

The command was

```text
python3 d6_k6_same_z0.py --workers 1 \
  --output d6_k6_same_z0_report.json
```

On macOS 14.5 arm64 with Python 3.11.15 it took 9.03 seconds wall time.
It was intentionally run serially while an unrelated full tetrad campaign
used the CPU.  The report records full per-graph decisions, the two exhaustive
certificates, input-index and dependency hashes, command, platform, Python
version, and exact quantifier/optional-zero semantics.

The known realizable 18-point construction passes all 32 of its required K6
seeds.

## Independent verification

`verify_d6_k6_same_z0.py` imports neither the production evaluator nor any K6
engine.  It separately reconstructs:

- graph validation, K6 enumeration, defect masks, eligible `Z0` subsets, and
  Lorentz components;
- both generic block-support orientations and the lightlike alternative;
- exact inertia through SymPy integer characteristic polynomials and Sturm
  root counts, instead of production rational congruence elimination;
- every allowed-mask matching by direct Hall-subset enumeration, instead of
  production augmenting-path search;
- every light-ray coloring and the actual-support domains, pairwise
  intersection constraints, shared `Z0` assignments, and structural-rank
  tests through a separately written DFS.

The serial command

```text
python3 verify_d6_k6_same_z0.py \
  --output d6_k6_same_z0_verification.json
```

recomputed all 990 graphs in 10.90 seconds and returned `PASS`.  It matched
the complete two-index rejection set, eight core decision fields for every
graph, both first impossible seeds, every per-`Z0` failure partition, and
every full cross-classification.  Its independent 18-point positive control
also passes all 32 K6 seeds.

## Tests and artifact hashes

Run the controls with

```text
python3 -m unittest -v test_d6_k6_same_z0.py
python3 -m py_compile \
  d6_k6_same_z0.py verify_d6_k6_same_z0.py test_d6_k6_same_z0.py
```

The five unit controls check frozen totals and source boundaries, every
certificate partition, the complementary existential witnesses, absence of
production imports in the verifier, the frozen independent verification, and
the realizable 18-point control.  All five pass.

```text
d6_k6_same_z0.py
  6d6dd948a245d5d8245f6ee10ce58c725aba29b9d3a9b325c8be11b59d519b5f
d6_k6_same_z0_report.json
  cea11ab95f24aa5fe79f06cc4d1a719b775ac5647551e0fb737b85835c443ca1
verify_d6_k6_same_z0.py
  4db05ff06910a11395e4bbdd2fd13161d6c19f58f81cc99851b8c37558d7000f
d6_k6_same_z0_verification.json
  9d1f06c41827193f3a324bf2d08f450cf816567cbae727e21a6258afd68a9fbf
test_d6_k6_same_z0.py
  67c404ca7f3ca40ef3400a459126afcba1f9ffa38b8f7828117702b0ec1afa3b
```
