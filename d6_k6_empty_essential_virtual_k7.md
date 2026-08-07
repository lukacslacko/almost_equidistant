# K6 empty-essential profile and exact virtual-K7 closure

## Result

This package adds two exact layers after the verified 805-graph K6-only
residue of `d6_k6_empty_support_budget_report.json`.

1. The empty-essential profiler checks all 25,354 required K6 seeds.  It finds
   an explicit all-nonempty Hall witness for 25,243 seeds.  The remaining 111
   seeds, in 49 graphs, have no all-nonempty branch in the current necessary
   relaxation and therefore require one or two empty K6 defect vectors.
2. The virtual-K7 layer promotes every possible empty vertex at those 111
   seeds.  All 291 promoted branches are eliminated by the exact generic K7
   cover/rank quantifier.  Consequently all 111 seeds are impossible and all
   49 graphs are rejected.

The ordered K6-only residue is reduced from 805 to **756** graphs, with hash

```text
2cfedbc83f6ff01b6e386440fb7ea066b274d22ffbc52371c964e99d17cdcb90.
```

The existing K7 class has 155 graphs and is untouched by this layer.  Thus,
relative to the 960-graph current boundary that combines 155 K7 graphs and
805 K6-only graphs, the new combined boundary has **911 graphs**.

## Empty-essential quantifier

Fix a required unit K6 `Q` and a common normal-zero set `Z0`.  In every
bipartite Lorentz component, the inherited exact Hall system considers:

- both generic sign orientations;
- the all-nonempty choice on the negative side;
- each allowed single empty negative-side singleton; and
- the lightlike alternative.

Nonbipartite components remain ignored conservatively.  A finite DP combines
the component choices while retaining at most two empty vertices.  Two state
spaces are kept:

- the inherited cardinality-only budget; and
- the sharper opposite-apex budget, which also requires a two-empty pair to
  be a candidate nonedge.

The second condition is exact in the branch.  In centered K6 coordinates an
empty defect has `u=0` and is one of

```text
r_+ = +sqrt(7/12)e,       r_- = -sqrt(7/12)e.
```

Two distinct empty vertices are therefore `r_+` and `r_-`, whose squared
distance is `7/3`, not 1.  This use of a candidate nonedge is only permissive:
no candidate nonedge is globally required to be non-unit.

The profiler stops a routine seed as soon as it records an all-nonempty
witness.  If no such witness exists, it exhausts every `Z0`, orientation,
component alternative, and global state.  The exact partition is

| seed class | count |
|---|---:|
| explicit all-nonempty witness | 25,243 |
| empty-essential | 111 |
| impossible from the opposite-pair rule alone | 0 |
| total K6 seeds | 25,354 |

The 111 exceptional seeds occur in 49 graphs.  Their state profile has:

- 291 distinct `(graph, K6 seed, apex)` promotion keys;
- 174 candidate opposite-apex pair keys;
- 104 seeds allowing both one- and two-empty relaxed states;
- 7 seeds allowing only one-empty states;
- 92 seeds with one particular empty vertex in every relaxed state; and
- 19 seeds with no individually mandatory vertex.

The known realizable 18-point configuration supplies a positive control: all
32 of its K6 seeds retain an explicit all-nonempty witness.

## Virtual-K7 branch theorem

If an outside point `x` has empty actual defect relative to `Q`, then every
pair `xq`, `q in Q`, has distance 1.  Adding only missing `x-Q` pairs as
required edges is therefore sound in this branch, and `Q union {x}` is an
actual regular unit K7.  No arbitrary virtual point is introduced; `x` is an
existing point of the candidate configuration.

For each of the 291 possible apices the production layer reconstructs the
augmented 19-vertex graph and exhausts the generic K7 zero-factor cover
quantifier:

- the zero-factor set is an eligible vertex cover of the exact `L` graph;
- orthogonal dimension gives size at most 7;
- cover sizes 4--7 fail the proved direct K7 dimension caps; and
- every cover of size at most 3 is checked by the exact support, subspace,
  component-inertia, positive-definite clique, perpendicular-degree,
  rational basis-kernel, and saturating-mask rules.

The cover audit is

| quantity | count |
|---|---:|
| promotions | 291 |
| promotions with no eligible cover | 259 |
| eligible covers over all promotions | 2,386 |
| direct-cap covers of size 4--7 | 2,237 |
| exactly checked covers of size 0--3 | 149 |
| passing covers | 0 |

The eligible-cover size histogram is

```text
size 1:   2       size 4: 353       size 7: 535
size 2:  24       size 5: 628
size 3: 123       size 6: 721
```

There are no size-zero eligible covers in these 291 target branches.  Failure
counts below overlap because a small cover can have several independently
valid obstructions:

```text
zero-factor support   54       PD clique              145
subspace K            10       perpendicular degree   119
component B          111       basis kernel            149
                               saturating mask          149
```

Every relaxed state at every empty-essential seed contains at least one of
these eliminated apices.  Hence all 111 seed conjunctions are impossible.
The double-apex test is unnecessary here: an impossible singleton promotion
already eliminates every pair state containing that vertex.

## Independent checking and controls

`verify_d6_k6_empty_essential_profile.py` imports neither the production
profiler nor the parent empty-support producer.  It independently rebuilds
the Hall alternatives and state DP and exactly compares all 25,354 seed
classifications.

`verify_d6_k6_empty_essential_virtual_k7.py` imports neither production
layer.  It independently:

1. reconstructs the ordered 805 input;
2. extracts the branch list from the independently verified profile;
3. adds each apex star;
4. reconstructs the K7 defect system;
5. brute-enumerates all eligible covers, including strict supercovers;
6. replays and compares all 149 small-cover certificates; and
7. recomputes the 111 seed and 49 graph conjunctions.

The generic K7 leaf kernel is shared and hash-pinned.  Its complete 16-test
suite passes.  The conjunction additionally has:

- a realizable regular K7 positive control, whose empty cover survives;
- a K8 negative control;
- an explicit all-nonempty-state preservation control;
- an opposite-apex required-edge exclusion control;
- an augmentation test showing that only the apex star is added;
- an independent strict-supercover enumeration control; and
- explicit report-hash tamper rejection.

All 9 package controls and the 25 inherited K7/K6 kernel controls pass.

## Commands and artifact hashes

Production and verification commands:

```text
python3 d6_k6_empty_essential_profile.py \
  --workers 11 \
  --output d6_k6_empty_essential_profile.json

python3 verify_d6_k6_empty_essential_profile.py \
  --report d6_k6_empty_essential_profile.json \
  --expected-report-sha256 \
    205c0877982c8de8ad6e5a685e33eff9ee05236c9907e60d488a593792d90781 \
  --output d6_k6_empty_essential_profile_verification.json \
  --workers 11

python3 d6_k6_empty_essential_virtual_k7.py \
  --workers 11 \
  --output d6_k6_empty_essential_virtual_k7.json

python3 verify_d6_k6_empty_essential_virtual_k7.py \
  --report d6_k6_empty_essential_virtual_k7.json \
  --expected-report-sha256 \
    63207da634a56e8fa46f18f87e1ed8b86007c62b51c49c5584eb80403f2b50ba \
  --output d6_k6_empty_essential_virtual_k7_verification.json \
  --workers 11

python3 -m unittest -v \
  test_d6_k6_empty_essential_virtual_k7.py \
  test_d6_k7_rank_reference.py \
  test_d6_k6_empty_support_budget.py
```

SHA-256 boundaries:

```text
dab7603b11946635a29c3c37daf44a895db342b2cdddb17c3699f87c182f507a  d6_k6_empty_essential_profile.py
205c0877982c8de8ad6e5a685e33eff9ee05236c9907e60d488a593792d90781  d6_k6_empty_essential_profile.json
8eeea41c5ebd04125ad2a33539757229472b2709a8771f56cc6af394c8cf3c2e  verify_d6_k6_empty_essential_profile.py
c8f73e8a9bed6d7a102ab168eb3842312a916b23dbd577cceb6d089c8832f719  d6_k6_empty_essential_profile_verification.json

dbfd44cef783be6849255504518697279d3dfc8dabf1aab72b221697f596a6c5  d6_k6_empty_essential_virtual_k7.py
63207da634a56e8fa46f18f87e1ed8b86007c62b51c49c5584eb80403f2b50ba  d6_k6_empty_essential_virtual_k7.json
83c04d7b6e035f63c9e1d38867dc4be2594397709b318cb8379a4ba42b443608  verify_d6_k6_empty_essential_virtual_k7.py
d64c36856e1288c8f64ea7918bc76640bde260b9b731e881f9c277c1a1e43485  d6_k6_empty_essential_virtual_k7_verification.json
e94e2fd92ad03a921b26440827b8365cb4740d9a58399379860f7106fec37db0  d6_k7_rank_reference.py
```

All mathematical decisions use exact combinatorics and rational arithmetic.
There is no floating-point trust assumption in these two layers.  A survivor
would mean only filter non-rejection; in the reported target quantifier there
are no surviving promotion branches.
