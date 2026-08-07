# Strongest K6 same-`Z0` conjunction

This exact incremental layer starts from the 625 survivors of the frozen
repeated-two-support-arm campaign.  It reuses the already proved same-`Z0`
principle, but replaces its old bipartite block test by the complete current
tight-Hall, singleton-fan, and repeated-arm kernel.  The source-bound discovery
target is the two-graph marginal rejection

```text
1442098, 1495055.
```

The ordered rejection-list stable hash is

```text
908325e53918685e3732180266c3c8408c545d5c77859e8ffd7c997ea926ec36
```

and the expected ordered 623-graph residue hash is

```text
b31aeac00b91d0d843ea51909c64f2c4ca3da49e5d2d33792312a45aa9f25e79
```

These counts remain a discovery result until the frozen producer and
independent verifier have completed their official full-corpus runs.

## Necessary common quantifier

Fix a required unit `K6` seed.  Every outside point has a geometrically
determined Lorentz factor.  In any putative realization there is one actual
set

```text
Z0 = {x : the Lorentz factor of x is zero}.
```

It is not sound to choose one `Z0` for a bipartite-component argument and a
different `Z0` for a nonbipartite-component argument.  For a fixed enumerated
`Z0`, let

```text
B_tight(Z0) = the complete current bipartite-component relaxation passes;
N_light(Z0) = the nonbipartite two-light-ray actual-support relaxation passes.
```

A realization necessarily supplies

```text
exists Z0: B_tight(Z0) and N_light(Z0).             (1)
```

The producer enumerates every eligible `Z0` of size at most six and rejects a
required K6 seed only when every choice fails (1).  One impossible required
seed rejects the graph.

### The `B_tight` subsystem

For every bipartite component of the Lorentz graph after deleting `Z0`, the
current frozen kernel exhausts:

- both generic sign orientations and the separate lightlike alternative;
- every inherited exact inertia, zero-forcing, PSD-Z, arbitrary-subset Hall,
  and tight-Hall branch;
- the possible empty negative singleton and the global two-empty-point rule;
- the exact singleton-fan and repeated-two-support-arm collision rules after
  every tight coordinate deletion.

Nonbipartite components are deliberately ignored by this subsystem.  Thus a
realization supplies a witness to `B_tight(Z0)`, but a passing branch is only a
necessary relaxation.

### The `N_light` subsystem

Every nonbipartite component after deleting `Z0` lies wholly on one of the two
Lorentz light rays.  The frozen two-light-ray system exhausts the
symmetry-reduced assignments of these components to the rays.  In each ray
bin it requires structural Hall rank and searches actual defect supports.

If `D_x` is the allowed defect mask, the search ranges over actual subsets

```text
empty != S_x subseteq D_x
```

for a nonzero light-ray point.  For `x in Z0`, it ranges over all subsets of
`D_x` of size at least three; that one support is shared between the two ray
bins.  Within either bin, two orthogonal nonzero defect vectors cannot have
supports intersecting in exactly one coordinate, and every assigned family
must pass the exact Hall condition.  Intersections of size at least two are
accepted as possibly cancellable, so this too is a relaxation.

This explicitly preserves optional-zero semantics.  An allowed coordinate
may vanish, and a candidate nonedge is never required to be genuinely
non-unit.

### A further sound relaxation

For one fixed `Z0`, production asks only that *some* witness pass `B_tight`
and that *some* (possibly different) actual-support witness pass `N_light`.
An actual realization would give compatible common data, so allowing the two
witnesses to differ can create survivors but cannot create a false rejection.
The new result therefore uses only the common geometric `Z0`, not an
unproved coupling of the remaining existential choices.

## Exhaustive certificates for the two discovery rejections

The first impossible seed of graph `1442098` is

```text
[0,2,12,13,15,16]    mask 110597.
```

It is the tenth required K6 seed in deterministic order.  Its seven eligible
zero-factor vertices give all
`sum_(k=0)^6 binom(7,k)=127` permitted choices.  Their complete partition is

| status of one `Z0` | count |
|---|---:|
| `B_tight` passes, `N_light` fails | 2 |
| `N_light` passes, `B_tight` fails | 46 |
| both fail | 79 |
| both pass | 0 |

Equivalently, 125 choices fail the tight bipartite subsystem first and the
remaining two fail the nonbipartite actual-support subsystem.

The first impossible seed of graph `1495055` is

```text
[3,4,11,14,16,18]    mask 346136.
```

It is the thirtieth required K6 seed.  Its 127 choices partition as

| status of one `Z0` | count |
|---|---:|
| `B_tight` passes, `N_light` fails | 1 |
| `N_light` passes, `B_tight` fails | 68 |
| both fail | 58 |
| both pass | 0 |

Thus in both graphs each constituent subsystem is feasible by itself, while
no common `Z0` passes both.  The certificate archive stores every one of the
254 choices, the first failure reason, and the independent cross-classification.

## Discovery census

The complete single-process discovery replay considered the frozen 625
graphs and preserved all 32 required K6 seeds of the known realizable
18-point construction.  Before graph-level early exit it recorded:

| quantity | count |
|---|---:|
| enumerated `Z0` choices | 19,403 |
| choices failing `B_tight` | 378 |
| choices passing `B_tight` | 19,025 |
| bipartite components processed | 136,958 |
| nonbipartite components deferred to `N_light` | 12,116 |
| `B_tight`-passing choices failing `N_light` | 120 |
| choices passing the conjunction | 18,905 |

The discovery-only source and output hashes are

```text
probe_d6_k6_tight_same_z0_conjunction.py
  43f2b367c9e3929d0c3c00524ee21c86ac8c15c2f3e6349c4770e2e5a684e757

/private/tmp/d6_k6_tight_same_z0_conjunction_pilot.json
  10f0039a22be485880fe401a08a36ea3f1ade42475d508864c55cb61669779a1
```

## Independent verification design

`verify_d6_k6_tight_same_z0.py` imports neither the producer nor the
discovery probe.  It combines two previously frozen independent
implementations on the identical independently enumerated `Z0`:

- the repeated-arm/tight-Hall verifier, using SymPy characteristic
  polynomials and Sturm inertia, a separately implemented simultaneous
  zero-forcing solver, brute-force Hall subsets, and 64-container tight-set
  enumeration;
- the older same-`Z0` verifier's separately transcribed two-light-ray
  coloring and actual-support DFS.

The checker reconstructs the 625-input boundary from the repeated-arm
report, recomputes every graph decision, and compares every rejected seed's
full `Z0` failure list and cross-classification.  Focused controls already
show exact agreement on both proposed rejections and all 254 of their `Z0`
rows.  The official verifier must additionally replay the entire 625-graph
corpus and the independent 32-seed positive control.

## Frozen parent boundary and commands

The producer pins the repeated-arm source, report, deterministic compressed
certificate archive, independent verification, and the old independently
verified same-`Z0` theorem and support implementations.  It also checks the
SHA-256 of the decompressed repeated-arm certificate payload.  The frozen
input is exactly the parent's ordered 625 survivors, with hash

```text
04d9475217e294efed0c6eac34a7fad88616c90ad7d37a2525776a713e6dd0e1
```

After committing the producer, verifier, tests, and this note as a source
boundary, run

```text
python3 -m unittest -v test_d6_k6_tight_same_z0.py

python3 d6_k6_tight_same_z0.py --workers 8

python3 verify_d6_k6_tight_same_z0.py --workers 8
```

The five focused controls cover the frozen 625 boundary, both fixed
rejections and their complete `Z0` partitions, row-for-row agreement with the
independent implementation, optional-zero actual-support domains, import
independence, and both versions of the 32-seed positive control.

A survivor is only a graph not rejected by these necessary conditions.  This
incremental layer does not settle the remaining K6 class or dimension six.
