# K6 tight-Hall equality support propagation

This exact layer starts from the 756 K6-only graphs in
`d6_current_residue_manifest_v5.json`.  Its source-bound production rejects
107 graphs and leaves 649.  The producer, independent verifier and focused
tests were frozen before the official run at commit

```text
f65178ed3f1e5ebfed0c91e3eabe91ac8b1d5b5c
```

## Tight equality lemma

Fix a K6 seed, a possible zero-Lorentz-factor set `Z0`, one bipartite
component of the disjoint-defect Lorentz graph, one generic sign orientation,
and (on the negative side) the possible identity of its one empty defect
vector.  The inherited Hall system supplies a finite family of mutually
orthogonal actual defect subspans.  For any selected subfamily `J`, write

```text
r_J = sum of its proved rank lower bounds,
C_J = union of its allowed seed-coordinate masks.
```

Every actual selected subspan is contained in the coordinate space
`R^{C_J}`.  If

```text
r_J = |C_J|,
```

then the direct orthogonal sum of the actual subspans has dimension at least
`|C_J|` and at most `|C_J|`.  It therefore equals `R^{C_J}`.  This remains
true even though every allowed support is only an upper bound and an actual
vector may use a strict subset: strict supports cannot lower the dimension
below the proved rank bounds.

Every other Hall span in the same Lorentz component is orthogonal to the
selected spans.  Its actual vectors consequently vanish on every coordinate
in `C_J`.  Production enumerates all tight unions by the inherited exact
64-state union DP.  The checker instead enumerates all 64 coordinate
containers and independently maximizes one subset choice in each original
span.  If that maximum equals the container size, Hall validity forces the
actual union to equal the container.

## Consequences after deleting a tight coordinate set

Let `F` be a connected required-unit block in one generic Lorentz side, and
replace each allowed defect mask by

```text
D'_x = D_x minus C_J.
```

Every vertex of `F` has a nonzero defect vector except for the one explicitly
branched empty negative singleton, whose whole singleton block is omitted.
The following failures are exact.

1. If `D'_x` is empty, the required nonzero Gram diagonal has no coordinate
   support.
2. If `xy` is a required edge of `F` and `D'_x` and `D'_y` are disjoint, its
   Gram entry cannot be nonzero.  On one generic Lorentz line,
   `u_x dot u_y = <ell_x,ell_y>_L/6`, and both Lorentz multipliers are
   nonzero, so this entry really is nonzero.
3. Two vertices in `F` cannot both be forced onto the same one-coordinate
   support.

For the last statement, put `u=a e_i` and let `r=z/c` be the common
projective Lorentz ratio on that side.  The K6 identities give `c=a+1` and

```text
(r^2-1)(a+1)^2 + 6a^2 - 6 = 0

  = (a+1)((r^2+5)a + (r^2-7)).
```

The root `a=-1` has `c=z=0`; it is the seed vertex `q_i`, and in any event
has zero Lorentz factor rather than a generic nonzero factor.  The only other
possible root is

```text
a = (7-r^2)/(r^2+5).
```

Thus there is at most one distinct non-seed point on the given seed
coordinate and projective Lorentz line.  Two such vertices collide.  This
argument never assigns a distance to a candidate nonedge.

## Complete quantifier

For every required K6 seed, production enumerates every eligible `Z0` of
size at most six whose allowed masks pass exact Hall matching.  For every
bipartite Lorentz component it checks:

- both generic Lorentz sign orientations;
- no empty defect, and every zero-capable negative singleton in turn as the
  unique empty defect on that side;
- every inherited arbitrary-subset/PSD-Z/inertia/zero-forcing Hall choice;
- every tight-union support consequence above;
- the separate lightlike alternative.

Generic choices are combined across components with the global fact that
there are only two empty K6 defect points.  If two occur, they are the
opposite common K6 neighbours at squared distance `7/3`, so their candidate
pair may not be a required unit edge.  Nonbipartite Lorentz components are
ignored conservatively.  A graph is rejected only when one required K6 seed
has no surviving relaxed branch.

Candidate nonedges remain unconstrained throughout.  In particular, an
allowed defect coordinate can be zero and no candidate nonedge is required
to be genuinely non-unit.

## Discovery census

The archived-witness pilot first tested the explicit witness stored for each
of all 23,601 K6 seeds:

| quantity | value |
|---|---:|
| input graphs / seeds | 756 / 23,601 |
| bipartite Lorentz components | 165,504 |
| generic orientation states newly failed | 1,210 |
| archived witnesses eliminated | 637 |
| graphs having at least one eliminated archived witness | 133 |

The subsequent exhaustive discovery run used every `Z0`, generic empty
state and lightlike alternative.  It checked 20,990 seeds before graph-level
early exit, considered 29,578 `Z0` choices, and found 107 impossible first
seeds, one in each of 107 graphs.  The ordered rejection-list stable hash is

```text
6d46d6db66236d294f3f5a61e1cb696492a2aad6b137bdfee044df260424e69d
```

The ordered 649-graph residue stable hash is

```text
6460baec9c378c1da3ba79f9ba6645e688158f2d9910bc6b88f4578c6e28e8d4
```

The raw exhaustive discovery report remains outside Git at
`/private/tmp/d6_k6_tight_support_exhaustive_full.json`; its byte SHA-256 is

```text
63d07cb6f756c24c73f23f499a61227df4028ac4b5282c250a3bf10e139b04bd
```

## Controls and commands

The focused tests cover tight equality with strict actual supports permitted,
loss of a required Gram entry, the exact singleton roots, fixed corpus
rejection `32033`, all 32 K6 seeds of the known realizable 18-point control,
the checker's independent container enumeration, and import independence.

```text
python3 -m unittest -v test_d6_k6_tight_hall_support.py

python3 d6_k6_tight_hall_support.py --workers 11

python3 verify_d6_k6_tight_hall_support.py --workers 11
```

The verifier imports neither the producer nor either discovery probe.  It
reconstructs inherited ranks through the prior SymPy/Sturm inertia and
simultaneous zero-forcing checker and compares every graph decision, first
impossible seed, and exhaustive rejected-seed `Z0` coverage.

The official 11-worker production completed in 12.1683 seconds and the
independent 11-worker verification completed in 31.7948 seconds.  The seven
focused controls passed in 8.170 seconds.  Artifact hashes are:

```text
d6_k6_tight_hall_support_report.json
  73b7fdd5b641d52a7eb3f2daaeea9d8e94057913e38f2bde06ef62691d01e864

d6_k6_tight_hall_support_certificates.json
  a3b29e28db706de5b1d63ba4f82e2b527a53d6b0ac479d098c3bb2cccb66198f

d6_k6_tight_hall_support_certificates.json.gz
  cbd35f5a9ba466df275474fc0fde384f0db8153e5f7d3e24f95b38e542c73594

d6_k6_tight_hall_support_verification.json
  5fa924e47f50b23a0fb6289d1beba7a6f445ea599d2078a8a8925fd23802523c
```

Git stores the deterministic `gzip -n` archive rather than the 65 MiB raw
JSON.  Decompress it to the raw filename above before rerunning the verifier;
the verifier and report both check the raw JSON hash.

Survival is only filter non-rejection, not a realization.  This layer alone
does not settle the K6 class or dimension six.
