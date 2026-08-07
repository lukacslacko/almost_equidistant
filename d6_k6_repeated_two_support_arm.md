# K6 repeated two-support arm collision

This incremental exact layer starts from the frozen 634 survivors of
`d6_k6_singleton_fan.py`.  It adds a point-reconstruction consequence after
a tight Hall coordinate deletion.  The complete discovery census rejects 9
more graphs and leaves 625.  The earlier tight-Hall and singleton-fan results
remain separate frozen milestones.

## Forbidden reduced-support motif

Fix one active block on one side of a generic bipartite K6 Lorentz component.
All Lorentz factors in the block are nonzero multiples of one non-lightlike
direction.  Suppose a tight Hall equality leaves three distinct points with
reduced allowed supports

```text
x: {i},
y: {i,j},
v: {i,j},
```

and `xy` and `xv` are required unit edges.  The distance between `y` and `v`
is not used.

The masks above are upper supports, not assertions that every listed
coordinate is nonzero.  The active-block and empty-singleton branching in the
parent kernel guarantee that `u(x)` is nonzero.  Write

```text
r     = z/c,
R     = r^2,
kappa = (1-R)/6,
w     = u/c.
```

The singleton anchor has `c != 0`: the alternative `c=0` gives
`u=-e_i,z=0`, a zero-Lorentz-factor seed collision, not an active generic
point.  Since all nonzero Lorentz factors on this side are scalar multiples
of the anchor factor, every arm also has nonzero `c` and the same signed ratio
`r`.  Dividing the K6 diagonal and pair identities by `c^2` and `c_x c_y`
gives

```text
||w||^2 = kappa + (1-sum w)^2,              (1)
w(a) dot w(b) = kappa                       (2)
```

whenever `ab` is a required same-side edge.

For `w(x)=tau e_i`, the exact singleton formula gives

```text
tau = (7-R)/12.
```

The value `R=7` makes the nonempty singleton coordinate zero; equivalently,
an anchor--arm edge would require `0=kappa=-1`.  Exclude that already
impossible case.  For an arm write

```text
w(y) = alpha e_i + a e_j.
```

Here `a=0` is still allowed.  Equation (2) for the anchor--arm edge fixes

```text
alpha = kappa/tau = 2(1-R)/(7-R).
```

Put

```text
beta  = 1-alpha = (R+5)/(7-R),
delta = kappa-alpha^2.
```

The arm diagonal (1) is

```text
a^2 = delta + (beta-a)^2.
```

The squares in `a` cancel, and `beta != 0` for real `R>=0`, `R!=7`.  Thus

```text
a = (delta+beta^2)/(2 beta) = (R+5)/12.      (3)
```

This resolves the optional-zero branch rather than assuming it away:
`a>0`, so the second coordinate is forced to be active.  More importantly,
equation (3) uniquely determines the whole normalized defect vector of an
arm.

## Finite reconstruction and collision

Normalization loses no scale information because `c=1+sum u`.  For either
arm,

```text
1 = c(1-sum w),

1-sum w = beta-a = (R+5)^2 / (12(7-R)) != 0,

c = 12(7-R)/(R+5)^2.                         (4)
```

Hence no finite-denominator exceptional branch remains.  Equations (3)--(4)
give identical `w`, `c`, `u=cw`, and `z=rc` for `y` and `v`.  K6 defect
coordinates and `z` reconstruct the Euclidean point uniquely, so `y=v`,
contrary to required point distinctness.

The two arms' mutual distance was never used.  In the level-19 candidate
corpus it is in fact forced to be an edge by `alpha(G)<=2`: both arms are
candidate-nonadjacent to the seed vertex `q_i`.  The implementation does not
need that corpus-specific shortcut and checks only the two anchor edges.

## Scope and optional-zero audit

The rule is installed only inside the inherited generic-orientation support
check.

- `R=1` is lightlike and remains in the separate inherited lightlike branch.
  The new rule does not prune that branch.
- Vertices in `Z0` have zero Lorentz factor and are not members of the active
  generic block to which the rule is applied.
- The one explicitly chosen empty negative singleton is omitted from the
  active block before support propagation.
- The anchor cannot use the empty subset of `{i}` because its block is
  explicitly nonempty.  Each arm's `i` coordinate is fixed by its required
  edge to the anchor; its optional `j` coordinate is forced nonzero by (3).
- If either arm retains a third allowed coordinate, its diagonal has two free
  arm coordinates and the cancellation no longer reconstructs a unique
  point.  The code does not fire.
- No candidate nonedge is assigned a distance, including the arm pair.

## Frozen parent and complete quantifier

The producer pins the singleton-fan source, report, deterministic compressed
certificate archive, and independent verification.  It also verifies the
raw SHA-256 of the decompressed parent archive and selects exactly the
parent's ordered 634 survivors.  The inherited solver still quantifies over:

- every required K6 seed;
- every eligible zero-factor set `Z0`;
- both generic sign orientations and each eligible empty negative singleton;
- all inherited exact inertia, zero-forcing, Hall, and tight-equality choices;
- the separate lightlike alternative;
- the global two-empty-point condition.

The independent verifier imports neither the producer nor the discovery
probe.  It builds on the frozen independent SymPy/Sturm, simultaneous
zero-forcing, and 64-container tight-Hall checker and adds a separately
written ordered-triple version of the collision predicate.

## Discovery census

The complete 634-residue pilot preserves all 32 K6 seeds of the known
realizable 18-point control and rejects exactly

```text
552851, 1335614, 1488658, 1904567, 1938285,
2639869, 3043357, 3165458, 3386554
```

The incremental rejection hash is

```text
9daf9767df740b642483c3046019043d9114fa9cabf010a126c51fe293234639
```

and the ordered 625-graph residue hash is

```text
04d9475217e294efed0c6eac34a7fad88616c90ad7d37a2525776a713e6dd0e1
```

Together with the frozen 107 tight-Hall and 15 singleton-fan rejections, the
131-element cumulative rejection set has stable hash

```text
c7fe84fc5b95683602fc58b9d89d7ea9427c11f6b04b193eb7da5e3d880504ac
```

Requiring the arm--arm edge gives the same nine hits on this corpus, as
predicted by `alpha(G)<=2`, but the frozen theorem uses only the two necessary
anchor edges.

## Commands and controls

Commit the producer, verifier, tests, and this note before running the
source-bound campaign.  Then run

```text
python3 -m unittest -v test_d6_k6_repeated_two_support_arm.py

python3 d6_k6_repeated_two_support_arm.py --workers 11

python3 verify_d6_k6_repeated_two_support_arm.py --workers 11
```

The focused controls cover the exact rational reconstruction, `R=7` and
lightlike scope, a synthetic motif whose arm pair is deliberately a candidate
nonedge, third-coordinate and different-second-coordinate no-fire cases, a
missing-anchor-edge no-fire case, fixed marginal witness `552851`, the frozen
634-graph input boundary, all 32 positive-control K6 seeds, and independent
import/transcription boundaries.

The official source-bound run used commit
`8ba9e4621eb246e8af7399dd18ba2039762f4ac4`.  Production with eight workers
completed in 10.4474 seconds; the independent eight-worker replay completed
in 34.1685 seconds.  It recomputed all 634 graph decisions and all 32 K6
seeds of the realizable 18-point control and returned `PASS`.

```text
d6_k6_repeated_two_support_arm_report.json
  ea8d9438d062b92857a1057b950e76b8ae08bbe1bf97b4327130fd4058171ebe

d6_k6_repeated_two_support_arm_certificates.json
  fbaebe286acce8b0c3603f89f05a5316010b51dceabfeface7e73ff55ce90ed2

d6_k6_repeated_two_support_arm_certificates.json.gz
  604ec86b943ec909d177acc1fd36ef623c3e3bff67b4a8a5e5ffa601fc202af5

d6_k6_repeated_two_support_arm_verification.json
  3abf48bab092671a9dff03c2dc47c0063a75515cd219e1d5a4bfa9a5fa797b5b
```

Git stores the deterministic `gzip -n` archive rather than the 3.7 MiB raw
JSON.  Decompress it to the raw filename before rerunning the independent
verifier; the report and verifier bind the raw hash above.

Survival is filter non-rejection, not realization.  This layer does not by
itself settle the K6 class or dimension six.
