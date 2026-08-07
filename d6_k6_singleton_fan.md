# K6 singleton--two-arm algebra

This incremental exact layer starts from the frozen 649 survivors of
`d6_k6_tight_hall_support.py`.  It adds one within-span algebraic consequence
after a tight Hall coordinate deletion.  The discovery census rejects 15
more graphs, leaving 634.  The parent 107 rejections remain a separate frozen
milestone.  The producer, verifier and tests were frozen before the official
run at commit

```text
c751ce5201596dfc6a20938ddd786f3423b0b2d1
```

## Forbidden reduced-support fan

Fix one generic side of a bipartite K6 Lorentz component.  Every Lorentz
factor on the side is a nonzero multiple of one non-lightlike direction.
Suppose a tight Hall equality has forced three distinct points to have
reduced allowed supports

```text
x: {i},
y: {i,j},
w: {i,k},

j != k,
```

and all three pairs are required unit edges.  The frozen support propagation
already guarantees that every point in the active block has a nonzero defect
vector.

Write the common projective Lorentz ratio as

```text
r = z/c,    R = r^2.
```

For the singleton point `x`, put `u_x=t e_i`.  The K6 diagonal identity gives

```text
c_x = 12/(R+5),
t   = (7-R)/(R+5).
```

The value `R=7` would give `t=0`, contrary to the forced nonempty singleton.
Because `x` has nonzero `c_x`, every nonzero Lorentz multiple on this side
also has nonzero `c`.

The required edge `xy` and the K6 pair identity imply

```text
t u_yi = ((1-R)/6)c_x c_y,

u_yi = alpha c_y,
alpha = 2(1-R)/(7-R).
```

The same formula holds for `w`.  If either arm used only coordinate `i`, it
would be the unique non-seed singleton point on this Lorentz line and would
collide with `x`, by the frozen exact singleton-root argument.  Thus the
second coordinate of each arm is genuinely nonzero.  Since `j != k`, the
only common coordinate of the two arms is `i`.  Their required edge now gives

```text
alpha^2 c_y c_w = ((1-R)/6)c_y c_w.
```

Canceling the nonzero `c` factors and clearing the nonzero denominator yields

```text
(R-1)(R+5)^2 = 0.
```

The generic direction is non-lightlike, so `R != 1`.  The remaining root is
`R=-5`, impossible because `R=r^2` for a real Lorentz ratio.  Hence the fan
cannot occur.

All three distances used by this proof are required unit edges.  No candidate
nonedge is assigned a distance.  The code asserts the fan only when the tight
deletion leaves the exact upper masks `{i}`, `{i,j}`, `{i,k}`.  If an arm has
a third allowed coordinate, optional-zero semantics prevent the rule from
firing.

## Complete quantifier and frozen parent

The new producer pins the source, report, certificate archive and independent
verification of the 107-hit parent.  It selects exactly the parent's ordered
649 survivors, then replays the complete inherited quantifier:

- every required K6 seed;
- every eligible zero-Lorentz-factor set `Z0`;
- both generic sign orientations and every possible empty negative
  singleton;
- all inherited Hall choices and tight equalities;
- the separate lightlike branch;
- the global two-empty-point condition.

Only the generic tight-support consequence is extended by the fan lemma.
The independent verifier imports neither this producer nor the discovery
probe.  It extends the prior independent SymPy/Sturm and simultaneous
zero-forcing checker with a separately written index-pair fan predicate.

## Discovery result

The complete 649-residue discovery pilot preserves all 32 K6 seeds of the
known realizable 18-point control and rejects exactly these 15 graphs:

```text
58458, 1345949, 1612547, 2261989, 2293477,
2297966, 2297969, 3339806, 3341086, 3357805,
3358042, 3401800, 3668326, 3942891, 3967473
```

Together with the frozen 107, the 122-element cumulative rejection list has
stable hash

```text
53ac1854ed59878df7e637b2797358f9165ad34b5799e69a8a98680af529390f
```

The ordered 634-graph residue has stable hash

```text
bafeef8b2f16419bbbfd1dbd8756081be78d14d0dc4187a1d954738a7edb9d88
```

The incremental production report records only the 15 new rejections, so its
incremental rejection hash differs from the cumulative hash above.

## Commands and controls

```text
python3 -m unittest -v test_d6_k6_singleton_fan.py

python3 d6_k6_singleton_fan.py --workers 11

python3 verify_d6_k6_singleton_fan.py --workers 11
```

The controls include the cleared exact `R` equation, a synthetic forbidden
fan, two optional-zero/no-false-fan examples, fixed marginal witness `58458`,
the frozen 649-graph input boundary, all 32 positive-control K6 seeds, and the
independent import/transcription boundary.

The official 11-worker production completed in 15.1348 seconds.  The
independent 11-worker verification completed in 46.5262 seconds, and all six
focused tests passed in 11.231 seconds.  Artifact hashes are:

```text
d6_k6_singleton_fan_report.json
  bb11b82e8783108e45ec97a8e92957037bf219dd475b71fc18faee50318a4e2d

d6_k6_singleton_fan_certificates.json
  b2fbf39d443679fee52d4984f3a54e4302acf416760087ed75e0953eb384d7bf

d6_k6_singleton_fan_certificates.json.gz
  fcce6c4e7758e4ffef107342a81d9e541d696d378592aad92701611280e5f467

d6_k6_singleton_fan_verification.json
  77011cf280b250aff59d43008af10525c936f853172040905c0cfc5324491820
```

Git stores the deterministic `gzip -n` archive rather than the 6.9 MiB raw
JSON.  Decompress it to the raw filename before rerunning the verifier; the
report and verifier check the raw hash above.

Survival is filter non-rejection, not realization.  This layer does not by
itself settle the K6 class or dimension six.
