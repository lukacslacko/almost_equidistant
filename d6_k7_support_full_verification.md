# Independent full verification of K7 support propagation

## Result

`verify_d6_k7_support_full.py` independently replayed all 1,536 exact
support-propagation rejections in the archived 17,764-graph K7 rank residue.
For each rejected graph it used the recorded failing K7 seed, found every
cover that passes the frozen baseline rank reference, and exhausted every
**labeled** actual support family of every such cover.  All 35,035 families
failed one of the two elementary propagation conditions:

- 22,807 produced an empty propagated support mask;
- 12,228 left a required edge whose endpoint support masks were disjoint.

Thus no family survived.  The replay covered 2,354 baseline-passing covers
and recorded 260,690 singleton-coordinate deletions.  The complete
machine-readable report, including per-graph and per-cover counts, is
`d6_k7_support_full_verification_report.json`.

This is exact integer and bit-mask computation.  It uses no floating-point
comparisons, numerical rank decisions, random choices, or transcendental
functions.

## Independent decision boundary

The verifier deliberately does **not** import
`d6_k7_support_propagation.py`.  It separately implements:

1. seven-coordinate support domains;
2. labeled Cartesian enumeration (including both labelings when two
   vertices have equal domains);
3. pairwise singleton-intersection exclusion;
4. Hall feasibility through a fresh augmenting-path matching routine;
5. singleton-overlap fixed-point propagation; and
6. the empty-mask and disjoint-required-edge decisions.

It uses the frozen `d6_k7_rank_reference.py` only to classify which covers
pass all pre-existing baseline rank filters.  The checker independently
constructs the K7 defect masks, the auxiliary graph `L`, all eligible covers
of size at most three, and the induced required-edge graph, and cross-checks
the seed and cover construction against the frozen reference.

The mathematical implication behind propagation is elementary.  Once an
actual zero-factor support `S_z` is fixed, orthogonality to that vector rules
out a nonzero vector whose current support intersects `S_z` in exactly one
coordinate: a one-term dot product of nonzero factors cannot vanish.  That
coordinate can therefore be deleted, and the argument may be iterated.  An
empty support is impossible.  Likewise, two nonzero-factor vertices joined
by a required unit edge cannot have disjoint actual supports.  The verifier
enumerates a safe superset of the algebraically realizable support families,
so eliminating all enumerated families is a sound rejection.

Checking one failing required K7 seed suffices: every realization of the
whole graph would also realize that required simplex.  Covers rejected by the
frozen baseline remain the responsibility of that separately verified rank
layer; this package exhausts every cover not already rejected there.

## Archive and population bindings

Before the replay, the verifier checked all of the following:

- exactly 17,764 archive rows, in the same order and with the same unique
  indices as `.runs/d6_k7_rank_survivors.json`;
- exactly 1,536 refined rejections and no baseline rejection;
- every scalar, rejection list, and column total in the production report;
- all rejected rows have a nonzero recorded failing seed, and no survivor
  does;
- the production report's input, frozen-reference, production-source,
  decision-log, and checkpoint hashes; and
- the generated input's binding to the independently verified full K7 rank
  decision archive.

Principal SHA-256 values are:

```text
independent verification report
69e07f23f615dd0ce12f04b6c36929d8512e24e979ff69bf9593b172156c3787

independent verifier source
c6b9dce47618bd9a4247eb5a9f093682265447729e1868c4dc0677a24f51c628

compressed support decision archive
3c35b228c4768b06f881dc9046329196b683d390f4763f83816c6a3464ecf938

uncompressed support decisions
117e057a79135fd48cac316ee4f1e08e11a2be84e00d9d10393bc17d4a518a31

production support report
5138d207543967b393baae9b14b47b49a276273cbc4468ee843856891917691e

17,764-graph input
7a0a350142930e4e13830ea44c4fd217f25aeb07bbd6286e51e316f72712c4f9

frozen rank reference
e94e2fd92ad03a921b26440827b8365cb4740d9a58399379860f7106fec37db0

production support source
578103f70d2906fbb3d449b7b7b5c1b6c86e5ae17afb0c341602a543bea69f0c
```

## Reproduction commands

The seven independent-kernel controls were run as:

```sh
D6_ALLOW_BACKGROUND_TEST_ONLY=1 \
  /Users/lukacs/claude/opengauss/venv/bin/python3 \
  -m unittest -v test_verify_d6_k7_support_full.py
```

The complete 11-worker replay was run as:

```sh
caffeinate -dimsu \
  /Users/lukacs/claude/opengauss/venv/bin/python3 \
  verify_d6_k7_support_full.py \
  --workers 11 \
  --progress-every 64 \
  --report d6_k7_support_full_verification_report.json
```

It completed with status `PASS` in 21.561094042 seconds of internal wall time
on the machine recorded in the JSON report.  The command is deterministic;
`--limit N --workers 1` provides a bounded control without changing the full
archive-binding checks.

## Trust scope

The remaining trust assumptions are Python's exact integer, container, gzip,
CSV, JSON, and SHA-256 implementations; the correctness of the documented
support lemmas; and the frozen baseline rank reference for classification of
pre-existing cover failures.  The independently replayed propagation kernel
does not depend on IEEE-754 or `libm` behavior.
