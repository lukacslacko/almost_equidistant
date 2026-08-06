# Bounded exact K7 support-capacity pilot

## Scope and outcome

This is a deliberately bounded, one-core pilot of the exact support-type
multiplicity theorem in `d6_k7_exact_support_multiplicity.md`.  It does not
modify the frozen rank, labeled-support propagation, sparse-value, strict-H,
degree-one Schur, degree-four tetrad, or pattern-954 packages.

The target population is the independently verified 258-graph complement of
the tetrad/pattern-954 union.  The pilot selected 32 graphs deterministically,
stratified equally across the four cover-mechanism signatures present in the
sample.  The result was:

```text
selected graphs                                      32
K7 seeds                                             49
eligible covers                                  17,570
current tetrad-passing covers                         74
labeled Z support families checked                   788
families failing singleton propagation               744
families reaching the capacity CSP                    44
capacity-feasible families                            44
capacity-infeasible or unresolved families             0
pre-capacity-passing covers                           41
capacity-passing covers                               41
joint pre-capacity-failing covers                     33
exact joint graph rejections                           7
incremental support-capacity graph rejections          0
unresolved graphs or families                          0
```

Thus the new multiplicity-cap relaxation itself had **zero hits** on the
bounded sample.  Every one of the 41 cover witnesses checked independently
uses the propagated masks themselves as actual supports; those identity
choices already respect all type capacities.  In accordance with the pilot
stop rule, no full 258-graph support-capacity run was launched.

The combination of already certified layers nevertheless exposed seven new
exact graph rejections:

```text
3958422  3945531  3958119  3333850
2592565  2591766  3957913
```

For at least one K7 seed of each graph, all covers surviving the strict-H,
degree-one, and tetrad certificates fail the frozen labeled-support/sparse-
value quantifiers.  These are cross-layer intersection rejections, not new
multiplicity-cap rejections.  The independent checker exhaustively replayed
all ten recorded failing seeds and every one of their remaining labeled
families.

## Exact capacity relaxation

Fix a K7 seed.  For a labeled zero-factor family, let its exact supports be
`Z_1,...,Z_r`, and let singleton propagation give masks `E_x` for the
nonzero-factor vertices.  The new finite relaxation assigns

```text
empty != S_x subseteq E_x
```

and imposes only:

1. `S_x` and `S_y` intersect when `xy` is a required edge;
2. the number of all fixed and assigned points of exact type `S` is at most
   `|S|`.

Candidate nonedges impose no condition.  Omitting all further coordinate
equations makes the CSP weaker than geometric realizability, so an exhaustive
infeasibility is a sound obstruction and a feasible assignment is only a
nonclaim.

The kernel uses capacity-aware backtracking.  At every node, an exact
bipartite b-matching test assigns the remaining vertices to support types with
their residual capacities while forgetting edges between unassigned
vertices.  This is a Hall relaxation and therefore a sound prune.  The exact
edge constraints are restored by backtracking.  A node cap can return only
`UNRESOLVED`; it can never generate a rejection.  On this pilot every search
found a witness without backtracking: 523 nodes, 523 flow checks, zero flow
prunes, maximum depth 12.

## Independent boundary

`verify_d6_k7_support_capacity_pilot.py` imports neither the pilot runner nor
the new capacity kernel.  It:

- reconstructs the deterministic 32-graph sample from the frozen 258-graph
  manifest;
- binds all inherited decision/certificate archives by SHA-256;
- reconstructs current-passing covers;
- uses `verify_d6_k7_support_full.py` and
  `verify_d6_k7_sparse_value_full.py`, the previously independent support and
  sparse-value kernels, to exhaust all ten failing seeds;
- directly checks each of the 41 recorded capacity witnesses from support
  containment, required-edge intersection, and exact type counts.

It reports `PASS`: seven joint exact rejections, 25 bounded survivors, zero
capacity-only hits, and zero unresolved cases.  A survivor remains only a
failure of this incomplete screen to reject.

The realizable 18-point lower-bound graph has no K7 seed, so this
K7-conditioned control is explicitly `NOT_APPLICABLE_NO_K7`.  Synthetic
positive and negative cover-level controls, plus 120 deterministic random
instances cross-checked against full Cartesian enumeration, test the capacity
kernel.  Candidate-nonedge optionality has its own independent control.

## Reproduction

```sh
/Users/lukacs/claude/opengauss/venv/bin/python3 \
  -m unittest -v \
  test_d6_k7_support_capacity.py \
  test_verify_d6_k7_support_capacity_pilot.py

/Users/lukacs/claude/opengauss/venv/bin/python3 \
  run_d6_k7_support_capacity_pilot.py \
  --limit 32 \
  --node-limit 200000 \
  --output d6_k7_support_capacity_pilot_report.json

/Users/lukacs/claude/opengauss/venv/bin/python3 \
  verify_d6_k7_support_capacity_pilot.py \
  --report-sha256 <printed report SHA-256> \
  --output d6_k7_support_capacity_pilot_verification.json
```

All mathematical decisions are finite integer/bit-mask computations.  No
floating-point rank test, numerical optimizer, transcendental function, or
candidate-nonedge distance assumption enters this layer.  Elapsed times are
the only floating-point values in the report.
