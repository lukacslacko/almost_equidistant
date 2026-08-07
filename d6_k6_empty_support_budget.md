# K6 empty-defect support budget

This exact layer is designed for the independently checked 822-graph residue
after the K6 arbitrary-subset Hall layer.  Its production runner and separate
checker are source-bound and ready for review, but the full run must not begin
until this source boundary has been committed and audited.  The bounded pilot
finds one new marginal witness, corpus index `202556`; no theorem-level
coverage is claimed until the full production run and independent replay are
complete.

## The missing distinctness condition

Fix a required unit K6 and use the exact coordinates

```text
u_i(x) = ||x-q_i||^2-1,
c_x    = 1+sum_i u_i(x),
z_x^2  = c_x^2+6-6||u(x)||^2.
```

Allowed defect masks are upper bounds.  Every allowed coordinate may vanish.
If the actual defect vector is empty, then

```text
u(x)=0, c_x=1, z_x=+sqrt(7) or -sqrt(7).
```

These are the two common unit neighbours of the K6.  Distinctness therefore
permits at most two empty defect vectors globally and at most one of either
sign.

The earlier K6 support layer forces nonempty supports for zero Lorentz factors
and lightlike components.  In a generic bipartite Lorentz component it bounds
the rank of each connected required-edge block in the two side Gram matrices.
A singleton block on the negative Lorentz side has sound rank lower bound
zero, however, so the previous Hall DP can silently allow any number of these
vertices to have `u=0`.  This layer closes precisely that gap.

## Component theorem

In a generic component write the factors on the two sides as

```text
ell_x=lambda_x ell_A,
ell_y=lambda_y ell_B,
```

where the two Lorentz norms are nonzero and have opposite signs.

On the positive side,

```text
||u_x||^2 = 1+kappa lambda_x^2 > 0,
```

so no vertex has empty defect.  Equivalently, the positive singleton rank
lower bound is already one.

On the negative side a singleton block may have zero diagonal.  Because its
Gram matrix is positive semidefinite, a zero diagonal forces its whole row to
zero, consistently giving `u_x=0`.  A connected block of order at least two
cannot contain such a vertex: every vertex has an incident required edge and
the generic Gram entry on that edge is nonzero.

All negative-side Lorentz factors in this connected component share one
projective line.  If `u_x=0`, then `c_x=1` fixes the scalar `lambda_x`, and
hence fixes `z_x` and the reconstructed point.  Two empty singletons on this
same side would therefore collide.  Thus:

```text
at most one negative-side singleton block per generic component may be empty.
```

For each generic orientation the implementation enumerates exact alternatives:

1. no zero-rank singleton is empty, so each is replaced by a rank-one
   nonempty-vector choice in the existing Hall DP;
2. exactly one chosen negative singleton is empty, contributes zero rank and
   zero support, and all the others are forced rank one.

It retains the two generic orientations and the lightlike alternative.  The
minimum number of empty vertices needed by each bipartite Lorentz component is
therefore `0`, `1`, or impossible.  For one fixed zero-factor set `Z0`, the
sum of these component minima must be at most two.  Treating signs as freely
assignable to two different components is conservative: an actual solution
must fit within this relaxation.  Non-bipartite components are ignored by the
new budget, which can only weaken rejection and avoids double-counting their
older independently checked constraints.

Candidate nonedges are never forced nonunit.  The only zero/nonzero decision
added here is about the *actual* defect vector, and every branch is exhausted.

## Controls, production runner, and bounded pilot

`test_d6_k6_empty_support_budget.py` checks:

* a positive singleton has rank one while a negative singleton may have rank
  zero;
* one empty singleton on one negative projective line is permitted;
* two empty singletons on the same line are forbidden by collision;
* two different components can consume the two global empty points, while
  three cannot;
* the known realizable 18-point construction passes all 32 K6 seeds;
* fixed corpus witness `202556` is rejected at seed
  `[3,8,11,12,13,18]`.
* checkpoints reject changes to the production source, semantic configuration,
  dependency hashes, input hash, or ordered completed prefix;
* the checker has no import of the production kernel.

The bounded pilot selects five stable-hash records from each quartile of the
822 parent residue ordered by K6-seed count, then adds witness `202556` if it
is absent.  On the current 21-record selection, only `202556` is rejected.
All 247 eligible/matchable `Z0` rows in its failing seed contain a bipartite
component for which both generic orientations remain impossible even after
the unique allowed empty singleton, and the lightlike alternative also fails.
Thus the pilot hit is local-component marginal coverage, not an artifact of
adding empty costs from unrelated components.

The production package hard-pins the ordered 822-index hash

```text
cf94aac41eba28904ee6254f9c50996e73d38dcc15316fb6a411b29ff5cd05d8
```

and the parent source, report, certificate, and independent-verification
hashes.  An initial run uses at most 11 ordered worker processes.  Every 32
completed records it atomically replaces a checkpoint containing the entire
ordered prefix.  `--resume` accepts that checkpoint only when its source,
semantic configuration, dependencies, input hash, and prefix ordering all
match.  A worker exception or interruption atomically writes a distinct
`INFRA_ABORT` artifact; that status explicitly carries no rejection claim.

Rejection certificates omit bulky internal Hall-failure traces.  They retain
every exhausted `Z0`, every component reached before the exact terminal, all
three component alternatives, and the actual support masks.  The independent
checker reconstructs the alternatives with its own Hall DP and compares every
archived row exactly.  It also recomputes all 822 graphs in a fresh process
pool, validates the ordered post-filter residue, and separately checks all 32
K6 seeds of the realizable 18-point control.  Verification refuses to start
without an explicit SHA-256 of the production report.

Source-audit commands are

```text
python3 -m unittest -v test_d6_k6_empty_support_budget.py
python3 d6_k6_empty_support_budget.py --pilot \
  --output /private/tmp/d6_k6_empty_support_budget_pilot_report.json
```

After source review and a commit, the production and resume commands are

```text
python3 d6_k6_empty_support_budget.py --full --workers 11 \
  --checkpoint-every 32 \
  --output d6_k6_empty_support_budget_report.json \
  --certificates d6_k6_empty_support_budget_certificates.json \
  --checkpoint d6_k6_empty_support_budget_checkpoint.json \
  --abort d6_k6_empty_support_budget_abort.json

python3 d6_k6_empty_support_budget.py --full --workers 11 \
  --checkpoint-every 32 --resume \
  --output d6_k6_empty_support_budget_report.json \
  --certificates d6_k6_empty_support_budget_certificates.json \
  --checkpoint d6_k6_empty_support_budget_checkpoint.json \
  --abort d6_k6_empty_support_budget_abort.json
```

Once the report exists, record its SHA-256 explicitly and run

```text
python3 verify_d6_k6_empty_support_budget.py \
  --report d6_k6_empty_support_budget_report.json \
  --certificates d6_k6_empty_support_budget_certificates.json \
  --expected-report-sha256 <exact-64-hex-report-sha256> \
  --workers 11 \
  --output d6_k6_empty_support_budget_verification.json
```

The production runner refuses to run if its kernel or the checker differs
from the corresponding Git blob at `HEAD`.  Thus pilot execution is possible
before review, while theorem-level execution is not.

## Virtual-K7 follow-up

An empty-defect vertex is an actual common unit neighbour of the K6, so its
six missing seed incidences may be added as required unit edges and the seed
promoted to a virtual K7.  For a pair of empty vertices, the vertices must be
the two opposite common neighbours and their mutual pair must be a candidate
nonedge.  Branching over the empty set, a singleton, or an eligible pair can
therefore feed K7 exact kernels on augmented graphs.  This is promising but
not part of the present layer: each K7 kernel must first be audited to ensure
that it uses only required-edge semantics and does not assume that its input
is a minimal or canonically enumerated corpus graph.
