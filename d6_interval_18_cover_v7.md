# Certified interval search on the v7 dense 18-deletion cover

## Scope

`run_d6_interval_18_cover_v7.py` turns the final 263-parent v7 residue into a
small, deterministic obstruction campaign.  For every parent it retains all
18-vertex deletions with the maximum number of required unit edges.  It then
uses greedy set cover, maximizing the number of newly covered parents and
breaking every tie by the **lexicographically largest** canonical class ID.
The frozen result is:

```text
v7 parents                                      263
selected canonical 18-deletion classes          179
selected classes from K7 parents                  12
selected classes from K6-only parents             167
greedy gains: 5x1, 4x6, 3x10, 2x42, 1x120       263
non-positive actionable backups                      2
production classes (179 base + 2 backups)           181
```

The selector is importable as `build_dense_deletion_cover(parents,
deletions)`.  It returns the 179 engine records and a complete provenance
object containing every parent maximum, every eligible dense class, the full
greedy chronology, exact coverage lists, stable hashes, and known-positive
annotations.  This can also be reused by heuristic screens; those screens do
not inherit any certified conclusion from the interval campaign.

The sole selected known-positive class is the only base-cover class for
parents `K6_only:364827` and `K6_only:3335955`.  A sound interval engine can
never reject those parents through that class.  Production therefore retains
the class as a full-parameter positive control and adds a non-positive dense
backup for each blocked parent.  The frozen deterministic rule maximizes the
backup's number of active v7 parents and then takes the lexicographically
largest class ID.  It chooses:

```text
parent 364827 -> ordinal 126,  class u18-f6d7f4c977df6677ea2ceb5682a0686e7f6255ebad113646aaca3eb8ce46cb43
parent 3335955 -> ordinal 5673, class u18-d3ce53599fbf37b60b6b5bd7f30468280eca97105df1044f097518392fd936d0
```

Both backups have 104 required edges and only two-circle placement orders.
At the production four-order cap ordinal 126 has four generated orders and
ordinal 5673 has two.  Neither is marked standard-18-compatible.  Ordinal
5673 also densely covers current parents `3335944` and `3335952`.

The selected class-index and class-ID roots are respectively

```text
ef4ac7e73ea924d3d9c6382f590dc73cdf2e65b357e86ddc69719372f717d6fe
c3b74892ca75a407fa3b4455546acab5899c851493bdf9cfcd34cb3dcc4827ee
```

The sorted 181-class production-index root is

```text
8c624c969b517bb5663dc32bdf577ab3addbe8a108458b2803443c17f294d8e6
```

## Mathematical semantics

Only edges in a selected deletion record are required to have distance one.
A missing edge is unconstrained and may also have distance one.  Therefore,
if a deletion class is certified `KILLED`, every v7 parent in that class's
coverage list is impossible.  This is ordinary subgraph monotonicity and
does not use the deleted vertex or any candidate nonedge as a distance
constraint.

The search status meanings are unchanged:

- `KILLED`: one placement order was certified empty on every theta slice;
- `ABORT`: a kernel call hit its node cap;
- `UNRESOLVED`: survivor boxes remained or no eligible order existed;
- `INFRA_ERROR`: an execution failure occurred.

Only `KILLED` is a rejection.  The other three statuses make no mathematical
claim.

## Circle orders and the one explicit exception

The 18-vertex inputs are not zero-circle cases.  The wrapper deliberately
sets `zero_circle_only=False` and `include_bulk_order=True`, so position zero
is not accidentally dropped.  With four generated orders per class, the
profile is:

```text
minimum one-circle classes                        123
minimum two-circle classes                         57
no eligible placement order                         1
all generated orders: 1 circle / 2 circles     249 / 120
```

The sole no-order class is deletion-manifest ordinal `7259`, class ID
`u18-d47fb2b0a2f8d5f03f8bc4d37ca035690fd7a33a2c1d60bbc6be8fd912a83416`.
It is the selected dense deletion covering K6-only parent `3950509`.  The
runner will checkpoint it as `UNRESOLVED` with zero kernel calls; it is never
silently omitted.

## Positive controls

Deletion ordinal `441`, class ID
`u18-81842b7afe827273043aa730eaec7a1fda421f033ad25e8b4a126f4fddd456f6`,
is a known standard-18-compatible support.  It covers current parents
`364827` and `3335955`, so those parents cannot be rejected through this
particular deletion.  The wrapper marks it as a positive control, exercises
it before the campaign, and raises an error if a production result ever calls
it `KILLED`.

The preflight controls also include:

- the flexible circle-family realization (`K7` plus eleven generic circle
  points), checked on every production slice;
- a `K8` negative control, checked as `KILLED` on every slice;
- the full standard 18-point unit graph;
- both independently exactified nonstandard 18-point unit graphs.

The exact positive controls are shallow kernel checks and must return an
overall status other than `KILLED`.  Their report and verification artifacts
are hash-pinned.

## Production discipline

The final v7 manifest and independent verification JSON are pinned by their
immutable SHA-256 roots.  Run the tests and commit this wrapper before
launching.  At runtime the wrapper checks that itself, the shared runner, the
Python driver, the C kernel, and the interval library all equal their blobs in
`HEAD` and in the newest commit touching any of those five sources.  That
source commit, the branch name, and the five hashes form the stable source
identity.  Unrelated later commits and unrelated staged/untracked files are
deliberately excluded from the campaign hash, so they cannot break checkpoint
resume.  Any executable-source dirt or mixed executable-source boundary fails
closed.  Nondeterministic control wall times are likewise excluded from the
hashed control witnesses; their statuses, orders, calls, nodes, and survivor
counters remain pinned.  Production must run on branch `codex/dimension6`.

Focused tests do not invoke the interval kernel:

```text
python3 -m unittest -v test_run_d6_interval_18_cover_v7.py
```

After the source commit exists, a restartable full pilot can be launched as
follows (record the exact outer command):

```text
caffeinate -dimsu python3 run_d6_interval_18_cover_v7.py \
  --workers 11 --orders 4 --cap 100000 --slices 24 \
  --checkpoint-dir .runs/d6_interval_18_cover_v7_cap100000 \
  --output .runs/d6_interval_18_cover_v7_cap100000_report.json \
  --decisions .runs/d6_interval_18_cover_v7_cap100000_decisions.tsv \
  --checkpoint-every 1 --progress-every 1 \
  --outer-launch-command '<the exact full command above>'
```

The checkpoint campaign hash includes the complete parent boundary, deletion
boundary, greedy cover, placement-order profile, search parameters, source
hashes, binary hash, controls, environment, Git commit, and launch command.
Changing any mathematical parameter requires a new checkpoint directory.
