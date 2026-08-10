# Steering for the local Codex campaign

This file is addressed to the **local Codex 5.6 Sol High session running on the user's Mac**.
Read it completely after every pull that changes it.

## 1. How this collaboration works

There are two separate model sessions:

1. **You, the local Codex agent**, have the checkout, compiler, local files, processes, and the Mac CPU. You implement, test, run, checkpoint, and publish commits.
2. **A separate ChatGPT GPT-5.6 Pro session** acts as mathematical reviewer and research steering partner. It does not see your live terminal, uncommitted files, running processes, or private Codex conversation.

The `main` branch is the shared state and message bus. Work directly on
`main`; the former `codex/dimension6` branch has been retired. The user
mediates the loop:

1. ChatGPT reviews a pushed branch and may add or revise steering in a commit.
2. The user tells you to pull.
3. You pull with a fast-forward, read this file plus `STATUS.md`, the relevant plan, and the incoming diff.
4. You do the local implementation and computation.
5. You commit and push auditable results, update the handoff in `STATUS.md`, and tell the user the exact commit hash.
6. The user pings ChatGPT, which reviews that exact Git state and supplies the next steering step.

Do not assume that either model session remembers the other session's reasoning. Put every durable fact, command, proof obligation, result, and unresolved issue in the repository. User instructions override this file.

## 2. Current project state and immediate interpretation

The branch currently reports:

- a completed computer-assisted campaign for `f(5) = 16`;
- a verified level-19 candidate corpus for dimension 6;
- 3,055,474 of 3,971,787 dimension-6 level-19 candidates certified by the current search;
- 916,313 deferred cases for which simple node-cap grinding has measured cost on the order of weeks, not hours.

Preserve all completed work. Do not restart the campaign or delete provenance. At the same time, treat claims in prose as claims to be audited: a publishable theorem needs a reproducible manifest, complete coverage check, explicit trust assumptions, controls, and preferably an independent checker.

The present strategic recommendation is:

> **Do not spend the next phase merely raising caps on all 916,313 deferred graphs. First measure theory-driven filters and reusable obstructions that can eliminate many candidates at once.**

Already-running jobs may be allowed to finish only when they are checkpoint-safe and the user wants them to continue. Do not launch another exhaustive pass until the profiling milestone below has reported coverage and projected cost.

## 3. Mathematical semantics that must never be changed

For an abstract unit-distance graph `G`:

- every edge of `G` must have Euclidean distance exactly 1;
- a nonedge is unconstrained and may also have distance 1;
- all represented vertices must be distinct;
- no arbitrary epsilon-separation may be imposed;
- numerical optimization, a nonzero residual, failure to find a solution, or exhaustion of a heuristic budget is not a proof;
- an unresolved or aborted branch must be reported, never silently discarded;
- a candidate can be rejected only by a sound combinatorial, interval, exact algebraic, exact SDP/SOS, or independently checkable certificate.

When using a candidate's complement as an allowed non-unit support, remember that allowed entries may still be zero. Do not accidentally require every candidate nonedge to be genuinely non-unit.

## 4. Proof-engineering requirements

Preserve the `d=4` and `d=5` reproduction packages and keep their controls passing.

For every theorem-level campaign:

- separate search from certificate checking as far as practical;
- maintain an independent completeness verifier over the candidate corpus and result manifests;
- record corpus hashes, source hashes, git commit, compiler and flags, OS/architecture, Python version, worker count, exact command lines, timestamps, and random seeds;
- make long jobs restartable from atomic checkpoints;
- log `KILLED`, `SURVIVOR`, `UNRESOLVED`, `ABORT`, and infrastructure failures distinctly;
- validate both negative and positive controls on the exact final kernel build;
- retain small witness data sufficient to diagnose any rejected or unresolved case;
- do not describe search-management changes as new mathematical pruning rules unless they really are;
- state the IEEE-754 and transcendental-function assumptions explicitly. Do not silently inherit an 8-ulp `libm` assumption on a new platform without either evidence, conservative padding, or a correctly rounded alternative;
- use numerical least squares only for triage and discovery. Reconstruct and certify any apparent realization exactly or by rigorous intervals.

Before a long run, commit the code that will run. After a long run, do not modify logs in place; add a manifest or derived summary whose inputs are hashed.

## 5. Theory-first filters for dimension 6

The next implementation phase should evaluate the following filters on the complete 3,971,787-graph level-19 corpus. Keep each filter independently testable and report both individual and cumulative coverage. A filter may be used for rejection only after its proof has been documented.

### 5.1 Clique-link bounds

If `Q` is a unit `K_m` in `R^d`, all common unit neighbours of `Q` lie on a sphere in an affine space of dimension `d-m+1`, with radius squared

```text
(m + 1) / (2m).
```

The common-neighbour set remains almost equidistant. In `R^6`, the following safe bounds are available:

```text
unit K2 link: at most 16   (use f(5) = 16)
unit K3 link: at most 12   (use f(4) = 12)
unit K4 link: at most 10   (use f(3) = 10)
unit K5 link: at most 4    (circle argument below)
unit K6 link: at most 2
unit K7 link: 0
```

For `K5`, the normalized inner product for a unit chord is `1/6`. The corresponding angular step is an irrational multiple of pi: if it were rational, `2 cos(theta) = 1/3` would be a rational algebraic integer, which it is not. Hence the finite unit-distance graph on that circle has maximum degree 2 and no cycles, so it is a union of paths; independence number at most 2 then gives at most 4 vertices.

A sharper `K4` link bound of 9 may follow from the uniqueness and non-cosphericity of the 10-point `R^3` extremal configuration, but do not encode 9 until that theorem and deduction are cited or reproved in the repository.

Implement these as graph-only common-neighbour tests. Measure how many of the 916,313 deferred graphs they reject, and also run them over all already-certified graphs as a consistency check.

### 5.2 Low-rank and inertia obstruction

For points `p_1, ..., p_n in R^d`, define

```text
M_ii = -1
M_ij = ||p_i - p_j||^2 - 1,  i != j.
```

Then

```text
rank(M) <= d + 2,
M has at most one positive eigenvalue,
-1/2 C M C is positive semidefinite of rank <= d,
C = I - J/n.
```

In dimension 6 at `n=19`:

```text
rank(M) <= 8,
nullity(M) >= 11.
```

All required unit edges force corresponding off-diagonal entries to zero. Candidate nonedges are only *allowed* nonzero positions; their values may also be zero.

Develop exact or certifiable filters around:

- fixed-diagonal sparse minimum rank;
- vanishing 9-by-9 minors;
- star-complement / Schur-complement reconstruction from an 8-vertex core;
- exact inertia tests after algebraic reconstruction;
- small induced or non-induced support patterns whose fixed-zero affine matrix space cannot contain rank-at-most-8 matrices with the required Euclidean inertia.

Floating-point rank is only a screen. A rejection needs exact rational elimination, interval certification, a checkable polynomial certificate, or another rigorous argument. Standard zero-forcing/minimum-rank theorems often assume an exact nonzero graph pattern; use them carefully because this problem has optional zeros.

### 5.3 `K7` simplex defect-set reduction

Stratify the level-19 candidates by whether they contain a unit `K7`. For each inequivalent `K7`, place a centered regular unit simplex `q_1, ..., q_7` and define

```text
u_i(x) = ||x - q_i||^2 - 1,
s(x)   = sum_i u_i(x).
```

The full-dimensional simplex identities in `R^6` are

```text
x = -sum_i u_i(x) q_i,
7 sum_i u_i(x)^2 - (s(x) + 1)^2 = 7,
||x-y||^2 - 1 = ((s(x)+1)(s(y)+1))/7 - sum_i u_i(x)u_i(y).
```

Define the actual defect set

```text
S_x = { i : u_i(x) != 0 }.
```

For a candidate graph, `S_x` is a nonempty subset of the seed vertices that are *allowed* to be non-unit from `x`; it need not equal the whole candidate non-neighbour set.

Useful exact consequences:

- `S_x` is nonempty, or the seed plus `x` would be a unit `K8`;
- if `S_x` and `S_y` overlap, then `x` and `y` must be unit, by the almost-equidistant condition with the shared seed vertex;
- each seed vertex has at most 7 actual non-unit neighbours, so for 12 outside vertices, `sum_x |S_x| <= 49`;
- a one-defect point has the nonduplicate solution `u_i = 4/3` (the other root `-1` is the existing seed vertex);
- there is at most one one-defect point of each type, and at most two one-defect points in total, since three distinct one-defect types are pairwise non-unit;
- a two-defect point on coordinates `i,j` satisfies

```text
3 u_i^2 + 3 u_j^2 - u_i u_j - u_i - u_j - 4 = 0;
```

- a fixed two-defect type has multiplicity at most two.

First build a purely discrete SAT/CSP layer assigning possible defect sets subject to candidate adjacency, overlap-forces-unit, degree capacities, and multiplicity rules. Only surviving support assignments should proceed to the rational quadratic and bilinear equations.

Do not assume every level-19 candidate contains a `K7`; `R(3,7)` does not force that at 19 vertices. Keep the `K6`-only class separate. For `K6` seeds, an analogous coordinate system has one remaining normal coordinate and should be derived explicitly rather than guessed.

### 5.4 Reusable obstruction library

This is currently the most promising way to turn millions of expensive searches into fast graph containment checks.

A valid obstruction is a smaller unit-edge graph `F` that has no realization by distinct points in `R^6` using only the unit constraints recorded in `F`. If `F` is a subgraph of a candidate unit graph `G`, then `G` is impossible. The obstruction certificate must not rely on omitted nonedges having any prescribed distance.

Recommended pilot:

1. Select stratified samples from already-certified and deferred level-19 graphs.
2. From certified graphs, minimize the certificate by deleting vertices and then edges while rechecking non-realizability, to find small cores.
3. Canonically label and deduplicate the cores.
4. Independently re-certify every core with generous budgets or a second method.
5. Index subgraph containment with bitsets/canonical invariants and measure coverage over all 3,971,787 candidates.
6. Report the coverage curve by obstruction size and number of library entries.

Begin with known exact obstructions (`K8`, the appropriate dimension lemma such as `K_{1,3,3,3}`, and any existing certified small controls). Do not launch a massive core-mining campaign until a pilot estimates whether the library gives at least an order-of-magnitude reduction.

### 5.5 Simplex-overlap reflection propagation

Two unit `K7` simplices sharing a `K6` facet have their two apices related by reflection across that facet. Build the `K7`-clique overlap complex of each candidate:

- quotient seed cliques by automorphisms;
- propagate exact coordinates through `K6` overlaps;
- inspect cycles of reflections for forced collisions or wrong distances;
- combine with the link bounds for vertices reached through smaller shared faces.

This can produce compact graph-specific proofs and may also yield small reusable obstructions.

## 6. Required profiling milestone before another exhaustive campaign

Produce a committed, machine-readable report over the full dimension-6 level-19 corpus containing at least:

- total candidate count and corpus hash;
- deferred-set count and hash;
- clique-number distribution and number containing `K7`;
- rejection counts for each clique-link rule, individually and cumulatively;
- feasibility counts for the `K7` defect-set CSP, grouped by seed orbit;
- a pilot low-rank/star-complement screening summary, clearly separating exact rejections from numerical suspicions;
- obstruction-library pilot size, certified core sizes, containment coverage, and verification method;
- representative engine costs for the residue after each filter;
- a projected core-hour range based on measured samples, not optimism;
- exact commands and a resume procedure.

Call out any filter that rejected a known realizable positive control: that is a bug and blocks use of the filter.

Update `STATUS.md` with a concise current summary and a clearly marked handoff section. Avoid an indefinitely growing stream of predictions. Move detailed chronological diagnostics to a dated log file when needed.

## 7. Mac execution discipline

Use the Mac for local CPU work, but keep the computation subordinate to reproducibility:

- default to at most `logical CPU count - 1` workers, or fewer if memory bandwidth or thermals reduce throughput;
- benchmark worker counts rather than assuming more processes are faster;
- prevent sleep for an intentional long run with a visible, recorded mechanism such as `caffeinate`, but do not leave unidentified orphan jobs;
- write PID, command, source hash, and checkpoint directory into the run manifest;
- make every output unit atomic and resumable;
- keep generated binaries and transient logs out of Git unless they are intentional proof artifacts;
- compress large immutable logs and publish their hashes;
- after a crash or reboot, verify checkpoint integrity before resuming;
- do not run two exhaustive campaigns that oversubscribe the same CPU;
- do not promise periodic updates or an ETA that the agent cannot actually deliver. Report measured throughput, completed work, and uncertainty.

## 8. Handoff format for the remote review session

Before asking the user to ping ChatGPT, push a focused commit and put this block near the top of `STATUS.md`:

```text
REMOTE REVIEW HANDOFF
Branch: main
Commit: <full SHA>
Goal of this milestone: ...
Commands run: ...
Machine/compiler: ...
Corpus and result hashes: ...
Tests and controls: ...
Exact mathematical conclusions: ...
Heuristic/numerical observations only: ...
Unresolved graphs or cases: ...
Known trust assumptions: ...
Files the reviewer should read first: ...
Specific questions for the reviewer: ...
Recommended next local action: ...
```

Never ask the reviewer to infer a result from an uncommitted working tree or a live process. The commit hash is the review boundary.

## 9. Immediate next action after pulling this file

1. Read the incoming commit and current `STATUS.md` without discarding local work.
2. Preserve or checkpoint any running process; do not relaunch the 916,313-case cap-grinding approach unchanged.
3. Commit any currently uncommitted sound engine changes separately.
4. Implement the dimension-6 profiling milestone, beginning with the graph-only clique-link filters and the `K7`/`K6` clique distribution.
5. In parallel only if cleanly isolated, run a small obstruction-library pilot on a stratified sample.
6. Push a milestone commit with the handoff block above. The user will then ask the remote ChatGPT session to review it.
