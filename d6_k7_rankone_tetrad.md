# Exact rank-one Schur tetrads for the dimension-six K7 residue

## Outcome and production decision

The near-saturating-clique layer reaches the main class missed by the
saturating-clique positive-polynomial dual.  On a deterministic 12-graph
sample from the current 12,839-graph K7 residue, the degree-four tetrad cone
eliminated the unique surviving cover of **all 12 graphs**, producing 12
exact marginal graph rejections.  On the older frozen sample it eliminated
20 of 26 covers having no saturating clique.

The production bundle deliberately implements degree four only.  An
exploratory staged degree-six run (tetrads plus triangle signs) added zero
cover rejections on the six degree-four-resistant frozen-sample covers while
increasing one-process time from about 14 seconds to about 58 seconds.  The
go/no-go rule is therefore:

- go: commit, benchmark, and run the degree-four layer, which exceeded both
  pilot thresholds (at least one graph rejection and at least 25% cover
  coverage);
- no-go: do not run degree six at scale unless a new pilot shows positive
  marginal coverage on degree-four survivors.

## Rank-one Schur identities

Let `K` be the normalized Gram matrix on the nonzero-factor set `N`, with

```text
K_xx > 1,
K_xy in {0,1} for x != y,
rank(K) <= U.
```

Take any required clique `C` of size `U-1`.  Its principal block is

```text
K[C,C] = J + diag(e_i),  e_i>0,
x_i=1/e_i>0,             T=1+sum_i x_i,
H=K[C,C]^-1=diag(x_i)-x x^T/T.
```

For a remainder vertex `y`, let `A_y` be its zero-one neighbour mask on
`C`, and let `w(A)=sum_{i in A}x_i`.  For distinct remainder vertices define

```text
g_yz = T(K_yz - p_y^T H p_z)
     = K_yz T - T w(A_y intersect A_z) + w(A_y)w(A_z).
```

The Schur complement

```text
S=K[R,R]-K[R,C] H K[C,R]
```

is positive semidefinite of rank at most one.  Thus `S=t t^T`; after scaling
`a_y=sqrt(T)t_y`, every known off-diagonal numerator is

```text
g_yz=a_y a_z.
```

Four distinct remainder vertices consequently satisfy the exact tetrads

```text
g_ij g_kl - g_ik g_jl = 0,
g_ij g_kl - g_il g_jk = 0.
```

Triples also satisfy the non-strict sign condition

```text
g_ij g_ik g_jk = a_i^2 a_j^2 a_k^2 >= 0.
```

Only the tetrad equalities enter production.

## Exact positive-polynomial certificate

Let `E_j(x)=0` be all nonzero tetrads for one near-saturating clique.  The
degree-four locator searches for rational multipliers `q_j(x)` such that

```text
P(x)=sum_j q_j(x)E_j(x)
```

has total degree at most four, is nonzero, and every coefficient of `P` is
nonnegative.  Since every `x_i>0`, this gives `P(x)>0`; the tetrads instead
give `P(x)=0`.  The cover is impossible.

HiGHS only locates a point in a rational polyhedral cone.  The locator rounds
candidate multipliers to rationals and accepts only after recomputing every
coefficient over arbitrary-precision `Fraction`.  The independent checker
reconstructs the Schur numerators and tetrads without importing either the
locator or production runner, then expands every serialized identity again.
No floating-point result, failed LP, or absent certificate is a rejection.

Candidate nonedges are not assigned non-unit distances.  Their normalized
zero entries arise only after fixing a zero-factor cover, exactly as in the
frozen K7 rank proof.

## Pilot and controls

### Current-residue pilot

`d6_k7_rankone_tetrad_selection.py` reproducibly extracts the first 12 rows
ordered by `(cover count,index)` satisfying all of:

```text
full degree-one status = SURVIVOR,
dual-passing covers = 1,
strict-H-passing saturating cliques = 0.
```

The compact sample pins the 17,764-graph rank input, the 12,941-graph prior
selection, and the completed degree-one decisions.  Exact result:

| quantity | count |
|---|---:|
| current-residue graphs | 12 |
| prior-passing covers | 12 |
| covers with a `K_(U-1)` | 12 |
| degree-four tetrad certificates | 12 |
| marginal graph rejections | **12** |
| near cliques tested | 12 |
| one-process wall time | 4.13 s |

The rejected indices are:

```text
60920 62178 62184 247398 247400 276919 549237
1466007 1567987 1568002 3726247 3963381
```

`d6_k7_rankone_tetrad_verify.py` checks all 12 identities, the deterministic
sample selection, the source hashes, and the fact that each certificate
eliminates the unique prior-passing cover.

### Frozen-sample control

The nine graphs in the old 84-graph degree-one sample that had any relevant
near/saturating clique give:

| quantity | count |
|---|---:|
| no-saturating-clique covers | 26 |
| covers with a `K_(U-1)` | 26 |
| degree-four tetrad cover rejections | **20** |
| graph rejections relative to that old pipeline | 1 |
| one-process wall time | 13.69 s |

All 20 identities in `d6_k7_rankone_tetrad_frozen_control_report.json` pass
the independent checker.  This control is coverage evidence, not an update
to the current residue count.

Pilot reproduction:

```text
python3 d6_k7_rankone_tetrad_selection.py
python3 d6_k7_rankone_tetrad_pilot.py \
  --indices 1466007 1567987 1568002 3726247 3963381 \
            60920 62178 62184 247398 247400 276919 549237 \
  --output-degree 4 --output d6_k7_rankone_tetrad_pilot_report.json
python3 d6_k7_rankone_tetrad_verify.py
python3 -m unittest -v d6_k7_rankone_tetrad_test.py
```

## Full production workflow

The selection builder freezes the 12,839 survivors of the independently
verified full degree-one run:

```text
python3 build_d6_k7_rankone_tetrad_full_selection.py
```

Before any launch, audit the frozen dependency pins in the runner and full
verifier, run the unit tests, commit the exact source, and use that commit as
the run boundary.

The required 64-graph timing/control benchmark is:

```text
caffeinate -dimsu python3 run_d6_k7_rankone_tetrad_full.py \
  --selection-report-sha256 \
    86e4f19a203bca111a47924ae05461ada5b21b6341be59d38b1428bebe1f9479 \
  --workers 11 --max-inflight 22 --limit 64 \
  --checkpoint-dir .runs/d6_k7_rankone_tetrad_benchmark64 \
  --decisions-archive /tmp/d6_k7_rankone_tetrad_benchmark64.tsv.gz \
  --certificate-archive /tmp/d6_k7_rankone_tetrad_benchmark64.jsonl.gz \
  --checkpoint-copy /tmp/d6_k7_rankone_tetrad_benchmark64_checkpoint.json \
  --report d6_k7_rankone_tetrad_full_benchmark.json \
  --pid-file .runs/d6_k7_rankone_tetrad_benchmark64.pid.json
```

The benchmark must have zero `INFRA_ERROR` rows, preserve the 12/12 and
20/26 controls, pass the independent verifier on all 64 rows, and provide a
measured projected full wall/core-hour range before the full launch.

The full command is the same without `--limit 64` and with the production
artifact paths.  Resume uses the identical command plus `--resume`; every
graph has an atomic JSON checkpoint keyed by ordinal and graph index.

After the production report is hashed, run:

```text
python3 verify_d6_k7_rankone_tetrad_full.py \
  --report-sha256 <exact production report SHA256> \
  --workers 11 --max-inflight 22
```

The independent verifier also has atomic per-graph checkpoints and supports
the identical `--resume` discipline.

## Trust and scope

- Exact certificates trust Python arbitrary-precision integer and `Fraction`
  arithmetic, plus the frozen K7 rank and strict-H layers.
- The full verifier consumes and rechecks every prior degree-one cover
  certificate needed by a selected graph; it does not trust aggregate counts
  alone.
- HiGHS, BLAS, IEEE-754, and `libm` are discovery/search components only and
  are absent from the independent identity checker.
- `SURVIVOR`, an LP failure, and an infrastructure failure make no
  non-realizability claim.
- The present pilot proves 12 additional sampled graphs impossible.  The
  full residue count changes only after the complete production run and
  independent full verification finish with zero errors.
