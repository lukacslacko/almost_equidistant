# Bounded special-inverse profile on the 189-graph K7 residue

## Outcome

The exact low-dimensional affine-`H` route has no target on the current K7
residue.  Reconstructing the full joint-cover boundary from committed,
hash-pinned inputs gives 357 jointly surviving covers.  Of these, 268 have no
saturating clique.  The other 89 covers contain 118 saturating cliques, but
every affine system is far from unique:

| saturating rank | symmetric variables | equation rank | free dimension | cliques |
|---:|---:|---:|---:|---:|
| 6 | 21 | 10 | 11 | 54 |
| 7 | 28 | 10 | 18 | 64 |

Thus there are **zero unique systems**, **zero one-parameter systems**, and
zero graph or cover rejections from the requested direct special-inverse
test.  A univariate gcd/Sturm layer is inapplicable.

This is useful negative evidence: further work on this route would require a
genuinely multivariate positivity or elimination argument, not a small
extension of the exact affine screen.

## Exact reconstruction boundary

`probe_d6_k7_saturating_inverse.py` reads no local per-graph checkpoint.  It
starts from:

```text
d6_k7_joint_support_full_report.json
  b3fd2752a371f8909b79b50b8e2b6135f93423d0adf7286aee26d4a8c0643326
d6_k7_joint_support_full_decisions.tsv.gz
  b6040e7796e3d7df511a4e460ad71082e63c074ac2cf2713c19d8f3b772bcc90
```

and calls the same hash-pinned baseline, strict-`H`, degree-one/tetrad
certificate, labeled-support propagation, and sparse-value implementations
used at the committed boundary.  It independently re-enumerates every K7
seed and eligible cover.  For each of the 189 survivor graphs, its current
cover count and joint-pre-capacity-passing cover count are required to equal
the corresponding committed decision row.  The ordered 189-index list has
SHA-256

```text
a399d45dfb06ffccd4c1c6f0f447c2623bc9faada47212b3b37766a134efe034
```

The reconstructed aggregate accounting is:

```text
graphs                                                   189
K7 seeds                                                 274
eligible covers                                       99,352
current-layer passing covers                             572
joint-support failing among those                        215
jointly surviving covers                                 357

joint covers with 0 saturating cliques                   268
joint covers with 1 saturating clique                     65
joint covers with 2 saturating cliques                    19
joint covers with 3 saturating cliques                     5
total saturating cliques                                  118
graphs with at least one saturating clique                 20
graphs with none                                          169
```

A separate one-worker reconstruction produced the identical exact profile
and the identical ordered six-system algebra sample.  It took 146.49 seconds;
the 11-worker run took 29.51 seconds.  Runtime is not part of any decision.

## Special-inverse characterization checked

For a saturating clique `C`, put

```text
H = K[C,C]^-1,
r = H 1,
delta = 1 - 1^T r.
```

If `H=(J+diag(e_i))^-1` with every `e_i>0`, then

```text
delta > 0,
H_ij = -r_i r_j / delta                 for i != j,
delta H_ii = r_i(1-r_i).
```

The probe reduces the outside-pair equations over `fractions.Fraction` using
an implementation local to the probe.  If a system is unique, it checks all
original affine equations again, checks the strict diagonal, off-diagonal,
row-sum, total-sum, outside-diagonal, and distinct-column quadratic-form
inequalities, applies exact Sylvester tests, checks the displayed polynomial
identities, and independently inverts `H` to verify off-diagonal entries one
and diagonal entries strictly above one.  No actual residue system reaches
that unique branch.

The positive control uses

```text
J + diag(2,3,5).
```

Its exact inverse passes the identities with

```text
delta = 30/61,       r = (15/61, 10/61, 6/61).
```

Forcing its six symmetric entries gives a feasible synthetic unique affine
system.  A rational off-diagonal perturbation is detected by nine independent
identity/inverse failures.  These are controls of the special-inverse kernel,
not realizable 19-point controls.

## Linearized quadratic profile

Writing `x_i=1/e_i`, `T=1+sum_i x_i`, and `w(A)=sum_{i in A}x_i`, every pair
of outside basis-neighbour masks gives the integral quadratic

```text
T w(A intersect B) - w(A)w(B) - target*T = 0,
target in {0,1}.
```

Linearizing in all monomials of degree at most two gives no hidden cheap rank
collapse:

| `u` | monomials through degree 2 | exact row rank | formal free dimension | systems |
|---:|---:|---:|---:|---:|
| 6 | 28 | 10 | 18 | 54 |
| 7 | 36 | 10 | 26 | 64 |

## Strictly bounded modular elimination pilot

Only the first three deterministically ordered systems for each of `u=6`
and `u=7` were screened.  SymPy used graded lexicographic order modulo each
of `101`, `103`, and `107`, with a hard eight-second wall timeout per child.
The 18 jobs ran in short-lived child processes, at most six concurrently.

All 18 children exited normally.  Their PIDs were

```text
86777 86778 86779 86780 86781 86782
86784 86785 86787 86788 86789 86790
86832 86833 86834 86835 86836 86837
```

They ran between `2026-08-06T21:24:17.337799Z` and
`2026-08-06T21:24:18.674140Z`; every recorded stop mode is `EXIT`.  Individual
algebra time was 0.27--0.43 seconds.  The modular bases had:

```text
basis polynomials          18--21
maximum total degree            3
maximum terms per basis element 16--27
total terms               214--314
```

None of the 18 modular bases was `[1]`.  Per the frozen rule, no exact-Q
Groebner job was launched and no rejection was credited.  A second serial run
with a one-second cap obtained the same nonunit result on all 18 jobs.

The small modular bases show that the computations themselves are cheap, but
their positive-dimensional, nonunit output supplies no obstruction.  Given
that 169 of 189 graphs have no jointly surviving saturating-clique cover at
all, indiscriminate continuation of this elimination route is not currently
recommended.  A future continuation should first identify an additional
sound equation or positivity certificate that can turn these
positive-dimensional ideals into contradictions.

## Reproduction

The full command was:

```sh
/Users/lukacs/claude/opengauss/venv/bin/python3 \
  probe_d6_k7_saturating_inverse.py \
  --workers 11 \
  --groebner-timeout 8 \
  --output /private/tmp/d6_k7_saturating_inverse_report.json
```

The generated report and source hashes for this run are:

```text
/private/tmp/d6_k7_saturating_inverse_report.json
  66577457ea0ecdd98d4632f9de1b79e519af26d350d049ffc477e7f58d0df74f
probe_d6_k7_saturating_inverse.py
  e556219295d223cdfc52d777795ecc1d92242db98b88698d5d3d504618eed2a7
```

The temporary report contains the exact identity, masks, targets, term-growth
statistics, PID, start time, finish time, and stop mode for every bounded
algebra job.  It assigns theorem credit only to exact rational checks; modular
nonunit output, timeout, and infrastructure failure are explicit nonclaims.
