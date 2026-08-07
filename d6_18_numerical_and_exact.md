# Dimension-6 18-deletion screen and two exact nonstandard 18-sets

## What was screened

`run_d6_18_numerical_screen.py` screened all 12,712 unrooted isomorphism
classes in `d6_residue_18_deletions.json`.  It compiled
`d6_18_lm_screen.c`, which reuses the existing `lm5.c` kernel at dimension 6,
and ran six deterministic Levenberg--Marquardt starts per class in 199 atomic
chunks of at most 64 classes.  Eleven independent worker processes kept the
11 usable logical CPU cores occupied.  The completed invocation took 152.3
seconds (83.45 classes/second including orchestration).

The exact command was

```text
python3 run_d6_18_numerical_screen.py \
  --workers 11 --chunk-size 64 --restarts 6 --seedbase 600180001
```

The same command resumes the hash-addressed checkpoint directory
`d6_18_numerical_checkpoints/3eaad0fa1ff43461`.  Each chunk is validated and
atomically renamed.  The report pins the corpus and verification hashes, both
C sources, the orchestration source, compiler version and flags, binary hash,
seed formula, thresholds, worker count, exact chunk hashes, and raw result
hash.

This screen is heuristic only.  It rejects no graph.  Candidate nonedges were
unconstrained, and a failure to reach a small residual is not evidence of
non-realizability.

## Numerical results

All 14 classes already known to embed into the standard 18 coordinates were
included as explicit seeded positive controls.  All 14 passed with residual
about `5e-30`, minimum separation 1, and numerical rigidity ranks 86 or 87
matching exact mod-13 rank lower bounds.

The random-start sweep found exactly two distinct near-zero endpoints:

| deletion class | random residual | minimum separation | numerical rank |
|---:|---:|---:|---:|
| 10261 | `2.35e-30` | `1/2` | 87 |
| 10887 | `1.43e-29` | `1/2` | 87 |

Both supports were already among the 14 standard-compatible controls, but
these endpoints are not the standard realization: the standard set has
minimum separation 1.  No noncontrol class reached residual `1e-20`, `1e-8`,
or even `1e-2`.  The smallest noncontrol residual was about `0.04094`; the
median over all classes was about `0.58894`.  These figures only prioritize
future work.

The numerical report retains the vertex-major binary64 coordinates for both
near-zero endpoints.  `verify_d6_18_numerical_screen.py` independently checks
the complete raw table, hashes, checkpoint coverage, seeded controls, summary
counts, and the coordinate residuals.  It still makes zero mathematical
rejections.

## Exact reconstruction

The two numerical squared-distance matrices round, with maximum error below
`1.1e-14`, to the following exact spectra:

| type | exact squared-distance multiset | unit pairs |
|---|---|---:|
| class 10261 | `1/4` (2), `1/2` (8), `7/8` (1), `1` (110), `3/2` (1), `2` (31) | 110 |
| class 10887 | `1/4` (1), `1/2` (5), `1` (111), `3/2` (1), `2` (35) | 111 |
| standard 18 | `1` (112), `3/2` (1), `2` (40) | 112 |

`reconstruct_d6_18_exact.py` records both full rational distance matrices.  If
`D` is either matrix and `C=I-J/18`, it computes

```text
B = -1/2 C D C.
```

For class 10261, the principal basis `J=[0,1,2,3,4,6]` has determinant
`3/1024`; for class 10887, `J=[0,1,2,3,4,5]` has determinant `1/432`.  Every
leading principal minor is positive, and exact rational arithmetic verifies

```text
B = B[:,J] inverse(B[J,J]) B[J,:].
```

Thus `B` is positive semidefinite of rank exactly 6.  It reproduces `D`, all
off-diagonal distances are positive, every required support edge has squared
distance 1, and an exhaustive exact check finds a unit pair in every triple.
The matrices therefore certify two genuine 18-point almost-equidistant sets in
`R^6`.  Their 110 and 111 unit-pair counts distinguish them from each other and
from the 112-edge standard set, proving at least three nonisometric 18-point
configurations.

## The switching construction behind them

Let the standard base consist of the 16 odd-parity vectors

```text
q in {+1,-1}^5 / (2 sqrt(2)),   ||q||^2 = 5/8,
```

in a five-flat, and let

```text
a_+ = (0,...,0,sqrt(3/8)),   a_- = -a_+.
```

Choose a base point `q`, choose an apex sign, and replace `q` by

```text
x_q = (a_sign - q)/2.
```

For every unchanged base point `p`, a direct inner-product calculation gives

```text
||x_q-p||^2 = 1     if ||q-p||^2 = 1,
||x_q-p||^2 = 1/2   if ||q-p||^2 = 2.
```

Also, `x_q` has squared distance `1/4` from its assigned apex and 1 from the
opposite apex.  One replacement is exactly the 111-edge type.

The distance-2 graph on the 16 base points is the Clebsch graph, with exact
strongly regular parameters `(16,5,0,2)`.  Two replacements give a valid set
exactly when the original base vertices are Clebsch-adjacent (their standard
squared distance is 2) and the assigned apex signs are opposite.  Their mutual
squared distance is then `7/8`, giving the 110-edge type.  The conditions are
forced:

- same-sign replacements and their assigned apex form a bad triple;
- Clebsch-nonadjacent base vertices have two common Clebsch neighbours, each
  non-unit from both replacements;
- every pair of replacements is non-unit, so three replacements themselves
  form a bad triple.

Sign flips in an even number of coordinates act transitively on the odd base
vectors.  After composing coordinate permutations with the even sign flips
that return a chosen base vertex, its stabilizer acts as `S_5` and is
transitive on distance-2 neighbours; swapping the apices swaps signs.  Hence
this switching operation has exactly the 0-, 1-, and 2-replacement types up to
symmetry.  This classifies the switching family only, not all 18-point
configurations.

## Exact global nonextension

Both new configurations are dead ends for growing back to 19, and this is
proved globally rather than only for the rooted parent attachments.

For any added point `x`, its non-unit neighbours in the 18-set must form a unit
clique; otherwise `x` and two of those neighbours make a bad triple.  Extend
that clique to a maximal unit clique `K`, so `x` must be unit from every vertex
outside `K`.

Write centered points as `q_i`, choose the six-point principal basis above,
write `x=sum_j y_j q_j`, and let `s=||x||^2`.  For every forced unit neighbour
`i`, the sphere equation is the exact linear equation

```text
-2 sum_j B[i,j] y_j + s = 1 - B[i,i].
```

The exact checker enumerates every maximal unit clique and solves these
rational systems:

| type | maximal cliques | clique sizes | linearly consistent | norm failures |
|---|---:|---|---:|---|
| class 10261 | 92 | 70 of size 5, 22 of size 6 | 2 | both have `s-y^T B[J,J]y=1/3` |
| class 10887 | 102 | 75 of size 5, 27 of size 6 | 17 | sixteen have discrepancy `1/3`, one has `10/27` |

Every other linear system is inconsistent.  None passes the quadratic norm
identity, so no point at all satisfies the necessary unit spheres.  In
particular, neither exact 18-set can extend to an almost-equidistant 19-set.

## Artifacts and hashes

```text
d6_18_numerical_screen_report.json
  0a3e1f0288b5eec4ecec4357d1ba1b75e91b1e4cd2a2d258273357c742a581d8
d6_18_numerical_screen_results.tsv.gz
  41807f0cec494344453c939dfaa3530fef47404e6c298b391d6e2e670b0bc95c
d6_18_numerical_screen_verification.json
  da4a876eef7e793c90d0ececbdf823839f899ee3548377382f53aa0139b81a7e
d6_18_numerical_checkpoints_3eaad0fa1ff43461.tar.gz
  3d3bfbdcd71ccbf6d654684934eb472b8f56438f02a5fecb82dd7c12dc240c61
d6_18_exact_reconstruction.json
  69fa1c0be5f46e285773f13cc05270219d606673703fd16b3e36a9551cb9a990
d6_18_exact_reconstruction_verification.json
  e86f8b538f4f78c8a0ace71a13d2e029be054e2135d19ed710faf4fc21aca887
```

Verification commands:

```text
python3 verify_d6_18_numerical_screen.py
python3 verify_d6_18_exact.py
python3 -m unittest -v \
  test_d6_18_numerical_screen.py \
  test_d6_18_exact_reconstruction.py
```
