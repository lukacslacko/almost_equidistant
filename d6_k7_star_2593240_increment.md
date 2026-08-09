# Production exact K7 rank-one-star increment for graph 2593240

## Claim and boundary

This package promotes the frozen exploratory contradiction for graph
`2593240` to a source-bound, independently checked one-graph increment.  Its
input is the ordered 15-graph output of
`d6_k7_schur_3949382_increment_report.json`; its output is that list with
`2593240` removed.

The producer is `build_d6_k7_star_2593240_increment.py`.  It imports the
frozen exploratory probe, checks its exact stable report hash, validates the
preceding exact increment and its verifier, and refuses an official launch
unless all seven proof/checking source files equal blobs in the launch commit.
The checker is `verify_d6_k7_star_2593240_increment.py`.  It imports neither
the producer nor the probe and independently reconstructs both seed
quantifiers, both exact Schur patterns, their isomorphism, and the algebra.

No floating-point step enters the rejection.  Candidate nonedges remain
optional.  A normalized-Gram entry is set to zero only when two independently
propagated support supersets are disjoint; a one entry comes only from a
required unit edge.

## Exhaustive two-seed quantifier

The graph has exactly two required K7 seeds:

```text
[2,4,7,10,13,15,18]
[6,7,8,9,12,15,18].
```

The independent checker enumerates 502 raw eligible zero-factor covers for
each seed.  The pinned support layers leave current covers `[0,8]` for the
first seed and `[0,32]` for the second.  Covers `8` and `32` respectively have
no passing support family.  Cover `0` has exactly one family for each seed,
with propagated support supersets

```text
[100,127,5,98,25,19,26,108,11,20,38,72]
[127,98,21,10,28,97,25,102,13,18,35,68].
```

The exact local K6 masks and remainder targets in lexicographic pair order are

```text
masks   = [43,30,54,58,46,11]
targets = 011111111110111

masks   = [57,29,43,53,45,11]
targets = 111110111111101.
```

The checker reconstructs all 102 displayed basis/remainder entries and aborts
if any optional entry is not forced by disjoint support supersets.  It also
checks the exact isomorphism from the first system to the second:

```text
remainder permutation = [2,1,3,0,4,5]
K6-coordinate map     = [1,0,2,3,4,5].
```

Only one infeasible required K7 seed is logically needed to reject the graph,
but this package quantifies and rejects both.

## Rank-one Schur contradiction

For a required K6 block, write

```text
M = J + diag(e_0,...,e_5),  e_i>0,
x_i=1/e_i>0,                T=1+sum_i x_i.
```

The block is positive definite.  Since the normalized Gram rank is at most
seven, its Schur complement is PSD of rank at most one.  If `A_i` is the
K6-neighbour mask of remainder vertex `i`, its known scaled off-diagonal entry
is

```text
g_ij = T k_ij - T w(A_i intersect A_j) + w(A_i)w(A_j).
```

Rank one requires `g_ij=h_i h_j`.  Hence the 30 four-distinct tetrads vanish.
Writing the positive inverse diagonal weights as `(a,b,c,d,e,f)`, one tetrad
is

```text
-f(a+1)T(c-e)=0,
```

so `e=c`.

The producer's frozen proof uses a two-stage grevlex ideal-membership
calculation.  The independent verifier instead performs one lexicographic
elimination after substituting `e=c`.  It obtains exact triangular
consequences including

```text
(c-1)^2(c+1)(2c+1)(3c+1)=0,
(c-f)(2c+1)=0,
```

followed by linear triangular relations for `d`, `b`, and `a`.  Strict
positivity forces the unique tetrad point

```text
(a,b,c,d,e,f) = (2/3,5,1,1,1,1).
```

At this point every four-distinct tetrad is zero, but the nonzero known Schur
entries are

```text
g_01=-8/3, g_02=8, g_03=-8/3, g_04=-8/3, g_05=-28/3,
```

and all leaf-leaf entries vanish.  In particular `g_01,g_02` are nonzero but
`g_12=0`.  A factorization `g_ij=h_i h_j` would make `h_0,h_1,h_2` nonzero,
forcing `g_12` nonzero.  This is the exact rank-one contradiction.

The positive point is retained as a control: it passes every old tetrad, so a
checker that accidentally omits the three-index support-closure condition
will fail the package test.

## Freeze and official run

First commit the four production source files without generating artifacts.
Then run serially:

```text
python3 -m unittest -v test_d6_k7_star_2593240_increment.py
python3 build_d6_k7_star_2593240_increment.py
python3 verify_d6_k7_star_2593240_increment.py \
  --report-sha256 <printed-report-sha256>
```

The report and verification JSON files are the union-ready evidence.  Commit
them without changing the frozen source boundary between producer and checker
runs.
