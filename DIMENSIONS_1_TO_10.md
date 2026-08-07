# Almost-equidistant sets in dimensions 1–10

_Last checked: 2026-08-06._

A finite set `X ⊂ R^d` is **almost equidistant** if every three points of
`X` contain two points at distance exactly `1`.  Write `f(d)` for the
maximum possible size.

This is a living tracker.  It deliberately distinguishes established
literature from the newer computer-assisted results and constructions in this
repository.

## Current table

| `d` | lower bound | upper bound | status and provenance |
|---:|---:|---:|---|
| 1 | 4 | 4 | elementary |
| 2 | 7 | 7 | classical; see BPSSV |
| 3 | 10 | 10 | Györey; independently treated in BPSSV |
| 4 | **12** | **12** | computer-assisted result in this repository; published literature before this work had `12 ≤ f(4) ≤ 13` |
| 5 | **16** | **16** | computer-assisted result in this repository; published literature before this work had `16 ≤ f(5) ≤ 20` |
| 6 | 18 | 26 | BPSSV; this repository has an active attempt to eliminate every 19-point candidate |
| 7 | 20 | 34 | BPSSV |
| 8 | 24 | 40 | lower bound BPSSV; upper bound from `R(3,10) ≤ 41` |
| 9 | 24 | 49 | lower bound BPSSV; upper bound from `R(3,11) ≤ 50` |
| 10 | **26** | 58 | 26-point construction proved below; upper bound from `R(3,12) ≤ 59` |

The bold entries in dimensions 4, 5, and 10 are newer project results rather
than claims imported from the older literature.

## Repository results in dimensions 4 and 5

### Dimension 4

[`f4_equals_12.pdf`](f4_equals_12.pdf) and its source
[`f4_equals_12.tex`](f4_equals_12.tex) give a computer-assisted proof of

```text
f(4) = 12.
```

The upper-bound computation enumerates the 59 abstract 13-vertex candidate
graphs and certifies that none has a realization by 13 distinct points in
`R^4`.  The interval engine, controls, candidate file, and reproduction driver
are included in the repository.

### Dimension 5

[`f5_equals_16.pdf`](f5_equals_16.pdf) and its source
[`f5_equals_16.tex`](f5_equals_16.tex) give a computer-assisted proof of

```text
f(5) = 16.
```

The campaign certifies all 21,827 minimal abstract candidates on 17–20
vertices as non-realizable in `R^5`:

```text
n = 17: 12,654
n = 18:  8,825
n = 19:    340
n = 20:      8
```

The 16-point half-cube/Clebsch lower bound is checked separately in integer
arithmetic.  See [`README.md`](README.md) for reproduction commands.

Both results are computer-assisted manuscripts.  They should remain labelled
as such until independently reproduced and subjected to normal mathematical
peer review.

## Why Ramsey numbers give upper bounds

Let `G` be the unit-distance graph of an almost-equidistant set in `R^d`.
Then

- `alpha(G) ≤ 2`, because an independent triple would contain no unit pair;
- `G` contains no `K_(d+2)`, because at most `d+1` equidistant points fit in
  `R^d`.

Equivalently, the complement of `G` is triangle-free and has no independent
set of size `d+2`.  Hence

```text
f(d) ≤ R(3,d+2) - 1.
```

For dimensions 6 and 7, BPSSV obtain the stronger specialized upper bounds 26
and 34.  For dimensions 8–10, the table uses the current Ramsey bounds

```text
R(3,10) ≤ 41,
R(3,11) ≤ 50,
R(3,12) ≤ 59.
```

## A 26-point construction in `R^10`

This section gives a self-contained proof that

```text
f(10) ≥ 26.
```

It combines a 16-point spherical almost-equidistant set in one copy of `R^5`
with a 10-point spherical almost-equidistant set in an orthogonal copy of
`R^5`.

### Orthogonal-join lemma

Suppose `X ⊂ R^a` and `Y ⊂ R^b` are almost equidistant, every point of `X`
has squared norm `r^2`, every point of `Y` has squared norm `s^2`, and

```text
r^2 + s^2 = 1.
```

Embed them in orthogonal coordinate subspaces of `R^(a+b)`.  Every cross pair
then has squared distance

```text
||(x,0) - (0,y)||^2 = ||x||^2 + ||y||^2 = 1.
```

Therefore their union is almost equidistant.

### The 16-point component of squared radius `5/8`

Let

```text
X = { eps/(2 sqrt(2)) : eps ∈ {±1}^5 and product_i eps_i = 1 }.
```

There are 16 points, each of squared norm `5/8`.  Two distinct sign vectors
have Hamming distance 2 or 4, so their squared Euclidean distance after
scaling is respectively 1 or 2.

The non-unit pairs cannot contain a triangle.  After multiplying coordinates
by one member of a hypothetical triangle, the other two members would both
have four minus signs.  Two distinct 4-subsets of a 5-set differ in only two
positions, not four.  Thus every three points of `X` contain a unit pair.

### The 10-point Petersen component of squared radius `3/8`

Represent the Petersen graph `P` by the ten 2-subsets of `{1,2,3,4,5}`, with
adjacency meaning disjointness.  Let `A` be its adjacency matrix and `J` the
10-by-10 all-ones matrix.  Define

```text
B = (4 I + 2 A - J)/3.
```

The Petersen graph is strongly regular and satisfies

```text
A^2 = 2 I - A + J.
```

Its adjacency spectrum is

```text
3^1, 1^5, (-2)^4.
```

Consequently `B` is positive semidefinite with spectrum

```text
2^5, 0^5.
```

It is therefore the Gram matrix of ten unit vectors
`v_1,...,v_10 ∈ R^5`.  Its diagonal entries are 1, while for distinct
vertices

```text
<v_i,v_j> =  1/3   if ij is an edge of P,
<v_i,v_j> = -1/3   otherwise.
```

Set

```text
y_i = sqrt(3/8) v_i.
```

Each `y_i` has squared norm `3/8`, and

```text
||y_i-y_j||^2 = 1/2   on Petersen edges,
||y_i-y_j||^2 = 1     on Petersen nonedges.
```

The Petersen graph is triangle-free—three pairwise-disjoint 2-subsets of a
5-set do not exist—so every three of the `y_i` contain a unit pair.  Thus
`Y={y_1,...,y_10}` is almost equidistant.

### Joining the components

Place `X` and `Y` in orthogonal copies of `R^5` inside `R^10`.  Their squared
radii add to one:

```text
5/8 + 3/8 = 1.
```

Every cross pair is at unit distance, while each component is internally
almost equidistant.  The union has `16+10=26` distinct points, proving

```text
f(10) ≥ 26.
```

This construction was derived during this project on 2026-08-06 and has
been machine-verified in exact rational arithmetic at the Gram level
([`verify_lower_bound_26_d10.py`](verify_lower_bound_26_d10.py): both
component Grams are PSD of rank 5, all 26 points distinct, every triple
contains a unit pair).  This file
makes no priority claim: a dedicated literature check should precede any
claim that the bound is new.  The Petersen component is closely related to the
standard ten equiangular lines of common angle `1/3` in `R^5`.

## Dimension 6 campaign

The known 18-point construction gives `f(6) ≥ 18`.  Since the property is
hereditary, proving that no 19-point set exists would immediately establish
`f(6)=18`.

### Provisional upper-bound improvement (2026-08-07): `f(6) ≤ 20`, conditional

[`f6_le_20.pdf`](f6_le_20.pdf) proves `f(6) ≤ 20` **conditional on the
level-19 certificate ledger** stated there: every one of the 3,971,787
level-19 candidates outside an explicit 644-graph residue is certified
non-realizable (3,855,094 by exact graph filters that have been
independently reimplemented with exact agreement; 2,451 by this
repository's interval certificates alone; 113,598 by the codex track's
later exact layers, whose independent reimplementation is still in
progress — this last class is why the result is labelled conditional and
the table above still shows 26).  The argument is a frontier/extension
computation: any 21-point set would force a chain of candidate
restrictions into the residue, and the complete level-21 frontier over
the residue is empty.  The same computation reduces `f(6) ≤ 19` to the
non-realizability of 34 explicit 20-vertex graphs
([`frontier20_true.txt`](frontier20_true.txt)), currently under interval
certification.  Computer-assisted, unpublished, not peer-reviewed.

The repository contains a verified corpus of 3,971,787 minimal abstract
19-vertex candidates and an active certified-realizability campaign.  Partial
coverage is research progress, but **does not improve the formal upper bound**
until every 19-point candidate is rigorously eliminated.  See
[`STATUS.md`](STATUS.md) and [`d6_plan.md`](d6_plan.md) for the current state.

## Update rules for this tracker

When changing the table:

1. give a proof or stable source for every lower-bound construction;
2. give the exact Ramsey or geometric theorem used for every upper bound;
3. mark computer-assisted, unpublished, and non-peer-reviewed claims clearly;
4. do not turn a partial candidate computation into an upper bound;
5. record the date and provenance of any project-only construction.

## References

1. M. Balko, A. Pór, M. Scheucher, K. Swanepoel, P. Valtr,
   “Almost-Equidistant Sets,” *Graphs and Combinatorics* 36 (2020), 729–754.
   DOI: <https://doi.org/10.1007/s00373-020-02149-w>;
   arXiv: <https://arxiv.org/abs/1706.06375>.
2. V. Angeltveit, “R(3,10) ≤ 41,” *Electronic Journal of Combinatorics*
   32(4) (2025), P4.30. DOI: <https://doi.org/10.37236/12936>.
3. S. Radziszowski, “Small Ramsey Numbers,” *Electronic Journal of
   Combinatorics*, Dynamic Survey DS1, revision of 2026-04-24.
   DOI: <https://doi.org/10.37236/21>.
4. Local result notes: [`f4_equals_12.pdf`](f4_equals_12.pdf) and
   [`f5_equals_16.pdf`](f5_equals_16.pdf).
