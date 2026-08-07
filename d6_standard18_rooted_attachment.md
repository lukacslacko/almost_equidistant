# Clebsch scope and the rooted standard-18 attachment audit

## What is Clebsch here?

Write the 16 base points, in coordinates scaled by `2 sqrt(2)`, as the odd
sign vectors in `{+-1}^5`.  Two base vectors differ in either two or four
coordinates.  Hamming distance two has scaled squared Euclidean distance `8`
and is therefore unit distance; Hamming distance four has scaled squared
distance `16` and is non-unit.

Thus the **non-unit graph** on the base is the 5-regular Clebsch graph `C`,
and the required/full unit graph on the base is its 10-regular complement.
This convention matters because both complementary graphs are sometimes
called “the Clebsch graph.”  The two apices are unit from all 16 base points
and non-unit from each other, so the non-unit graph of the full standard
18-set is

```text
C disjoint-union K2.
```

The Clebsch graph has strongly regular parameters `(16,5,0,2)`: adjacent
vertices have no common neighbour and nonadjacent vertices have two.  It is
therefore triangle-free and maximal triangle-free.  If an abstract
18-vertex required-unit support `G` embeds into the standard coordinates,
its complement must contain `C disjoint-union K2`.  The complement of `G`
is triangle-free.  After choosing the two pole vertices, its restriction to
the remaining 16 vertices is consequently a triangle-free supergraph of the
maximal triangle-free graph `C`, and hence equals `C`.  This is why
`build_d6_residue_18_deletions.py` recognizes an induced Clebsch graph in
the **candidate complement**.

This is a statement about the abstract support and one known coordinate
embedding.  A candidate nonedge is unconstrained and may become unit in a
different realization.

## What is and is not unique in dimension five

The local computer-assisted result `f(5)=16` proves a cardinality upper bound:
it eliminates all minimal abstract supports on 17 through 20 vertices and
verifies one 16-point lower-bound construction.  It does not classify all
16-point almost-equidistant sets or all embeddings of their supports; see
`f5_equals_16.tex`, especially lines 43--48 and 74--94.

There is a separate and narrower uniqueness theorem for **two-distance
sets**.  Lisoněk classified the maximum two-distance sets through dimension
seven; in dimension five the maximum has 16 points and is the Clebsch
half-cube, unique up to similarity.  Nozaki and Shinohara explicitly restate
the dimension-five entry as `DS_5(2)=16` with the optimal set isomorphic to
the Clebsch construction.  These results do not imply that every extremal
almost-equidistant set is a two-distance set.

References:

- P. Lisoněk, [*New Maximal Two-Distance Sets*](https://doi.org/10.1006/JCTA.1997.2749), Journal of Combinatorial Theory, Series A 77 (1997), 318--338.
- H. Nozaki and M. Shinohara, [*On a generalization of distance sets*](https://arxiv.org/abs/0906.0199), Journal of Combinatorial Theory, Series A 117 (2010), 810--826; see Section 4.3.
- M. Balko, A. Pór, M. Scheucher, K. Swanepoel and P. Valtr, [*Almost-equidistant sets*](https://arxiv.org/abs/1706.06375), Graphs and Combinatorics 36 (2020), 729--754.  Its uniqueness claims concern the then-known extremal cases in dimensions two and three, not dimension five.

## Exact rooted attachment result

The deletion manifest contains 14 unique supports, 39 deletion occurrences,
and 16 parents compatible with the standard 18 coordinates.  Three supports
have two pole-pair witnesses, giving 42 rooted occurrence/witness
representatives.

For every one of the 42 representatives, the missing vertex has, after exact
transport through an occurrence isomorphism,

```text
11 or 12 neighbours in the half-cube base, and exactly one apex neighbour.
```

Every 11 base vertices affinely span the five-flat.  If a point `x=(y,t)` is
equidistant from those base neighbours, subtracting sphere equations forces
`y=0`.  In the scaled coordinates, equality with the sphere about the one
apex then forces

```text
t = -sqrt(3)/3  for the +sqrt(3) apex,
t = +sqrt(3)/3  for the -sqrt(3) apex.
```

Its common scaled squared distance is `5+t^2=16/3`, whereas a unit distance
has scaled square `8`.  Hence the required attachment is impossible at every
represented standard-coordinate embedding orbit.  This is exact arithmetic
in `Q(sqrt(3))`, with no numerical tolerance.

This still gives **zero parent rejections by itself**.  A compatible support
could have another noncongruent realization.  The certificate is intended to
remove the known standard component inside a complete algebraic or interval
classification of the rooted support.

## Rigidity split and suggested use

At the standard coordinates, the rigidity matrices of the 14 compatible
supports have exact ranks

```text
rank 87 (infinitesimally rigid): 9 classes
rank 86 (one extra infinitesimal flex): 5 classes
```

Rank 87 is certified by reduction modulo 13 with `sqrt(3)=4`, together with
the `108-21=87` rigid-motion upper bound.  For each rank-86 class, the report
also archives a nonzero exact infinitesimal flex which vanishes on an
affinely spanning unit `K6`; this gives rank at most 86, while the modular
calculation gives rank at least 86.

Fourteen of the 16 compatible parents have at least one rank-87 compatible
deletion.  The only parents whose compatible deletions all have rank 86 are
corpus indices `3950119` and `3950509`.

A focused next step is therefore:

1. anchor a `K6` and certify all realization components of a rank-87
   18-deletion support;
2. discard the exact standard root with the rooted `16/3 != 8` certificate;
3. test the missing attachment on every other isolated root;
4. for the two rank-86-only parents, first try a different deletion, or keep
   the missing attachment equations active while resolving the one
   infinitesimal direction.

Infinitesimal rigidity is only local.  It is not a global uniqueness theorem.

## Reproduction

```sh
python3 d6_standard18_rooted_attachment.py
python3 verify_d6_standard18_rooted_attachment.py
python3 -m unittest -v test_d6_standard18_rooted_attachment.py
```

Pinned/generated SHA-256 values:

```text
d6_residue_18_deletions.json                       9d08f9ec579434c9601ad3908bd2ae51f10e09e9d06a7f38fbc25794a634b732
d6_residue_18_deletions_verification.json          50ab5d88441b2869e60a05f3e1485db0b25008c3a816960a20a3bb21c1c4defd
d6_standard18_rooted_attachment.py                 ec06eaf6a09a754597d9a66b35925da35c28314dafe200d186f65bc04287cd71
d6_standard18_rooted_attachment_report.json        28391ec07b07d299d4c3fc75fa7c8f94663f08e25a29168ae681d105cd1140ce
d6_standard18_rooted_attachment_verification.json  67ac89d11330454972fe37339c27631116b5c702ab0b3f1298c3706ae882d14d
```
