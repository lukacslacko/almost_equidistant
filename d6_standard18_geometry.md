# Exact facts about the standard 18-point construction in R6

The construction consists of the 16 odd half-cube points in a five-flat,

```text
{s in {+1,-1}^5 : s has odd parity} / sqrt(8),
```

and the two points at height `+sqrt(3/8)` and `-sqrt(3/8)` orthogonal to that
flat.  This note concerns those coordinates and their full actual unit graph.
It is deliberately not a classification theorem.

## Exact audited facts

The full unit graph has 112 edges: 80 inside the half-cube and 32 from the
two poles to all base vertices.  The poles are not unit from one another.
Its maximal cliques are 80 copies of `K5` and 32 copies of `K6`; it contains
no `K7`.  The base unit graph has clique number five.

The 112-edge framework is infinitesimally rigid in `R6`.  Its rigidity matrix
has exact rank

```text
6*18 - binom(7,2) = 87,
```

certified by a nonzero modular minor over `F_13` with `sqrt(3)=4`.  The base
framework has exact rank 65 in `R5`.  These ranks prove local rigidity of the
displayed frameworks; they do not prove global uniqueness.

## Exact nonextension of these coordinates

Suppose a distinct nineteenth point `x` extended this particular full
18-point set.  Its non-unit neighbours among the 16 base vertices must form a
unit clique: otherwise `x` and two mutually non-unit base vertices would be a
bad triple.  The base clique number is five, so `x` is unit from at least 11
base vertices.

Every 11 of the 16 base points affinely span the five-flat (all 4,368 subsets
were checked by exact modular rank).  Subtracting the unit-sphere equations
for those 11 points therefore forces the horizontal projection of `x` to be
zero.  Each base point has squared norm `5/8`; one sphere equation then gives

```text
x_6^2 = 1 - 5/8 = 3/8.
```

Thus `x` equals one of the two existing poles, contradicting distinctness.
So the displayed standard 18-point coordinates admit no nineteenth point.

## What f(5)=16 does and does not prove

The completed dimension-five campaign proves the cardinality statement
`f(5)=16`: it supplies the half-cube construction and certifies that no
17-point almost-equidistant set exists.  It does **not** classify all
16-point extremizers.  In particular it does not prove that every
16-point set is congruent to the half-cube, nor that every 18-point set in
`R6` contains that half-cube in a five-flat.

Accordingly, the exact failure to extend the displayed 18-point component is
not by itself a rejection of a 19-vertex abstract support.  To make that
deduction one would need either a global classification of every realization
of the relevant 18-support, or a direct rooted grow-back obstruction valid on
all of its realization space.

## Artifacts

```text
d6_standard18_geometry_report.json
  97904946b382eb1e6cc325083e1c1ab20e5e8c426fcbaaed1b79357646c846ef
d6_standard18_geometry_verification.json
  12e3e1df3919dd465a54b1c470f524a3907db2a19110a34876137eeedeb926d6
```

The independent checker reconstructs the coordinates, graph, clique data,
affine ranks, rigidity ranks, and nonextension calculation without importing
the producer.
