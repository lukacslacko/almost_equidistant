# Exact triple-petal obstruction for K7 graph 3936435

## Result

The source-bound package

```text
build_d6_k7_triple_petals_3936435_increment.py
verify_d6_k7_triple_petals_3936435_increment.py
test_d6_k7_triple_petals_3936435_increment.py
```

certifies that current K7-residue graph **3936435** is not realizable by
distinct points in `R^6` with every candidate edge at unit distance.  Removing
that graph changes the ordered K7 residue from 12 graphs to 11; its new ordered
index hash is

```text
380936282d98f2e561c71680a04d04401fabc032e8a00cebf12ae4d18967a874
```

The proof is exact arithmetic in `Q(sqrt(7))`.  It uses no floating-point
decision, optimizer failure, tolerance, or transcendental function.  A
candidate nonedge remains unconstrained and may still be a unit pair.

## Frozen parent boundary and quantifier

The input is the independently verified v8 current residue.  Its K7 class is
the ordered list

```text
316173 2581209 3648882 3729907 3935560 3936310
3936435 3945490 3945555 3945557 3945564 3947605
```

with ordered hash

```text
e545ea5a8fc9691ed98e885919e149f84dc7e8b4796d0fe9e6cb3a519df6fd8c
```

For graph `3936435`, fix the required K7 seed

```text
[1,3,6,8,11,13,17],        seed mask 141642.
```

The outside order is

```text
[0,2,4,5,7,9,10,12,14,15,16,18].
```

The producer reconstructs the seed defects, disjoint-defect Lorentz graph,
zero-factor eligibility, and all 255 raw eligible covers.  The independent
checker enumerates the same covers in combination order and then sorts them;
their ordered-list hash is

```text
cf605db4a5fed03da735e569687d1e84d2bdde24c45d95a9fc4b61490e1a491d.
```

The independently checked one/two-free-star boundary has already eliminated
253 of these covers.  Exactly two current covers remain:

```text
Z = 0,       one propagated support family;
Z = 2048,    one propagated support family.
```

For `Z=2048`, outside local vertex 11 is the sole zero-factor vertex and its
actual support is mask `63`.  It is not one of the five vertices used below.
Consequently every vector in the new obstruction has a nonzero normalization
factor in both covers.  The propagated masks of the five vertices are the
same in both families.

This accounts for all quantifiers: a geometric realization supplies one of
the 255 eligible covers and one of its actual-support branches.  The parent
certificate eliminates 253 covers; the two remaining one-family covers are
both eliminated here.  Since every realization of the graph must realize
this required K7 seed, one impossible seed rejects the graph.

## Normalized defect identity

For the fixed centered regular unit K7, write

```text
u_i(x) = ||x-q_i||^2-1,
c_x    = 1+sum_i u_i(x).
```

For a nonzero-factor vertex (`c_x != 0`) normalize

```text
w(x) = sqrt(7) u(x) / c_x.
```

A required unit edge `xy` gives

```text
w(x).w(y) = 1.                                      (1)
```

The diagonal simplex identity gives

```text
||w(x)||^2 = 1 + (sum_i w_i(x)-sqrt(7))^2.          (2)
```

The support of `w(x)` equals the support of `u(x)` and is only known to be a
subset of its propagated mask.  The proof uses this in the safe direction:
coordinates outside a propagated mask vanish; coordinates allowed by the
mask are never required nonzero merely because they are allowed.

## The five-vertex support pattern

In outside-local labels, use hubs `A,D` and petals `B,C,E`:

```text
vertex       local    global    propagated support superset
A              0        0       {2,3,4,5}       mask 60
B              1        2       {1,3}           mask 10
C              3        5       {1,2}           mask  6
D              6       10       {0,2,3,5}       mask 45
E              9       15       {1,5}           mask 34
```

All ten pairs among these five vertices are required graph edges.  Thus they
form a required K5 and every pair has normalized inner product one by (1).
No candidate nonedge supplies an equation.

The three petals intersect pairwise only in center coordinate 1.  Write their
center coordinates as `x,y,z`.  Their three required edges give

```text
xy = xz = yz = 1.
```

Therefore

```text
x=y=z=epsilon,       epsilon in {-1,+1}.             (3)
```

Each petal has one remaining possible coordinate, on leaves 3, 2, and 5
respectively.  Equation (2) gives the same leaf value `lambda` for all three:

```text
lambda = (sqrt(7)-epsilon)/2.                         (4)
```

Each hub meets each petal only at the corresponding leaf.  The six hub-petal
edge equations therefore force the hub values on leaves 2, 3, and 5 all to
equal `1/lambda`.  The two hub masks intersect in exactly those three leaves,
so the required hub-hub edge gives

```text
3/lambda^2 = 1,       hence lambda^2=3.               (5)
```

But (4) gives

```text
lambda^2 = 2 - epsilon sqrt(7)/2.
```

Equality with 3 would imply `epsilon sqrt(7)=-2`, and squaring would imply
`7=4`.  This contradicts both choices of `epsilon` and completes the proof.

Notice that extra hub coordinates 4 and 0 are harmless: the two hubs do not
share them.  If an allowed coordinate happens to vanish, every displayed
deduction remains valid.  Required edges themselves force every coordinate
that is divided by to be nonzero.

## Independent checking

The verifier imports neither the producer nor a discovery module.  It:

1. hash-checks the v8 manifest/checker and parent star report/checker;
2. reconstructs the graph, seed, defects, Lorentz edges, and eligibility;
3. re-enumerates all raw covers in a different order;
4. verifies the two archived current covers and independently remaps each
   nonzero-local propagated-mask list;
5. reconstructs the required K5 and all mask intersections;
6. checks the algebra by symbolic polynomial reduction modulo `r^2-7`, rather
   than using the producer's `Q7` pair-arithmetic class;
7. checks the known realizable 18-point control (which has zero K7 seeds) and
   an optional-hub-nonedge control.

The parent independent checker boundary covers 19,932 eligible covers and 88
then-current passing support families.  This package does not silently rerun
or replace that established certificate; it pins its report and independent
verification by SHA-256 and reconstructs the complete 255-cover quantifier of
the one seed used here.

## Source freeze and reproduction

The theorem report may be generated only after all four package sources are
committed on `codex/dimension6`.  The producer checks that every working
source equals its blob in the launch commit and rejects tracked or staged
worktree dirt.  Unrelated untracked exploratory files are recorded in the
launch porcelain hash but do not alter the committed source boundary.

After the source-only commit:

```bash
python3 -m unittest -v test_d6_k7_triple_petals_3936435_increment.py

python3 build_d6_k7_triple_petals_3936435_increment.py

python3 verify_d6_k7_triple_petals_3936435_increment.py \
  --report-sha256 <sha256-of-production-report>
```

The verifier output should be committed separately from the source freeze.
The final trust boundary consists of exact Python integer/Fraction arithmetic,
Sympy rational polynomial reduction in the independent checker, Git blob and
SHA-256 bindings, and the four pinned parent artifacts.  IEEE-754 and `libm`
assumptions do not enter the rejection.

## Scope

This package removes one graph from the current K7 class.  It does not assert
that any survivor is realizable and does not settle the remaining K7 class,
the K6-only class, or the dimension-six extremal value.
