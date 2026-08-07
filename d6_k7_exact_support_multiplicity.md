# Exact defect-type multiplicity after a K7 seed

This follow-on note records a general exact support theorem suggested by the
one- and two-defect value layer.  It is not yet implemented and does not
change the frozen sparse-value sample report.

## The theorem

Fix a regular unit K7 seed in `R^6`.  Let `S` be a nonempty subset of its
seven vertices, with `|S|=k`.  Among all outside points, at most `k` distinct
points can have actual defect support exactly `S`.

Here "actual" is essential: a zero defect coordinate means that the
corresponding seed distance really is one.  Candidate nonedges remain
optional and are not assigned a distance.

## Proof

Any two points `x,y` of exact type `S` are both non-unit from every seed
vertex in `S`.  Choosing one such seed vertex and applying the
almost-equidistant condition to that triple shows that `x,y` must be a unit
pair.  Thus all points of type `S` are pairwise unit.

For `k<=6`, these points are also common unit neighbours of the complementary
seed clique `K_(7-k)`.  The clique-link sphere lies in an affine space of
dimension `k` and has radius squared

```text
R_k^2 = (8-k)/(2(7-k)).
```

Suppose there were `k+1` points of exact type `S`.  Being pairwise unit, they
would form a regular unit `k`-simplex and span the whole `k`-dimensional link
space.  Its circumradius squared is

```text
r_k^2 = k/(2(k+1)).
```

A full-dimensional simplex has a unique circumsphere in that affine space,
so it would require `R_k^2=r_k^2`.  But after putting the two fractions over
a common positive denominator, the numerator of their difference is

```text
(8-k)(k+1)-k(7-k)=8.
```

They are never equal.  Hence there are at most `k` points.

For `k=7`, points of exact type `S` are still pairwise unit, and a unit clique
in `R^6` has at most seven vertices.  This again gives the bound `k`.

The cases `k=1,2` recover the existing one-defect uniqueness and fixed
two-defect multiplicity-at-most-two rules without separate coordinate
elimination.

## Proposed exact CSP extension

For a fixed labeled zero-factor support family, singleton orthogonality
propagation gives a nonzero-factor mask `E_x` containing the unknown actual
support `S_x`.  A next finite CSP can use the domains

```text
empty != S_x subseteq E_x
```

and impose simultaneously:

1. for every exact type `S`, its multiplicity among the already fixed
   zero-factor supports and assigned nonzero-factor supports is at most
   `|S|`;
2. actual supports of a required edge inside the nonzero-factor set
   intersect, since their normalized dot product is one;
3. the existing Mobius cycle, one-defect path, branch, and parallel-type
   rules whenever an assigned support has size one or two.

A cheap first stage would enumerate only vertices with `|E_x|<=3`.  Each
three-bit mask has only seven nonempty submasks.  Vertices with larger masks
can be left as unconstrained escape variables, preserving the sound
direction.  If that pilot is useful, capacity-aware backtracking or a flow
relaxation can extend it before any full-residue campaign.

