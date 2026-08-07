# The 18-vertex deletion corpus of the d=6 residue

Deleting any vertex from a realizable 19-point support leaves a realizable
18-point support.  This makes induced deletions a natural source of smaller,
reusable obstructions and, in the opposite direction, realizations that can
be tested for one-vertex grow-back.  This package constructs the deletion
corpus but does not reject a graph geometrically.

## Exact corpus

The source is the independently verified 960-graph v4 residue.  All 19
deletions of every parent are retained, including the deleted vertex's exact
required-unit attachment mask.

```text
19-vertex parents                                      960
  K7-containing                                        155
  K6-only                                               805
labeled deletion occurrences                        18,240
unique unrooted 18-vertex classes                   12,712
  occurring below a K7 parent                        2,288
  occurring below a K6-only parent                  10,486
  overlap between those two class sets                  62

unique classes containing K7                         1,616
unique classes containing K6 but no K7              11,096
unique classes containing no K6                          0
```

Unrooted isomorphism is used only to share work.  A grow-back certificate is
rooted data: the same 18-vertex graph can occur with inequivalent attachments
of the nineteenth vertex.  The manifest therefore stores every parent index,
deleted label, attachment mask, and attachment degree rather than replacing
the 18,240 occurrences by 12,712 bare labels.

## Matches to the standard 18-point coordinates

Let `U` be the full unit graph of the 16 odd half-cube vertices and two
orthogonal poles.  Exact containment testing gives:

```text
unrooted deletion classes whose required edges embed in U       14
corresponding labeled deletion occurrences                       39
corresponding 19-vertex parents                                  16
```

These 14 classes are exact positive 18-point supports: their required edges
are unit in the displayed coordinates.  This says nothing yet about the
deleted vertex.  In particular, none of the 16 parents is declared realizable
or impossible until its stored attachment is tested against every relevant
embedding/realization of the deletion.

There are also 1,960 labeled edge-minimal almost-equidistant supports carried
by the standard coordinates, in seven unrooted types.  Exactly 26 deletion
occurrences from 14 parents equal one of those minimal supports.  The broader
14/39/16 containment figures above are the ones relevant to optional
nonedges, because a candidate nonedge is allowed to be unit.

## Why the Clebsch graph occurs

Write the 16 base points as odd-parity vectors in `{+1,-1}^5`, scaled by
`1/sqrt(8)`.  Two such vectors differ in an even number of coordinates:
Hamming distance two gives squared Euclidean distance one, while Hamming
distance four gives squared distance two.  Thus the **non-unit** graph of the
base has degree five.  After translating the odd parity coset to
`F_2^4`, its five difference vectors are

```text
e1, e2, e3, e4, e1+e2+e3+e4,
```

which is the standard Cayley model of the Clebsch graph.  Equivalently it is
strongly regular with parameters `(16,5,0,2)`.  The base unit graph is the
degree-ten complement of the Clebsch graph, not the Clebsch graph itself.

This also makes the standard-coordinate containment test exact.  The
complement of `U` is a Clebsch component plus the edge joining the two poles.
For a candidate support `G`, its complement is triangle-free.  Once a
candidate nonedge pair is chosen as the poles, every Clebsch edge must occur
among the remaining candidate nonedges.  The Clebsch graph is maximal
triangle-free: each of its nonadjacent pairs has two common neighbours, so
adding any edge creates a triangle.  Consequently the remaining 16-vertex
complement must equal the Clebsch graph.  Conversely, that equality and the
pole nonedge give an explicit embedding of every required edge of `G` into
`U`.  The manifest stores the full permutation witness.

This recognition is conditional on the displayed standard coordinates.  It
does not assert that every 16- or 18-point almost-equidistant set has this
form.

## Verification and next use

```text
d6_residue_18_deletions.json
  9d08f9ec579434c9601ad3908bd2ae51f10e09e9d06a7f38fbc25794a634b732
d6_residue_18_deletions_verification.json
  50ab5d88441b2869e60a05f3e1485db0b25008c3a816960a20a3bb21c1c4defd
```

The independent checker imports neither the builder nor a geometric
production kernel.  It reconstructs all 18,240 deletions, repeats canonical
labelling, checks all attachment masks and explicit standard embeddings, and
rebuilds the seven minimal-support types.

Immediate uses are:

1. certify small 18-vertex nonrealizable cores once and reject every parent
   occurrence containing one;
2. numerically triage the 12,712 classes for likely flexible/rigid
   realizations, without using optimization failures as proofs;
3. for exact 18-realizations, solve the stored one-point sphere-intersection
   problem for each rooted attachment;
4. mine repeated failed attachments as smaller rooted algebraic or interval
   obstructions.
