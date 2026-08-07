# Correlated pinned-product one-free-edge obstruction

## Scope

This layer extends the frozen singleton one-free-edge theorem without
changing it.  Let `xy` be a required unit edge in a propagated nonzero-factor
`K7` branch.  Assume that each endpoint has exactly one allowed coordinate
not pinned by a nonbipartite singleton-intersection coordinate component,
and that it is the same coordinate `i` at both endpoints.  Unlike the first
lemma, the two propagated masks may also overlap in pinned coordinates.

For a mask of size `k`, pinned-sign sum `P`, and last component `z`, the exact
diagonal identity still gives

```text
z = F(k,P)
  = (9-k+P^2-2 sqrt(7) P) / (2(sqrt(7)-P)).
```

If `S` is the sum of the products of endpoint signs over all pinned shared
coordinates, the required edge equation becomes

```text
F(k,P) F(l,Q) + S = 1.                           (1)
```

The branch is impossible if no assignment of the pinned component signs
satisfies (1).

## Why the signs are correlated

For coordinate `j`, an edge in its singleton-intersection coordinate graph
gives

```text
w_j(a) w_j(b) = 1.
```

An odd cycle pins every value in its connected component to a sign.  Along
each edge the two signs are equal, so the entire nonbipartite component has
one common unknown sign.  Therefore sign variables are keyed by

```text
(coordinate, connected component).
```

If both target endpoints lie in the same pinned component for a shared
coordinate, that coordinate contributes `+1` to `S`.  If they lie in
different components, the contribution may be either sign.  The checker
enumerates the component variables, rather than treating `P`, `Q`, and `S`
independently.  This retains every correlation justified by the required
singleton-intersection edges while deliberately ignoring other possible
global relations.  Ignoring extra relations only enlarges the domain and is
safe for rejection.

There are at most six pinned coordinates at each endpoint, hence at most 12
local component variables and 4096 assignments.  Every calculation uses
`Fraction` pairs for `a+b sqrt(7)`.

## Complete integer-product audit

Because `S` is an integer, (1) requires `F(k,P)F(l,Q)` to be a rational
integer.  Exhausting all 28 one-free types gives the following complete list,
up to exchanging the two factors:

```text
F(1,0)   F(4,+/-3) =  4
F(3,-2)  F(3,2)    = -1
F(3,0)   F(4,+/-3) =  3
F(4,+/-3)F(4,+/-3) =  7
F(4,+/-3)F(5,0)    =  2
F(4,+/-3)F(7,0)    =  1
F(5,-2)  F(5,-2)   =  4
F(5,2)   F(5,2)    =  4
F(5,-2)  F(5,2)    = -4
```

The size-one/four product `4` and the opposite extreme size-three product
`-1` are easy to miss; both are included in the source test.  Products not
listed here have a nonzero `sqrt(7)` coefficient or are nonintegral
rationals, and therefore cannot equal `1-S`.

If all pinned signs at the two endpoints are first allowed to vary
independently, while respecting their occurrence in `P`, `Q`, and `S`, the
only locally compatible mask-size/intersection patterns have sizes four and
five and either one or three pinned shared coordinates.  They require

```text
P = +/-3,  Q = 0,  S = -1
```

up to exchanging endpoints.  Actual component identifications can remove
even these possibilities.  The production checker does not rely on this
coarse corollary: it directly enumerates the exact component variables for
each candidate edge.

## Certificate and semantics

`d6_k7_correlated_one_free_edge.py` records:

- the target required edge, both masks, and their common free coordinate;
- a nonbipartite component/conflict certificate for every pinned coordinate
  at each endpoint;
- the number of distinct component-sign variables and the complete number
  of assignments exhausted.

Verification reconstructs the coordinate graphs, validates every component,
reconstructs shared component identities, and repeats the exhaustive exact
equation check.  The layer requires at least one pinned shared coordinate,
so its certificates are disjoint in form from the frozen singleton-overlap
lemma.  Candidate nonedges remain unconstrained and allowed mask entries may
be zero until the exact equations force them nonzero.

## Tests

```sh
python3 -m unittest -v test_d6_k7_correlated_one_free_edge.py
```

Controls cover the complete integer-product table, an incompatible branch,
JSON round-trip and tampering, the surviving independent-sign four/five
pattern, and the same pattern becoming impossible when a component bridge
forces its pinned shared product to `+1`.
