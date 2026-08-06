# Audit of the K6 actual-support intersection refinement

This note documents the exact scope of `d6_k6_support_reference.py`.  It is a
prototype layered on top of the completed odd-component/two-light-ray filter;
it does not alter that milestone's production artifacts.

## Necessary condition

Fix a required K6 seed and one possible zero-factor set `Z0`.  After deleting
`Z0` from the disjoint-defect required-edge graph `L`, each connected
non-bipartite component has its nonzero Lorentz factors on one of the two
lightlike projective lines.  Components assigned the same line, together with
`Z0`, give an orthonormal family of nonzero `u` vectors in `R^6`.

For each relevant vertex choose its actual support

```text
S_x = {i : u_i(x) != 0} subseteq D_x.
```

The prototype exhaustively imposes:

1. `|S_x|>=3` for `x in Z0`.  Here `ell_x=0` gives
   `||u(x)||^2=1` and `sum u_i(x)=-1`; support at most two forces
   `u(x)=-e_i`, which duplicates a seed point.
2. `S_x` is nonempty for a nonzero light-ray vertex, since
   `||u(x)||^2=1`.
3. Within either orthonormal bin, `|S_x intersect S_y|` is zero or at least
   two.  If the intersection were exactly `{i}`, then the dot product would
   be the single nonzero term `u_i(x)u_i(y)`, so it could not vanish.
4. Each bin's family of *actual* supports admits a matching into the six seed
   coordinates.  This is a necessary structural-rank condition for a linearly
   independent, hence orthonormal, vector family.

The same actual support is chosen for a zero-factor vertex in both bins.
Choosing separate supports in the two bins would be unsound.  Conversely, no
intersection condition is imposed between vertices on opposite lightlike
lines: their Lorentz product need not vanish.

There is no hidden coefficient assumption in item 3.  Intersection of size at
least two is merely retained as a necessary possibility; the prototype does
not claim that suitable coefficients cancel.  Thus passing is not a
realization claim, while exhaustive failure is a sound obstruction.

Candidate graph nonedges remain unconstrained.  Every `D_x` is only an allowed
support, and the exhaustive search permits every subset with the required
minimum cardinality.

## Quantifiers and safe optimization

For every K6 seed the implementation exhausts, in order:

1. every eligible `Z0` of size at most six;
2. every two-coloring of the non-bipartite components of `L-Z0`, modulo only
   global exchange of the two color names;
3. every joint actual-support assignment satisfying the conditions above.

A graph is rejected if any one K6 seed exhausts all choices.

Vertices in an initially bipartite connected component of `L` are omitted
from the `Z0` candidates.  This is safe because deleting vertices preserves
bipartiteness, while removing such a vertex from a proposed `Z0` only removes
one vector from both constrained bins.  The optimization retains the whole of
every initially non-bipartite component, including tree attachments and
articulation vertices away from its odd-cycle core.

The six-dimensional capacity after fixing `Z0` is already encoded: either bin
contains `Z0` and has at most six vectors, equivalently at most `6-|Z0|`
nonzero light-ray vectors.  Unlike the K7 rank refinement, vectors in the two
opposite K6 light-ray bins are not mutually orthogonal, so no additional total
rank cap across both colors follows.

## Controls and measurement

`test_d6_k6_support.py` checks:

- the realizable 18-point configuration over all 32 K6 seeds;
- the exact one-coordinate noncancellation rule;
- conservative acceptance when two coordinates overlap;
- allowed-mask matching being strictly weaker than actual-support feasibility;
- a case that would pass if `Z0` supports were chosen independently per bin
  but correctly fails with a shared support;
- all eight fixed real-corpus rejection witnesses.

The complete input is the 1,106-graph residue reconstructed by
`extract_d6_k6_support_residue.py` from the hashed full decision TSV and
candidate corpus.  The 12-worker run checked 36,221 K6 seeds and rejects eight
graphs (ten impossible seeds), leaving 1,098.  It considered 36,538 `Z0`
subsets, 38,249 symmetry-reduced colorings, and 103,879 actual-support DFS
nodes.  No floating-point arithmetic participates in any decision.
