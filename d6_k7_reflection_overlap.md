# Exact propagation through the K7 overlap complex

## Result of the pilot

The required-unit-edge graph determines exact rational coordinates throughout
every connected component of its graph of `K7` cliques joined across `K6`
facets.  A coordinate conflict, forced collision, or wrong distance on a
required edge is therefore a sound obstruction to a realization by distinct
points in `R6`.

The deterministic 1,024-graph pilot on the current 12,941-graph exact `K7`
residue found no new rejection.  The structure needed by the filter was very
rare:

- 106 sampled graphs had more than one required `K7`, but only graph 3,784,405
  had two `K7`s sharing a `K6`;
- that graph had three `K7`s in total, one overlap edge, and one nontrivial
  overlap component containing two cliques and eight vertices;
- its forced reflection was consistent, so all 1,024 sampled graphs survived.

This low hit rate makes a substantial extension of this filter unattractive.
A full pass is nevertheless warranted after the current CPU campaign because
the pilot's exact propagation took 0.083 seconds on one core; completing the
12,941-graph accounting should cost only seconds, not a material compute run.

## Facet-reflection lemma

Let `p_0,...,p_6` be seven distinct pairwise unit points in `R6`.  They are the
vertices of a full-dimensional regular simplex.  Fix an index `i` and let

```text
F_i = {p_j : j != i},             h_i = (1/6) sum_{j != i} p_j.
```

The six points of `F_i` are affinely independent and span a hyperplane.  Their
circumcenter in that hyperplane is `h_i`, and their squared circumradius is
`5/12`.  Subtracting the equations `||x-p_j||^2=1` shows that the orthogonal
projection of any common unit neighbour `x` of `F_i` is `h_i`.  Its squared
normal displacement is

```text
1 - 5/12 = 7/12.
```

There are consequently exactly two common unit neighbours.  One is `p_i` and
the other is its reflection through the facet hyperplane,

```text
p_i' = 2 h_i - p_i = (1/3) sum_{j != i} p_j - p_i.       (1)
```

Thus, if two graph `K7`s share a `K6` facet, their differently labelled
apices must occupy the two positions in (1): injectivity excludes using the
same position.  This conclusion uses only required unit edges and distinctness.

## Rational affine coordinates

Choose one `K7` in an overlap component as root and give its vertices the
affine coordinates `e_0,...,e_6`.  Formula (1) maps rational affine
coordinates to rational affine coordinates and preserves coordinate sum one.
It therefore propagates an exact coordinate for every vertex reached through
successive `K6` overlaps.

The metric also has a rational form.  If `a` and `b` are two affine coordinate
vectors and `delta=a-b`, then `sum_i delta_i=0`.  Since every off-diagonal
root-simplex squared distance is one,

```text
||sum_i delta_i p_i||^2 = (1/2) sum_i delta_i^2.          (2)
```

One way to derive (2) is to expand
`sum_ij delta_i delta_j ||p_i-p_j||^2`: its value is
`-sum_i delta_i^2`, while the zero-sum condition makes the same expansion
equal to `-2 ||sum_i delta_i p_i||^2`.

All propagation and distance decisions can hence be checked with rational
integer arithmetic.  No floating-point tolerance enters the filter.

## Sound contradictions and certificates

For each overlap component, the implementation breadth-first propagates from
one root `K7`.  Each successful trace record names two required `K7`s, their
six-vertex intersection, and the exchanged apex labels.  A terminal record is
accepted only in one of three cases:

1. **coordinate conflict:** a previously reached graph vertex is assigned a
   different rational coordinate by another overlap path;
2. **collision:** a newly derived vertex has the same rational coordinate as
   a different previously reached graph vertex;
3. **wrong required distance:** equation (2) gives a value different from one
   for an edge present in the required-unit graph.

The root simplex determines an affine frame up to Euclidean isometry, and the
facet lemma makes every trace step compulsory.  Each terminal condition is
therefore impossible in an injective unit-edge realization.  Conversely,
survival is only a necessary-condition result: it does not assert that the
graph is realizable.

Candidate nonedges are never tested.  In particular, the checker does not
require a candidate nonedge to have nonunit distance.  A certificate remains
valid after adding more required edges, so it is also a reusable unit-edge
subgraph obstruction.

`d6_k7_reflection_verify.py` independently replays the certificate format.  It
does not import the generator and separately checks the graph, every clique
and overlap, formula (1), the trace invariants, and the terminal rational
contradiction.

## Relation to the earlier reflection rule

The original `profile_d6.c` rule fixes one `K7` and looks at all outside
vertices adjacent to exactly six seed vertices.  Such a vertex creates a
second `K7` across the corresponding facet.  It rejects two vertices of the
same type (a collision) and a required edge between vertices of different
types (their squared distance is `16/9`).

The present filter contains that one-step argument and then repeats it from
every derived `K7`.  It additionally checks consistency around overlap cycles
and all required edges among coordinates in the same propagated component.
The current residue has already passed the one-step rule, explaining why
nontrivial overlap components are scarce.

## Controls and reproducibility

The nine controls include:

- a single `K7`, an honest reflected pair, and a realized two-step chain as
  positive controls;
- a required edge between opposite apices, three extensions of one facet, and
  an inconsistent overlap cycle as negative controls;
- direct checks of the exact distances `7/3` and `16/9`;
- acceptance of an unconstrained nonedge and rejection of a tampered
  certificate;
- replay of every negative certificate by both the generator-side checker and
  the independent checker.

Run them with

```text
python3 -m unittest -v d6_k7_reflection_overlap_test.py
```

The pilot command was

```text
python3 d6_k7_reflection_overlap.py \
  --sample-size 1024 \
  --output d6_k7_reflection_overlap_sample.json
```

The input and deterministic selection are pinned in the report:

```text
.runs/d6_k7_rank_survivors.json
  SHA-256 7a0a350142930e4e13830ea44c4fd217f25aeb07bbd6286e51e316f72712c4f9
d6_k7_positive_dual_selection.json
  SHA-256 814c80651746ce05048f72f0cfa7e49a921554abaf167bd5c8cf4063a34d35c0
```

Frozen pilot source hashes are:

```text
d6_k7_reflection_overlap.py
  SHA-256 13e6b0a4d6780550d02d9a5c370c5c042684060d8696191ed98776a14a2f9582
d6_k7_reflection_verify.py
  SHA-256 86ae0104bda81ab34deaf0f8db2e87136d76b9594c8026604735a6bf8c67016f
d6_k7_reflection_overlap_test.py
  SHA-256 a5c33d29b32fd8f9d00f68a141671bfbbe61591193d0b14d648afe5db3b59271
d6_k7_reflection_overlap_sample.json
  SHA-256 92c62866f60807cbdd33d67d0ef9ec2cd1ebf755a54062d96378e8117cc7db69
```

The proof decisions trust Python arbitrary-precision integers, `Fraction`,
JSON parsing, and the audited generator/checker logic.  Platform timing fields
are informational and are not certificate inputs.
