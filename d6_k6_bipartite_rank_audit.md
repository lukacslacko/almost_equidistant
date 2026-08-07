# Exact K6 bipartite Lorentz-component rank refinement

This note proves and profiles the filter in
`d6_k6_bipartite_rank_reference.py`.  It is layered on the pure K6 Lorentz
filter and the joint actual-support refinement.  All decisions use finite
bit masks and integer zero forcing; there is no floating-point or numerical
rank test.

## Statement

Fix a required K6 seed, its allowed defect masks `D_x`, and the graph `L` on
the thirteen outside vertices: `xy` is an edge of `L` exactly when `xy` is a
required unit edge and `D_x` and `D_y` are disjoint.  For one possible actual
zero-factor set

```text
Z0 = {x : ell_x=(c_x,z_x)=0},
```

consider a connected bipartite component of `L-Z0`, with sides `A,B`.  If
`Zf(H)` is the ordinary zero-forcing number of a graph `H`, every realization
in `R^6` necessarily satisfies

```text
|Z0| + (|A|-Zf(G[A])) + (|B|-Zf(G[B])) <= 6.       (1)
```

The condition is imposed separately on every bipartite component.  Spans
belonging to different components need not be orthogonal and are not added.

## Proof

The K6 coordinates give, for every outside pair,

```text
M_xy = <ell_x,ell_y>_L/6 - u(x) dot u(y),
<ell_x,ell_y>_L = c_x c_y-z_x z_y.
```

A required edge has `M_xy=0`.  If `xy` is a candidate graph nonedge, then
`D_x` and `D_y` are disjoint: a common seed coordinate would make that seed
vertex, `x`, and `y` an independent triple, contrary to `alpha(G)<=2`.
Consequently `u(x) dot u(y)=0` for every candidate nonedge.  This deduction
uses only *allowed* supports.  It never assumes that an allowed coordinate is
actually nonzero or that a candidate nonedge has non-unit distance.

All factors in `L-Z0` are nonzero.  Along an `L` edge their Lorentz product
is zero.  In the nondegenerate two-dimensional Lorentz plane, taking the
orthogonal line is an involution.  Starting from one vertex of a connected
bipartite component therefore gives one projective line on `A` and its
orthogonal line on `B`.

There are two cases.

### Generic direction

If the first line is not lightlike, neither is its orthogonal line.  Write

```text
ell_x=lambda_x ell_A  (x in A),
ell_y=lambda_y ell_B  (y in B),
```

where every scalar is nonzero, both Lorentz self-products are nonzero, and
`<ell_A,ell_B>_L=0`.

For two vertices of `A`, a required edge gives

```text
u(x) dot u(y)
  = lambda_x lambda_y <ell_A,ell_A>_L/6 != 0,
```

while a candidate nonedge gives dot product zero by disjoint allowed masks.
Thus the Gram matrix of the `u(x)`, `x in A`, has off-diagonal graph exactly
`G[A]`.  The identical statement holds on `B`.  The ordinary zero-forcing
minimum-rank theorem yields

```text
dim span u(A) >= |A|-Zf(G[A]),
dim span u(B) >= |B|-Zf(G[B]).
```

For a cross pair in `A x B`, a required edge has zero dot product because
the Lorentz directions are orthogonal; a candidate nonedge has zero dot
product because its allowed masks are disjoint.  Hence the two spans are
orthogonal.

For `z in Z0`, the diagonal identity gives `||u(z)||=1`.  A required pair
with any other vertex is orthogonal by the unit equation, while a candidate
nonedge is orthogonal by disjoint masks.  Therefore the `u(z)` are an
orthonormal family, orthogonal to both component-side spans.  These three
mutually orthogonal subspaces lie in `R^6`, proving (1).

### Lightlike direction

The Lorentz orthogonal complement of a lightlike line is the same line.
Connectivity then puts every factor in `A union B` on one light ray, so all
their pairwise Lorentz products vanish.  Required edges again give zero
defect-vector dot products and candidate nonedges give zero by disjoint
masks.  The diagonal identity now gives `||u(x)||=1` for every component
vertex.  Together with `Z0`, the whole component is an orthonormal family:

```text
|Z0|+|A|+|B| <= 6.
```

This is stronger than (1), so (1) is valid in both Lorentz cases.  An
isolated component vertex is harmless: its one-vertex graph has zero-forcing
number one and contributes zero to (1).

## Exact quantifiers and combination with earlier filters

A zero factor has norm-one defect vector with coordinate sum `-1`.
Distinctness from the six seed vertices forces its actual support to have at
least three coordinates.  Hence `Z0` is a subset of the vertices with
`|D_x|>=3`; orthonormality gives `|Z0|<=6`, and the allowed masks of `Z0`
must admit a matching into the six coordinates.

The implementation enumerates **every** such eligible subset.  It does not
use the old optimization that omitted zero candidates in initially
bipartite components.  That optimization was safe while bipartite components
were ignored, but is not safe here: deleting a zero vertex can split or
change the exact component whose two side ranks occur in (1).

For each `Z0` passing (1), every non-bipartite component of `L-Z0` is assigned
to one of the two light rays, modulo only global interchange of their names.
The two old allowed-mask matchings and the joint actual-support intersection
CSP are then applied unchanged.  A seed passes if at least one complete
choice survives.  A graph is rejected only if at least one of its required
K6 seeds exhausts every choice.

## Controls

`test_d6_k6_bipartite_rank.py` checks:

- exact zero-forcing numbers of an isolate, paths, a cycle, and a clique;
- bipartite and odd-cycle component classification;
- a fixed two-side kernel example requiring seven dimensions;
- a regression in which the only passing choice puts a zero in an initially
  bipartite component;
- the known realizable 18-point graph over all 32 of its K6 seeds;
- the fixed complete-corpus rejection witness at index `461363`.

`verify_d6_k6_bipartite_rank.py --full` checks all bound source/input hashes,
re-sums every aggregate from the per-graph decisions, and independently
recomputes all 1,098 graph decisions.  A separate serial production run was
also compared structurally with the 11-worker report after removing only
runtime/worker metadata; every decision and exact counter agreed.

## Complete profile

The exact input extractor removes the eight actual-support rejections from
the hashed 1,106-graph pure-Lorentz residue, leaving all 1,098 current K6
survivors.  On macOS 14.5, arm64, Python 3.11.15:

| quantity | exact value |
|---|---:|
| input graphs | 1,098 |
| K6 seeds checked | 35,954 |
| all-eligible `Z0` subsets checked | 36,515 |
| `Z0` subsets failing (1) | 31 |
| bipartite components checked | 247,793 |
| graphs rejected | 1 |
| graphs surviving | 1,097 |
| 11-worker wall time | 0.640 s |
| serial wall time | 4.536 s |

The new rejection is corpus index `461363`, with one impossible K6 seed
`[0,2,5,8,17,18]`.  Its eligible zero candidates are outside vertices
`4,10,15`, so all eight possible `Z0` subsets are explicitly checked.  The
same bipartite component survives every deletion because none of those three
vertices lies in it.  At `Z0=empty`, its sides have sizes four and six,
zero-forcing numbers one and two, and hence rank lower bounds three and four:

```text
0 + 3 + 4 = 7 > 6.
```

Adding any eligible zero only increases the left side, so all eight choices
fail exactly.

Commands:

```text
python3 extract_d6_k6_bipartite_rank_input.py
taskpolicy -a python3 d6_k6_bipartite_rank_reference.py \
  d6_k6_bipartite_rank_input.json --workers 11 \
  --output d6_k6_bipartite_rank_report.json
taskpolicy -a python3 d6_k6_bipartite_rank_reference.py \
  d6_k6_bipartite_rank_input.json --workers 1 \
  --output /tmp/d6_k6_bipartite_rank_serial_report.json
python3 -m unittest -v test_d6_k6_bipartite_rank.py
python3 extract_d6_k6_bipartite_rank_input.py --check
python3 verify_d6_k6_bipartite_rank.py --full
```

Tracked artifact hashes:

```text
d6_k6_bipartite_rank_input.json
  dfee979e6ca2759e80cf803c3054c8954d98543875d467644b91ee434c22e845
d6_k6_bipartite_rank_reference.py
  e760fef74c423a89f1e1f9d08909fa1223e0389900cf3d45b0140ee296dea922
d6_k6_bipartite_rank_report.json
  ede04a71a7f9b6901bbc1e4c198a117797bfd90d4d1366912ca9efd1ca559bdc
```

The report itself binds the earlier pure-Lorentz/support artifacts and both
imported dependency source hashes back to the complete 3,971,787-graph
corpus SHA-256.
