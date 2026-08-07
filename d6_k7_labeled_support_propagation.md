# Exact labeled-support propagation after a K7 cover

This note proves the refinement implemented independently in
`d6_k7_support_propagation.py`.  The frozen K7 rank reference is not changed.

## Setup

Fix a required K7 seed and an eligible zero-factor cover `Z`, with `|Z|<=3`.
Write `N` for the remaining outside vertices.  For every outside vertex `x`,
`D_x` is the seven-bit candidate defect mask and `S_x subseteq D_x` is its
actual (nonzero-coordinate) support.

The simplex identity for a required unit pair is

```text
u(x) dot u(y) = c_x c_y.
```

Here `c_z=0` for `z in Z`, while `c_x!=0` for `x in N`.  The existing support
lemma gives `|S_z|>=3`; the frozen support solver's matching and
no-singleton-intersection rules are necessary for the labeled family
`(S_z:z in Z)`.

## Singleton deletion

For every `z in Z` and `x in N`, the vectors `u(z)` and `u(x)` are orthogonal.
If `zx` is a required edge this is the displayed identity.  If it is a
candidate nonedge, then `D_z` and `D_x` are disjoint: otherwise a seed vertex
in their intersection would form an independent triple with `z,x`.

Suppose a current over-approximation `E_x` to `S_x` has

```text
E_x intersect S_z = {i}.
```

If `i` belonged to `S_x`, the dot product would contain exactly the one
nonzero term `u_i(x)u_i(z)` and could not vanish.  Therefore `i` is absent
from `S_x` and may be deleted from `E_x`.  Repeating this deletion to a fixed
point is sound because every intermediate mask remains a superset of the
unknown actual support.  The simultaneous iteration used by the checker is
canonical and terminates after at most seven coordinate deletions per vertex.

An empty propagated mask is impossible for `x in N`.  Moreover, if `xy` is a
required edge inside `N`, then the propagated masks cannot be disjoint:
disjoint actual supports would give `u(x) dot u(y)=0`, contradicting
`c_x c_y!=0` in the unit-pair identity.

## Sharpened rank caps

For a fixed labeled support family, let `E_x` be the propagated masks and set

```text
nu_N = term-rank(E_x : x in N),
nu_T = term-rank((S_z : z in Z), (E_x : x in N)).
```

Every actual support is contained in its corresponding mask, so these are
upper bounds on the relevant column ranks.  Consequently the existing bounds
remain valid with the sharper values

```text
U_K = min(|N|, nu_N, nu_T-|Z|, 7-|Z|),
U_B = min(|N|, U_K+1, 8-|Z|).
```

The checker reruns, at these assignment-specific caps, the exact K zero
forcing bound, component/inertia B bound, positive-definite clique rule,
perpendicular-degree rule, rational basis-kernel rule, and saturating-clique
Schur-mask rules.  A cover is rejected only if **every** labeled support
family is exhausted.  A seed is rejected only if every eligible cap-three
cover is rejected, and a graph is rejected if one required K7 seed is thus
impossible.  All quantifiers therefore have the same sound direction as the
frozen existential-cover reference.

## Trust boundary

All operations are finite bit-mask enumeration, exact integer matching,
exact rational elimination inherited from the frozen reference, or standard
zero-forcing/component arguments.  Floating-point arithmetic is used only
for elapsed-time reporting and plays no role in a decision.

The campaign writer fsyncs decision prefixes and records their committed row
count in an atomically replaced checkpoint manifest.  Resume requires exact
agreement of the input, propagation source, frozen reference, graph count,
and decision schema hashes; an uncommitted or torn TSV tail is discarded.
The report also records the exact command, Python/platform provenance, oracle
hash, and final decision/checkpoint hashes.

## Exact 512-graph sample profile

On the frozen deterministic K7 sample (SHA-256 `bf699f36...`), a clean
11-worker run reproduced all 426 baseline enhanced-rank rejections and then
rejected six additional graphs:

```text
265354  2215178  2646332  2743446  3331579  3558001
```

The 561 baseline-passing covers became 326 passing covers: exhaustive labeled
support propagation rejected 235.  It checked 45,660 labeled families.  Of
the failed assignments, 30,282 produced an empty propagated mask and 15,052
made the propagated masks of a required N edge disjoint.  No assignment in
this sample first failed at one of the recomputed rank/Schur screens, although
that code path has a focused positive/negative control.

The independently developed strict-H sample layer rejected two graphs,
`2215178` and `3331579`; both lie in the six above.  Thus this support layer
adds four sample decisions beyond strict H and leaves 80 of the original 512
graphs unresolved by the combined baseline/support test.  These are exact
sample statements, not yet full-residue coverage claims.
