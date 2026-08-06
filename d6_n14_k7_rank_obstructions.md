# Generic `K7` support/rank audit of the 89 selected `n=14` patterns

## Outcome

The dimension-independent part of the new `K7` support/rank theory certifies
no additional member of the 89 LM-selected `n=14` obstruction hypotheses.
The exact classification therefore remains:

- 78 patterns rejected by the existing eligible-cover bound;
- 3 further patterns (`303`, `304`, and `443`) rejected by the previously
  proved `K5` link bound;
- 81 exact patterns in the combined union;
- 8 patterns still heuristic: `317`, `367`, `368`, `803`, `905`, `936`,
  `952`, and `954`.

This is a useful negative result rather than a failed search.  The audit
exhausted every eligible zero-factor cover of every required `K7` in every
selected pattern, including covers of sizes four through seven.  It did
**not** use the level-19-only theorem `|Z|<=3`.

The generated machine-readable record is
`d6_n14_k7_rank_obstructions.json`, produced by
`analyze_d6_n14_k7_rank_obstructions.py`.

## Exact rules audited

Fix a required seed `Q=K7` in an `alpha(G)<=2` pattern and let `T` be its
seven outside vertices.  For `x in T`, let `D_x` be the set of seed vertices
which are allowed to be non-unit from `x`.  If `Z` is the zero-factor set in
the simplex Schur-complement equations, then:

1. `Z` is contained in `{x: |D_x|>=3}` and is a vertex cover of the graph
   `L` whose edges are required outside pairs `xy` with `D_x` and `D_y`
   disjoint.
2. The vectors indexed by `Z` form an orthonormal family in seven
   coordinates.  Their actual supports `S_x subseteq D_x` consequently
   satisfy `|S_x|>=3`, no pair of supports intersects in exactly one
   coordinate, and the support family has a matching saturating `Z`.
3. When `|Z|=7`, the seven-by-seven orthogonal support matrix also obeys the
   square row rules: a row has support size one or at least three, and two
   row supports cannot intersect in exactly one column.
4. Put `N=T\Z`.  On `N`, the required-edge graph is exactly the intersection
   graph of the allowed masks.  If `nu_N` and `nu_T` are the corresponding
   support term ranks, the normalized Gram matrices obey

   ```text
   U_K = min(7-|Z|, nu_N, nu_T-|Z|),
   rank(K) <= U_K,
   rank(B=K-J) <= min(U_K+1, 8-|Z|).
   ```

5. Exact ordinary zero forcing supplies lower bounds for both matrix ranks.
   For `B`, the stronger component/inertia bound uses that at most one
   nontrivial component block can be indefinite and every other connected
   positive-semidefinite Z-matrix block has nullity at most one.
6. Every required clique principal submatrix of `K` is
   `J+diag(r_x^2)` and is positive definite.  Thus
   `omega(G[N])<=U_K`.  If `F` is the complement of `G[N]`, the independent
   Gram vectors in the `F`-neighbourhood of a vertex lie in its perpendicular
   space, so `Delta(F)<=U_K-1` whenever `N` is nonempty.

All support choices and all covers of sizes zero through seven are finite
bit-mask enumerations.  A feasible support pattern is only a necessary
condition and is never reported as a realization.  A rejection is reported
only when one seed has no cover, or every eligible cover for that seed fails
at least one necessary exact condition.

## Why the certificates survive non-induced containment

The selected `n=14` graph is used as a required-edge pattern; its nonedges
are not prescribed to be non-unit.  Let `H` be any `alpha<=2` edge
supergraph of a pattern `G`.

- The seed `K7` and every edge of `L_G` persist in `H`, while every allowed
  mask can only shrink: `D_x(H) subseteq D_x(G)`.
- Hence an eligible cover for `H` is also one of the covers enumerated for
  `G`.  Support-CSP infeasibility over the larger pattern domains cannot be
  repaired by shrinking those domains.
- If a pattern nonedge `xy` is upgraded to an edge in `H`, then
  `D_x(G)` and `D_y(G)` were disjoint (otherwise `x,y` and a common seed
  non-neighbour would be an independent triple).  The pair is therefore an
  edge of `L_H`, so `x` and `y` cannot both lie in `N`.  It follows that
  `H[N]=G[N]` for the cover under consideration.
- The term ranks `nu_N` and `nu_T` can only decrease after mask shrinkage.
  The zero-forcing graphs, clique number, and complementary maximum degree
  on `N` stay the same while their rank upper bounds can only decrease.

Thus failure of one of the enumerated necessary conditions for `G` implies
failure for every edge supergraph `H`.  The exact 81-pattern union is valid
for the non-induced containment hits recorded by the obstruction pilot.

## Enumeration results

The audit checked 89 patterns, 238 inequivalent-or-not-yet-quotiented `K7`
seeds, and 1,710 eligible covers.  Their size distribution was:

| `|Z|` | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| covers | 48 | 244 | 502 | 526 | 292 | 84 | 13 | 1 |

At cover level, the support CSP rejected 143 covers and the new
perpendicular-degree rule rejected 30.  The cross-orthogonal zero-forcing
rank test, positive-definite clique test, ordinary `B` test, and exact
component/inertia `B` test rejected no individual cover.  Neither the 143
support failures nor the 30 degree failures eliminated every cover for a
previously surviving seed, so all graph-level incremental counts were zero.

The eight still-heuristic patterns had the following cover totals:

| pattern | seeds | covers | support-failing covers | conclusion |
|---:|---:|---:|---:|:---|
| 317 | 4 | 40 | 0 | survives audited conditions |
| 367 | 4 | 44 | 0 | survives audited conditions |
| 368 | 4 | 48 | 0 | survives audited conditions |
| 803 | 5 | 40 | 0 | survives audited conditions |
| 905 | 2 | 32 | 0 | survives audited conditions |
| 936 | 2 | 32 | 0 | survives audited conditions |
| 952 | 1 | 8 | 0 | survives audited conditions |
| 954 | 1 | 32 | 14 | at least one cover still survives |

## Frequent patterns and fixed-sample coverage

Patterns `569` and `571` are already rigorous obstructions by the bounded
eligible-cover certificate; no LM or numerical assertion is needed:

| pattern | exact method | known hits in fixed 800-target sample |
|---:|:---|---:|
| 569 | bounded cover | 497 |
| 571 | bounded cover | 569 |

The union of all 81 exact patterns hits 619 of the fixed 800 stratified
targets:

| stratum | exact-union hits / sample |
|:---|---:|
| certified `K7`, preexisting exact rejected | 160 / 160 |
| certified `K7`, preexisting exact residue | 160 / 160 |
| deferred `K7`, preexisting exact rejected | 160 / 160 |
| deferred `K7`, preexisting exact residue | 139 / 160 |
| deferred `K6`-only residue | 0 / 160 |

The last zero is structural: every selected pattern contains a required
`K7`, so it cannot occur in a `K6`-only graph.  Because the sample was
deliberately stratified, `619/800` is a fixed-sample coverage statement, not
an unweighted population estimate.

## Reproduction and controls

The reported file was generated with 12 processes:

```sh
python3 analyze_d6_n14_k7_rank_obstructions.py \
  --workers 12 \
  --output d6_n14_k7_rank_obstructions.json
```

A one-process rerun was compared after deleting only worker count and timing
fields.  Every mathematical result, certificate summary, counter, and
coverage field agreed:

```sh
python3 analyze_d6_n14_k7_rank_obstructions.py \
  --workers 1 \
  --output /tmp/d6_n14_k7_rank_obstructions_w1.json

jq -s \
  'map(del(.workers,.runtime) |
       .per_pattern_reference |= map(del(.elapsed_seconds))) |
   .[0] == .[1]' \
  d6_n14_k7_rank_obstructions.json \
  /tmp/d6_n14_k7_rank_obstructions_w1.json
```

The comparison returned `true`.  The 12-process run took 0.833 wall seconds;
the sum of per-pattern CPU times was 2.862 seconds.  A GPU is not useful for
this tiny, branch-heavy exact bit-mask enumeration.

Input and artifact SHA-256 hashes:

```text
83059e1c0f689949ec281b47577457ba8b5988c11884cfa0c10c58dd124073f4  analyze_d6_n14_k7_rank_obstructions.py
fa77bf042cd8bcdddbe4425d654bea9136026cb7a0ee92c02f6e7dc8892c42f9  d6_n14_k7_rank_obstructions.json
72bce307675f47c68858c4d4029d1e76a4aa2c31def229852d9fc376a4149464  d6_obstruction_pilot.json
8f4b434f01a5ce01478d24ae9452ead48cb1f3db13b36651e8af259cc763d86c  d6_obstruction_pilot_incremental.json
0e3d74c081b272731848d655da0c68cfba09435e39a2fd2107c32c2bd3b378f0  aeq_d6_n14.txt
3de1d74affbb1ff68264ac51650fc57c3d8ddce55b7601b3635b106366145504  d6_k7_rank_reference.py
```
