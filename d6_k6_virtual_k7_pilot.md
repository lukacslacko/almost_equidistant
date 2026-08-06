# Virtual K7 promotion from an empty K6 defect

## Exact branch theorem

Fix a required unit `K6`

```text
Q = {q_1,...,q_6}
```

in a hypothetical realization in `R^6`.  Center it in its five-dimensional
span `W`, choose a unit normal `e`, and retain the exact K6 coordinates

```text
u_i(x) = ||x-q_i||^2-1,
c_x    = 1+sum_i u_i(x),
z_x^2  = c_x^2+6-6||u(x)||^2,
x      = -sum_i u_i(x)q_i + (z_x/sqrt(12))e.
```

If the *actual* defect vector is empty, `u(x)=0`, then

```text
c_x=1,       z_x=+sqrt(7) or -sqrt(7),
x=r_+ or r_-,       r_+-r_- = 2sqrt(7/12)e.
```

Thus the only two possible empty-defect points are

```text
r_+ = +sqrt(7/12)e,       r_- = -sqrt(7/12)e.
```

They are unit from every vertex of `Q`, because

```text
||r_+-q_i||^2 = ||r_--q_i||^2 = 7/12+5/12 = 1.
```

Distinctness gives the complete empty-branch classification:

1. there are at most two empty-defect outside vertices;
2. two distinct ones must be the opposite pair `r_+,r_-`;
3. their squared distance is `7/3`, so their candidate pair must be a
   nonedge (a candidate edge is required unit); and
4. relative to the promoted K7 using either apex, the opposite apex is an
   actual one-defect point of value `7/3-1=4/3` in the new apex coordinate.

For an outside vertex `x`, let `G[Q -> x]` be the graph obtained by adding
all missing pairs `xq_i` as required edges and changing no other pair.  In the
branch `u_Q(x)=0`, every added pair is genuinely unit, so the same geometric
realization realizes this augmented required-edge graph.  Moreover
`Q union {x}` is an actual regular unit K7.  Therefore:

> Any sound K7 impossibility proof for `G[Q -> x]` rules out
> `u_Q(x)=0`, even if some second outside vertex is also empty.

This is a branch implication, not an assumption that an arbitrary virtual
apex belongs to the point set.  It avoids the unsound use warned about in the
K6 coordinate note: the seventh coordinate is used only after branching on
an existing outside point and adding only distances forced by that branch.

If `B_Q` is the set of outside vertices whose promoted graph survives the
chosen necessary K7 filters, every realization must have

```text
E_Q = {x outside Q : u_Q(x)=0}
```

in one of the following deliberately relaxed classes:

```text
E_Q = empty;
E_Q = {x},                 x in B_Q;
E_Q = {x,y},               x,y in B_Q and xy is a candidate nonedge.
```

Passing a filter is only necessary, so this list may contain false branches.
If `B_Q` is empty, the virtual-K7 argument proves only that all 13 K6 defect
vectors are nonempty.  It does **not** reject the K6 seed or graph.  A graph
rejection would additionally have to eliminate the all-nonempty branch and
every surviving singleton/pair branch with one common exact quantifier.

## Audit of committed K7 code

The audit distinguished generic mathematical kernels from full-campaign
drivers and archives.  Adding required edges preserves `alpha(G)<=2`, and
all kernels used by this pilot reconstruct the K7 seed, allowed defect masks,
eligible covers, and induced graphs from the supplied augmented adjacency.
They do not look up the corpus index or assume that the graph was a minimal
canonical representative.

| committed code | arbitrary 19-vertex augmented graph? | pilot use / caveat |
|---|---|---|
| `d6_reference_filters.py` | yes | Generic clique-Hall and eligible-cover seed routines.  They are superseded by the stronger rank reference and were not separately counted. |
| `d6_k7_rank_reference.py` | yes, for a 19-vertex graph with `alpha<=2` | **Used.**  `seed_instance`, all nonminimal eligible covers, exact support/rank bounds, zero forcing, component inertia, clique, degree, basis-kernel, and saturating-mask rules have no index/archive dependency.  Its direct cover-cap transitions use the level-19 fact that a K7 has 12 outside vertices, which holds here. |
| `d6_k7_support_propagation.py` | yes at the core-function level | **Used.**  `labeled_support_families` and `analyze_support_assignment` are generic.  The CLI's sample defaults and pinned reference hash are campaign plumbing, not a mathematical precondition. |
| `d6_k7_propagated_link_caps.py` | yes | **Used.**  It consumes only fixed actual Z supports and propagated N masks. |
| `d6_k7_small_support_value.py` | yes | **Used.**  Its one/two-defect CSP is local and archive-free. |
| `d6_k7_small_support_matching.py` | yes | Generic, but proved redundant when `alpha<=2` and required-edge support intersection are already present; omitted. |
| `d6_k7_support_capacity.py` | yes | Generic finite relaxation with an explicit `UNRESOLVED` status.  It had zero marginal coverage in its committed campaign and was outside this small pilot. |
| `d6_k7_special_h_reference.py` | yes through `analyze_graph` / exact clique assessment | Generic exact Schur certificates; a numerical locator may miss a certificate but accepted conclusions are rationally checked.  Excluded from the requested graph/support/sparse-value pilot. |
| `d6_k7_positive_polynomial_dual.py` and `d6_k7_dual_degree2_pilot.py` | yes through their core analyzers | Generic rationally verified polynomial certificates, with numerical LP used only as locator.  Excluded here. |
| `d6_k7_rankone_tetrad_pilot.py` | yes through its graph analyzer | The kernel reconstructs its cover/rank preconditions, but its full runner and selections are archive-bound.  Excluded here. |
| `d6_k7_ranktwo_pentad.py` | yes per cover after reconstructing the rank-two precondition | The polynomial kernel is generic; committed target/full reports apply only to their hash-pinned covers and cannot be transferred by index.  Excluded here. |
| `d6_k7_arbitrary_basis_psd_dual.py` | certificate kernel only | A certificate rules out one proposed basis; cover rejection requires exhausting every admissible basis.  The bounded pilot does not supply that quantifier. |
| `d6_k7_reflection_overlap.py` | yes through `analyze_graph(adj)` | The reflection kernel is graph-generic and exactly verifies its witness.  Its CLI input/selection loader is hash-bound.  Excluded here because this pilot examines the promoted seed's defect system only. |
| `d6_k7_h_mps_triage.py` | no theorem conclusion | Numerical triage only; never a rejection source. |

The full runners, selection builders, conjunction builders, result manifests,
decision archives, and certificate archives remain tied to their recorded
corpus and cover identities.  In particular, an old certificate associated
with the same original corpus index is not automatically a certificate for
the augmented graph: the added edges can change K7 seeds, defect masks, cover
sets, and Schur systems.  Such layers can be reused only by rerunning their
generic kernel and reconstructing the complete quantifier on the new graph.

## Read-only 20-graph pilot

The input is the exact ordered 822-graph complement of the nine committed K6
arbitrary-subset Hall rejections.  Its stable ordered-index hash is

```text
cf94aac41eba28904ee6254f9c50996e73d38dcc15316fb6a411b29ff5cd05d8.
```

The deterministic sample divides the residue, sorted by `(K6 count,
SHA256(decimal corpus index))`, into four equal-size K6-count strata and
takes the five smallest hash priorities in each stratum.  The 20 indices are

```text
2668615, 2668925, 3335382, 2669331, 3926657,
2621024, 3359366, 3331610, 3950119, 3743629,
3910757, 765267, 2147597, 490738, 2593611,
1511757, 2620376, 3668326, 1122794, 1090729.
```

Their ordered list hash is

```text
a45245b9a7fe3b73de131698f1f42d6f68c712e68078a2055c10d2b0d24fd37e.
```

The sample contains 640 required K6 seeds.  The probe examined all 13 outside
vertices for each seed, hence all 8,320 singleton promotions.  It analyzes
the intended promoted K7 seed.  For each promotion it:

1. enumerates every eligible zero-factor cover, including nonminimal covers;
2. uses the exact level-19 direct cap to discard cover sizes 4--7;
3. applies every generic enhanced rank/reference rule to covers of size at
   most three;
4. exhausts every labeled actual Z-support family on baseline-passing covers;
5. applies singleton-overlap propagation, the generic rank rules again, and
   propagated clique-link caps; and
6. applies the exact one/two-defect sparse-value CSP.

The exact sequential branch partition is

| promotion status | count | share |
|---|---:|---:|
| generic graph/rank/cover eliminated | 8,055 | 96.815% |
| labeled-support propagation eliminated | 47 | 0.565% |
| propagated link-cap eliminated | 0 | 0% |
| sparse-value eliminated | 82 | 0.986% |
| survivor | 136 | 1.635% |
| total | 8,320 | 100% |

Thus 8,184 of 8,320 possible empty-defect promotions, **98.365%**, are
eliminated exactly in this sample.  The generic first stage breaks down as

```text
7,503 promotions with no eligible zero-factor cover;
   40 with eligible covers only at directly impossible sizes 4--7;
  512 with small eligible covers but no enhanced-rank passing cover.
```

The 753 baseline-passing covers generated 18,775 labeled support families.
Only 308 families survived propagation; none failed the propagated link caps;
172 of those 308 failed the sparse-value CSP and 136 supplied the first
surviving branch witness.  Surviving searches short-circuit after that first
witness, while every eliminated promotion exhausts its full cover/support
space.

At the K6-seed level:

```text
531 / 640 seeds force every K6 defect vector nonempty;
634 / 640 seeds allow no candidate-nonedge pair of surviving apices;
  6 / 640 seeds retain at least one possible opposite-apex pair;
 38       candidate-nonedge pairs survive in total.
```

Three sampled graphs have no surviving empty-defect promotion at any K6 seed:

```text
2147597, 490738, 1122794.
```

This means a realization of any of those graphs would have to use the
all-nonempty branch at every K6 seed.  It is not a graph rejection.  Across
the whole report `original_graph_rejections_claimed` is explicitly zero.

The stable hash of all `(graph index,K6 seed mask,apex)` keys is

```text
914356bd118fa45fe1d5e1cafaf33fd576251518645dcc8babc3285608316c20.
```

The corresponding hash restricted to the 8,184 eliminated keys is

```text
d587fe540e86bbcb495aa586803b4647174c249268f18bf917913cf84fec310f.
```

## Commands, controls, and artifacts

The 11-worker pilot command was

```text
/Users/lukacs/claude/opengauss/venv/bin/python3 \
  probe_d6_k6_virtual_k7.py \
  --workers 11 \
  --output d6_k6_virtual_k7_pilot_report.json
```

It completed in 3.275 measured wall seconds on
`macOS-14.5-arm64-arm-64bit` with Python 3.11.15.  A second 11-worker run to
`/private/tmp/d6_k6_virtual_k7_pilot_rerun.json` completed in 3.191 seconds.
After removing runtime metadata and each graph's elapsed-time field, the two
reports were byte-canonically equal with stable content hash

```text
50c1a4d7d49593652d0d5d702d27618c72191ce6bd93429f174840af4912e6df.
```

The four inherited kernel suites passed all 43 controls:

```text
/Users/lukacs/claude/opengauss/venv/bin/python3 -m unittest -v \
  test_d6_k7_rank_reference.py \
  test_d6_k7_support_propagation.py \
  test_d6_k7_propagated_link_caps.py \
  test_d6_k7_small_support_value.py
```

Byte-code compilation and `git diff --check` for the probe also pass.

```text
probe_d6_k6_virtual_k7.py
  SHA-256 531f4d123f249de03438e411f03fb7c70a591f3398b2a162b92158a31fd655b0

d6_k6_virtual_k7_pilot_report.json
  bytes 4,036,266
  SHA-256 9f43436fddd5b1aa3511318251cdd992c0621d76807200166ffd49e7fb872148
```

The report also binds the exact parent report, adjacency input, and four K7
kernel sources by SHA-256.

## Interpretation and next exact conjunction

This pilot shows that virtual-K7 promotion is not merely a rare corner case:
generic K7 cover constraints rule out almost every candidate empty apex at
very low cost.  The natural exact use is not to run it as a standalone graph
filter.  Instead, feed each seed's surviving set `B_Q` into the K6
empty-support budget/Hall quantifier:

- a negative-side singleton may be chosen empty only if its absolute vertex
  lies in `B_Q`;
- at most two chosen empties exist globally;
- if two are chosen, their vertices must be a candidate-nonedge pair; and
- every all-nonempty, singleton, and eligible pair alternative must still be
  tested.

That conjunction would let the existing K6 Hall system eliminate the
all-nonempty branch while the K7 promotion eliminates most escape vertices.
It should first be piloted on the same 20 graphs, then source-frozen and
independently checked before any full 822-graph theorem claim.

All decisions here use integer/rational graph, support, zero-forcing, inertia,
Schur-mask, link-cap, and sparse-value logic.  Floating point appears only in
elapsed-time metadata.  Candidate nonedges remain unconstrained and allowed
defect coordinates may vanish.  A surviving branch is not a realization, and
an unexamined or surviving branch is never converted into a rejection.
