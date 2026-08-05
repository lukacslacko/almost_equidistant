# Erratum: BPSSV Table 3, d=6 column (arXiv:1706.06375, Graphs Combin. 36 (2020))

**Claim.** In Table 3 of Balko-Pór-Scheucher-Swanepoel-Valtr,
*Almost-equidistant sets*, the d=6 column (numbers of all abstract
almost-equidistant graphs in R^6 on n vertices) is incorrect for n >= 11.
The published values fail to apply the K_{1,3,3,3} exclusion (their own
Lemma 11, even-d case): the correct class forbids K_8 and K_{1,3,3,3} as
subgraphs in addition to independent 3-sets.

| n  | published (d=6) | correct | difference |
|----|-----------------|---------|------------|
| 11 | 103,333         | 103,194 | 139        |
| 12 | 1,217,849       | 1,210,392 | 7,457    |
| 13 | 19,170,728      | 18,682,440 | 488,288 |

(Values agree for n <= 10, where K_{1,3,3,3} cannot yet occur in the
class or occurs too rarely to matter; first divergence is at n = 11,
the smallest order with room for the 10-vertex forbidden pattern plus
an extra vertex.)

**Evidence.**

1. Independent from-scratch enumeration (`enumaeq6.c`: vertex growth with
   minimum-degree normalization, incremental checks, isomorphism
   rejection) gives the corrected values; the same program's d=5 variant
   (`enumaeq5.c`) reproduces every published d=5 entry exactly
   (n = 4..13), validating the framework.
2. A second, structurally different route at n = 11: enumerate ALL
   105,071 graphs with triangle-free complement (count matches the
   published last column of Table 3 and the known number of
   triangle-free graphs on 11 vertices), then filter by K_8 and
   K_{1,3,3,3}: 105,071 - 1,738 (K_8) - 139 (K_{1,3,3,3}) = 103,194.
3. Definition-level brute force (`verify_k1333.py`, pure itertools over
   apex + triple partitions, no bitmask tricks) confirms each of the
   1,738 + 139 rejections at n = 11: kept = 103,194 exactly.
4. Smoking gun: 105,071 - 1,738 = 103,333 -- the published value equals
   the count with ONLY the K_8 filter applied, i.e. the published d=6
   Table 3 pipeline never rejected K_{1,3,3,3}.

**What is NOT affected.**

- BPSSV's Table 2 (minimal abstract almost-equidistant graphs) is
  correct: our pipeline reproduces the d=6 minimal counts at every
  order we checked (n <= 13: ..., 54, 130, 339; and n = 17: 105,238 of
  164,796 maximal triangle-free graphs, rejecting 9,268 for
  K_{1,3,3,3}) — their minimal-graph pipeline (triangleramsey +
  forbidden-subgraph filter) clearly did apply the filter.
- All theorems of the paper: the bounds rest on the minimal-graph
  enumerations and on Lemma 11 itself, not on the Table 3 d=6 totals.
  (A too-large candidate class would in any case only weaken upper
  bounds, not invalidate them.)

The d=5 column of Table 3 is correct (reproduced exactly at every order
n = 4..13 by our enumeration).
