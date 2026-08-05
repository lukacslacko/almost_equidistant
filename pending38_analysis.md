# The 38 pending level-17 graphs: structure, commonalities, and realizability assessment

Data file: `aeq_d5_n17_pending38.txt` (17 adjacency bitmasks per graph,
indices into `aeq_d5_n17.txt`). Status: these are the last 38 of the
12654 minimal abstract almost-equidistant graphs on 17 vertices not yet
*certified* non-realizable in R^5; the certified interval searches on
them are running. Everything below is about what they look like and how
likely a realization is.

## What they have in common (measured)

| property | value across all 38 |
|---|---|
| vertices / edges | 17 / 86-95 (most 90-93); degrees only 10, 11, 12 |
| K_6 count | few: 1-10 (median ~4-5) — sparse in rigid seed cliques |
| K_5 count | high: 57-98 |
| max common nbhd of a non-adjacent pair | 6-9 (never >= 10) |
| best LM residual, distinct points | 0.409 - 1.672 (nowhere near 0 = 1e-30 scale of realizable controls) |
| best LM residual, coincidences allowed | identical to the distinct-point value for every single graph |
| complements | near-regular maximal triangle-free graphs, degrees 4/5/6 (mostly 5-regular-ish) — the "quintic" Ramsey(3,7)-flavoured family |

Reading of the table:

1. **They are the barely-overdetermined middle of the family.** With
   90-93 edges against 70 degrees of freedom (85 coordinates minus the
   15-dimensional isometry group), and degrees pinned to 10-12, these
   are the most "regular" candidates - no vertex is strongly
   constrained, no non-adjacent pair has a large common neighbourhood.
   The engine's two coincidence-discard rules need 5 (resp. 6) common
   placed neighbours of a non-adjacent pair; with only 6-9 common
   neighbours in total (`maxcn` column), those rules activate only late
   in a placement order and with no slack. That is *why* these graphs
   are the slow tail for the interval search: wide intermediate boxes
   survive longer before any rule can bite.

2. **Few K_6's = few rigid anchors.** The bulk of the 12654 died in
   milliseconds because a K_6 seed plus a cascade of >= 5-neighbour
   sphere placements keeps every box tight. The pending graphs have 1-10
   K_6's and their placement orders are forced through a circle stage
   (every one of the 38 is a one-circle graph), where boxes are
   parametrized by an angle and everything downstream inherits the
   angle's width.

3. **The complements cluster.** All 38 complements are maximal
   triangle-free graphs with degree sequences built from 4/5/6, most
   nearly 5-regular (largest families: 14x deg-5 + 3x deg-6, etc.).
   Near-regular quintic maximal triangle-free graphs on 17 vertices are
   exactly the "extremal-Ramsey-like" zone; several pending indices come
   in consecutive runs (11377-11405, 12519-12595) because
   triangleramsey emits structurally similar graphs together. In short:
   the 38 are one structural family with small variations, not 38
   unrelated hard cases.

## What point constellations are still "considered" by the search

The unresolved regions are NOT neighbourhoods of near-solutions. Three
independent observations:

1. **Numerics, global:** 300-restart Levenberg-Marquardt per graph
   (both random and clique-seeded initializations) never gets below
   residual 0.409 on any of the 38 - against 1e-30 for realizable
   control graphs. For a 90-ish-edge system, residual 0.4-1.7 means
   several edges are off by ~5-10% of a unit length at the best
   attainable configuration: qualitatively far, not "almost there".

2. **Numerics, degenerate:** for every one of the 38, the best residual
   with coincidences allowed equals the distinct-point best - and for
   the hardest graph (11402) we additionally contracted each of its 44
   non-adjacent pairs (u=v identified) and re-optimized: best residual
   1.034. So there is no exact degenerate (collapsed) configuration
   either - unlike the d=4 campaign's channel, which did contain exact
   coincident families.

3. **Interval geometry:** direct probes of a "hot" arc show that
   splitting it once makes all sub-arcs die in ~65 search nodes. A
   region containing (or even approaching) a true solution cannot die
   at all - killed means *certified empty*. The hot arcs are therefore
   width artifacts of interval evaluation (dependency blow-up at a
   particular box width), not geometric features. The still-running
   ladder is bookkeeping, not suspense about new point constellations.

**Assessment:** the probability that one of the 38 is realizable is,
on this evidence, negligible. Every quantitative signal - global
numerics, degenerate-branch numerics, and the behaviour of the certified
search itself - points to non-realizability, matching the other 12616
graphs of the level. The remaining work is to convert that into
certificates, which only the interval engine can do; no numerical
observation is treated as proof, which is exactly why the grind
continues.

(This file will be superseded by the final result note when the last
certificates land; the survivor-box dump from a live hot core and
LM restarts seeded from those boxes are appended below when available.)

## Appendix: survivor-box dump from the live hot core (completed)

The debug kernel (env-gated dump of any parked survivor leaf) was run on
the hottest known sub-arc of the hardest family's graph 11402
(decomposition 0, theta in [78.28, 78.75] degrees, coarse floor, 3M-node
budget). Result: **ABORT at the node budget with zero unresolved cells
and zero survivor leaves** — in three million nodes the search parked
nothing at all: no candidate constellation ever reached even the
"unresolved, report later" state. The hot arc is an expensive kill-tree,
not a reservoir of candidate configurations. Consequently there are no
box midpoints to seed LM from — the closest thing to "surviving
structure" is the empty set; the standard sweeps (clique-seeded and
random, all far from solvability) remain the operative numerics.
