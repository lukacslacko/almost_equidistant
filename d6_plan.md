# f(6) campaign plan (queued behind f(5) completion)

Published state (BPSSV 2020): 18 <= f(6) <= 26. Lower bound: the
Larman-Rogers 16-point set in a 5-dim subspace plus two apex points at
+-sqrt(3/8) e6 (18 points, verified the same way as d=5). Upper bound:
no minimal Ramsey(3,8) graph on 27 vertices is abstract (BPSSV via
McKay's lists), i.e. zero candidates at level 27.

Abstract almost-equidistant graph in R^6 (even d rule):
complement triangle-free + no K_8 + no K_{1,3,3,3} subgraph.

Candidate counts (BPSSV Table 2, d=6 column):
  n=19: 3,971,787 minimal graphs; n=20..25: nonzero but unknown counts;
  n=26: unknown; n=27: 0.

## Why bottom-up this time

Top-down (as done for f(5)) stalls at the enumeration step: level 26
needs all maximal triangle-free graphs on 26 vertices with alpha<=7
(Ramsey(3,8)-maximal graphs) — plausibly 10^10..10^13 graphs. Not
enumerable here.

Bottom-up decides everything in ONE level: killing all 3,971,787
candidates on 19 vertices proves f(6) <= 18, hence f(6) = 18 (any
realizable >19-point set contains a 19-point subset). Extrapolating the
f(5) phase-A rate (0.5 ms/graph at n=17, d=5; expect ~5-20x at n=19,
d=6), the bulk pass is ~80-300 core-hours ≈ 3-10 h wall on this machine;
the channel tail is the wildcard (the f(5) tail was ~0.5% of graphs).
If a realizable graph exists at 19, the numerical sweep should find it
(f(6) >= 19: bigger news) and the campaign target moves up.

No known construction beats 18 in R^6: the triangle-free strongly
regular graphs (the natural sources, cf. Clebsch for d=5) skip from 16
vertices (Clebsch) to 50 (Hoffman-Singleton) — nothing lands in R^6
range. So f(6) = 18 is the expected outcome, but the computation decides.

## Steps

1. Generate + verify candidates at n=19: triangleramsey 19 (13.7M maximal
   TF graphs, seconds) filtered by alpha<=7 (no K_8 in complement) and a
   new K_{1,3,3,3} checker; count must equal 3,971,787. Independent
   validation: extend enumaeq5.c to d=6 rules and check the published
   totals for n <= 12 (Table 3: 7, 14, 38, 107, 409, 1888, 12064,
   103333, 1217849, ...) and minimal counts (2,3,4,6,10,15,29,54,130,
   339,1052,...); also set-compare at some small n.
2. Port the engine: ckernel6.c (D=6, MAXN 28, K_7/K_6 seeds — unit
   6-simplex enclosures, sphere stage >= 6 placed neighbours, circle
   stage exactly 5, Rule 1 with 6-subsets, Rule 2 with 7 common
   neighbours / 6x6 interval determinant). Controls: K_8 killed;
   K_{2,2,2,2,2,2} (6-dim cross-polytope) survives; apex+cross killed;
   the 18-point LR graph survives; K_{1,3,3,3} killed (it IS the
   dimension lemma for even d - a strong new control unavailable at
   d=5 since K_{3,3,3} has no K_5/K_6 seed, whereas K_{1,3,3,3} has
   omega=4... also no seed. Fallback controls as at d=5).
3. LM sweep (lm5.c generalizes trivially; D as a parameter) over all
   3.97M graphs, prioritized by edge count ascending; flag anything
   with residual -> 0.
4. Bulk kill campaign at n=19 with the phase A/B/finisher pipeline
   from f(5) (reproduce5/finish17 generalized).
5. If complete: f(6) = 18; write f6_equals_18.pdf; push. If survivors:
   investigate realizability; the result may instead be f(6) >= 19.

## Reusable machinery from f(5)

- multicode filter (filter_mtf.c): parametrize forbidden subgraphs.
- bulk/racing/slicing drivers: parametrize dimension and level.
- The channel diagnosis + finisher (96-slice coarse scan, hot-arc
  drilling with cap ladder): expect the same phenomenon around the
  18-point degenerate collapses.
