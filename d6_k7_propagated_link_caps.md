# Clique-link caps after labeled K7 support propagation

This note proves the independent graph/mask rule implemented in
`d6_k7_propagated_link_caps.py`.  It does not modify the frozen sparse-value
package or production runner.

## Rule

Fix a required K7 seed.  For a subset `S` of its coordinates, `|S|=k`, every
outside point whose actual defect support is contained in `S` is unit from
all `7-k` seed vertices outside `S`.  It is therefore a common unit neighbour
of the complementary seed `K_(7-k)`.

The `k` seed vertices in `S` are also distinct common unit neighbours of that
complementary clique.  The exact clique-link bounds in `R^6` consequently
give the following caps on *outside* points with support contained in `S`:

| `k` | complementary clique | full link bound | outside cap |
|---:|---:|---:|---:|
| 1 | K6 | 2 | 1 |
| 2 | K5 | 4 | 2 |
| 3 | K4 | 10 | 7 |
| 4 | K3 | 12 | 8 |
| 5 | K2 | 16 | 11 |

Now fix a labeled zero-factor support family.  Its supports are actual.  For
each nonzero-factor point, singleton orthogonality propagation produces a
mask `E_x` containing its unknown actual support.  If

```text
E_x subseteq S,
```

then the actual support is necessarily contained in `S`, so the point may be
counted without choosing its support values.  The checker adds all fixed
zero-factor supports contained in `S`, adds all propagated nonzero masks
contained in `S`, and rejects if the displayed outside cap is exceeded.

The quantifier direction is sound: points with masks not contained in `S`
are ignored even if their eventual actual supports might be contained in
`S`.  Thus the count is a lower bound on the actual link population.  A
rejection uses only required unit distances and distinctness.  Candidate
nonedges remain unconstrained and may still have unit distance.

All 119 nonempty coordinate subsets of sizes one through five are checked in
fixed size/numeric order.  The implementation uses only exact seven-bit
masks and integer counts.

## Sample assessment

The exact single-process profiler replayed the 61 survivors of the frozen
sparse-value sample.  It reached 70 propagation- and sparse-value-passing
labeled families across 70 seed checks.  Every one passed the propagated
link caps, so this layer added **zero** graph or family rejections on that
sample.  This is a precise negative pilot result, not evidence that the rule
is redundant on the much larger full residue.  Its cost is only 119 small
mask counts per surviving family, so a full-pass measurement remains cheap
if it is integrated after independent review.
