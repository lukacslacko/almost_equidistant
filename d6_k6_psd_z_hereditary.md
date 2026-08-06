# Hereditary PSD--Z support Hall: exact K6 pilot

This pilot adds a principal-submatrix consequence to the exact K6 PSD
Z-matrix layer.  On the 861 graphs surviving `d6_k6_psd_zmatrix.py`, it
rejects 30 graphs and leaves 831.  Every rejection is already obtained by
the hereditary bipartite system alone; none needs the frozen non-bipartite
system as a final conjunction.  The known realizable 18-point control passes
all 32 of its required K6 seeds.

The implementation is `probe_d6_k6_psd_z_hereditary.py`.  This is an exact
pilot, not yet a production certificate package or an independent checker.

## Principal-submatrix lemma

Let `A` be a real symmetric positive-semidefinite Z-matrix whose graph of
strictly negative off-diagonal entries is connected.  Then every proper
principal submatrix of `A` is positive definite.

If `A` is positive definite this follows immediately.  Otherwise choose
`rho > max_i A_ii` and put

```text
B = rho I - A.
```

The matrix `B` is nonnegative and irreducible.  Since `A` is singular and
positive semidefinite, `rho` is the Perron root of `B`.  Perron--Frobenius
therefore makes `ker A` one-dimensional and spanned by a vector with every
coordinate strictly positive.  If a proper principal submatrix `A[T,T]`
were singular, a nonzero vector `x` supported on `T` would satisfy
`x^T A x=0` after extension by zero.  Positive semidefiniteness implies
`Ax=0`, contradicting the full support of every nonzero kernel vector.

For a PSD--Z-applicable connected K6 side-graph component `F_i`, diagonal
congruence and, in the positive bipartite case, signature switching identify
its Gram block `Q_i` with such an `A`.  Congruence respects principal
submatrices.  Consequently

```text
rank Q_i[T,T] = |T|                 for nonempty proper T subset F_i.
```

For `T=F_i` the pilot retains the strongest old componentwise lower bound

```text
r_i = max(ordinary zero-forcing, exact inertia, PSD--Z n-1).
```

It does not pretend that a full singular block has rank `|F_i|`.  A
positive-sign non-bipartite `F_i` is not PSD--Z-applicable; for that block the
pilot permits only the empty choice or the full block with its existing
lower bound.

## Strongest additive Hall formulation used

Fix one eligible `Z0`, one bipartite Lorentz component, and one of its two
generic sign orientations.  Decompose each of the two side graphs into
connected components.  Choose independently:

* at most one subset `T_i` from each connected side component `F_i`;
* any subset of the singleton `Z0` blocks.

For a PSD--Z-applicable `F_i`, every `T_i` is allowed and contributes
`|T_i|` when proper or `r_i` when full.  For a nonapplicable `F_i`, only the
full choice contributes `r_i`.  The chosen spans from different `F_i`, from
opposite Lorentz sides, and from `Z0` are mutually Gram-orthogonal.  Thus
every such selection must obey

```text
sum_i rank_lower(T_i) + number_of_chosen_Z0
    <= |union of the selected vertices' allowed seed-coordinate masks|.   (H)
```

Only one subset is selected from a given connected span.  In particular,
the ranks of two overlapping (or disjoint) proper subsets of the same `F_i`
are never added; their union is simply another single choice.  Selecting
every full component recovers the previous whole-side block Hall test, so
this system genuinely refines the frozen parent test.

The actual support of a defect vector may be any subset of its candidate
allowed mask.  Replacing actual supports by allowed masks enlarges the right
side of (H), hence makes the test weaker.  Failure of (H) therefore never
requires a candidate nonedge to be genuinely non-unit or an allowed
coordinate to be nonzero.

The implementation compresses the Cartesian product exactly.  After each
connected span it keeps, for each of the 64 possible coordinate unions, the
largest attainable rank.  A smaller rank with the same union can never
create a later violation that the larger one cannot create.  This dominance
argument is the only compression; there is no numerical or heuristic prune.
The generic systems remain alternative to the separate lightlike case.

## Controls

The probe runs exact integer controls before reading the corpus:

* the singular path Laplacian has determinant zero and every proper
  principal determinant positive;
* a positive-definite connected Z-matrix has all tested principal
  determinants positive;
* signature switching the bipartite path preserves those determinants;
* a reducible PSD Z-matrix supplies the expected counterexample, showing why
  connectedness is needed;
* a synthetic support system has an old full-block inequality `2 <= 2` but a
  hereditary proper-pair violation `2 > 1`;
* the known realizable 18-point graph passes all 32 required K6 seeds with
  zero impossible seeds.

## Full 861-graph pilot

The parent boundary is the ordered 861-index complement of the exact PSD
Z-matrix rejection set, with stable hash

```text
09ebce17d2b72fa6514fc6a8a938376626d373161d2e1b4dd8d31c5e00c18db5.
```

Commands:

```text
/Users/lukacs/claude/opengauss/venv/bin/python3 \
  probe_d6_k6_psd_z_hereditary.py --workers 11 \
  --output /private/tmp/d6_k6_psd_z_hereditary_full.json

/Users/lukacs/claude/opengauss/venv/bin/python3 \
  probe_d6_k6_psd_z_hereditary.py --workers 1 \
  --output /private/tmp/d6_k6_psd_z_hereditary_full_serial.json
```

The 11-worker run took 2.56 seconds; the serial run took 13.15 seconds.  The
reports are byte-different only in their runtime records and are exactly
equal after removing the `runtime` object.  The stable hash of that common
content is

```text
1397cb3d354ebebc7f0b0c96eab1a64dbf63e661f6f9557fd353fb75f07a2508.
```

The parallel raw report has SHA-256
`7ddaed302ac94c718626c61801b5906f61e326830432ec59493be0c5066d3def`.
It is 44,685,877 bytes because it retains an exhaustive `Z0` diagnostic for
the first impossible seed of every rejected graph, so it remains in
`/private/tmp` rather than being proposed as a Git artifact.

| quantity | exact pilot value |
|---|---:|
| pinned input / rejected / surviving | 861 / 30 / 831 |
| standalone hereditary rejections | 30 |
| same-`Z0` conjunction-only rejections | 0 |
| K6 seeds checked / impossible | 26,616 / 30 |
| `Z0` choices considered / matchable | 31,303 / 31,303 |
| bipartite Lorentz components checked | 205,528 |
| parent-failing components | 4,262 |
| parent-passing components newly failed | 227 |
| generic orientations checked | 411,056 |
| parent-passing orientations newly failed | 556 |
| exact dominance-DP transitions | 3,061,385 |

The 30 rejected indices are

```text
95496, 565234, 940626, 967444, 967512, 1069852, 1167032, 1314598,
1719587, 1719589, 1784375, 1948930, 2096710, 2237945, 2292297,
2382072, 2569222, 2569269, 2637179, 2820201, 2862937, 2862950,
2862961, 2870953, 2902883, 3005616, 3328693, 3555822, 3682937,
3734545.
```

Their stable list hash is

```text
bf256b846be361287784260dfd77b2f50afad63731d9bb81adbe889eab4d82ee.
```

For example, graph `95496` has an impossible K6 seed
`[0,1,3,6,14,18]`.  Its parent system had four common passing `Z0` choices;
the hereditary system has none.  At `Z0=[]`, one orientation selects a
rank-three full nonapplicable component on vertices `[2,12,16]` and an
independent proper pair `[4,10]`.  Their allowed masks use only four seed
coordinates, while their mutually orthogonal spans require rank five.  The
opposite orientation has the symmetric rank-two proper-pair/rank-three full
obstruction.  The lightlike alternative also fails.

## Source-bound production package

The pilot has now been translated into a source-only production boundary:

* `d6_k6_psd_z_hereditary.py` writes a deterministic report, a compact
  exhaustive first-seed certificate archive, and an atomic ordered-prefix
  checkpoint.  Resume requires both the exact input hash and production
  source hash to match.
* `verify_d6_k6_psd_z_hereditary.py` imports neither the production evaluator,
  this probe, nor a shared new kernel.  It reconstructs Lorentz components,
  side components, bipartiteness, sign applicability and subset Hall; uses
  the frozen independent SymPy/Sturm inertia and simultaneous zero-forcing
  implementations; brute-forces one alternative per connected span; and
  validates every archived Hall witness, including that no span label is
  selected twice.
* `test_d6_k6_psd_z_hereditary.py` contains seven source-bound controls.  In
  particular, two overlapping independent subsets of one three-vector span
  would naively give the false inequality `4 > 3`; the correct one-choice
  formulation passes.  Treating the same data as two genuinely distinct
  orthogonal spans correctly fails.

The cheap source-bound tests include the known positive graph and an
independent replay of fixed pilot rejection `95496`.  The final 861-graph
production and independent verification are intentionally deferred until
these sources are committed, so their source hashes form the run boundary.
