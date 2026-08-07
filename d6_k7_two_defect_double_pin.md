# Exact K7 two-defect double-pin conjunction

## Certified result

The exact double-pin obstruction rejects **130** of the 155 K7-containing
graphs in the v5 dimension-six residue and leaves 25 unresolved.  The full
seed/cover/support conjunction has been replayed by an independently
transcribed checker, which returns `PASS`.

This is not a proof of `f(6)=18`.  The 25 K7 survivors and the separate
K6-only residue remain unresolved unless another independently verified
layer rejects them.

The separately packaged odd-cycle/full-support generalization rejects one
further survivor, 3919831, leaving the exact 24-graph K7 list and hash in
`d6_k7_full_pin_odd_cycle.md`.  The two packages are combined by set union:
the full-pin increment is explicitly conditioned on the independently PASS
double-pin base, so no graph is double-counted.

## The local obstruction

Fix a regular unit K7 seed.  For an outside point put

```text
u_i(x) = ||x-q_i||^2-1,
t_x    = 1+sum_i u_i(x).
```

On the nonzero-factor set `N={x:t_x!=0}`, normalize

```text
w_i(x) = sqrt(7) u_i(x)/t_x.
```

If `xy` is a required unit edge, the simplex bilinear identity gives

```text
w(x).w(y)=1.                                      (1)
```

Suppose three vertices form a required triangle and every pair of their
propagated allowed support masks intersects in the same singleton coordinate
`i`.  Actual supports are contained in those masks, so each dot product in
(1) is a single coordinate product.  Writing the three values as `a,b,c`
gives

```text
ab=ac=bc=1.
```

Therefore `a=b=c` and `a^2=1`: the triangle pins all three `i`-coordinates
to a common value in `{+1,-1}`.

Now let a vertex `x` have propagated mask `{i,j}`.  If `x` lies in a pinning
triangle at `i` and another pinning triangle at `j`, both coordinates are
nonzero, so its actual support is exactly `{i,j}` and both normalized values
belong to `{+1,-1}`.  But the normalized point equation for an exact
two-defect point is

```text
(w_i-sqrt(7))(w_j-sqrt(7))=3.                     (2)
```

No sign pair satisfies (2).  If the two signs sum to zero, its expanded
left-minus-right equation equals `3`; otherwise (2) would make `sqrt(7)`
rational.  This is the double-pin contradiction.

Only required edges are used.  Candidate nonedges remain unconstrained and
may also have distance one.  The propagated masks are supersets of actual
supports; singleton intersections become exact only because (1) forces the
corresponding product to be nonzero.

## Complete exact accounting

The production builder reconstructs every inherited cover gate on the 155
embedded adjacencies and then exhausts every labeled zero-factor support
family:

```text
K7 graphs                                             155
K7 seeds                                              237
eligible covers                                    87,270
inherited-current covers                              535
labeled zero-factor support families               19,716
families reaching the double-pin test                  555
double-pin-infeasible families                         336
double-pin-passing families                            219

current covers infeasible after conjunction            419
current covers still passing                           116
graphs rejected                                        130
graphs surviving                                        25
```

Of the 535 current covers, 215 were already support/sparse-value infeasible;
the double-pin rule eliminates a further 204 covers.  A graph is rejected
only when one required K7 seed loses every current cover.  No cap, timeout,
floating-point rank decision, or numerical optimization enters the result.

The ordered 25-survivor list is

```text
226183 316173 423661 424226 2581209 2592657 2593240
3595554 3624785 3648882 3729907 3785980 3888410 3919831
3935560 3936176 3936177 3936310 3936435 3945490 3945555
3945557 3945564 3947605 3949382
```

Its stable JSON SHA-256 is

```text
16872c94cce864ca17689714c05cb3a49e2cf6cb658d808c270d047c2794ac04
```

## Relation to the interval campaign

The first two interval-certified K7 graphs, 379078 and 2280137, exposed the
pattern.  The exact filter rejects 21 of the 24 independently interval-killed
K7 indices at cap 100,000 and adds 109 further rejections.  The interval-only
indices are 423661, 424226, and 3936176.  Thus the interval result remains
useful independent corroboration and contributes three graphs beyond this
algebraic layer when the two verified rejection sets are unioned.

## Independent verification and trust boundary

`verify_d6_k7_double_pin_conjunction.py` imports neither the production
builder nor the double-pin locator.  It independently:

1. reconstructs all K7 seeds and all size-at-most-seven eligible covers;
2. replays the inherited baseline, strict-H, degree-one, and tetrad gates;
3. enumerates labeled zero-factor supports with the existing independent
   support implementation;
4. recomputes propagation and sparse-value constraints with the existing
   independent checker;
5. locates every pinning triangle and double pin with separately written
   loops; and
6. compares every graph, seed, cover, family count, certificate, passing
   witness, rejection, and survivor against the production report.

It reproduces all 155 graph records, all 336 local certificates, and the
ordered 130/25 partition.  The proof trusts the displayed Euclidean algebra,
Python arbitrary-precision integer/bit-mask/container operations, the
hash-pinned inherited exact kernels, and SHA-256 binding.  It does not trust
IEEE-754, BLAS, `libm`, or failed numerical search.

## Reproduction

Run the focused controls:

```sh
python3 -m unittest -v \
  test_d6_k7_two_defect_double_pin.py \
  test_d6_k7_double_pin_conjunction.py
```

Build the theorem artifact with eleven workers:

```sh
python3 build_d6_k7_double_pin_conjunction.py \
  --workers 11 \
  --output d6_k7_double_pin_conjunction_report.json
```

Bind the independent replay to the printed report hash:

```sh
python3 verify_d6_k7_double_pin_conjunction.py \
  --report d6_k7_double_pin_conjunction_report.json \
  --report-sha256 <printed-report-sha256> \
  --workers 11 \
  --output d6_k7_double_pin_conjunction_verification.json
```

The official source-bound run used commit
`daac9c54a31e76f1687a457413a349ea25eb4823`.  Eight-worker production took
22.8246 seconds; the independent eight-worker replay took 21.3231 seconds and
checked 155 graphs, 237 K7 seeds, 87,270 eligible covers, 19,716 labeled
support families, and all 336 local certificates.

```text
d6_k7_double_pin_conjunction_report.json
  7c90c9a518243de4095f4ec4394bb7b1d22c5889844e90cd7fe4d895f22ac873
d6_k7_double_pin_conjunction_verification.json
  890568a80b8363de84997db86ab5271fee13f707544f91145d15f4976512a33b
```
