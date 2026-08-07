# Singleton one-free-edge obstruction

## Exact statement

Fix a regular unit `K7` in `R^6`, a zero-factor cover, and a propagated
support family on the nonzero-factor vertices.  Suppose a required unit edge
`xy` has the following properties:

1. the propagated masks of `x` and `y` intersect in the singleton `{i}`;
2. every allowed coordinate of `x` except `i` is pinned to a sign by a
   nonbipartite singleton-intersection coordinate component;
3. every allowed coordinate of `y` except `i` is pinned in the same sense;
4. coordinate `i` itself is not pinned at either endpoint.

Then the branch is impossible.  This is an exact required-edge-only
obstruction.  A propagated allowed entry may still be zero until an equation
proves otherwise, and no candidate nonedge is assigned any distance.

## Diagonal equation with one free component

On the nonzero-factor set normalize

```text
w_i(x) = sqrt(7) u_i(x) / t_x.
```

The diagonal simplex identity is

```text
Q - T^2 + 2 sqrt(7) T - 8 = 0,                 (1)
T = sum_i w_i,  Q = sum_i w_i^2.
```

Let a propagated mask have size `k`.  If `k-1` coordinates have been pinned
to signs whose sum is `P`, and the remaining coordinate is `z`, substituting
`T=P+z` and `Q=k-1+z^2` into (1) gives

```text
F(k,P) = z
       = (9-k+P^2-2 sqrt(7) P) / (2(sqrt(7)-P)).       (2)
```

Here `1 <= k <= 7` and

```text
P in {-(k-1), -(k-3), ..., k-3, k-1}.
```

The denominator is nonzero because `P` is an integer.  The numerator cannot
vanish in `Q(sqrt(7))` for these values, so (2) also proves that the last
allowed coordinate is actually nonzero.  This is why optional-zero support
semantics do not weaken the argument.  The case `k=1` is included: its empty
sum is `P=0` and `F(1,0)=4/sqrt(7)`.

Rationalizing (2) gives a convenient exact representation:

```text
F(k,P) = [P(P^2-k-5) + (9-k-P^2)sqrt(7)] / [2(7-P^2)].  (3)
```

## Finite reciprocal classification

Because the masks of `x` and `y` overlap only at `i`, their required unit
edge equation is

```text
F(k,P) F(l,R) = 1.                               (4)
```

Exhausting the 28 admissible `(k,P)` types exactly in `Q(sqrt(7))` gives only
the following ordered solutions to (4):

```text
(k,P)=(4,-3), (l,R)=(7,0)
(k,P)=(4, 3), (l,R)=(7,0)
(k,P)=(7, 0), (l,R)=(4,-3)
(k,P)=(7, 0), (l,R)=(4, 3)
```

This finite check uses integer/Fraction arithmetic, not floating point, in
`test_d6_k7_one_free_edge.py`.  It multiplies `a+b sqrt(7)` and
`c+d sqrt(7)` by comparing the rational and irrational coefficients
separately.

The apparent exception cannot occur in the hypothesis: two subsets of a
seven-set with sizes four and seven have intersection size at least four,
whereas the edge masks must intersect in one coordinate.  Equivalently,
singleton intersection implies `k+l <= 8`, and the exact table has no
reciprocal pair in that range.  Thus (4) is always contradictory.

## Certificate

`d6_k7_one_free_edge.py` records:

- the required edge and its singleton shared coordinate;
- both propagated masks;
- for each endpoint and each nonshared allowed coordinate, the vertex set
  and same-color conflict edge of a nonbipartite coordinate component.

The verifier reconstructs every coordinate graph from required edges and
propagated masks, rechecks each nonbipartite component certificate, and
repeats the exact finite reciprocal test.  The production locator is
deterministic.

The source package is deliberately separate from the already frozen
double-pin and full-support-pin layers.  Corpus-level incremental coverage
must be measured only after this source, its independent verifier, and the
upstream exact artifacts are frozen to commit hashes.

## Tests

```sh
python3 -m unittest -v test_d6_k7_one_free_edge.py
```

The tests cover the complete exact reciprocal table, the combinatorial
elimination of its sole exception, nonvanishing of every `F(k,P)`, a positive
synthetic certificate and JSON round trip, a missing-pin negative control,
and certificate tampering.
