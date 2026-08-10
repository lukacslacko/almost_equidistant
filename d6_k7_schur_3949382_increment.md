# Production exact K7 Schur increment for graph 3949382

## Claim and boundary

This package promotes the frozen exploratory Schur contradiction for graph
`3949382` to a source-bound, independently checked one-graph increment.  Its
input is the ordered 16-graph survivor list in
`d6_k7_one_two_star_increment_report.json`; its output is that list with
`3949382` removed.

The producer is `build_d6_k7_schur_3949382_increment.py`.  It imports the
frozen exploratory probe, checks its exact stable report hash, and refuses an
official launch unless the producer, checker, tests, note, and frozen probe
package all equal blobs in the launch commit.  The checker is
`verify_d6_k7_schur_3949382_increment.py`.  It imports neither the producer nor
the probe and reconstructs the certificate independently.

The rejection has no floating-point step.  Candidate nonedges remain optional.
An entry is set to zero only when two independently propagated support
supersets are disjoint; a one entry is supplied only by a required unit edge.

## Exhaustive seed quantifier

The required K7 seed is

```text
[2,5,7,12,14,16,18].
```

The checker reconstructs all 128 eligible raw zero-factor covers.  The pinned,
independently checked upstream layers leave covers `0`, `2048`, and `3072`.
Covers `2048` and `3072` have no passing support family.  Cover `0` has one,
with propagated support supersets

```text
[84,36,127,66,34,33,81,88,40,29,82,61].
```

Thus the new certificate has exactly one branch to eliminate.  Every
realization of the graph would realize this required seed, so infeasibility of
this seed rejects the graph; it is unnecessary to reject the graph's other K7
seeds.

The exact seven-vertex clique among the 12 outside vertices is locally

```text
[0,2,6,7,9,10,11],
```

and the five remainder vertices have clique-neighbour masks

```text
83,47,98,86,90.
```

In the ten-pair order `(0,1),(0,2),(1,2),(0,3),(1,3),(2,3),(0,4),
(1,4),(2,4),(3,4)`, the exact targets are `0111011011`.  The checker
reconstructs all 45 displayed entries and aborts if any optional entry is not
forced by disjoint support supersets.

## Independent exact elimination

The clique principal block is positive definite and saturates the normalized
Gram rank bound seven.  Its Schur complement therefore vanishes.  With
strictly positive diagonal parameters `a,b,c,d,e,f,g`, the ten equations first
give `a=c=d=t`.  Writing the four distinct reduced equations as `E0,...,E3`,
the independent checker verifies coefficientwise

```text
(f+2t+1) A0 = (f+t) E0 + (f+2t) E1,
(b+t-1) B0 = (b-2)(E0+E1) + t(E2+E3),

A0 = b(f+t)-f-t^2-2t,
B0 = bf+bt+b-ft-2f-gt-4t-2.
```

Since `f,t>0`, `A0=0` and

```text
b = 1 + t(t+1)/(f+t) > 1.
```

Consequently `b+t-1>0`, so the second identity gives `B0=0`.  Solving
`A0=B0=E0=0` for `b,g,e` and substituting into `E2` gives

```text
(f+2t+1)^2 [t^2-t-f(t+1)] / [t(f+t)] = 0.
```

All omitted factors are strictly positive, hence

```text
f = t(t-1)/(t+1).
```

Because `f>0` and `t>0`, this forces `t>1`; the reconstructed value

```text
g = -(t-1)(t+1)/(2t^2)
```

is then negative, contradicting `g>0`.  As a strictness control,
`(a,b,c,d,e,f,g)=(1,3,1,1,2,0,0)` satisfies all ten equations but lies on the
boundary rather than in the positive orthant.

## Freeze and official run

The four production source files must first be committed.  Then run serially:

```text
python3 -m unittest -v test_d6_k7_schur_3949382_increment.py
python3 build_d6_k7_schur_3949382_increment.py
python3 verify_d6_k7_schur_3949382_increment.py \
  --report-sha256 <printed-report-sha256>
```

The generated report and verification JSON files are the union-ready evidence.
They must be committed without changing the frozen source boundary between the
producer and checker runs.
