# Two exact nonlinear special-H obstruction templates

This note records two exact mask obstructions for a saturating required
six-clique in the K7 normalized-Gram reduction.  They were found by using the
Apple-MPS run only to rank targets; the proofs below are rational identities
and positivity arguments independent of that numerical run.

## Common setup

Let `C` be a required clique whose order six equals the rank upper bound for
the normalized Gram matrix `K`.  Its principal block is

```text
M = J + diag(e_i),    e_i>0,
H = M^-1.
```

Put `x_i=1/e_i>0`, `X=sum_i x_i`, and `T=1+X`.  For a basis-neighbour mask
`S subseteq C`, write `w(S)=sum_{i in S} x_i`.  Sherman--Morrison gives

```text
K_ST = w(S intersect T) - w(S)w(T)/T.
```

Thus a required edge between outside masks `S,T` gives

```text
T w(S intersect T) - w(S)w(T) = T,
```

whereas an outside candidate nonedge gives the same left side equal to zero.
Every division below is by an explicitly positive quantity.

## Type I: one outside nonedge

Up to relabelling the basis coordinates, the five outside masks are

```text
A={c,e,f},       B={b,c,d,e,f},       C={d,e,f},
D={a,b,e},       E={c,d,e}.
```

All outside pairs are required edges except `BD`, which is a candidate
nonedge in the fixed nonzero-factor K pattern.

Because `A,C,E` are subsets of `B`, their three edge equations with `B` are

```text
w(A)(1+a)=T,   w(C)(1+a)=T,   w(E)(1+a)=T.
```

They imply `c=d=f`.  Reuse `c` for their common positive weight and put

```text
P=2c+e.
```

Then

```text
T=P(1+a).                                           (1)
```

The intersections of each of `A,C,E` with `D` are the singleton `{e}`.
Their edge equations give

```text
b=a(e-2)-1.                                         (2)
```

Comparing (1) with

```text
T=1+a+b+3c+e
```

and using (2) yields

```text
a=c/(2c+1).                                         (3)
```

The pairwise intersections among `A,C,E` have weight `c+e`.  Their edge
equation, together with (1), is

```text
(1+a)(c+e-1)=2c+e.
```

Using (3) reduces this to

```text
ce=c^2+4c+1,                                        (4)
b=c^2/(2c+1).                                       (5)
```

Finally, `B intersect D={b,e}`.  The required nonedge equation `K_BD=0`
would be

```text
P(b+e)=(1+a)(P-1)(e-1).                             (6)
```

Substituting (3)--(5) into the left side minus the right side of (6) gives

```text
(3c+1)/c,
```

which is strictly positive because `c>0`.  This is a contradiction.

## Type II: a three-nonedge star

The complementary template has masks

```text
A={a,c,e,f},     B={b,c,d,e,f},       C={a,d,e,f},
D={b,e},         E={a,c,d,e}.
```

Here `B` is nonadjacent in the fixed K pattern to `A,C,E`; every other
outside pair is a required edge.

For each `Y in {A,C,E}`, the `BY` nonedge equation is

```text
(1+a)w(Y)=Ta.
```

Consequently `c=d=f`.  Put

```text
P=w(A)=w(C)=w(E)=a+2c+e.
```

Then

```text
T=P(1+a)/a.                                         (7)
```

Since `D subset B`, the `BD` edge equation gives

```text
w(D)=T/(1+a)=P/a.
```

The intersections of `D` with `A,C,E` are all `{e}`.  Their edge equations
give

```text
P=(1+a)(e-1).                                       (8)
```

Writing `w(D)=b+e=P/a` in the total-weight identity and comparing with (7)
gives

```text
e=c+1.                                              (9)
```

Equations (8) and (9) now say

```text
a+2c+1=ac.                                          (10)
```

The pairwise intersections among `A,C,E` have weight `a+c+e`.  Their edge
equation is

```text
(1+a)(a+c+e-1)=aP.
```

Using (9) and `P=a+3c+1` reduces this to

```text
c(2-a)=0.
```

Positivity gives `a=2`, but substituting this into (10) gives `3=0`.  This is
the required contradiction.

## Combinatorial detector and scope

`d6_k7_special_h_templates.py` recognizes each pattern invariantly under
basis-coordinate and outside-vertex relabelling.  For Type I it checks the
unique outside nonedge, mask sizes `3,5,3,3,3`, the three-mask triangle, and
the common singleton intersections with `D`.  For Type II it checks the
three-nonedge star, mask sizes `4,5,4,2,4`, and the analogous triangle after
adding the basis coordinate omitted by `B`.

These are cover-specific obstructions.  A graph is rejected only if, for at
least one required K7 seed, every eligible cover is eliminated.  On sample
graph 2592168, four of the eight covers surviving the exact strict-affine-H
layer are eliminated: each of those four covers has two saturating cliques,
one of Type I and one of Type II.  The other four covers survive, so the
template layer makes no whole-graph rejection on the 84-graph sample.  This
distinction between failing cliques and failing covers is checked in
`d6_k7_special_h_templates_report.json`.  Full selected-residue counts are
recorded separately in the machine-readable report; no claim in this note
depends on the heuristic MPS residuals.
