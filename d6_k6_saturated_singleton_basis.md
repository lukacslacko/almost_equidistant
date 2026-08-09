# Saturated singleton K6 ray bases

This note proves the exact necessary condition implemented in
`d6_k6_saturated_singleton_basis.py` and independently reimplemented in
`verify_d6_k6_saturated_singleton_basis.py`.  Its frozen input is the ordered
251-graph residue certified by the K6 opposite-ray rank conjunction package.
The rule rejects exactly two further graphs in the source-bound census,
indices `3138618` and `3673988`, leaving 249.

All arithmetic in the new rule is exact.  There is no numerical optimization,
floating-point rank decision, tolerance, or transcendental-function
assumption.  A candidate nonedge remains geometrically unconstrained and may
still have unit distance.  An allowed defect coordinate may vanish.

## 1. K6 coordinates and the quantifiers

Fix a centered regular unit `K6` in `R^6`.  For an outside point `x`, use the
six defect coordinates and the scaled normal coordinate of the parent
package:

```text
u_i(x) = ||x-q_i||^2-1,
c_x    = 1+sum_i u_i(x),
z_x    = sqrt(12) times the normal coordinate.
```

They obey

```text
z_x^2 = c_x^2+6-6||u_x||^2,                         (1)

||x-y||^2-1
  = (c_x c_y-z_x z_y)/6-u_x.u_y.                    (2)
```

Let `D_x` be the coordinates at which `u_x` is *allowed* to be nonzero by the
candidate graph, and let `S_x=supp(u_x)` be its actual support.  The only
support assertion is

```text
S_x subseteq D_x.
```

The parent package exhausts all admissible choices of its zero Lorentz set
`Z0`, light-ray coloring, and actual supports.  The new test is called only
at a complete actual-support leaf, after every frozen parent leaf condition
has passed.  If the new condition fails, only that leaf is pruned; the parent
depth-first search continues through every other actual-support choice.

This placement is essential.  It would be invalid to find one saturated
support witness, reject it, and ignore alternative support assignments.

## 2. A six-singleton bin is the coordinate basis

Suppose `Z0` is empty and one nonzero light-ray bin has six vectors, each with
a singleton actual support.  Vectors in one ray bin are orthonormal.  Two of
their singleton supports therefore cannot be equal, so the six supports use
the six coordinates bijectively.  Relabel the bin as `x_1,...,x_6`, with

```text
S_{x_i}={i}.
```

A nonzero light-ray vector has `z_x=sigma c_x`, where the common ray sign is
`sigma in {+1,-1}`.  Equation (1) gives `||u_x||=1`.  Its single nonzero
coordinate is consequently `+1` or `-1`.  The value `-1` would give
`c_x=0`, hence also `z_x=0`, contrary to membership in a nonzero light-ray
bin.  Therefore

```text
u_{x_i}=e_i,       c_{x_i}=2,       z_{x_i}=2 sigma. (3)
```

## 3. Reconstructing every remaining support

For any remaining outside vertex `y`, let `N_y` be its required graph
neighbours among `x_1,...,x_6`, represented on the corresponding seed
coordinates, and write `m=|N_y|`.

If `x_i y` is a required unit edge, equations (2) and (3) give

```text
(u_y)_i = (c_y-sigma z_y)/3 =: t_y.                 (4)
```

If `x_i y` is a candidate nonedge, then `(u_y)_i=0`.  Indeed, every graph in
the corpus has independence number at most two.  Thus candidate-nonadjacent
vertices `x_i,y` cannot share a seed nonneighbour: such a seed vertex would
form an independent triple with them.  Hence `D_{x_i}` and `D_y` are
disjoint.  Coordinate `i` belongs to the actual support of `x_i` and thus to
`D_{x_i}`, so it cannot belong to `D_y`.

Consequently, if `m>0`,

```text
u_y=t_y 1_{N_y},       S_y=N_y.                     (5)
```

The implementation rejects the leaf if this reconstructed support is not a
subset of `D_y`.  If `y` belongs to the other assigned light bin, it also
requires the complete leaf's assigned actual support to equal `N_y`.  For a
point not assigned by this particular support DFS, no equality to a guessed
support is imposed.

Since `c_y=1+m t_y` and, by (4),
`sigma z_y=c_y-3t_y`, substituting (5) in (1) cancels `m` and gives

```text
3t_y^2-2t_y-2=0.
```

Thus there are exactly two possibilities,

```text
t_y=(1+epsilon_y sqrt(7))/3,
epsilon_y in {+1,-1}.                               (6)
```

If `m=0`, then `u_y=0`, `c_y=1`, and (1) gives

```text
z_y=delta_y sigma sqrt(7),
delta_y in {+1,-1}.                                 (7)
```

This is the pair of empty-support K6 apices, expressed relative to the basis
ray sign.

## 4. Exact pair constraints

Take two vertices `y,z` with nonempty basis neighbourhoods.  Put

```text
m=|N_y|,  n=|N_z|,  k=|N_y intersect N_z|,
t=t_y,    s=t_z.
```

Equations (2), (4), and (5) reduce their distance defect to

```text
M_yz = (t+s)/2 + ((m+n-3)/2-k) t s.                 (8)
```

For equal root signs in (6), the irrational and rational coefficients in
`M_yz=0` cannot both vanish.  For opposite signs,

```text
t+s=2/3,       ts=-2/3,
```

so a required unit edge is possible exactly when

```text
m+n-2k=4       and       epsilon_y=-epsilon_z.      (9)
```

Substituting (6) and (7) in (2) gives the remaining required-edge rules:

```text
empty--nonempty:  the nonempty degree is 4 and delta=epsilon;
empty--empty:     impossible.                       (10)
```

No equation is imposed for a candidate nonedge.  This is deliberate: the
candidate graph records required unit edges, not required non-unit pairs.

Finally, two vertices with the same basis neighbourhood and the same root
sign have identical `u,c,z`, hence are the same point.  Distinctness forbids
that collision.  The rule includes the equal-empty-neighbourhood case, so it
also forbids duplicating one of the two empty apices.

## 5. Finite exact CSP and soundness

There are seven outside vertices beyond the six-vector basis.  Assign one
binary root sign to each of them.  The producer exhausts all at most `2^7`
assignments and tests (9), (10), and the collision rule with integer masks and
counts.  The independent verifier instead constructs equal/opposite parity
constraints and checks them by graph propagation, providing a differently
organized transcription of the same exact condition.

Any geometric realization supplies its true `Z0`, ray coloring, actual
supports, and root signs.  The parent enumeration reaches that complete
support leaf.  Sections 2--4 show that its signs satisfy every CSP constraint.
Therefore an empty sign CSP proves that particular leaf impossible.  If all
parent alternatives fail, the seeded candidate is impossible; if one K6 seed
is impossible, the whole candidate graph is impossible.

The two new graph rejections are often witnessed even earlier inside this
extension by a reconstructed nonempty support not contained in the allowed
seed-defect mask.  That remains the same theorem: equation (5) reconstructs
the actual support before any root signs need be tried.

## 6. Frozen source and verification boundary

The production loader pins the exact hashes of the parent producer, report,
certificate archive, and independent `PASS` verification.  It selects the
ordered 251 indices from that report and reconstructs only those records from
the parent input loader.  A theorem run additionally requires all four new
package sources to be committed on `codex/dimension6`, with no tracked or
staged worktree dirt, and records their commit blobs and launch status.

The source-bound discovery census found:

```text
input graphs                                      251
saturated-singleton support leaves checked       291
such leaves rejected by the new rule             291
newly rejected graphs                              2
ordered residue                                  249
```

The production run refuses to emit a completed report if any of those counts
or the ordered rejection/residue hashes differ.  The independent verifier
then reloads the report and compressed certificates, independently replays
all 251 graphs, compares every graph decision and rejected seed certificate,
and reruns all 32 K6 seeds of the known realizable 18-point control.

After the four source files have been committed, the intended commands are:

```bash
python3 -m unittest -v test_d6_k6_saturated_singleton_basis.py
python3 d6_k6_saturated_singleton_basis.py --workers 7
python3 verify_d6_k6_saturated_singleton_basis.py --workers 7
```

The worker count may be changed without changing the exact enumeration.
Producer and verifier artifacts are generated only from the committed source
boundary; discovery output is not a theorem artifact.

## 7. Scope

The package certifies two additional non-realizable level-19 candidate unit
graphs under the existing K6 reduction.  It does not assert realizability of
the 249 survivors, does not settle the K6-only class, and does not by itself
prove the dimension-six extremal value.
