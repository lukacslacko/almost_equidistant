# Five-singleton K6 fixed-remainder obstruction

This source-only discovery package adds one exact necessary condition to the
certified 249-graph K6-only residue:

- producer: `probe_d6_k6_five_singleton_fixed_remainder.py`;
- independent checker: `verify_probe_d6_k6_five_singleton_fixed_remainder.py`;
- focused controls: `test_probe_d6_k6_five_singleton_fixed_remainder.py`.

It rejects nine further graphs and leaves 240.  No numerical optimization,
floating-point comparison, or unproved nonedge condition is used.

## Coordinates and quantifier

Fix a centered unit `K6`, and use the parent package's coordinates

```text
u_i(y)=||y-q_i||^2-1,
c_y=1+sum_i u_i(y),
z_y=sqrt(12) times the normal coordinate.
```

They satisfy

```text
z_y^2=c_y^2+6-6||u_y||^2,                           (1)

||y-w||^2-1=(c_y c_w-z_y z_w)/6-u_y.u_w.            (2)
```

The rule runs only at a complete actual-support leaf after the frozen parent
tests.  Thus it prunes one fully quantified support alternative.  A graph is
rejected only when the parent search exhausts every zero-factor choice,
light-ray coloring, and actual-support assignment for one of its K6 seeds.

For every outside point, the candidate mask `D_y` is only an upper bound on
the actual support `S_y`.  Optional coordinates may vanish.  A candidate
nonedge may still have unit distance.

## Five singleton basis points

Suppose `Z0` is empty and five vertices in one nonzero light-ray bin have
distinct singleton actual supports.  Relabel them as `x_1,...,x_5`, leaving
seed coordinate 6 unused.  Same-ray orthonormality gives

```text
u(x_i)=e_i,  c_{x_i}=2,  z_{x_i}=2 sigma,
sigma in {+1,-1}.                                   (3)
```

For another point `y`, let `N_y` be its required neighbours among the five
basis points, represented as a five-coordinate mask, and let `m=|N_y|`.
Equation (2) with each required edge `x_i y` gives one common value

```text
u_i(y)=(c_y-sigma z_y)/3=:t_y  for i in N_y.        (4)
```

When `x_i y` is a candidate nonedge, coordinate `i` is zero: otherwise seed
vertex `q_i`, basis point `x_i`, and `y` would be an independent triple.  The
input corpus has independence number at most two.  The implementation checks
this local precondition and treats the extension as not applicable rather
than silently forcing an uncontrolled coordinate.

## Fixed remainder points

A remainder point is *fixed* when coordinate 6 is forced to zero, either
because `6` is absent from `D_y` or because the complete support leaf assigns
a support not containing 6.  Then

```text
u_y=t_y 1_{N_y},
c_y=1+m t_y,
sigma z_y=1+(m-3)t_y.
```

Substitution in (1) cancels `m`:

```text
3t_y^2-2t_y-2=0,
t_y=(1+epsilon_y sqrt(7))/3,
epsilon_y in {+1,-1}.                               (5)
```

Both roots are nonzero.  Consequently the actual support of a fixed point is
exactly `N_y`.  The code checks `N_y subseteq D_y`, and checks equality with
the actual support only when the parent leaf assigned that point.  Points
whose sixth coordinate is not fixed are ignored completely; this is a sound
relaxation and keeps the trusted theorem path finite.

## Exact required-edge and collision rules

For two fixed points `y,w`, put

```text
m=|N_y|, n=|N_w|, k=|N_y intersect N_w|,
t=t_y, s=t_w.
```

Equation (2) reduces a required unit edge to

```text
t+s+(m+n-3-2k)ts=0.                                 (6)
```

If their signs in (5) are opposite, then `t+s=2/3` and `ts=-2/3`.
If the signs agree, irrational and rational parts cannot both vanish.
Therefore (6) holds exactly when

```text
|N_y symmetric_difference N_w|=4
and epsilon_y != epsilon_w.                         (7)
```

The producer checks (7) with integer masks.  The independent checker instead
constructs the two SymPy radicals and substitutes them directly in (6).

There is no equation for a candidate nonedge.  Distinctness is separate and
mandatory: two fixed points collide exactly when their neighbourhood masks
and signs are both equal.  The sign enumeration rejects that equality even
when the abstract pair is a candidate nonedge.

## Complete accounting and independent replay

The producer records every labeled K6 seed it actually checks (no orbit
quotient), every complete parent-passing support leaf reaching the hook,
every invoked five-basis extension, and every fixed sign outcome.  Across the
249 graphs it checks 7,477 labeled K6 seeds.  Each per-seed row stores its
exact seed mask, feasibility decision, counter deltas, and extension slice.
Every applicable extension stores all of its branch outcomes, not a sample.

The producer census is

```text
complete parent leaves reaching the hook                         8,785
  rejected first by the older six-basis rule                       278
  Z0 nonempty                                                        33
  Z0 empty                                                        8,474
five-basis extensions                                            1,039
fixed sign branches                                             16,760
  fixed support not allowed                                      2,064
  required pair has wrong symmetric difference                  14,300
  required pair has equal signs                                    396
  assigned-support failures / collisions                             0
passing five-basis branches                                           0
```

The independently ordered support DFS reaches 8,770 leaves, 1,024 five-basis
extensions, and 17,304 sign branches.  The work totals differ only on
surviving graphs because each existential parent DFS can find a different
first witness.  It confirms the identical ordered graph decisions, all
per-graph seed counts, and every production branch-coverage invariant.

Both engines run the known realizable 18-point configuration through the
same final theorem path.  All 32 labeled K6 seeds pass, and the independently
generated positive-seed accounting hashes are identical.

## Result

The new exact rejected indices are

```text
428414, 430135, 1688505, 2669917, 2802776,
3655397, 3931365, 3950363, 3968088.
```

The hashes are

```text
ordered nine-index rejection hash
b7de35407737cc884b7c227e880654c8d738f75889678d3f0d842e5f3696a452

ordered 240-index residue hash
7b9c2e26e2ac178ab4e179dc0623b20f746bbd8eba85afb5c9c2c2013df8956e

producer graph-row hash
f5297d9136b20cd9e307dbbfb8778798cd6d09cf5513dc174398a3e05300df02

independent graph-row hash
51a3fb891230cf26808deaf634d2ac823f23c251f93ebec0560423944936719d

positive 32-seed row hash (both engines)
c340fd60403f18ef8f0f9a0353c4edfa342d637a7581d6466a310ec06ccebb9c
```

The current source-only commands write their large exhaustive transcripts to
`/private/tmp`:

```bash
python3 -m unittest -v \
  test_probe_d6_k6_five_singleton_fixed_remainder.py
python3 probe_d6_k6_five_singleton_fixed_remainder.py --workers 1
python3 verify_probe_d6_k6_five_singleton_fixed_remainder.py --workers 1
```

The broader free-coordinate equations are not part of this trusted rule.  If
`v_y=u_6(y)` is free, the correct conic is

```text
3t_y^2-2(1+v_y)t_y+2v_y^2-2=0,
```

and required edges to known points are linear over `Q(sqrt(7))`.  That is a
possible future extension, but it is unnecessary for these nine exact
rejections and is intentionally excluded from the production/checker code.
