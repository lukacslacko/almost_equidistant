# Exact completion of the K7 basis-kernel diagonal test

This note completes the PSD/rank audit of the saturating-clique
basis-kernel rule.  It is an auxiliary design note and does not modify the K7
implementation.

Fix a candidate cover and a required clique `C` of size `U_K`.  Let
`R=N\C`, let `P=K[C,R]`, and put

```text
K_R = D+A,
```

where `A` is the zero-one adjacency matrix of `G[R]` and `D` is diagonal with
entries strictly greater than one.  Since `C` is a basis,

```text
ker(K_R)=ker(P).
```

Let the columns of the rational matrix `L` form a basis of `ker(P)`.  The
compatibility equations are

```text
(D+A)L=0.
```

Row by row, if `L_i` is nonzero then `-(AL)_i` must be a common scalar
multiple `d_i L_i`, uniquely fixing a rational `d_i>1`.  If `L_i=0`, then
`(AL)_i` must vanish.

## Complete PSD/kernel criterion with partially free diagonal

Let

```text
S={i : L_i != 0},   T=R\S,
H=diag(d_i:i in S)+A[S,S].
```

Then a choice of the still-free diagonal entries `d_i>1`, `i in T`, exists
such that `D+A` is positive semidefinite and has kernel exactly `ker(P)` if
and only if

```text
H is positive semidefinite,
ker(H)=column_space(L[S,:]).
```

Equivalently, after compatibility has passed,

```text
H >= 0,              rank(H)=|S|-dim ker(P).
```

### Necessity

The first condition holds because `H` is a principal submatrix of the
positive-semidefinite Gram matrix `K_R`.  Compatibility gives
`column_space(L[S,:]) subseteq ker(H)`.  If `H` had an additional kernel
vector `w`, positive semidefiniteness of the full block matrix would force
`A[T,S]w=0`: a vector with zero quadratic form in a PSD matrix lies in its
kernel.  Thus `(w,0)` would be an additional kernel vector of `K_R`, contrary
to `ker(K_R)=ker(P)`.

### Sufficiency

The zero rows of `L` and compatibility already give

```text
A[T,S] L[S,:]=0.
```

Under the displayed kernel equality, this says the columns of `A[S,T]` lie
in the range of `H`.  Let `H^+` be its rational generalized inverse on that
range.  The generalized Schur complement is

```text
diag(d_T)+A[T,T]-A[T,S] H^+ A[S,T].
```

Choose every free `d_i` sufficiently large (and hence greater than one) so
this matrix is positive definite.  The completed block matrix is then PSD
and its kernel is exactly `(column_space(L[S,:]),0)=ker(P)`.

Hence checking only the fully constrained case is unnecessary: exact PSD and
rank of the smaller determined principal block give the complete answer for
every compatibility instance.  PSD can be checked by rational `LDL^T` or all
principal minors; rank is exact rational elimination.

## Sample observation

An independent exact-Fraction prototype on the 145 survivors of the current
K7 joint layer found 41 graphs rejected by basis-kernel compatibility.  The
complete PSD/kernel criterion rejected no additional cover and no additional
graph.  All 41 compatibility rejections were already contained in the 59
graphs rejected by the stronger saturating-clique Schur mask rules documented
in `d6_k7_saturating_clique_schur.md`; those mask rules reject 18 further
sample graphs.
