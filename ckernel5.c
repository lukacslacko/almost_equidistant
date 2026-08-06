/* Certified branch-and-prune kernel deciding unit-distance realizability
 * in R^5 — the d=5 generalization of ckernel.c (which certified f(4)=12;
 * algorithm documented in Appendix B of f4_equals_12.pdf, dimension-generic
 * parts unchanged).  Differences from the d=4 kernel, all forced by D=5:
 *   - points are 5-dim interval vectors; adjacency masks are uint32 (n<=22);
 *   - a "sphere" placement needs >=5 placed neighbours (4 difference
 *     equations, 5 unknowns, 1-dim nullspace -> quadratic in t);
 *   - a "circle" placement has exactly 4 placed neighbours (3x3 interval
 *     Gram system for the circumcentre, 2-dim orthogonal complement);
 *   - Krawczyk contraction is 5x5 on the best-conditioned 5-subset;
 *   - coincidence Rule 1 (root separation) uses 5-subsets of common
 *     neighbours; Rule 2 (bisector hyperplane) needs 6 common placed
 *     neighbours with certifiably nonzero affine volume (5x5 interval det).
 * A vertex with exactly 3 placed neighbours would need a 2-sphere stage;
 * the driver only emits orders that avoid this (the kernel aborts loudly
 * if it ever happens, exactly like the <3 case in the d=4 kernel).
 * Interval arithmetic with outward rounding via nextafter (IEEE basic ops
 * are correctly rounded); cos/sin get 8-ulp padding.
 * Compiled automatically by cdriver5.py:
 *   POSIX:   cc -O2 -shared -o ckernel5.so ckernel5.c -lm
 *   Windows: clang -O2 -shared -fuse-ld=lld -o ckernel5.dll ckernel5.c
 *        or  cl /O2 /LD ckernel5.c /Fe:ckernel5.dll (from a VS x64 prompt)
 */
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifdef _WIN32
#define EXPORT __declspec(dllexport)
#else
#define EXPORT
#endif

#define D 5
#define MAXN 22
#define TWO_PI 6.283185307179586476925286766559

typedef struct { double a, b; } IV;

static inline double rlo(double x){ return nextafter(x, -INFINITY); }
static inline double rhi(double x){ return nextafter(x,  INFINITY); }
static inline IV iv(double a, double b){ IV r={a,b}; return r; }
static inline IV ivp(double x){ IV r={x,x}; return r; }
static inline IV iadd(IV x, IV y){ return iv(rlo(x.a+y.a), rhi(x.b+y.b)); }
static inline IV isub(IV x, IV y){ return iv(rlo(x.a-y.b), rhi(x.b-y.a)); }
static inline IV ineg(IV x){ return iv(-x.b, -x.a); }
static inline IV imul(IV x, IV y){
    double p1=x.a*y.a, p2=x.a*y.b, p3=x.b*y.a, p4=x.b*y.b;
    double mn=fmin(fmin(p1,p2),fmin(p3,p4)), mx=fmax(fmax(p1,p2),fmax(p3,p4));
    return iv(rlo(mn), rhi(mx));
}
static inline int contains0(IV x){ return x.a <= 0 && 0 <= x.b; }
static inline IV idivi(IV x, IV y, int *ok){
    if (contains0(y)) { *ok = 0; return ivp(0); }
    double p1=x.a/y.a, p2=x.a/y.b, p3=x.b/y.a, p4=x.b/y.b;
    double mn=fmin(fmin(p1,p2),fmin(p3,p4)), mx=fmax(fmax(p1,p2),fmax(p3,p4));
    *ok = 1; return iv(rlo(mn), rhi(mx));
}
static inline IV isq(IV x){
    if (x.a >= 0) return iv(rlo(x.a*x.a), rhi(x.b*x.b));
    if (x.b <= 0) return iv(rlo(x.b*x.b), rhi(x.a*x.a));
    double m = fmax(-x.a, x.b); return iv(0.0, rhi(m*m));
}
static inline IV isqrtI(IV x, int *ok){
    if (x.b < 0){ *ok=0; return ivp(0); }
    double a = x.a > 0 ? x.a : 0.0;
    *ok=1; return iv(rlo(sqrt(a)), rhi(sqrt(x.b)));
}
static inline double wid(IV x){ return x.b - x.a; }
static inline double mid(IV x){ return 0.5*(x.a+x.b); }
static inline double mig(IV x){
    if (contains0(x)) return 0.0;
    double aa=fabs(x.a), bb=fabs(x.b); return aa<bb?aa:bb;
}
static double pad8lo(double x){ for(int i=0;i<8;i++) x=nextafter(x,-INFINITY); return x; }
static double pad8hi(double x){ for(int i=0;i<8;i++) x=nextafter(x, INFINITY); return x; }
static IV icosI(IV x){
    if (wid(x) >= TWO_PI) return iv(-1.0, 1.0);
    double ca=cos(x.a), cb=cos(x.b);
    double lo = fmin(pad8lo(ca), pad8lo(cb));
    double hi = fmax(pad8hi(ca), pad8hi(cb));
    long k0 = (long)floor(x.a / M_PI), k1 = (long)floor(x.b / M_PI);
    for (long k=k0; k<=k1; k++){
        double t = k * M_PI;
        if (x.a <= t && t <= x.b){ if (k % 2 == 0) hi = 1.0; else lo = -1.0; }
    }
    if (lo < -1.0) lo = -1.0; if (hi > 1.0) hi = 1.0;
    return iv(lo, hi);
}
static IV isinI(IV x){ return icosI(isub(iv(M_PI/2, M_PI/2), x)); }

typedef IV V5[D];
static void vsubv(const V5 u, const V5 v, V5 r){ for(int i=0;i<D;i++) r[i]=isub(u[i],v[i]); }
static void vaddv(const V5 u, const V5 v, V5 r){ for(int i=0;i<D;i++) r[i]=iadd(u[i],v[i]); }
static void vscalev(IV t, const V5 u, V5 r){ for(int i=0;i<D;i++) r[i]=imul(t,u[i]); }
static IV vdotv(const V5 u, const V5 v){
    IV s=ivp(0); for(int i=0;i<D;i++) s=iadd(s, imul(u[i],v[i])); return s;
}
static IV vnorm2v(const V5 u){
    IV s=ivp(0); for(int i=0;i<D;i++) s=iadd(s, isq(u[i])); return s;
}
static int overlap5(const V5 x, const V5 y){
    for(int i=0;i<D;i++) if (x[i].a > y[i].b || y[i].a > x[i].b) return 0;
    return 1;
}

/* ---- generic interval Gauss-Jordan: m rows, nuns unknowns (<=5) ---- */
static int gauss(IV A[D][D], IV b[D], int m, int nuns, IV xp[D], IV nullv[D], int *nnull){
    int usedr[D]={0}, usedc[D]={0};
    int piv_r[D], piv_c[D];
    for (int step=0; step<m; step++){
        double best=0.0; int br=-1, bc=-1;
        for (int i=0;i<m;i++){ if (usedr[i]) continue;
            for (int j=0;j<nuns;j++){ if (usedc[j]) continue;
                double g = mig(A[i][j]);
                if (g > best){ best=g; br=i; bc=j; } } }
        if (br < 0 || best == 0.0) return 0;
        usedr[br]=1; usedc[bc]=1; piv_r[step]=br; piv_c[step]=bc;
        for (int i=0;i<m;i++){
            if (i==br) continue;
            int ok; IV f = idivi(A[i][bc], A[br][bc], &ok);
            if (!ok) return 0;
            for (int j=0;j<nuns;j++) A[i][j] = isub(A[i][j], imul(f, A[br][j]));
            b[i] = isub(b[i], imul(f, b[br]));
            A[i][bc] = ivp(0);
        }
    }
    int freec[D], nf=0;
    for (int j=0;j<nuns;j++) if (!usedc[j]) freec[nf++]=j;
    /* particular solution: free vars = 0 */
    for (int j=0;j<nuns;j++) xp[j]=ivp(0);
    for (int s=m-1;s>=0;s--){
        int i=piv_r[s], j=piv_c[s];
        IV acc = b[i];
        for (int jj=0;jj<nuns;jj++) if (jj!=j) acc = isub(acc, imul(A[i][jj], xp[jj]));
        int ok; xp[j] = idivi(acc, A[i][j], &ok);
        if (!ok) return 0;
    }
    *nnull = nf;
    if (nf == 1){
        for (int j=0;j<nuns;j++) nullv[j]=ivp(0);
        nullv[freec[0]] = ivp(1.0);
        for (int s=m-1;s>=0;s--){
            int i=piv_r[s], j=piv_c[s];
            IV acc = ivp(0);
            for (int jj=0;jj<nuns;jj++) if (jj!=j) acc = isub(acc, imul(A[i][jj], nullv[jj]));
            int ok; nullv[j] = idivi(acc, A[i][j], &ok);
            if (!ok) return 0;
        }
    }
    return 1;
}

/* ---- problem state ---- */
typedef struct {
    int n;
    uint32_t adj[MAXN];
    int seed_len, seed[D+1];
    int order_len, order[MAXN];
    double seedc[D+1][D][2];
    double theta_min, max_width;
    int64_t max_nodes;
    double th0_lo, th0_hi;
    /* runtime */
    int64_t nodes, unresolved;
    int aborted, th0_used;
    /* per-cell node quota: a circle cell whose subtree exceeds the quota
       is abandoned and split (management only; never affects kill
       soundness). cellq_limit==0 means no quota active. */
    int64_t cellq_start, cellq_limit;
    int cellq_over;
} Prob;

#define CELL_QUOTA 4000

typedef struct { V5 p[MAXN]; uint32_t placedmask; int porder[MAXN]; int np; } State;

static void dfs(Prob *P, State *st, int oi, double cellw, int has_cellw);

/* place via 5 spheres: cands out; returns -1 illcond, else count (0..2);
   seg: xp,n,t for each cand */
typedef struct { V5 x; V5 xp, nv; IV t; } Cand;
static int place_sphere_c(const V5 nb[D], Cand out[2]){
    IV A[D][D], b[D]; V5 d;
    for (int r=1;r<D;r++){
        for (int c=0;c<D;c++) A[r-1][c] = imul(ivp(2.0), isub(nb[r][c], nb[0][c]));
        b[r-1] = isub(vnorm2v(nb[r]), vnorm2v(nb[0]));
    }
    IV xp[D], nv[D]; int nn;
    if (!gauss(A, b, D-1, D, xp, nv, &nn)) return -1;
    if (nn != 1) return -1;
    vsubv(xp, nb[0], d);
    IV a = vnorm2v(nv), bq = imul(ivp(2.0), vdotv(nv, d)), c = isub(vnorm2v(d), ivp(1.0));
    if (!(a.a > 0)) return -1;
    IV disc = isub(isq(bq), imul(ivp(4.0), imul(a, c)));
    if (disc.b < 0) return 0;
    int ok; IV sq = isqrtI(disc, &ok);
    if (!ok) return 0;
    if (contains0(disc)){
        int ok2; IV t = idivi(iadd(ineg(bq), iv(-sq.b, sq.b)), imul(ivp(2.0), a), &ok2);
        if (!ok2) return -1;
        V5 tn; vscalev(t, nv, tn); vaddv(xp, tn, out[0].x);
        memcpy(out[0].xp, xp, sizeof(V5)); memcpy(out[0].nv, nv, sizeof(V5)); out[0].t = t;
        return 1;
    }
    for (int s=0;s<2;s++){
        int ok2; IV t = idivi(iadd(ineg(bq), s? ineg(sq): sq), imul(ivp(2.0), a), &ok2);
        if (!ok2) return -1;
        V5 tn; vscalev(t, nv, tn); vaddv(xp, tn, out[s].x);
        memcpy(out[s].xp, xp, sizeof(V5)); memcpy(out[s].nv, nv, sizeof(V5)); out[s].t = t;
    }
    return 2;
}

static int edge_ok(const V5 x, const V5 p){
    V5 d; vsubv(x, p, d);
    IV e = isub(vnorm2v(d), ivp(1.0));
    return contains0(e);
}

/* float 5x5 inverse; returns 0 on singular */
static int inv5(double M[D][D], double Y[D][D]){
    double a[D][2*D];
    for (int i=0;i<D;i++){ for (int j=0;j<D;j++){ a[i][j]=M[i][j]; a[i][j+D]=(i==j);} }
    for (int c=0;c<D;c++){
        int p=-1; double best=0;
        for (int r=c;r<D;r++){ double v=fabs(a[r][c]); if (v>best){best=v;p=r;} }
        if (p<0 || best<1e-14) return 0;
        if (p!=c){ for (int j=0;j<2*D;j++){ double t=a[c][j]; a[c][j]=a[p][j]; a[p][j]=t; } }
        double d=a[c][c];
        for (int j=0;j<2*D;j++) a[c][j]/=d;
        for (int r=0;r<D;r++){ if (r==c) continue; double f=a[r][c];
            for (int j=0;j<2*D;j++) a[r][j]-=f*a[c][j]; }
    }
    for (int i=0;i<D;i++) for (int j=0;j<D;j++) Y[i][j]=a[i][j+D];
    return 1;
}

/* choose a well-conditioned 5-subset of neighbour indices: greedy pivoted
   Gram-Schmidt on the midpoint difference vectors (O(k*D^2), deterministic;
   selection only affects performance, never soundness) */
static void best_quint(const V5 x, const V5 nb[], int k, int q[D]){
    double xm[D]; for (int i=0;i<D;i++) xm[i]=mid(x[i]);
    if (k==D){ for (int i=0;i<D;i++) q[i]=i; return; }
    double rows[MAXN][D]; int used[MAXN];
    for (int r=0;r<k;r++){ used[r]=0;
        for (int c=0;c<D;c++) rows[r][c] = 2*(xm[c]-mid(nb[r][c])); }
    int got=0;
    for (int s=0;s<D;s++){
        int best=-1; double bn=-1;
        for (int r=0;r<k;r++){ if (used[r]) continue;
            double n2=0; for (int c=0;c<D;c++) n2 += rows[r][c]*rows[r][c];
            if (n2 > bn){ bn=n2; best=r; } }
        if (best<0) break;
        used[best]=1; q[got++]=best;
        double nn = bn > 1e-300 ? bn : 1e-300;
        for (int r=0;r<k;r++){ if (used[r]) continue;
            double d=0; for (int c=0;c<D;c++) d += rows[r][c]*rows[best][c];
            d /= nn;
            for (int c=0;c<D;c++) rows[r][c] -= d*rows[best][c]; }
    }
    for (int r=0; got<D && r<k; r++)
        if (!used[r]){ used[r]=1; q[got++]=r; }
}

/* Krawczyk contraction on 5 chosen spheres; 1=ok (x updated), 0=empty */
static int krawczyk(V5 x, const V5 nb[], int k, int iters){
    int q[D]; best_quint(x, nb, k, q);
    const IV *ps[D]; for (int i=0;i<D;i++) ps[i]=nb[q[i]];
    for (int it=0; it<iters; it++){
        double mm[D]; for (int i=0;i<D;i++) mm[i]=mid(x[i]);
        double Jm[D][D], Y[D][D];
        for (int r=0;r<D;r++) for (int c=0;c<D;c++) Jm[r][c]=2*(mm[c]-mid(ps[r][c]));
        if (!inv5(Jm, Y)) return 1;
        IV m[D]; for (int i=0;i<D;i++) m[i]=ivp(mm[i]);
        IV Fm[D];
        for (int r=0;r<D;r++){
            V5 dd; for (int c=0;c<D;c++) dd[c]=isub(m[c], ps[r][c]);
            Fm[r] = isub(vnorm2v(dd), ivp(1.0));
        }
        IV K[D]; int shrunk=0;
        IV JX[D][D];
        for (int r=0;r<D;r++) for (int c=0;c<D;c++)
            JX[r][c] = imul(ivp(2.0), isub(x[c], ps[r][c]));
        IV dxm[D]; for (int i=0;i<D;i++) dxm[i]=isub(x[i], m[i]);
        for (int i=0;i<D;i++){
            IV s = m[i];
            for (int j=0;j<D;j++) s = isub(s, imul(ivp(Y[i][j]), Fm[j]));
            for (int j=0;j<D;j++){
                IV acc = (i==j)? ivp(1.0): ivp(0.0);
                for (int c=0;c<D;c++) acc = isub(acc, imul(ivp(Y[i][c]), JX[c][j]));
                s = iadd(s, imul(acc, dxm[j]));
            }
            K[i]=s;
        }
        for (int i=0;i<D;i++){
            double lo = fmax(x[i].a, K[i].a), hi = fmin(x[i].b, K[i].b);
            if (lo > hi) return 0;
            if (hi-lo < wid(x[i]) - 1e-15) shrunk=1;
            x[i]=iv(lo,hi);
        }
        if (!shrunk) break;
    }
    return 1;
}

/* interval determinants: 3x3, 4x4 (cofactor), 5x5 (cofactor over 4x4) */
static IV det3I(IV m[3][3]){
    IV t1 = imul(m[0][0], isub(imul(m[1][1],m[2][2]), imul(m[1][2],m[2][1])));
    IV t2 = imul(m[0][1], isub(imul(m[1][0],m[2][2]), imul(m[1][2],m[2][0])));
    IV t3 = imul(m[0][2], isub(imul(m[1][0],m[2][1]), imul(m[1][1],m[2][0])));
    return iadd(isub(t1, t2), t3);
}
static IV det4I_m(IV m[4][4]){
    IV tot = ivp(0);
    for (int j=0;j<4;j++){
        int cols[3], nc=0;
        for (int k2=0;k2<4;k2++) if (k2!=j) cols[nc++]=k2;
        IV s[3][3];
        for (int r=1;r<4;r++) for (int c=0;c<3;c++) s[r-1][c]=m[r][cols[c]];
        IV term = imul(m[0][j], det3I(s));
        if (j % 2) term = ineg(term);
        tot = iadd(tot, term);
    }
    return tot;
}
/* interval det of 5 difference vectors (rows) */
static IV det5I(V5 rows[D]){
    IV tot = ivp(0);
    for (int j=0;j<D;j++){
        int cols[4], nc=0;
        for (int k2=0;k2<D;k2++) if (k2!=j) cols[nc++]=k2;
        IV s[4][4];
        for (int r=1;r<D;r++) for (int c=0;c<4;c++) s[r-1][c]=rows[r][cols[c]];
        IV term = imul(rows[0][j], det4I_m(s));
        if (j % 2) term = ineg(term);
        tot = iadd(tot, term);
    }
    return tot;
}

/* greedy farthest-point ordering of placed vertices (by box midpoints):
   well-spread subsets give well-conditioned determinants and sphere
   systems; selection only affects which subsets are tried, never
   soundness */
static int spread_order(State *st, const int *common, int nc, int *ord){
    double mids[MAXN][D];
    int used[MAXN];
    for (int i=0;i<nc;i++){ used[i]=0;
        for (int c=0;c<D;c++) mids[i][c]=mid(st->p[common[i]][c]); }
    ord[0]=0; used[0]=1;
    for (int s=1;s<nc;s++){
        double bd=-1; int bi=-1;
        for (int i=0;i<nc;i++){
            if (used[i]) continue;
            double mind=1e300;
            for (int t=0;t<s;t++){
                double d2=0;
                for (int c=0;c<D;c++){
                    double dd=mids[i][c]-mids[ord[t]][c]; d2+=dd*dd; }
                if (d2<mind) mind=d2;
            }
            if (mind>bd){ bd=mind; bi=i; }
        }
        ord[s]=bi; used[bi]=1;
    }
    return nc;
}

/* Rule 2: 6 placed common neighbours of u,v with certifiably nonzero
   affine volume force x_u = x_v (no injective realization on branch) */
static int forced_coincident_c(Prob *P, State *st, int u, int v){
    int common[MAXN], nc=0;
    for (int wpi=0; wpi<st->np; wpi++){
        int w = st->porder[wpi];
        if (w==u || w==v) continue;
        if (((P->adj[v]>>w)&1) && ((P->adj[u]>>w)&1)) common[nc++]=w;
    }
    if (nc < D+1) return 0;
    int ord[MAXN];
    spread_order(st, common, nc, ord);
    int m = nc < 9 ? nc : 9;      /* 6-subsets of the 9 best-spread */
    int idx[D+1], tried = 0;
    for (idx[0]=0; idx[0]<m-5; idx[0]++)
    for (idx[1]=idx[0]+1; idx[1]<m-4; idx[1]++)
    for (idx[2]=idx[1]+1; idx[2]<m-3; idx[2]++)
    for (idx[3]=idx[2]+1; idx[3]<m-2; idx[3]++)
    for (idx[4]=idx[3]+1; idx[4]<m-1; idx[4]++)
    for (idx[5]=idx[4]+1; idx[5]<m; idx[5]++){
        if (++tried > 64) return 0;   /* rule only ever discards */
        V5 rows[D];
        for (int r=0;r<D;r++)
            vsubv(st->p[common[ord[idx[r+1]]]], st->p[common[ord[idx[0]]]], rows[r]);
        IV dt = det5I(rows);
        if (dt.a > 0 || dt.b < 0) return 1;
    }
    return 0;
}

/* noninjective placement filter (Rule 1 root separation + Rule 2) */
static int noninjective_c(Prob *P, State *st, int v, const V5 x){
    for (int upi=0; upi<st->np; upi++){
        int u = st->porder[upi];
        if ((P->adj[v]>>u)&1) continue;
        if (!overlap5(x, st->p[u])) continue;
        /* hyperplane rule with v's box in place */
        State tmp = *st;
        memcpy(tmp.p[v], x, sizeof(V5));
        tmp.placedmask |= (uint32_t)1<<v;
        tmp.porder[tmp.np++] = v;
        if (forced_coincident_c(P, &tmp, u, v)) return 1;
        /* root-separation rule */
        int common[MAXN], nc=0;
        for (int wpi=0; wpi<st->np; wpi++){
            int w = st->porder[wpi];
            if (w==u) continue;
            if (((P->adj[v]>>w)&1) && ((P->adj[u]>>w)&1)) common[nc++]=w;
        }
        if (nc < D) continue;
        int ordu[MAXN];
        spread_order(st, common, nc, ordu);
        int mu = nc < 8 ? nc : 8;   /* 5-subsets of the 8 best-spread */
        int idx[D], tried = 0;
        for (idx[0]=0; idx[0]<mu-4; idx[0]++)
        for (idx[1]=idx[0]+1; idx[1]<mu-3; idx[1]++)
        for (idx[2]=idx[1]+1; idx[2]<mu-2; idx[2]++)
        for (idx[3]=idx[2]+1; idx[3]<mu-1; idx[3]++)
        for (idx[4]=idx[3]+1; idx[4]<mu; idx[4]++){
            if (++tried > 56) goto next_u;      /* cap attempts */
            V5 nbv[D];
            for (int r=0;r<D;r++) memcpy(nbv[r], st->p[common[ordu[idx[r]]]], sizeof(V5));
            Cand cd[2];
            int r = place_sphere_c((const V5*)nbv, cd);
            if (r < 0) continue;
            if (r == 0) return 1;               /* v must solve: empty */
            if (r != 2) continue;
            if (overlap5(cd[0].x, cd[1].x)) continue;
            int d1 = overlap5(st->p[u], cd[0].x), d2 = overlap5(st->p[u], cd[1].x);
            int c1 = overlap5(x, cd[0].x), c2 = overlap5(x, cd[1].x);
            if (!d1 && !d2) return 1;
            if (!c1 && !c2) return 1;
            if ((d1 && !d2 && c1 && !c2) || (d2 && !d1 && c2 && !c1)) return 1;
            goto next_u;
        }
        next_u: ;
    }
    return 0;
}

static int local_sweep_c(Prob *P, State *st, int v0){
    int frontier[MAXN], nf=1; frontier[0]=v0;
    uint32_t seedm=0;
    for (int i=0;i<P->seed_len;i++) seedm |= (uint32_t)1<<P->seed[i];
    for (int round=0; round<2 && nf>0; round++){
        int nxt[MAXN], nn=0;
        for (int fi=0; fi<nf; fi++){
            int v = frontier[fi];
            for (int upi=0; upi<st->np; upi++){
                int u = st->porder[upi];
                if ((seedm>>u)&1) continue;
                if (!((P->adj[v]>>u)&1)) continue;
                double mw=0;
                for (int i=0;i<D;i++) mw=fmax(mw, wid(st->p[u][i]));
                if (mw < 1e-9) continue;
                V5 nbv[MAXN]; int k=0;
                int nbrids[MAXN];
                for (int wpi=0; wpi<st->np; wpi++){
                    int w = st->porder[wpi];
                    if ((P->adj[u]>>w)&1){
                        memcpy(nbv[k], st->p[w], sizeof(V5)); nbrids[k]=w; k++; } }
                if (k < D) continue;
                V5 xx; memcpy(xx, st->p[u], sizeof(V5));
                if (!krawczyk(xx, (const V5*)nbv, k, 1)) return 0;
                int chg=0;
                for (int i=0;i<D;i++)
                    if (wid(xx[i]) < wid(st->p[u][i]) - 1e-14) chg=1;
                memcpy(st->p[u], xx, sizeof(V5));
                if (chg && nn < MAXN) nxt[nn++]=u;
                for (int w2=0; w2<k; w2++)
                    if (!edge_ok(st->p[u], st->p[nbrids[w2]])) return 0;
            }
        }
        memcpy(frontier, nxt, nn*sizeof(int)); nf=nn;
    }
    return 1;
}

static int injective_or_unresolved_c(Prob *P, State *st){
    for (int upi=0; upi<st->np; upi++){
        int u = st->porder[upi];
        for (int vpi=0; vpi<st->np; vpi++){
            int v = st->porder[vpi];
            if (u >= v) continue;
            if ((P->adj[v]>>u)&1) continue;
            if (!overlap5(st->p[u], st->p[v])) continue;
            if (forced_coincident_c(P, st, u, v)) return 0;
            int common[MAXN], nc=0;
            for (int wpi=0; wpi<st->np; wpi++){
                int w = st->porder[wpi];
                if (w==u||w==v) continue;
                if (((P->adj[v]>>w)&1) && ((P->adj[u]>>w)&1)) common[nc++]=w;
            }
            if (nc < D) continue;
            int ordp[MAXN];
            spread_order(st, common, nc, ordp);
            int mp = nc < 8 ? nc : 8;
            int idx[D], tried = 0;
            for (idx[0]=0; idx[0]<mp-4; idx[0]++)
            for (idx[1]=idx[0]+1; idx[1]<mp-3; idx[1]++)
            for (idx[2]=idx[1]+1; idx[2]<mp-2; idx[2]++)
            for (idx[3]=idx[2]+1; idx[3]<mp-1; idx[3]++)
            for (idx[4]=idx[3]+1; idx[4]<mp; idx[4]++){
                if (++tried > 56) goto next_pair;   /* cap attempts */
                V5 nbv[D];
                for (int r=0;r<D;r++) memcpy(nbv[r], st->p[common[ordp[idx[r]]]], sizeof(V5));
                Cand cd[2];
                int r = place_sphere_c((const V5*)nbv, cd);
                if (r < 0) continue;
                if (r == 0) return 0;
                if (r != 2) continue;
                if (overlap5(cd[0].x, cd[1].x)) continue;
                int du1=overlap5(st->p[u], cd[0].x), du2=overlap5(st->p[u], cd[1].x);
                int dv1=overlap5(st->p[v], cd[0].x), dv2=overlap5(st->p[v], cd[1].x);
                if ((!du1&&!du2)||(!dv1&&!dv2)) return 0;
                if (du1==dv1 && du2==dv2 && du1+du2==1) return 0;
                goto next_pair;
            }
            next_pair: ;
        }
    }
    return 1;
}

static int final_sweep_c(Prob *P, State *stt){
    State s = *stt;
    uint32_t seedm=0;
    for (int i=0;i<P->seed_len;i++) seedm |= (uint32_t)1<<P->seed[i];
    for (int sweep=0; sweep<6; sweep++){
        int changed=0;
        for (int vpi=0; vpi<s.np; vpi++){
            int v = s.porder[vpi];
            if ((seedm>>v)&1) continue;
            double mw=0; for (int i=0;i<D;i++) mw=fmax(mw, wid(s.p[v][i]));
            if (mw < 1e-9) continue;
            V5 nbv[MAXN]; int k=0; int ids[MAXN];
            for (int wpi=0; wpi<s.np; wpi++){
                int w = s.porder[wpi];
                if ((P->adj[v]>>w)&1){
                    memcpy(nbv[k], s.p[w], sizeof(V5)); ids[k]=w; k++; } }
            if (k < D) continue;
            V5 xx; memcpy(xx, s.p[v], sizeof(V5));
            if (!krawczyk(xx, (const V5*)nbv, k, 2)) return 0;
            for (int i=0;i<D;i++)
                if (wid(xx[i]) < wid(s.p[v][i]) - 1e-15) changed=1;
            memcpy(s.p[v], xx, sizeof(V5));
            for (int w2=0; w2<k; w2++)
                if (!edge_ok(s.p[v], s.p[ids[w2]])) return 0;
        }
        if (!changed) break;
    }
    return injective_or_unresolved_c(P, &s);
}

/* ---- circle & segment stages ---- */
typedef struct { double t0, t1; } Cell;

/* circle through intersection of 4 unit spheres centred at p0..p3:
   circumcentre c solves the 3x3 Gram system, radius r = sqrt(1-R^2),
   u1,u2 = orthonormal basis of the orthogonal complement of the
   affine hull (dim 3) of {p0..p3}. Returns 1 ok, 0 certified empty,
   -1 ill-conditioned. */
static int circle_setup_c(const V5 p0, const V5 p1, const V5 p2, const V5 p3,
                          V5 c, IV *r, V5 u1, V5 u2){
    V5 dv[3];
    vsubv(p1, p0, dv[0]); vsubv(p2, p0, dv[1]); vsubv(p3, p0, dv[2]);
    IV A[D][D], b[D];
    for (int i=0;i<3;i++){
        for (int j=0;j<3;j++) A[i][j] = vdotv(dv[i], dv[j]);
        int okd; b[i] = idivi(vdotv(dv[i], dv[i]), ivp(2.0), &okd);
    }
    IV ab[D], nl[D]; int nn;
    if (!gauss(A, b, 3, 3, ab, nl, &nn)) return -1;
    V5 s; for (int i=0;i<D;i++) s[i]=ivp(0);
    for (int q2=0;q2<3;q2++){
        V5 t; vscalev(ab[q2], dv[q2], t); vaddv(s, t, s);
    }
    vaddv(p0, s, c);
    V5 cd; vsubv(c, p0, cd);
    IV R2 = vnorm2v(cd);
    IV r2 = isub(ivp(1.0), R2);
    if (r2.b < 0) return 0;                     /* certified empty */
    int ok; *r = isqrtI(r2, &ok);
    if (!ok) return 0;
    /* orthonormalize hull dirs, then complement */
    V5 ob[3]; int no=0;
    for (int q2=0; q2<3; q2++){
        V5 e; memcpy(e, dv[q2], sizeof(V5));
        for (int b2=0;b2<no;b2++){
            IV pr = vdotv(e, ob[b2]); V5 t; vscalev(pr, ob[b2], t); vsubv(e, t, e);
        }
        IV nrm = vnorm2v(e);
        if (!(nrm.a > 0) || mig(nrm) <= 1e-8) return -1;
        int ok2; IV s2 = isqrtI(nrm, &ok2);
        for (int i2=0;i2<D;i2++){ int ok3; ob[no][i2]=idivi(e[i2], s2, &ok3); if(!ok3) return -1; }
        no++;
    }
    int got=0;
    for (int k=0;k<D && got<2;k++){
        V5 e; for (int i2=0;i2<D;i2++) e[i2] = (i2==k)? ivp(1.0): ivp(0.0);
        for (int b2=0;b2<no;b2++){
            IV pr=vdotv(e, ob[b2]); V5 t; vscalev(pr, ob[b2], t); vsubv(e, t, e);
        }
        if (got==1){
            IV pr=vdotv(e, u1); V5 t; vscalev(pr, u1, t); vsubv(e, t, e);
        }
        IV nrm = vnorm2v(e);
        if (nrm.a > 0 && mig(nrm) > 1e-6){
            int ok2; IV s2 = isqrtI(nrm, &ok2);
            IV *dst = got==0 ? u1 : u2;
            for (int i2=0;i2<D;i2++){ int ok3; dst[i2]=idivi(e[i2], s2, &ok3); if(!ok3) return -1; }
            got++;
        }
    }
    if (got != 2) return -1;
    return 1;
}

static void segment_stage_c(Prob *P, State *st, int oi, int v,
                            int nbrids[], int nnb, Cand *seg,
                            double cellw, int has_cellw);

static void circle_stage_c(Prob *P, State *st, int oi, int v,
                           int nbrids[], int nnb, double outerw, int has_outer){
    /* find working quadruple */
    V5 c, u1, u2; IV r;
    int rquad[4] = {-1,-1,-1,-1};
    int found = 0;
    int idx[4];
    for (idx[0]=0; idx[0]<nnb-3 && !found; idx[0]++)
    for (idx[1]=idx[0]+1; idx[1]<nnb-2 && !found; idx[1]++)
    for (idx[2]=idx[1]+1; idx[2]<nnb-1 && !found; idx[2]++)
    for (idx[3]=idx[2]+1; idx[3]<nnb && !found; idx[3]++){
        int rc = circle_setup_c(st->p[nbrids[idx[0]]], st->p[nbrids[idx[1]]],
                                st->p[nbrids[idx[2]]], st->p[nbrids[idx[3]]],
                                c, &r, u1, u2);
        if (rc == 0) return;                     /* certified empty: kill */
        if (rc == 1){ rquad[0]=nbrids[idx[0]]; rquad[1]=nbrids[idx[1]];
                      rquad[2]=nbrids[idx[2]]; rquad[3]=nbrids[idx[3]]; found=1; }
    }
    if (!found){ P->unresolved++; return; }      /* degenerate: park */
    double floor_ = P->theta_min;
    if (has_outer && outerw/4.0 > floor_) floor_ = outerw/4.0;
    double lo0 = 0.0, hi0 = TWO_PI;
    if (!has_outer && !P->th0_used && (P->th0_lo != 0.0 || P->th0_hi != TWO_PI)){
        lo0 = P->th0_lo; hi0 = P->th0_hi; P->th0_used = 1;
    }
    int cap = 1<<20;
    Cell *stack = malloc(sizeof(Cell)*cap);
    int sp = 0;
    const int NINIT = 32;
    for (int i=0;i<NINIT;i++){
        stack[sp].t0 = lo0 + i*(hi0-lo0)/NINIT;
        stack[sp].t1 = lo0 + (i+1)*(hi0-lo0)/NINIT; sp++;
    }
    while (sp > 0){
        if (P->aborted || P->cellq_over) break;
        Cell cl = stack[--sp];
        /* own the per-cell quota if no outer stage holds one; cells at
           the width floor run quota-free (they cannot split further, so
           abandoning them would only park work that might complete) */
        int own = (P->cellq_limit == 0 && (cl.t1 - cl.t0) > floor_);
        if (own){
            P->cellq_limit = CELL_QUOTA;
            P->cellq_start = P->nodes;
            P->cellq_over = 0;
        }
        IV th = iv(cl.t0, cl.t1);
        V5 x, tc, ts, sum;
        IV cth = icosI(th), sth = isinI(th);
        IV rc_ = imul(r, cth), rs_ = imul(r, sth);
        vscalev(rc_, u1, tc); vscalev(rs_, u2, ts);
        vaddv(tc, ts, sum); vaddv(c, sum, x);
        /* certified mean-value form: for every real theta in the cell and
           every real configuration in the branch, coordinatewise
             x_i(theta) = x_i(tm) + x_i'(xi) (theta - tm)   for some xi in
           the cell; enclosing x_i(tm) and x_i' over the whole cell gives a
           second valid enclosure, O(w^2)-tight on narrow cells. Intersect
           with the direct form (both contain the true value). */
        {
            double tmv = 0.5*(cl.t0 + cl.t1);
            IV thm = iv(tmv, tmv);
            IV cm = icosI(thm), sm = isinI(thm);
            IV dth = isub(th, thm);
            int deadmv = 0;
            for (int i=0;i<D && !deadmv;i++){
                IV base = iadd(c[i], imul(r,
                              iadd(imul(cm, u1[i]), imul(sm, u2[i]))));
                IV deriv = imul(r, iadd(imul(ineg(sth), u1[i]),
                                        imul(cth, u2[i])));
                IV xmv = iadd(base, imul(deriv, dth));
                double lo = fmax(x[i].a, xmv.a), hi = fmin(x[i].b, xmv.b);
                if (lo > hi){ deadmv = 1; break; }
                x[i] = iv(lo, hi);
            }
            if (deadmv) goto cell_done;  /* enclosures disjoint: empty */
        }
        {
        int dead = 0;
        for (int w=0; w<nnb; w++){
            if (!edge_ok(x, st->p[nbrids[w]])){
                int isdef = (nbrids[w]==rquad[0]||nbrids[w]==rquad[1]
                          || nbrids[w]==rquad[2]||nbrids[w]==rquad[3]);
                if (isdef){ fprintf(stderr, "BUG circle defining violated\n"); exit(9); }
                dead = 1; break;
            }
        }
        if (dead) goto cell_done;
        int64_t before = P->unresolved;
        State st2 = *st;
        memcpy(st2.p[v], x, sizeof(V5));
        st2.placedmask |= (uint32_t)1<<v;
        st2.porder[st2.np++] = v;
        if (!local_sweep_c(P, &st2, v)) goto cell_done;
        dfs(P, &st2, oi+1, cl.t1-cl.t0, 1);
        if (own && P->cellq_over){
            /* subtree too expensive at this width: abandon and split
               (children re-explore; partial parks are discarded) */
            P->cellq_over = 0;
            if (cl.t1-cl.t0 > floor_ && sp+2 <= cap){
                P->unresolved = before;
                double tm = 0.5*(cl.t0+cl.t1);
                stack[sp].t0=cl.t0; stack[sp].t1=tm; sp++;
                stack[sp].t0=tm; stack[sp].t1=cl.t1; sp++;
            } else {
                P->unresolved = before + 1;   /* park at width floor */
            }
            goto cell_done;
        }
        if (P->unresolved > before){
            if (cl.t1-cl.t0 > floor_ && sp+2 <= cap){
                P->unresolved = before;
                double tm = 0.5*(cl.t0+cl.t1);
                stack[sp].t0=cl.t0; stack[sp].t1=tm; sp++;
                stack[sp].t0=tm; stack[sp].t1=cl.t1; sp++;
            }
        }
        }
      cell_done:
        if (own){ P->cellq_limit = 0; P->cellq_over = 0; }
    }
    free(stack);
}

static void segment_stage_c(Prob *P, State *st, int oi, int v,
                            int nbrids[], int nnb, Cand *seg,
                            double cellw, int has_cellw){
    double lo = seg->t.a, hi = seg->t.b;
    if (!(hi > lo)) return;
    double floor_ = 1e-7;
    if (has_cellw && cellw/4.0 > floor_) floor_ = cellw/4.0;
    int cap = 1<<18;
    Cell *stack = malloc(sizeof(Cell)*cap);
    int sp=0;
    const int NIN=16;
    for (int i=0;i<NIN;i++){
        stack[sp].t0 = lo + i*(hi-lo)/NIN;
        stack[sp].t1 = lo + (i+1)*(hi-lo)/NIN; sp++;
    }
    while (sp > 0){
        if (P->aborted || P->cellq_over) break;
        Cell cl = stack[--sp];
        IV tc = iv(cl.t0, cl.t1);
        V5 x, tn; vscalev(tc, seg->nv, tn); vaddv(seg->xp, tn, x);
        int dead=0;
        for (int w=0; w<nnb; w++)
            if (!edge_ok(x, st->p[nbrids[w]])){ dead=1; break; }
        if (dead) continue;
        if (noninjective_c(P, st, v, x)) continue;
        int64_t before = P->unresolved;
        State st2 = *st;
        memcpy(st2.p[v], x, sizeof(V5));
        st2.placedmask |= (uint32_t)1<<v;
        st2.porder[st2.np++] = v;
        if (!local_sweep_c(P, &st2, v)) continue;
        dfs(P, &st2, oi+1, cellw, has_cellw);
        if (P->unresolved > before){
            if (cl.t1-cl.t0 > floor_ && sp+2 <= cap){
                P->unresolved = before;
                double tm=0.5*(cl.t0+cl.t1);
                stack[sp].t0=cl.t0; stack[sp].t1=tm; sp++;
                stack[sp].t0=tm; stack[sp].t1=cl.t1; sp++;
            }
        }
    }
    free(stack);
}

static int dump_left = -1;   /* CK5_DUMP=N: print first N parked leaves */
static void dfs(Prob *P, State *st, int oi, double cellw, int has_cellw){
    if (P->aborted || P->cellq_over) return;
    if (++P->nodes > P->max_nodes){ P->aborted = 1; return; }
    if (P->cellq_limit && P->nodes - P->cellq_start > P->cellq_limit){
        P->cellq_over = 1; return;     /* unwind to the quota owner */
    }
    if (oi == P->order_len){
        if (final_sweep_c(P, st)){
            P->unresolved++;
            if (dump_left < 0){
                const char *e = getenv("CK5_DUMP");
                dump_left = e ? atoi(e) : 0;
            }
            if (dump_left > 0){
                dump_left--;
                fprintf(stderr, "SURVIVOR leaf:\n");
                for (int pi = 0; pi < st->np; pi++){
                    int v = st->porder[pi];
                    fprintf(stderr, "  v%d:", v);
                    for (int d = 0; d < D; d++)
                        fprintf(stderr, " [%.9f,%.9f]",
                                st->p[v][d].a, st->p[v][d].b);
                    fprintf(stderr, "\n");
                }
            }
        }
        return;
    }
    int v = P->order[oi];
    int nbrids[MAXN], nnb=0;
    for (int pi=0; pi<st->np; pi++){
        int u = st->porder[pi];
        if ((P->adj[v]>>u)&1) nbrids[nnb++]=u;
    }
    if (nnb >= D){
        /* try 5-subsets until well-conditioned (capped) */
        Cand cd[2]; int r=-1;
        int idx[D]; int tried=0;
        for (idx[0]=0; idx[0]<nnb-4 && r<0; idx[0]++)
        for (idx[1]=idx[0]+1; idx[1]<nnb-3 && r<0; idx[1]++)
        for (idx[2]=idx[1]+1; idx[2]<nnb-2 && r<0; idx[2]++)
        for (idx[3]=idx[2]+1; idx[3]<nnb-1 && r<0; idx[3]++)
        for (idx[4]=idx[3]+1; idx[4]<nnb && r<0; idx[4]++){
            if (++tried > 64) goto scan_done;
            V5 nbv[D];
            for (int t=0;t<D;t++) memcpy(nbv[t], st->p[nbrids[idx[t]]], sizeof(V5));
            r = place_sphere_c((const V5*)nbv, cd);
            if (r < 0) { r = -1; continue; }
            break;
        }
        scan_done: ;
        if (r < 0){
            if (has_cellw && cellw <= 64*P->theta_min)
                circle_stage_c(P, st, oi, v, nbrids, nnb, cellw, 1);
            else
                P->unresolved++;
            return;
        }
        for (int s=0; s<r; s++){
            V5 x; memcpy(x, cd[s].x, sizeof(V5));
            V5 nbv[MAXN];
            for (int t=0;t<nnb;t++) memcpy(nbv[t], st->p[nbrids[t]], sizeof(V5));
            if (!krawczyk(x, (const V5*)nbv, nnb, 3)) continue;
            double mw=0; for (int i=0;i<D;i++) mw=fmax(mw, wid(x[i]));
            if (r == 1 && has_cellw && cellw <= 8*P->theta_min && mw > 1e-7){
                /* merged roots (disc ~ 0): tangency channel; resolve along the
                   root segment so near-coincident vs mirror solutions separate */
                segment_stage_c(P, st, oi, v, nbrids, nnb, &cd[s], cellw, has_cellw);
                continue;
            }
            if (mw > P->max_width){
                if (has_cellw && cellw <= 8*P->theta_min)
                    segment_stage_c(P, st, oi, v, nbrids, nnb, &cd[s], cellw, has_cellw);
                else
                    P->unresolved++;
                continue;
            }
            if (noninjective_c(P, st, v, x)) continue;
            int okall=1;
            for (int w=0; w<nnb; w++)
                if (!edge_ok(x, st->p[nbrids[w]])){ okall=0; break; }
            if (!okall) continue;
            State st2 = *st;
            memcpy(st2.p[v], x, sizeof(V5));
            st2.placedmask |= (uint32_t)1<<v;
            st2.porder[st2.np++] = v;
            if (has_cellw && !local_sweep_c(P, &st2, v)) continue;
            dfs(P, &st2, oi+1, cellw, has_cellw);
        }
    } else if (nnb == D-1){
        circle_stage_c(P, st, oi, v, nbrids, nnb,
                       has_cellw? cellw: 0.0, has_cellw);
    } else {
        fprintf(stderr, "vertex %d has %d placed nbrs\n", v, nnb);
        exit(8);
    }
}

/* entry point */
EXPORT int64_t decide5_c(int n, uint32_t *adj, int seed_len, int *seed,
              int order_len, int *order, double *seedc /*6*5*2*/,
              double theta_min, double th0_lo, double th0_hi,
              int64_t max_nodes, int64_t *nodes_out, int64_t *unres_out){
    Prob P;
    memset(&P, 0, sizeof P);
    P.n = n;
    for (int i=0;i<n;i++) P.adj[i]=adj[i];
    P.seed_len=seed_len;
    for (int i=0;i<seed_len;i++) P.seed[i]=seed[i];
    P.order_len=order_len;
    for (int i=0;i<order_len;i++) P.order[i]=order[i];
    P.theta_min=theta_min; P.max_width=1.0; P.max_nodes=max_nodes;
    P.th0_lo=th0_lo; P.th0_hi=th0_hi;
    State st; memset(&st, 0, sizeof st);
    for (int i=0;i<seed_len;i++){
        int v = seed[i];
        for (int j=0;j<D;j++){
            st.p[v][j] = iv(seedc[(i*D+j)*2], seedc[(i*D+j)*2+1]);
        }
        st.placedmask |= (uint32_t)1<<v;
        st.porder[st.np++] = v;
    }
    dfs(&P, &st, 0, 0.0, 0);
    *nodes_out = P.nodes; *unres_out = P.unresolved;
    if (P.aborted) return 2;
    return P.unresolved == 0 ? 0 : 1;
}
