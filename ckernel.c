/* Certified branch-and-prune kernel deciding unit-distance realizability
 * in R^4 (algorithm documented in Appendix B of f4_equals_12.pdf).
 * Interval arithmetic with outward rounding via nextafter (IEEE basic ops
 * are correctly rounded); cos/sin get 8-ulp padding.
 * Compiled automatically by cdriver.py:
 *   POSIX:   cc -O2 -shared -o ckernel.so ckernel.c -lm
 *   Windows: cl /O2 /LD ckernel.c /Fe:ckernel.dll   (from a VS x64 prompt)
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

#define D 4
#define MAXN 13
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

typedef IV V4[D];
static void vsubv(const V4 u, const V4 v, V4 r){ for(int i=0;i<D;i++) r[i]=isub(u[i],v[i]); }
static void vaddv(const V4 u, const V4 v, V4 r){ for(int i=0;i<D;i++) r[i]=iadd(u[i],v[i]); }
static void vscalev(IV t, const V4 u, V4 r){ for(int i=0;i<D;i++) r[i]=imul(t,u[i]); }
static IV vdotv(const V4 u, const V4 v){
    IV s=ivp(0); for(int i=0;i<D;i++) s=iadd(s, imul(u[i],v[i])); return s;
}
static IV vnorm2v(const V4 u){
    IV s=ivp(0); for(int i=0;i<D;i++) s=iadd(s, isq(u[i])); return s;
}
static int overlap4(const V4 x, const V4 y){
    for(int i=0;i<D;i++) if (x[i].a > y[i].b || y[i].a > x[i].b) return 0;
    return 1;
}

/* ---- generic interval Gauss-Jordan: m rows, nuns unknowns (<=4) ---- */
static int gauss(IV A[4][4], IV b[4], int m, int nuns, IV xp[4], IV nullv[4], int *nnull){
    int usedr[4]={0,0,0,0}, usedc[4]={0,0,0,0};
    int piv_r[4], piv_c[4];
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
    int freec[4], nf=0;
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
    unsigned short adj[MAXN];
    int seed_len, seed[5];
    int order_len, order[MAXN];
    double seedc[5][D][2];
    double theta_min, max_width;
    int64_t max_nodes;
    double th0_lo, th0_hi;
    /* runtime */
    int64_t nodes, unresolved;
    int aborted, th0_used;
} Prob;

typedef struct { V4 p[MAXN]; unsigned short placedmask; int porder[MAXN]; int np; } State;

static void dfs(Prob *P, State *st, int oi, double cellw, int has_cellw);

/* place via 4 spheres: cands out; returns -1 illcond, else count (0..2);
   seg: xp,n,t for each cand */
typedef struct { V4 x; V4 xp, nv; IV t; } Cand;
static int place_sphere_c(const V4 nb[4], Cand out[2]){
    IV A[4][4], b[4]; V4 d;
    for (int r=1;r<4;r++){
        for (int c=0;c<D;c++) A[r-1][c] = imul(ivp(2.0), isub(nb[r][c], nb[0][c]));
        b[r-1] = isub(vnorm2v(nb[r]), vnorm2v(nb[0]));
    }
    IV xp[4], nv[4]; int nn;
    if (!gauss(A, b, 3, 4, xp, nv, &nn)) return -1;
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
        V4 tn; vscalev(t, nv, tn); vaddv(xp, tn, out[0].x);
        memcpy(out[0].xp, xp, sizeof(V4)); memcpy(out[0].nv, nv, sizeof(V4)); out[0].t = t;
        return 1;
    }
    for (int s=0;s<2;s++){
        int ok2; IV t = idivi(iadd(ineg(bq), s? ineg(sq): sq), imul(ivp(2.0), a), &ok2);
        if (!ok2) return -1;
        V4 tn; vscalev(t, nv, tn); vaddv(xp, tn, out[s].x);
        memcpy(out[s].xp, xp, sizeof(V4)); memcpy(out[s].nv, nv, sizeof(V4)); out[s].t = t;
    }
    return 2;
}

static int edge_ok(const V4 x, const V4 p){
    V4 d; vsubv(x, p, d);
    IV e = isub(vnorm2v(d), ivp(1.0));
    return contains0(e);
}

/* float 4x4 inverse; returns 0 on singular */
static int inv4(double M[4][4], double Y[4][4]){
    double a[4][8];
    for (int i=0;i<4;i++){ for (int j=0;j<4;j++){ a[i][j]=M[i][j]; a[i][j+4]=(i==j);} }
    for (int c=0;c<4;c++){
        int p=-1; double best=0;
        for (int r=c;r<4;r++){ double v=fabs(a[r][c]); if (v>best){best=v;p=r;} }
        if (p<0 || best<1e-14) return 0;
        if (p!=c){ for (int j=0;j<8;j++){ double t=a[c][j]; a[c][j]=a[p][j]; a[p][j]=t; } }
        double d=a[c][c];
        for (int j=0;j<8;j++) a[c][j]/=d;
        for (int r=0;r<4;r++){ if (r==c) continue; double f=a[r][c];
            for (int j=0;j<8;j++) a[r][j]-=f*a[c][j]; }
    }
    for (int i=0;i<4;i++) for (int j=0;j<4;j++) Y[i][j]=a[i][j+4];
    return 1;
}

/* choose best-conditioned quad of neighbour indices (float det) */
static void best_quad(const V4 x, const V4 nb[], int k, int q[4]){
    double xm[4]; for (int i=0;i<D;i++) xm[i]=mid(x[i]);
    if (k==4){ q[0]=0;q[1]=1;q[2]=2;q[3]=3; return; }
    double best=-1; int bq[4]={0,1,2,3};
    int idx[4];
    for (idx[0]=0; idx[0]<k-3; idx[0]++)
    for (idx[1]=idx[0]+1; idx[1]<k-2; idx[1]++)
    for (idx[2]=idx[1]+1; idx[2]<k-1; idx[2]++)
    for (idx[3]=idx[2]+1; idx[3]<k; idx[3]++){
        double J[4][4];
        for (int r=0;r<4;r++) for (int c=0;c<D;c++)
            J[r][c] = 2*(xm[c]-mid(nb[idx[r]][c]));
        /* det via LU-lite */
        double m[4][4]; memcpy(m,J,sizeof m); double det=1;
        int sing=0;
        for (int cc=0;cc<4;cc++){
            int p=-1; double bv=0;
            for (int r=cc;r<4;r++){ double v=fabs(m[r][cc]); if (v>bv){bv=v;p=r;} }
            if (p<0||bv<1e-300){ sing=1; break; }
            if (p!=cc){ for(int j=0;j<4;j++){double t=m[cc][j];m[cc][j]=m[p][j];m[p][j]=t;} det=-det; }
            det*=m[cc][cc];
            for (int r=cc+1;r<4;r++){ double f=m[r][cc]/m[cc][cc];
                for (int j=cc;j<4;j++) m[r][j]-=f*m[cc][j]; }
        }
        double ad = sing?0:fabs(det);
        if (ad > best){ best=ad; memcpy(bq, idx, sizeof bq); }
        if (best > 0.75) goto done;
    }
done:
    memcpy(q, bq, 4*sizeof(int));
}

/* Krawczyk contraction on 4 chosen spheres; 1=ok (x updated), 0=empty */
static int krawczyk(V4 x, const V4 nb[], int k, int iters){
    int q[4]; best_quad(x, nb, k, q);
    const IV *ps[4]; for (int i=0;i<4;i++) ps[i]=nb[q[i]];
    for (int it=0; it<iters; it++){
        double mm[4]; for (int i=0;i<D;i++) mm[i]=mid(x[i]);
        double Jm[4][4], Y[4][4];
        for (int r=0;r<4;r++) for (int c=0;c<D;c++) Jm[r][c]=2*(mm[c]-mid(ps[r][c]));
        if (!inv4(Jm, Y)) return 1;
        IV m[4]; for (int i=0;i<D;i++) m[i]=ivp(mm[i]);
        IV Fm[4];
        for (int r=0;r<4;r++){
            V4 dd; for (int c=0;c<D;c++) dd[c]=isub(m[c], ps[r][c]);
            Fm[r] = isub(vnorm2v(dd), ivp(1.0));
        }
        IV K[4]; int shrunk=0;
        IV JX[4][4];
        for (int r=0;r<4;r++) for (int c=0;c<D;c++)
            JX[r][c] = imul(ivp(2.0), isub(x[c], ps[r][c]));
        IV dxm[4]; for (int i=0;i<D;i++) dxm[i]=isub(x[i], m[i]);
        for (int i=0;i<4;i++){
            IV s = m[i];
            for (int j=0;j<4;j++) s = isub(s, imul(ivp(Y[i][j]), Fm[j]));
            for (int j=0;j<4;j++){
                IV acc = (i==j)? ivp(1.0): ivp(0.0);
                for (int c=0;c<4;c++) acc = isub(acc, imul(ivp(Y[i][c]), JX[c][j]));
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

/* interval det of 4 difference vectors */
static IV det4I(V4 rows[4]){
    /* cofactor expansion along first row */
    IV tot = ivp(0);
    for (int j=0;j<4;j++){
        int cols[3], nc=0;
        for (int k2=0;k2<4;k2++) if (k2!=j) cols[nc++]=k2;
        IV m[3][3];
        for (int r=1;r<4;r++) for (int c=0;c<3;c++) m[r-1][c]=rows[r][cols[c]];
        IV d3 = isub(imul(m[0][0], isub(imul(m[1][1],m[2][2]), imul(m[1][2],m[2][1]))),
               isub(imul(m[0][1], isub(imul(m[1][0],m[2][2]), imul(m[1][2],m[2][0]))),
                    imul(m[0][2], isub(imul(m[1][0],m[2][1]), imul(m[1][1],m[2][0])))));
        /* careful: d3 = a*(ei-fh) - b*(di-fg) + c*(dh-eg); rewrite properly */
        IV a_=m[0][0], b_=m[0][1], c_=m[0][2];
        IV t1 = imul(a_, isub(imul(m[1][1],m[2][2]), imul(m[1][2],m[2][1])));
        IV t2 = imul(b_, isub(imul(m[1][0],m[2][2]), imul(m[1][2],m[2][0])));
        IV t3 = imul(c_, isub(imul(m[1][0],m[2][1]), imul(m[1][1],m[2][0])));
        d3 = iadd(isub(t1, t2), t3);
        IV term = imul(rows[0][j], d3);
        if (j % 2) term = ineg(term);
        tot = iadd(tot, term);
    }
    return tot;
}

static int forced_coincident_c(Prob *P, State *st, int u, int v){
    int common[MAXN], nc=0;
    for (int wpi=0; wpi<st->np; wpi++){
        int w = st->porder[wpi];
        if (w==u || w==v) continue;
        if (((P->adj[v]>>w)&1) && ((P->adj[u]>>w)&1)) common[nc++]=w;
    }
    if (nc < 5) return 0;
    int idx[5];
    for (idx[0]=0; idx[0]<nc-4; idx[0]++)
    for (idx[1]=idx[0]+1; idx[1]<nc-3; idx[1]++)
    for (idx[2]=idx[1]+1; idx[2]<nc-2; idx[2]++)
    for (idx[3]=idx[2]+1; idx[3]<nc-1; idx[3]++)
    for (idx[4]=idx[3]+1; idx[4]<nc; idx[4]++){
        V4 rows[4];
        for (int r=0;r<4;r++)
            vsubv(st->p[common[idx[r+1]]], st->p[common[idx[0]]], rows[r]);
        IV dt = det4I(rows);
        if (dt.a > 0 || dt.b < 0) return 1;
    }
    return 0;
}

/* noninjective placement filter (mirror of Python) */
static int noninjective_c(Prob *P, State *st, int v, const V4 x){
    for (int upi=0; upi<st->np; upi++){
        int u = st->porder[upi];
        if ((P->adj[v]>>u)&1) continue;
        if (!overlap4(x, st->p[u])) continue;
        /* hyperplane rule with v's box in place */
        State tmp = *st;
        memcpy(tmp.p[v], x, sizeof(V4));
        tmp.placedmask |= (unsigned short)(1<<v);
        tmp.porder[tmp.np++] = v;
        if (forced_coincident_c(P, &tmp, u, v)) return 1;
        /* root-separation rule */
        int common[MAXN], nc=0;
        for (int wpi=0; wpi<st->np; wpi++){
            int w = st->porder[wpi];
            if (w==u) continue;
            if (((P->adj[v]>>w)&1) && ((P->adj[u]>>w)&1)) common[nc++]=w;
        }
        if (nc < 4) continue;
        int idx[4];
        for (idx[0]=0; idx[0]<nc-3; idx[0]++)
        for (idx[1]=idx[0]+1; idx[1]<nc-2; idx[1]++)
        for (idx[2]=idx[1]+1; idx[2]<nc-1; idx[2]++)
        for (idx[3]=idx[2]+1; idx[3]<nc; idx[3]++){
            const V4 nb[4] = {{{0}}};
            V4 nbv[4];
            for (int r=0;r<4;r++) memcpy(nbv[r], st->p[common[idx[r]]], sizeof(V4));
            Cand cd[2];
            int r = place_sphere_c((const V4*)nbv, cd);
            (void)nb;
            if (r < 0) continue;
            if (r == 0) return 1;               /* v must solve: empty */
            if (r != 2) continue;
            if (overlap4(cd[0].x, cd[1].x)) continue;
            int d1 = overlap4(st->p[u], cd[0].x), d2 = overlap4(st->p[u], cd[1].x);
            int c1 = overlap4(x, cd[0].x), c2 = overlap4(x, cd[1].x);
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
    unsigned short seedm=0;
    for (int i=0;i<P->seed_len;i++) seedm |= (unsigned short)(1<<P->seed[i]);
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
                V4 nbv[MAXN]; int k=0;
                int nbrids[MAXN];
                for (int wpi=0; wpi<st->np; wpi++){
                    int w = st->porder[wpi];
                    if ((P->adj[u]>>w)&1){
                        memcpy(nbv[k], st->p[w], sizeof(V4)); nbrids[k]=w; k++; } }
                if (k < 4) continue;
                V4 xx; memcpy(xx, st->p[u], sizeof(V4));
                if (!krawczyk(xx, (const V4*)nbv, k, 1)) return 0;
                int chg=0;
                for (int i=0;i<D;i++)
                    if (wid(xx[i]) < wid(st->p[u][i]) - 1e-14) chg=1;
                memcpy(st->p[u], xx, sizeof(V4));
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
            if (!overlap4(st->p[u], st->p[v])) continue;
            if (forced_coincident_c(P, st, u, v)) return 0;
            int common[MAXN], nc=0;
            for (int wpi=0; wpi<st->np; wpi++){
                int w = st->porder[wpi];
                if (w==u||w==v) continue;
                if (((P->adj[v]>>w)&1) && ((P->adj[u]>>w)&1)) common[nc++]=w;
            }
            if (nc < 4) continue;
            int idx[4];
            for (idx[0]=0; idx[0]<nc-3; idx[0]++)
            for (idx[1]=idx[0]+1; idx[1]<nc-2; idx[1]++)
            for (idx[2]=idx[1]+1; idx[2]<nc-1; idx[2]++)
            for (idx[3]=idx[2]+1; idx[3]<nc; idx[3]++){
                V4 nbv[4];
                for (int r=0;r<4;r++) memcpy(nbv[r], st->p[common[idx[r]]], sizeof(V4));
                Cand cd[2];
                int r = place_sphere_c((const V4*)nbv, cd);
                if (r < 0) continue;
                if (r == 0) return 0;
                if (r != 2) continue;
                if (overlap4(cd[0].x, cd[1].x)) continue;
                int du1=overlap4(st->p[u], cd[0].x), du2=overlap4(st->p[u], cd[1].x);
                int dv1=overlap4(st->p[v], cd[0].x), dv2=overlap4(st->p[v], cd[1].x);
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
    unsigned short seedm=0;
    for (int i=0;i<P->seed_len;i++) seedm |= (unsigned short)(1<<P->seed[i]);
    for (int sweep=0; sweep<6; sweep++){
        int changed=0;
        for (int vpi=0; vpi<s.np; vpi++){
            int v = s.porder[vpi];
            if ((seedm>>v)&1) continue;
            double mw=0; for (int i=0;i<D;i++) mw=fmax(mw, wid(s.p[v][i]));
            if (mw < 1e-9) continue;
            V4 nbv[MAXN]; int k=0; int ids[MAXN];
            for (int wpi=0; wpi<s.np; wpi++){
                int w = s.porder[wpi];
                if ((P->adj[v]>>w)&1){
                    memcpy(nbv[k], s.p[w], sizeof(V4)); ids[k]=w; k++; } }
            if (k < 4) continue;
            V4 xx; memcpy(xx, s.p[v], sizeof(V4));
            if (!krawczyk(xx, (const V4*)nbv, k, 2)) return 0;
            for (int i=0;i<D;i++)
                if (wid(xx[i]) < wid(s.p[v][i]) - 1e-15) changed=1;
            memcpy(s.p[v], xx, sizeof(V4));
            for (int w2=0; w2<k; w2++)
                if (!edge_ok(s.p[v], s.p[ids[w2]])) return 0;
        }
        if (!changed) break;
    }
    return injective_or_unresolved_c(P, &s);
}

/* ---- circle & segment stages ---- */
typedef struct { double t0, t1; } Cell;

static int circle_setup_c(const V4 p0, const V4 p1, const V4 p2,
                          V4 c, IV *r, V4 u1, V4 u2){
    V4 d1, d2;
    vsubv(p1, p0, d1); vsubv(p2, p0, d2);
    IV r11=vdotv(d1,d1), r12=vdotv(d1,d2), r22=vdotv(d2,d2);
    IV A[4][4], b[4];
    A[0][0]=r11; A[0][1]=r12; A[1][0]=r12; A[1][1]=r22;
    int okd;
    b[0]=idivi(r11, ivp(2.0), &okd); b[1]=idivi(r22, ivp(2.0), &okd);
    IV ab[4], nl[4]; int nn;
    if (!gauss(A, b, 2, 2, ab, nl, &nn)) return -1;
    V4 t1v, t2v;
    vscalev(ab[0], d1, t1v); vscalev(ab[1], d2, t2v);
    V4 s; vaddv(t1v, t2v, s); vaddv(p0, s, c);
    V4 cd; vsubv(c, p0, cd);
    IV R2 = vnorm2v(cd);
    IV r2 = isub(ivp(1.0), R2);
    if (r2.b < 0) return 0;                     /* certified empty */
    int ok; *r = isqrtI(r2, &ok);
    if (!ok) return 0;
    /* orthonormalize hull dirs, then complement */
    V4 ob[2]; int no=0;
    V4 cand[2]; memcpy(cand[0], d1, sizeof(V4)); memcpy(cand[1], d2, sizeof(V4));
    for (int q2=0; q2<2; q2++){
        V4 e; memcpy(e, cand[q2], sizeof(V4));
        for (int b2=0;b2<no;b2++){
            IV pr = vdotv(e, ob[b2]); V4 t; vscalev(pr, ob[b2], t); vsubv(e, t, e);
        }
        IV nrm = vnorm2v(e);
        if (!(nrm.a > 0) || mig(nrm) <= 1e-8) return -1;
        int ok2; IV s2 = isqrtI(nrm, &ok2);
        for (int i2=0;i2<D;i2++){ int ok3; ob[no][i2]=idivi(e[i2], s2, &ok3); if(!ok3) return -1; }
        no++;
    }
    int got=0;
    for (int k=0;k<D && got<2;k++){
        V4 e; for (int i2=0;i2<D;i2++) e[i2] = (i2==k)? ivp(1.0): ivp(0.0);
        for (int b2=0;b2<no;b2++){
            IV pr=vdotv(e, ob[b2]); V4 t; vscalev(pr, ob[b2], t); vsubv(e, t, e);
        }
        if (got==1){
            IV pr=vdotv(e, u1); V4 t; vscalev(pr, u1, t); vsubv(e, t, e);
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
    /* find working triple */
    V4 c, u1, u2; IV r;
    int rtri[3] = {-1,-1,-1};
    int found = 0;
    int idx[3];
    for (idx[0]=0; idx[0]<nnb-2 && !found; idx[0]++)
    for (idx[1]=idx[0]+1; idx[1]<nnb-1 && !found; idx[1]++)
    for (idx[2]=idx[1]+1; idx[2]<nnb && !found; idx[2]++){
        int rc = circle_setup_c(st->p[nbrids[idx[0]]], st->p[nbrids[idx[1]]],
                                st->p[nbrids[idx[2]]], c, &r, u1, u2);
        if (rc == 0) return;                     /* certified empty: kill */
        if (rc == 1){ rtri[0]=nbrids[idx[0]]; rtri[1]=nbrids[idx[1]];
                      rtri[2]=nbrids[idx[2]]; found=1; }
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
        if (P->aborted) break;
        Cell cl = stack[--sp];
        IV th = iv(cl.t0, cl.t1);
        V4 x, tc, ts, sum;
        IV rc_ = imul(r, icosI(th)), rs_ = imul(r, isinI(th));
        vscalev(rc_, u1, tc); vscalev(rs_, u2, ts);
        vaddv(tc, ts, sum); vaddv(c, sum, x);
        int dead = 0;
        for (int w=0; w<nnb; w++){
            if (!edge_ok(x, st->p[nbrids[w]])){
                int isdef = (nbrids[w]==rtri[0]||nbrids[w]==rtri[1]||nbrids[w]==rtri[2]);
                if (isdef){ fprintf(stderr, "BUG circle defining violated\n"); exit(9); }
                dead = 1; break;
            }
        }
        if (dead) continue;
        int64_t before = P->unresolved;
        State st2 = *st;
        memcpy(st2.p[v], x, sizeof(V4));
        st2.placedmask |= (unsigned short)(1<<v);
        st2.porder[st2.np++] = v;
        if (!local_sweep_c(P, &st2, v)) continue;
        dfs(P, &st2, oi+1, cl.t1-cl.t0, 1);
        if (P->unresolved > before){
            if (cl.t1-cl.t0 > floor_ && sp+2 <= cap){
                P->unresolved = before;
                double tm = 0.5*(cl.t0+cl.t1);
                stack[sp].t0=cl.t0; stack[sp].t1=tm; sp++;
                stack[sp].t0=tm; stack[sp].t1=cl.t1; sp++;
            }
        }
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
        if (P->aborted) break;
        Cell cl = stack[--sp];
        IV tc = iv(cl.t0, cl.t1);
        V4 x, tn; vscalev(tc, seg->nv, tn); vaddv(seg->xp, tn, x);
        int dead=0;
        for (int w=0; w<nnb; w++)
            if (!edge_ok(x, st->p[nbrids[w]])){ dead=1; break; }
        if (dead) continue;
        if (noninjective_c(P, st, v, x)) continue;
        int64_t before = P->unresolved;
        State st2 = *st;
        memcpy(st2.p[v], x, sizeof(V4));
        st2.placedmask |= (unsigned short)(1<<v);
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

static void dfs(Prob *P, State *st, int oi, double cellw, int has_cellw){
    if (P->aborted) return;
    if (++P->nodes > P->max_nodes){ P->aborted = 1; return; }
    if (oi == P->order_len){
        if (final_sweep_c(P, st)) P->unresolved++;
        return;
    }
    int v = P->order[oi];
    int nbrids[MAXN], nnb=0;
    for (int pi=0; pi<st->np; pi++){
        int u = st->porder[pi];
        if ((P->adj[v]>>u)&1) nbrids[nnb++]=u;
    }
    if (nnb >= 4){
        /* try 4-subsets until well-conditioned */
        Cand cd[2]; int r=-1;
        int idx[4];
        for (idx[0]=0; idx[0]<nnb-3 && r<0; idx[0]++)
        for (idx[1]=idx[0]+1; idx[1]<nnb-2 && r<0; idx[1]++)
        for (idx[2]=idx[1]+1; idx[2]<nnb-1 && r<0; idx[2]++)
        for (idx[3]=idx[2]+1; idx[3]<nnb && r<0; idx[3]++){
            V4 nbv[4];
            for (int t=0;t<4;t++) memcpy(nbv[t], st->p[nbrids[idx[t]]], sizeof(V4));
            r = place_sphere_c((const V4*)nbv, cd);
            if (r < 0) { r = -1; continue; }
            break;
        }
        if (r < 0){
            if (has_cellw && cellw <= 64*P->theta_min)
                circle_stage_c(P, st, oi, v, nbrids, nnb, cellw, 1);
            else
                P->unresolved++;
            return;
        }
        for (int s=0; s<r; s++){
            V4 x; memcpy(x, cd[s].x, sizeof(V4));
            V4 nbv[MAXN];
            for (int t=0;t<nnb;t++) memcpy(nbv[t], st->p[nbrids[t]], sizeof(V4));
            if (!krawczyk(x, (const V4*)nbv, nnb, 3)) continue;
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
            memcpy(st2.p[v], x, sizeof(V4));
            st2.placedmask |= (unsigned short)(1<<v);
            st2.porder[st2.np++] = v;
            if (has_cellw && !local_sweep_c(P, &st2, v)) continue;
            dfs(P, &st2, oi+1, cellw, has_cellw);
        }
    } else if (nnb == 3){
        circle_stage_c(P, st, oi, v, nbrids, nnb,
                       has_cellw? cellw: 0.0, has_cellw);
    } else {
        fprintf(stderr, "vertex %d has %d placed nbrs\n", v, nnb);
        exit(8);
    }
}

/* entry point */
EXPORT int64_t decide_c(int n, unsigned short *adj, int seed_len, int *seed,
              int order_len, int *order, double *seedc /*5*4*2*/,
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
        st.placedmask |= (unsigned short)(1<<v);
        st.porder[st.np++] = v;
    }
    dfs(&P, &st, 0, 0.0, 0);
    *nodes_out = P.nodes; *unres_out = P.unresolved;
    if (P.aborted) return 2;
    return P.unresolved == 0 ? 0 : 1;
}
