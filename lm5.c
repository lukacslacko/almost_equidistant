/* Numerical realizability sweep: Levenberg-Marquardt over unit-distance
 * constraints in R^5. For each graph (format: "n mask0 ... mask{n-1}" per
 * line), minimize  F(p) = sum_{ij in E} (|pi-pj|^2 - 1)^2,  p in R^{5n},
 * from many random starts. Reports per graph: best residual F, the minimum
 * pairwise distance at the best solution, and #near-solutions.
 *
 * A best F < 1e-20 with min pairwise distance > 1e-3 is a REALIZABLE
 * candidate (would imply f(5) >= n) and is flagged loudly. This program is
 * heuristic only - the certified engine decides; but it corroborates kills
 * (best residuals staying far from 0) and can detect realizability early.
 *
 * Usage: lm5 <graphfile> <start> <end> <restarts> <seedbase> [> results]
 * One output line per graph:
 *   idx n bestF mindist nsol edges
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#define D 5
#define MAXN 22
#define MAXE 260

typedef struct { int n, m; int ei[MAXE], ej[MAXE]; } Graph;

/* xorshift RNG (deterministic per seed) */
static uint64_t rng_s;
static inline uint64_t xr(void){
    rng_s ^= rng_s << 13; rng_s ^= rng_s >> 7; rng_s ^= rng_s << 17; return rng_s;
}
static inline double urand(void){ return (double)(xr() >> 11) / 9007199254740992.0; }
static double gauss_rand(void){
    double u1 = urand(), u2 = urand();
    if (u1 < 1e-300) u1 = 1e-300;
    return sqrt(-2.0 * log(u1)) * cos(6.283185307179586 * u2);
}

static double residual(const Graph *g, const double *p, double *r){
    double F = 0;
    for (int e = 0; e < g->m; e++){
        const double *a = p + g->ei[e]*D, *b = p + g->ej[e]*D;
        double s = 0;
        for (int d = 0; d < D; d++){ double t = a[d]-b[d]; s += t*t; }
        double v = s - 1.0;
        if (r) r[e] = v;
        F += v*v;
    }
    return F;
}

/* solve (A + lam*diag(A)) x = b in-place via Cholesky; A is NxN sym (packed full) */
static int chol_solve(int N, double *A, double lam, const double *b, double *x){
    static double L[110*110];
    for (int i = 0; i < N; i++)
        for (int j = 0; j <= i; j++){
            double v = A[i*N+j];
            if (i == j) v += lam * (A[i*N+i] > 1e-12 ? A[i*N+i] : 1e-12);
            L[i*N+j] = v;
        }
    for (int i = 0; i < N; i++){
        for (int j = 0; j <= i; j++){
            double s = L[i*N+j];
            for (int k = 0; k < j; k++) s -= L[i*N+k]*L[j*N+k];
            if (i == j){
                if (s <= 0) return 0;
                L[i*N+i] = sqrt(s);
            } else L[i*N+j] = s / L[j*N+j];
        }
    }
    static double y[110];
    for (int i = 0; i < N; i++){
        double s = b[i];
        for (int k = 0; k < i; k++) s -= L[i*N+k]*y[k];
        y[i] = s / L[i*N+i];
    }
    for (int i = N-1; i >= 0; i--){
        double s = y[i];
        for (int k = i+1; k < N; k++) s -= L[k*N+i]*x[k];
        x[i] = s / L[i*N+i];
    }
    return 1;
}

/* one LM run from p; returns final F */
static double lm_run(const Graph *g, double *p, int iters){
    int N = g->n * D;
    static double r[MAXE], JtJ[110*110], Jtr[110], delta[110], pn[110], rn[MAXE];
    double F = residual(g, p, r);
    double lam = 1e-3;
    for (int it = 0; it < iters; it++){
        /* build JtJ, Jtr: row e has 2(a-b) at i-block, -2(a-b) at j-block */
        memset(JtJ, 0, sizeof(double)*N*N);
        memset(Jtr, 0, sizeof(double)*N);
        for (int e = 0; e < g->m; e++){
            int i = g->ei[e], j = g->ej[e];
            double dv[D];
            for (int d = 0; d < D; d++) dv[d] = 2.0*(p[i*D+d]-p[j*D+d]);
            /* Jtr */
            for (int d = 0; d < D; d++){
                Jtr[i*D+d] += dv[d]*r[e];
                Jtr[j*D+d] -= dv[d]*r[e];
            }
            /* JtJ blocks: ii += dv dv^T, jj += dv dv^T, ij -= dv dv^T */
            for (int a = 0; a < D; a++)
                for (int b = 0; b < D; b++){
                    double v = dv[a]*dv[b];
                    JtJ[(i*D+a)*N + i*D+b] += v;
                    JtJ[(j*D+a)*N + j*D+b] += v;
                    JtJ[(i*D+a)*N + j*D+b] -= v;
                    JtJ[(j*D+a)*N + i*D+b] -= v;
                }
        }
        int accepted = 0;
        for (int tries = 0; tries < 8; tries++){
            if (chol_solve(N, JtJ, lam, Jtr, delta)){
                for (int k = 0; k < N; k++) pn[k] = p[k] - delta[k];
                double Fn = residual(g, pn, rn);
                if (Fn < F){
                    memcpy(p, pn, sizeof(double)*N);
                    memcpy(r, rn, sizeof(double)*g->m);
                    F = Fn;
                    lam = lam > 1e-12 ? lam * 0.3 : 1e-12;
                    accepted = 1;
                    break;
                }
            }
            lam *= 10.0;
            if (lam > 1e10) break;
        }
        if (!accepted) break;
        if (F < 1e-28) break;
    }
    return F;
}

/* greedy placement initializer: exact unit-simplex on a (max<=6)-clique,
   then vertices in max-placed-neighbour order, each solved by a small
   Gauss-Newton on its own 5 coordinates from a random kick. */
static const double SIMP[6][5] = {
    {0,0,0,0,0},
    {1,0,0,0,0},
    {0.5, 0.8660254037844386, 0, 0, 0},
    {0.5, 0.2886751345948129, 0.816496580927726, 0, 0},
    {0.5, 0.2886751345948129, 0.2041241452319315, 0.7905694150420949, 0},
    {0.5, 0.2886751345948129, 0.2041241452319315, 0.15811388300841897, 0.7745966692414834},
};
static void place_init(const uint32_t *adj, int n, double *p){
    /* greedy clique up to 6 */
    int clq[6], nc = 0;
    uint32_t cand = (1u << n) - 1;
    int best = 0, bd = -1;
    for (int v = 0; v < n; v++){ int d = 0; for (int u = 0; u < n; u++) d += (adj[v]>>u)&1;
        if (d > bd){ bd = d; best = v; } }
    clq[nc++] = best; cand &= adj[best];
    while (nc < 6 && cand){
        int bv = -1, bcnt = -1;
        uint32_t t = cand;
        while (t){
            int v = t & 1 ? 0 : 0; /* placeholder */
            v = __builtin_ctz(t); t &= t-1;
            int c2 = 0; uint32_t w = cand & adj[v];
            while (w){ w &= w-1; c2++; }
            if (c2 > bcnt){ bcnt = c2; bv = v; }
        }
        clq[nc++] = bv; cand &= adj[bv];
    }
    int placed[MAXN], np = 0;
    int isplaced[MAXN]; memset(isplaced, 0, sizeof isplaced);
    for (int i = 0; i < nc; i++){
        for (int d = 0; d < D; d++) p[clq[i]*D+d] = SIMP[i][d];
        isplaced[clq[i]] = 1; placed[np++] = clq[i];
    }
    while (np < n){
        int bv = -1, bcnt = -1;
        for (int v = 0; v < n; v++){
            if (isplaced[v]) continue;
            int c2 = 0;
            for (int u = 0; u < np; u++) c2 += (adj[v] >> placed[u]) & 1;
            if (c2 > bcnt){ bcnt = c2; bv = v; }
        }
        /* init near centroid of placed neighbours + random kick */
        double c0[D] = {0,0,0,0,0}; int k = 0;
        int nbr[MAXN]; int nk = 0;
        for (int u = 0; u < np; u++)
            if ((adj[bv] >> placed[u]) & 1){
                nbr[nk++] = placed[u];
                for (int d = 0; d < D; d++) c0[d] += p[placed[u]*D+d];
                k++;
            }
        double x[D];
        for (int d = 0; d < D; d++)
            x[d] = (k ? c0[d]/k : 0) + 0.6 * gauss_rand();
        /* Gauss-Newton on x against up to all placed neighbours */
        for (int it = 0; it < 25; it++){
            double JtJ[D*D], Jtr[D];
            memset(JtJ, 0, sizeof JtJ); memset(Jtr, 0, sizeof Jtr);
            for (int q = 0; q < nk; q++){
                double dv[D], s = 0;
                for (int d = 0; d < D; d++){ dv[d] = x[d]-p[nbr[q]*D+d]; s += dv[d]*dv[d]; }
                double r0 = s - 1.0;
                for (int a = 0; a < D; a++){
                    Jtr[a] += 2*dv[a]*r0;
                    for (int b = 0; b < D; b++) JtJ[a*D+b] += 4*dv[a]*dv[b];
                }
            }
            for (int a = 0; a < D; a++) JtJ[a*D+a] += 1e-6;
            /* tiny 5x5 solve (gauss elim) */
            double M[D][D+1];
            for (int a = 0; a < D; a++){ for (int b = 0; b < D; b++) M[a][b] = JtJ[a*D+b]; M[a][D] = -Jtr[a]; }
            for (int c2 = 0; c2 < D; c2++){
                int pv = c2;
                for (int r2 = c2+1; r2 < D; r2++) if (fabs(M[r2][c2]) > fabs(M[pv][c2])) pv = r2;
                if (fabs(M[pv][c2]) < 1e-14) break;
                if (pv != c2) for (int b = 0; b <= D; b++){ double t = M[c2][b]; M[c2][b] = M[pv][b]; M[pv][b] = t; }
                for (int r2 = 0; r2 < D; r2++){
                    if (r2 == c2) continue;
                    double f2 = M[r2][c2]/M[c2][c2];
                    for (int b = c2; b <= D; b++) M[r2][b] -= f2*M[c2][b];
                }
            }
            double step = 0;
            for (int a = 0; a < D; a++){
                double da = M[a][D]/M[a][a];
                x[a] += da; step += da*da;
            }
            if (step < 1e-24) break;
        }
        for (int d = 0; d < D; d++) p[bv*D+d] = x[d];
        isplaced[bv] = 1; placed[np++] = bv;
    }
}

int main(int argc, char **argv){
    if (argc < 5){
        fprintf(stderr, "usage: lm5 <graphfile> <start> <end> <restarts> [seedbase]\n");
        return 1;
    }
    FILE *f = fopen(argv[1], "r");
    if (!f){ perror("open"); return 1; }
    long start = atol(argv[2]), end = atol(argv[3]);
    int restarts = atoi(argv[4]);
    uint64_t seedbase = argc > 5 ? strtoull(argv[5], 0, 10) : 12345;
    char line[4096];
    long idx = -1;
    while (fgets(line, sizeof line, f)){
        idx++;
        if (idx < start) continue;
        if (idx >= end) break;
        Graph g; g.m = 0;
        char *tok = strtok(line, " \t\n");
        g.n = atoi(tok);
        uint32_t adj[MAXN];
        for (int i = 0; i < g.n; i++){
            tok = strtok(NULL, " \t\n");
            adj[i] = (uint32_t)strtoul(tok, 0, 10);
        }
        for (int i = 0; i < g.n; i++)
            for (int j = i+1; j < g.n; j++)
                if ((adj[i] >> j) & 1){ g.ei[g.m]=i; g.ej[g.m]=j; g.m++; }
        double bestF = 1e30;              /* over all runs (degenerate ok) */
        double bestFd = 1e30, bestFd_mind = 0;  /* over distinct solutions */
        int nsol = 0, nsold = 0;
        static double p[110];
        const double sigmas[3] = {0.35, 0.5, 0.75};
        for (int rs = 0; rs < restarts; rs++){
            rng_s = seedbase * 1000003ULL + (uint64_t)idx * 7919ULL + rs + 1;
            if (rs % 2 == 0){
                uint32_t adjl[MAXN];
                for (int i = 0; i < g.n; i++) adjl[i] = 0;
                for (int e = 0; e < g.m; e++){
                    adjl[g.ei[e]] |= 1u << g.ej[e];
                    adjl[g.ej[e]] |= 1u << g.ei[e];
                }
                place_init(adjl, g.n, p);
            } else {
                double sg = sigmas[(rs/2) % 3];
                for (int k = 0; k < g.n*D; k++) p[k] = sg * gauss_rand();
            }
            double F = lm_run(&g, p, 400);
            if (F < bestF) bestF = F;
            if (F < 1e-20) nsol++;
            /* distinctness of THIS run's endpoint */
            double mind = 1e30;
            for (int i = 0; i < g.n; i++)
                for (int j = i+1; j < g.n; j++){
                    double s = 0;
                    for (int d = 0; d < D; d++){
                        double t = p[i*D+d]-p[j*D+d]; s += t*t;
                    }
                    if (s < mind) mind = s;
                }
            mind = sqrt(mind);
            if (mind > 1e-3 && F < bestFd){ bestFd = F; bestFd_mind = mind; }
            if (mind > 1e-3 && F < 1e-20) nsold++;
            if (nsold >= 3) break;
        }
        printf("%ld %d %.3e %.3e %.3e %d %d %d\n", idx, g.n,
               bestF, bestFd, bestFd_mind, nsol, nsold, g.m);
        fflush(stdout);
        if (bestFd < 1e-20 && bestFd_mind > 1e-3)
            fprintf(stderr, "REALIZABLE-CANDIDATE idx=%ld n=%d F=%.3e mind=%.3e\n",
                    idx, g.n, bestFd, bestFd_mind);
    }
    fclose(f);
    return 0;
}
