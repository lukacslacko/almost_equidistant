/* ext6.c — one-vertex extension enumerator + frontier test for the d=6
 * almost-equidistant campaign.
 *
 * Corpus line format (as produced by filter_mtf6): "n m0 m1 ... m_{n-1}"
 * where m_i are adjacency bitmasks of the candidate graph G (edges = unit
 * distances); H = complement(G) is maximal triangle-free; G has no K_8 and
 * no K_{1,3,3,3} (BPSSV Lemma 11, even d).
 *
 * Extension parametrization (see STATUS.md 2026-08-07): every level-(n+1)
 * candidate C with C - w  >=  R (required-edge supergraph of corpus member
 * R on the old vertices) is determined by N = the non-neighborhood of the
 * new vertex w in C:
 *     D    = all edges of H0 = comp(R) inside N   (these flip to G-edges),
 *     C[V] = R + D,   w adjacent to exactly V \ N.
 * Validity is checked DIRECTLY on C: comp(C) triangle-free (theory: always
 * holds; counted as anomaly if not) and maximal, no K_8, no K_{1,3,3,3}.
 * |N| <= 7 is forced (N becomes a clique; K_8-free), and N must dominate
 * V \ N in H0 (part of maximality; used as a cheap prefilter).
 *
 * Modes:
 *   ext6 canontest <corpus> [nperm]     canonical-form invariance self-test
 *   ext6 canon <corpus>                 canonical corpus line per input line
 *   ext6 extend <corpus> <sel|-> [prov] extensions of selected 0-based lines
 *                                       (dedup by canonical form) -> stdout;
 *                                       provenance "parent Nmask" lines -> prov
 *   ext6 frontier <extfile> <resfile>   phase 2: for every graph C in extfile
 *                                       and every vertex u, greedy-complete
 *                                       comp(C-u) to maximal TF; its
 *                                       complement must be in resfile's
 *                                       canonical set, else C is killed.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#define MAXV 30
#define BSU64 ((MAXV*(MAXV-1)/2 + 63)/64)
typedef uint32_t u32;
typedef uint64_t u64;

#ifdef _MSC_VER
#include <intrin.h>
static inline int POP(u32 x){ return (int)__popcnt(x); }
static inline int CTZ(u32 x){ unsigned long i; _BitScanForward(&i, x); return (int)i; }
#else
static inline int POP(u32 x){ return __builtin_popcount(x); }
static inline int CTZ(u32 x){ return __builtin_ctz(x); }
#endif

/* ---- verified subgraph tests, verbatim from filter_mtf6.c (MAXV bumped) ---- */
static int is_triangle_free(const u32 *adj, int n){
    for (int v = 0; v < n; v++){
        u32 t = adj[v];
        while (t){
            int u = CTZ(t); t &= t - 1;
            if (adj[u] & adj[v] & (~0u << (u + 1))) return 0;
        }
    }
    return 1;
}
static int is_maximal_tf(const u32 *adj, int n){
    for (int u = 0; u < n; u++)
        for (int v = u + 1; v < n; v++)
            if (!((adj[u] >> v) & 1) && !(adj[u] & adj[v])) return 0;
    return 1;
}
static int has_clique(const u32 *adj, u32 cand, int size, int want){
    if (size == want) return 1;
    while (cand){
        int v = CTZ(cand); cand &= cand - 1;
        if (has_clique(adj, cand & adj[v], size + 1, want)) return 1;
    }
    return 0;
}
static int has_k1333(const u32 *adj, int n){
    int t2[3];
    for (t2[0] = 0; t2[0] < n; t2[0]++)
    for (t2[1] = t2[0]+1; t2[1] < n; t2[1]++)
    for (t2[2] = t2[1]+1; t2[2] < n; t2[2]++){
        u32 T2 = (1u<<t2[0])|(1u<<t2[1])|(1u<<t2[2]);
        u32 cn2 = adj[t2[0]] & adj[t2[1]] & adj[t2[2]] & ~T2;
        if (POP(cn2) < 7) continue;
        int lst[MAXV], k = 0; u32 s = cn2;
        while (s){ int i = CTZ(s); s &= s-1; lst[k++] = i; }
        for (int a = 0; a < k; a++)
        for (int b = a+1; b < k; b++)
        for (int c = b+1; c < k; c++){
            u32 T3 = (1u<<lst[a])|(1u<<lst[b])|(1u<<lst[c]);
            u32 cn3 = adj[lst[a]] & adj[lst[b]] & adj[lst[c]];
            u32 B = cn2 & cn3 & ~T3 & ~T2;
            if (POP(B) < 4) continue;
            u32 bb = B;
            while (bb){
                int ap = CTZ(bb); bb &= bb - 1;
                if (POP(B & adj[ap] & ~(1u << ap)) >= 3) return 1;
            }
        }
    }
    return 0;
}

/* ------------------------- graphs and parsing ------------------------- */
typedef struct { int n; u32 a[MAXV]; } G;

static int read_graph(FILE *f, G *g){
    if (fscanf(f, "%d", &g->n) != 1) return 0;
    if (g->n < 1 || g->n > MAXV){ fprintf(stderr, "bad n=%d\n", g->n); exit(1); }
    for (int i = 0; i < g->n; i++)
        if (fscanf(f, "%u", &g->a[i]) != 1){ fprintf(stderr, "short line\n"); exit(1); }
    return 1;
}
static void print_graph(FILE *f, const G *g){
    fprintf(f, "%d", g->n);
    for (int i = 0; i < g->n; i++) fprintf(f, " %u", g->a[i]);
    fprintf(f, "\n");
}
static void sanity_graph(const G *g, long lineno, const char *tag){
    u32 full = (g->n == 32) ? ~0u : ((1u << g->n) - 1);
    for (int i = 0; i < g->n; i++){
        if (g->a[i] & ~full){ fprintf(stderr, "%s line %ld: mask overflow\n", tag, lineno); exit(1); }
        if ((g->a[i] >> i) & 1){ fprintf(stderr, "%s line %ld: self-loop\n", tag, lineno); exit(1); }
        for (int j = 0; j < g->n; j++)
            if (((g->a[i] >> j) & 1) != ((g->a[j] >> i) & 1)){
                fprintf(stderr, "%s line %ld: asymmetric\n", tag, lineno); exit(1);
            }
    }
}

/* --------------------- canonical form (mini-nauty) ---------------------
 * Exact refinement: vertex signature = (color, sorted colors of neighbors),
 * classes split by full lexicographic comparison (no hashing anywhere).
 * Search: individualize each vertex of the first non-singleton class (class
 * order = color order, an isomorphism-invariant choice), refine, recurse.
 * Canonical form = lexicographically LARGEST upper-triangle bitstring over
 * all leaves; no pruning, so the leaf set — hence the max — is invariant. */
static const G *cg;
static u64 cbest[BSU64];
static int cbest_set;
static int cperm_best[MAXV];
static u64 cleaves;
/* node-invariant (partition trace) pruning: at each depth store the best
 * invariant vector [nc, class sizes in color order] seen; prune nodes with
 * a lex-smaller vector; a lex-greater vector invalidates the stored best
 * leaf and all deeper vectors. The trace is an isomorphism invariant and
 * ties are kept, so the surviving leaf set — hence the max — is invariant. */
static int tinv[MAXV+1][MAXV+2], tinv_len[MAXV+1], tinv_depth;
/* discovered automorphisms (leaf bitstring == current best leaf): used for
 * orbit pruning of branch choices under the pointwise stabilizer of the
 * individualized prefix. Orbit-equivalent siblings generate isomorphic
 * subtrees with identical leaf-bitstring sets, so skipping them is sound.
 * Automorphisms are properties of the graph, so trace-invalidations of the
 * best leaf do not invalidate stored generators. Each stored generator is
 * re-verified to preserve adjacency before use. */
#define MAXGENS 512
static int agen[MAXGENS][MAXV], nagen;
static int cprefix[MAXV+1];
static int cleafperm_valid;
static int cleafperm[MAXV];   /* permutation of the current best leaf */

static int uf_find(int *p, int x){ while (p[x] != x){ p[x] = p[p[x]]; x = p[x]; } return x; }

static int sigcmp(const int *x, int lx, const int *y, int ly){
    int m = lx < ly ? lx : ly;
    for (int i = 0; i < m; i++) if (x[i] != y[i]) return x[i] < y[i] ? -1 : 1;
    return lx == ly ? 0 : (lx < ly ? -1 : 1);
}
static int refine(int n, const u32 *adj, int *col){
    int sig[MAXV][MAXV+1], sl[MAXV], ord[MAXV], newcol[MAXV];
    for (;;){
        for (int v = 0; v < n; v++){
            int k = 0; sig[v][k++] = col[v];
            int tmp[MAXV], t = 0;
            u32 m = adj[v];
            while (m){ int u = CTZ(m); m &= m-1; tmp[t++] = col[u]; }
            for (int i = 1; i < t; i++){ int x = tmp[i], j = i;
                while (j > 0 && tmp[j-1] > x){ tmp[j] = tmp[j-1]; j--; } tmp[j] = x; }
            for (int i = 0; i < t; i++) sig[v][k++] = tmp[i];
            sl[v] = k;
        }
        for (int v = 0; v < n; v++) ord[v] = v;
        for (int i = 1; i < n; i++){ int x = ord[i], j = i;
            while (j > 0 && sigcmp(sig[ord[j-1]], sl[ord[j-1]], sig[x], sl[x]) > 0){
                ord[j] = ord[j-1]; j--; } ord[j] = x; }
        int nc = 0;
        for (int i = 0; i < n; i++){
            if (i > 0 && sigcmp(sig[ord[i-1]], sl[ord[i-1]], sig[ord[i]], sl[ord[i]]) != 0) nc++;
            newcol[ord[i]] = nc;
        }
        nc++;
        int changed = 0;
        for (int v = 0; v < n; v++) if (newcol[v] != col[v]) changed = 1;
        memcpy(col, newcol, n * sizeof(int));
        if (!changed || nc == n) return nc;
    }
}
static void cleaf(const int *col){
    int n = cg->n, perm[MAXV];
    for (int v = 0; v < n; v++) perm[col[v]] = v;
    u64 bs[BSU64]; memset(bs, 0, sizeof bs);
    int bit = 0;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++){
            if ((cg->a[perm[i]] >> perm[j]) & 1) bs[bit >> 6] |= 1ull << (bit & 63);
            bit++;
        }
    if (!cbest_set || memcmp(bs, cbest, sizeof bs) > 0){
        memcpy(cbest, bs, sizeof bs);
        memcpy(cperm_best, perm, sizeof perm);
        cbest_set = 1;
    } else if (memcmp(bs, cbest, sizeof bs) == 0 && nagen < MAXGENS){
        /* equal leaves: sigma maps cperm_best[i] -> perm[i]; verify + store */
        int sigma[MAXV], ident = 1;
        for (int i = 0; i < n; i++){
            sigma[cperm_best[i]] = perm[i];
            if (cperm_best[i] != perm[i]) ident = 0;
        }
        if (!ident){
            int okau = 1;
            for (int u = 0; u < n && okau; u++)
                for (int v = 0; v < n && okau; v++)
                    if (((cg->a[u] >> v) & 1) !=
                        ((cg->a[sigma[u]] >> sigma[v]) & 1)) okau = 0;
            if (okau){
                memcpy(agen[nagen], sigma, sizeof sigma);
                nagen++;
            }
        }
    }
}
static long g_line = -1;
static void csearch(int *col, int depth){
    int n = cg->n;
    int nc = refine(n, cg->a, col);
    if (++cleaves > 50000000ull){
        fprintf(stderr, "canon blowup at input line %ld (gens=%d)\n", g_line, nagen);
        exit(2);
    }
    /* node invariant: [nc, class sizes in color order] */
    int inv[MAXV+2], il = 0;
    inv[il++] = nc;
    for (int c = 0; c < nc; c++){
        int cnt = 0;
        for (int v = 0; v < n; v++) if (col[v] == c) cnt++;
        inv[il++] = cnt;
    }
    if (depth >= tinv_depth){
        memcpy(tinv[depth], inv, il * sizeof(int));
        tinv_len[depth] = il;
        tinv_depth = depth + 1;
    } else {
        int c = sigcmp(inv, il, tinv[depth], tinv_len[depth]);
        if (c < 0) return;                      /* worse trace: prune */
        if (c > 0){                             /* better: restart bests */
            memcpy(tinv[depth], inv, il * sizeof(int));
            tinv_len[depth] = il;
            tinv_depth = depth + 1;
            cbest_set = 0;
        }
    }
    if (nc == n){ cleaf(col); return; }
    int target = -1;
    for (int c = 0; c < nc && target < 0; c++){
        int cnt = 0;
        for (int v = 0; v < n; v++) if (col[v] == c) cnt++;
        if (cnt >= 2) target = c;
    }
    /* orbits of the subgroup of discovered automorphisms fixing the prefix */
    int uf[MAXV];
    for (int v = 0; v < n; v++) uf[v] = v;
    for (int gi = 0; gi < nagen; gi++){
        int fixes = 1;
        for (int d = 0; d < depth && fixes; d++)
            if (agen[gi][cprefix[d]] != cprefix[d]) fixes = 0;
        if (!fixes) continue;
        for (int v = 0; v < n; v++){
            int a = uf_find(uf, v), b = uf_find(uf, agen[gi][v]);
            if (a != b) uf[a] = b;
        }
    }
    int tried_roots[MAXV], ntried = 0;
    for (int v = 0; v < n; v++){
        if (col[v] != target) continue;
        int r = uf_find(uf, v), dup = 0;
        for (int i = 0; i < ntried && !dup; i++)
            if (tried_roots[i] == r) dup = 1;
        if (dup) continue;
        tried_roots[ntried++] = r;
        int col2[MAXV];
        memcpy(col2, col, n * sizeof(int));
        col2[v] = nc;            /* individualize v out of its class */
        cprefix[depth] = v;
        csearch(col2, depth + 1);
    }
}
static void canon_graph(const G *in, G *out){
    int col[MAXV];
    cg = in; cbest_set = 0; cleaves = 0; tinv_depth = 0; nagen = 0;
    for (int v = 0; v < in->n; v++) col[v] = 0;
    csearch(col, 0);
    out->n = in->n;
    for (int i = 0; i < in->n; i++) out->a[i] = 0;
    int bit = 0;
    for (int i = 0; i < in->n; i++)
        for (int j = i + 1; j < in->n; j++){
            if ((cbest[bit >> 6] >> (bit & 63)) & 1){
                out->a[i] |= 1u << j; out->a[j] |= 1u << i;
            }
            bit++;
        }
}

/* ----------------------------- hash set ----------------------------- */
typedef struct { G *arr; u32 *tab; size_t cap, count, mask; } Set;
static u64 ghash(const G *g){
    u64 h = 1469598103934665603ull;
    h = (h ^ (u64)g->n) * 1099511628211ull;
    for (int i = 0; i < g->n; i++) h = (h ^ g->a[i]) * 1099511628211ull;
    return h;
}
static int geq(const G *x, const G *y){
    if (x->n != y->n) return 0;
    return memcmp(x->a, y->a, x->n * sizeof(u32)) == 0;
}
static void set_init(Set *s, size_t cap){
    s->cap = cap; s->mask = cap - 1; s->count = 0;
    s->arr = malloc(cap * sizeof(G));
    s->tab = malloc(cap * sizeof(u32));
    memset(s->tab, 0xff, cap * sizeof(u32));
}
static int set_insert(Set *s, const G *g);
static void set_grow(Set *s){
    size_t nc = s->cap * 2;
    G *oa = s->arr; size_t ocount = s->count;
    Set ns; set_init(&ns, nc);
    for (size_t i = 0; i < ocount; i++) set_insert(&ns, &oa[i]);
    free(s->arr); free(s->tab);
    *s = ns;
}
static int set_insert(Set *s, const G *g){   /* 1 if new */
    if (s->count * 10 >= s->cap * 7) set_grow(s);
    u64 h = ghash(g);
    size_t i = h & s->mask;
    while (s->tab[i] != 0xffffffffu){
        if (geq(&s->arr[s->tab[i]], g)) return 0;
        i = (i + 1) & s->mask;
    }
    s->arr[s->count] = *g;
    s->tab[i] = (u32)s->count++;
    return 1;
}
static int set_has(const Set *s, const G *g){
    u64 h = ghash(g);
    size_t i = h & s->mask;
    while (s->tab[i] != 0xffffffffu){
        if (geq(&s->arr[s->tab[i]], g)) return 1;
        i = (i + 1) & s->mask;
    }
    return 0;
}

/* --------------------------- rng for tests --------------------------- */
static u64 rstate = 0x9e3779b97f4a7c15ull;
static u32 rnd(u32 m){ rstate ^= rstate << 13; rstate ^= rstate >> 7; rstate ^= rstate << 17;
    return (u32)(rstate % m); }

/* ------------------------------ modes ------------------------------ */
static int mode_canontest(const char *fn, int nperm){
    FILE *f = fopen(fn, "r");
    if (!f){ perror(fn); return 1; }
    G g, cref, gp, cp;
    long line = 0;
    u64 worst = 0; long worstline = -1;
    while (read_graph(f, &g)){
        sanity_graph(&g, line, "canontest");
        g_line = line;
        canon_graph(&g, &cref);
        if (cleaves > worst){ worst = cleaves; worstline = line; }
        for (int t = 0; t < nperm; t++){
            int perm[MAXV];
            for (int i = 0; i < g.n; i++) perm[i] = i;
            for (int i = g.n - 1; i > 0; i--){
                int j = (int)rnd((u32)(i + 1));
                int tmpv = perm[i]; perm[i] = perm[j]; perm[j] = tmpv;
            }
            gp.n = g.n;
            for (int i = 0; i < g.n; i++) gp.a[i] = 0;
            for (int i = 0; i < g.n; i++)
                for (int j = 0; j < g.n; j++)
                    if ((g.a[i] >> j) & 1) gp.a[perm[i]] |= 1u << perm[j];
            canon_graph(&gp, &cp);
            if (!geq(&cref, &cp)){
                fprintf(stderr, "CANON MISMATCH at line %ld perm %d\n", line, t);
                return 1;
            }
        }
        line++;
    }
    fclose(f);
    printf("canontest PASS: %ld graphs x %d perms (worst nodes %llu at line %ld)\n",
           line, nperm, (unsigned long long)worst, worstline);
    return 0;
}

static int mode_canon(const char *fn){
    FILE *f = fopen(fn, "r");
    if (!f){ perror(fn); return 1; }
    G g, c;
    long line = 0;
    while (read_graph(f, &g)){
        sanity_graph(&g, line++, "canon");
        canon_graph(&g, &c);
        print_graph(stdout, &c);
    }
    fclose(f);
    return 0;
}

/* enumerate subsets N (size 1..7) of [0,n) that dominate V\N in H0 */
static long ext_anomalies = 0;
static void extend_one(const G *R, long lineno, Set *seen, FILE *prov,
                       long *nvalid){
    int n = R->n;
    u32 full = (1u << n) - 1;
    u32 H0[MAXV];
    for (int i = 0; i < n; i++) H0[i] = full & ~R->a[i] & ~(1u << i);
    int idx[7];
    for (int k = 1; k <= 7; k++){
        /* combinations idx[0..k-1] increasing */
        for (int i = 0; i < k; i++) idx[i] = i;
        for (;;){
            u32 N = 0;
            for (int i = 0; i < k; i++) N |= 1u << idx[i];
            /* prefilter: N dominates V\N in H0 */
            int ok = 1;
            for (int u = 0; u < n && ok; u++)
                if (!((N >> u) & 1) && !(H0[u] & N)) ok = 0;
            if (ok){
                G C;
                C.n = n + 1;
                for (int i = 0; i < n; i++) C.a[i] = R->a[i];
                /* D = H0 edges inside N flip to G-edges */
                u32 t = N;
                while (t){
                    int u = CTZ(t); t &= t - 1;
                    u32 s = H0[u] & N & (~0u << (u + 1));
                    while (s){
                        int v = CTZ(s); s &= s - 1;
                        C.a[u] |= 1u << v; C.a[v] |= 1u << u;
                    }
                }
                C.a[n] = full & ~N;
                for (int i = 0; i < n; i++)
                    if (!((N >> i) & 1)) C.a[i] |= 1u << n;
                /* direct validity checks on C */
                u32 fulln1 = (1u << (n + 1)) - 1;
                u32 H[MAXV];
                for (int i = 0; i <= n; i++) H[i] = fulln1 & ~C.a[i] & ~(1u << i);
                if (!is_triangle_free(H, n + 1)){
                    ext_anomalies++;       /* theory says impossible */
                } else if (is_maximal_tf(H, n + 1)
                           && !has_clique(C.a, fulln1, 0, 8)
                           && !has_k1333(C.a, n + 1)){
                    (*nvalid)++;
                    G canon;
                    canon_graph(&C, &canon);
                    if (set_insert(seen, &canon)){
                        print_graph(stdout, &canon);
                        if (prov) fprintf(prov, "%ld %u\n", lineno, N);
                    }
                }
            }
            /* next combination */
            int i = k - 1;
            while (i >= 0 && idx[i] == n - k + i) i--;
            if (i < 0) break;
            idx[i]++;
            for (int j = i + 1; j < k; j++) idx[j] = idx[j-1] + 1;
        }
    }
}

static int mode_extend(const char *fn, const char *sel, const char *provfn){
    FILE *f = fopen(fn, "r");
    if (!f){ perror(fn); return 1; }
    char *want = NULL;
    long nwant = -1;
    if (strcmp(sel, "-") != 0){
        FILE *sf = fopen(sel, "r");
        if (!sf){ perror(sel); return 1; }
        long cap = 1 << 22;
        want = calloc(cap, 1);
        long v, mx = -1;
        while (fscanf(sf, "%ld", &v) == 1){
            if (v < 0 || v >= cap){ fprintf(stderr, "sel index %ld out of range\n", v); return 1; }
            want[v] = 1; if (v > mx) mx = v;
        }
        fclose(sf);
        nwant = mx;
    }
    FILE *prov = provfn ? fopen(provfn, "w") : NULL;
    Set seen; set_init(&seen, 1 << 20);
    G g;
    long line = 0, used = 0, nvalid = 0;
    while (read_graph(f, &g)){
        sanity_graph(&g, line, "extend");
        if (!want || want[line]){
            extend_one(&g, line, &seen, prov, &nvalid);
            used++;
        }
        line++;
        if (want && nwant >= 0 && line > nwant) break;
    }
    fclose(f);
    if (prov) fclose(prov);
    fprintf(stderr, "extend: %ld parents scanned, %ld used, %ld valid extensions "
            "(with multiplicity), %zu distinct up to iso, %ld TF-anomalies\n",
            line, used, nvalid, seen.count, ext_anomalies);
    if (ext_anomalies){ fprintf(stderr, "ANOMALY: complement-TF violated\n"); return 1; }
    return 0;
}

/* greedy completion of a TF graph to maximal TF: repeatedly add the first
 * addable pair (lexicographic), until none remains */
static void greedy_complete(u32 *H, int n){
    for (;;){
        int added = 0;
        for (int u = 0; u < n && !added; u++)
            for (int v = u + 1; v < n && !added; v++)
                if (!((H[u] >> v) & 1) && !(H[u] & H[v])){
                    H[u] |= 1u << v; H[v] |= 1u << u;
                    added = 1;
                }
        if (!added) return;
    }
}

static int mode_frontier(const char *extfn, const char *resfn){
    FILE *rf = fopen(resfn, "r");
    if (!rf){ perror(resfn); return 1; }
    Set S; set_init(&S, 1 << 12);
    G g, c;
    long rn = 0;
    int resn = -1;
    while (read_graph(rf, &g)){
        sanity_graph(&g, rn++, "residue");
        if (resn < 0) resn = g.n;
        else if (resn != g.n){ fprintf(stderr, "residue size mix\n"); return 1; }
        canon_graph(&g, &c);
        set_insert(&S, &c);
    }
    fclose(rf);
    fprintf(stderr, "frontier: residue set %ld graphs (%zu distinct canon) on %d vertices\n",
            rn, S.count, resn);
    FILE *ef = fopen(extfn, "r");
    if (!ef){ perror(extfn); return 1; }
    long line = 0, killed = 0, survive = 0;
    while (read_graph(ef, &g)){
        sanity_graph(&g, line, "ext");
        if (g.n != resn + 1){ fprintf(stderr, "ext size %d != residue+1\n", g.n); return 1; }
        int dead = -1;
        for (int u = 0; u < g.n && dead < 0; u++){
            /* delete u, relabel compactly */
            G d; d.n = g.n - 1;
            int map[MAXV], k = 0;
            for (int i = 0; i < g.n; i++) if (i != u) map[i] = k++;
            for (int i = 0; i < d.n; i++) d.a[i] = 0;
            for (int i = 0; i < g.n; i++){
                if (i == u) continue;
                for (int j = 0; j < g.n; j++){
                    if (j == u || !((g.a[i] >> j) & 1)) continue;
                    d.a[map[i]] |= 1u << map[j];
                }
            }
            u32 fulld = (1u << d.n) - 1;
            u32 H[MAXV];
            for (int i = 0; i < d.n; i++) H[i] = fulld & ~d.a[i] & ~(1u << i);
            greedy_complete(H, d.n);
            G c3;
            c3.n = d.n;
            for (int i = 0; i < d.n; i++) c3.a[i] = fulld & ~H[i] & ~(1u << i);
            G cc;
            canon_graph(&c3, &cc);
            if (!set_has(&S, &cc)) dead = u;
        }
        if (dead >= 0){ killed++; printf("KILLED %ld del=%d\n", line, dead); }
        else { survive++; printf("SURVIVES %ld ", line); print_graph(stdout, &g); }
        line++;
    }
    fclose(ef);
    fprintf(stderr, "frontier: %ld graphs, %ld killed by greedy re-completion, %ld survive\n",
            line, killed, survive);
    return 0;
}

/* escalate mode: like frontier, but enumerates ALL maximal TF completions
 * of every deletion (branching: first addable pair either included, or
 * excluded by adding a blocking cherry through each possible centre).
 * A graph survives only if every completion of every deletion is in the
 * residue set. */
static Set esc_S;
static int esc_all_in(u32 *H, int n, long *budget){
    /* recursively complete H; return 0 as soon as some completion's
       complement is outside the residue set */
    if (--(*budget) < 0){ fprintf(stderr, "escalate budget exceeded\n"); exit(2); }
    int fu = -1, fv = -1;
    for (int u = 0; u < n && fu < 0; u++)
        for (int v = u + 1; v < n; v++)
            if (!((H[u] >> v) & 1) && !(H[u] & H[v])){ fu = u; fv = v; break; }
    if (fu < 0){
        u32 fulln = (1u << n) - 1;
        G c3, cc;
        c3.n = n;
        for (int i = 0; i < n; i++) c3.a[i] = fulln & ~H[i] & ~(1u << i);
        canon_graph(&c3, &cc);
        return set_has(&esc_S, &cc);
    }
    u32 H2[MAXV];
    memcpy(H2, H, n * sizeof(u32));
    H2[fu] |= 1u << fv; H2[fv] |= 1u << fu;
    if (!esc_all_in(H2, n, budget)) return 0;
    for (int z = 0; z < n; z++){
        if (z == fu || z == fv) continue;
        memcpy(H2, H, n * sizeof(u32));
        int good = 1;
        int pairs[2][2] = {{fu, z}, {fv, z}};
        for (int pi = 0; pi < 2 && good; pi++){
            int x = pairs[pi][0], y = pairs[pi][1];
            if ((H2[x] >> y) & 1) continue;
            if (H2[x] & H2[y]){ good = 0; break; }
            H2[x] |= 1u << y; H2[y] |= 1u << x;
        }
        if (good && !esc_all_in(H2, n, budget)) return 0;
    }
    return 1;
}
static int mode_escalate(const char *extfn, const char *resfn){
    FILE *rf = fopen(resfn, "r");
    if (!rf){ perror(resfn); return 1; }
    set_init(&esc_S, 1 << 18);
    G g, c;
    long rn = 0;
    int resn = -1;
    while (read_graph(rf, &g)){
        sanity_graph(&g, rn++, "residue");
        if (resn < 0) resn = g.n;
        canon_graph(&g, &c);
        set_insert(&esc_S, &c);
    }
    fclose(rf);
    fprintf(stderr, "escalate: residue %ld graphs (%zu canon)\n", rn, esc_S.count);
    FILE *ef = fopen(extfn, "r");
    if (!ef){ perror(extfn); return 1; }
    long line = 0, killed = 0, survive = 0;
    while (read_graph(ef, &g)){
        sanity_graph(&g, line, "esc");
        int dead = 0;
        for (int u = 0; u < g.n && !dead; u++){
            G d; d.n = g.n - 1;
            int map[MAXV], k = 0;
            for (int i = 0; i < g.n; i++) if (i != u) map[i] = k++;
            for (int i = 0; i < d.n; i++) d.a[i] = 0;
            for (int i = 0; i < g.n; i++){
                if (i == u) continue;
                for (int j = 0; j < g.n; j++){
                    if (j == u || !((g.a[i] >> j) & 1)) continue;
                    d.a[map[i]] |= 1u << map[j];
                }
            }
            u32 fulld = (1u << d.n) - 1;
            u32 H[MAXV];
            for (int i = 0; i < d.n; i++) H[i] = fulld & ~d.a[i] & ~(1u << i);
            long budget = 5000000;
            if (!esc_all_in(H, d.n, &budget)) dead = 1;
        }
        if (dead){ killed++; printf("KILLED %ld\n", line); }
        else { survive++; printf("SURVIVES %ld ", line); print_graph(stdout, &g); }
        line++;
    }
    fclose(ef);
    fprintf(stderr, "escalate: %ld graphs, %ld killed, %ld survive\n",
            line, killed, survive);
    return 0;
}

int main(int argc, char **argv){
    if (argc < 3){
        fprintf(stderr, "usage: ext6 canontest|canon|extend|frontier|escalate ...\n");
        return 1;
    }
    if (!strcmp(argv[1], "canontest"))
        return mode_canontest(argv[2], argc > 3 ? atoi(argv[3]) : 8);
    if (!strcmp(argv[1], "canon"))
        return mode_canon(argv[2]);
    if (!strcmp(argv[1], "extend"))
        return mode_extend(argv[2], argc > 3 ? argv[3] : "-", argc > 4 ? argv[4] : NULL);
    if (!strcmp(argv[1], "frontier")){
        if (argc < 4){ fprintf(stderr, "frontier <extfile> <resfile>\n"); return 1; }
        return mode_frontier(argv[2], argv[3]);
    }
    if (!strcmp(argv[1], "escalate")){
        if (argc < 4){ fprintf(stderr, "escalate <extfile> <resfile>\n"); return 1; }
        return mode_escalate(argv[2], argv[3]);
    }
    fprintf(stderr, "unknown mode %s\n", argv[1]);
    return 1;
}
