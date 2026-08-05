/* Independent enumeration of ALL abstract almost-equidistant graphs in R^5
 * (BPSSV Section 3): graphs G with
 *   (1) no independent set of size 3 (complement triangle-free),
 *   (2) no K_7,
 *   (3) no K_{3,3,3} as a (not necessarily induced) subgraph.
 * Grown one vertex at a time (the new vertex has minimum degree, WLOG:
 * the three properties are hereditary under induced subgraphs, so every
 * witness arises by deleting a minimum-degree vertex; and deg >= n-6 since
 * the non-neighbourhood of a vertex is a clique, of size <= 6 by (2)).
 *
 * Expected totals (BPSSV Table 3, d=5), n = 4..13:
 *   7, 14, 38, 106, 402, 1817, 11132, 86053, 803299, 7623096
 * Expected minimal counts (BPSSV Table 2, d=5), n = 10..13: 25, 46, 106, 242
 * ("minimal": removing any edge creates an independent triple, i.e. edge uv
 *  removable iff N(u) u N(v) u {u,v} = V; minimal iff no removable edge.)
 *
 * Single file, no dependencies:  cc -O2 -o enumaeq5 enumaeq5.c
 * Usage: enumaeq5 [maxn]      (default maxn = 13)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

#ifdef _MSC_VER
#include <intrin.h>
static inline int POP(uint32_t x){ return (int)__popcnt(x); }
static inline int CTZ(uint32_t x){ unsigned long i; _BitScanForward(&i, x); return (int)i; }
#else
static inline int POP(uint32_t x){ return __builtin_popcount(x); }
static inline int CTZ(uint32_t x){ return __builtin_ctz(x); }
#endif

#define MAXV 15

typedef struct { uint32_t adj[MAXV]; } G;

/* ------------------------------ storage ------------------------------ */
static G *pool = NULL; static size_t pool_len = 0, pool_cap = 0;
static size_t pool_push(const G *g){
    if (pool_len == pool_cap){
        pool_cap = pool_cap ? pool_cap * 2 : 1u << 20;
        pool = (G*)realloc(pool, pool_cap * sizeof(G));
        if (!pool){ fprintf(stderr, "OOM pool\n"); exit(1); }
    }
    pool[pool_len] = *g; return pool_len++;
}

/* hash table: invariant hash -> chain of graph indices (same level) */
typedef struct HNode { uint64_t h; size_t gi; struct HNode *next; } HNode;
static HNode **ht = NULL; static size_t ht_size = 0;
static HNode *hnodes = NULL; static size_t hn_len = 0, hn_cap = 0;
static void ht_init(size_t size){
    ht_size = size; ht = (HNode**)calloc(ht_size, sizeof(HNode*));
    hn_len = 0;
    if (!ht){ fprintf(stderr, "OOM ht\n"); exit(1); }
}
static HNode *hnode_new(void){
    if (hn_len == hn_cap){
        fprintf(stderr, "OOM hnodes\n"); exit(1);
    }
    return &hnodes[hn_len++];
}

/* ------------------------- invariant + hashing ------------------------- */
/* per-vertex: (deg, #triangles at v, sorted neighbour degrees); the graph
   invariant is the sorted multiset of these; hashed with FNV-1a. */
static uint64_t invariant_hash(const G *g, int n){
    uint64_t items[MAXV];
    for (int v = 0; v < n; v++){
        uint32_t nb = g->adj[v];
        int deg = POP(nb);
        int tri = 0;
        uint32_t t = nb;
        while (t){
            int i = CTZ(t); t &= t - 1;
            tri += POP(g->adj[i] & nb & (~0u << (i + 1)));
        }
        int nd[MAXV], k = 0;
        uint32_t s = nb;
        while (s){ int i = CTZ(s); s &= s - 1; nd[k++] = POP(g->adj[i]); }
        /* insertion sort */
        for (int a = 1; a < k; a++){ int x = nd[a], b = a - 1;
            while (b >= 0 && nd[b] > x){ nd[b+1] = nd[b]; b--; } nd[b+1] = x; }
        uint64_t h = 1469598103934665603ULL;
        h = (h ^ (uint64_t)deg) * 1099511628211ULL;
        h = (h ^ (uint64_t)tri) * 1099511628211ULL;
        for (int a = 0; a < k; a++) h = (h ^ (uint64_t)nd[a]) * 1099511628211ULL;
        items[v] = h;
    }
    for (int a = 1; a < n; a++){ uint64_t x = items[a]; int b = a - 1;
        while (b >= 0 && items[b] > x){ items[b+1] = items[b]; b--; } items[b+1] = x; }
    uint64_t h = 14695981039346656037ULL;
    for (int a = 0; a < n; a++) h = (h ^ items[a]) * 1099511628211ULL;
    return h;
}

/* ------------------------- exact isomorphism ------------------------- */
static int iso_n;
static const G *iso_g1, *iso_g2;
static int iso_order[MAXV], iso_map[MAXV], iso_used[MAXV];
static int iso_d1[MAXV], iso_d2[MAXV];

static int iso_bt(int idx){
    if (idx == iso_n) return 1;
    int v = iso_order[idx];
    for (int w = 0; w < iso_n; w++){
        if (iso_used[w] || iso_d2[w] != iso_d1[v]) continue;
        int ok = 1;
        for (int t = 0; t < idx; t++){
            int u = iso_order[t];
            if (((iso_g1->adj[v] >> u) & 1) != ((iso_g2->adj[w] >> iso_map[u]) & 1)){ ok = 0; break; }
        }
        if (ok){
            iso_map[v] = w; iso_used[w] = 1;
            if (iso_bt(idx + 1)) return 1;
            iso_used[w] = 0;
        }
    }
    return 0;
}

static int isomorphic(const G *a, const G *b, int n){
    int da[MAXV], db[MAXV];
    long sa = 0, sb = 0;
    for (int i = 0; i < n; i++){ da[i] = POP(a->adj[i]); db[i] = POP(b->adj[i]); }
    int ca[MAXV+1] = {0}, cb[MAXV+1] = {0};
    for (int i = 0; i < n; i++){ ca[da[i]]++; cb[db[i]]++; sa += da[i]; sb += db[i]; }
    if (sa != sb) return 0;
    for (int i = 0; i <= n; i++) if (ca[i] != cb[i]) return 0;
    iso_n = n; iso_g1 = a; iso_g2 = b;
    memcpy(iso_d1, da, sizeof da); memcpy(iso_d2, db, sizeof db);
    for (int i = 0; i < n; i++){ iso_order[i] = i; iso_used[i] = 0; iso_map[i] = -1; }
    /* order by (degree) ascending for rarity */
    for (int x = 1; x < n; x++){ int o = iso_order[x], b2 = x - 1;
        while (b2 >= 0 && da[iso_order[b2]] > da[o]){ iso_order[b2+1] = iso_order[b2]; b2--; }
        iso_order[b2+1] = o; }
    return iso_bt(0);
}

/* --------------------------- full validators --------------------------- */
static int has_indep3(const G *g, int n){
    uint32_t full = (n == 32) ? 0xffffffffu : ((1u << n) - 1);
    for (int a = 0; a < n; a++)
        for (int b = a + 1; b < n; b++){
            if ((g->adj[a] >> b) & 1) continue;
            uint32_t m = ~(g->adj[a] | g->adj[b]) & full & ~(1u << a) & ~(1u << b);
            if (m) return 1;
        }
    return 0;
}
static int clique_ext(const G *g, uint32_t cand, int size, int want){
    if (size == want) return 1;
    while (cand){
        int v = CTZ(cand); cand &= cand - 1;
        if (clique_ext(g, cand & g->adj[v], size + 1, want)) return 1;
    }
    return 0;
}
static int has_k7(const G *g, int n){
    uint32_t full = (1u << n) - 1;
    return clique_ext(g, full, 0, 7);
}
/* K_{3,3,3} subgraph: disjoint triples T1,T2,T3, all 27 cross edges */
static int has_k333(const G *g, int n){
    int t2[3], t3[3];
    for (t2[0] = 0; t2[0] < n; t2[0]++)
    for (t2[1] = t2[0]+1; t2[1] < n; t2[1]++)
    for (t2[2] = t2[1]+1; t2[2] < n; t2[2]++){
        uint32_t T2 = (1u<<t2[0])|(1u<<t2[1])|(1u<<t2[2]);
        uint32_t cn2 = g->adj[t2[0]] & g->adj[t2[1]] & g->adj[t2[2]] & ~T2;
        if (POP(cn2) < 6) continue;
        /* choose T3 inside cn2 */
        int lst[MAXV], k = 0; uint32_t s = cn2;
        while (s){ int i = CTZ(s); s &= s-1; lst[k++] = i; }
        for (int a = 0; a < k; a++)
        for (int b = a+1; b < k; b++)
        for (int c = b+1; c < k; c++){
            t3[0]=lst[a]; t3[1]=lst[b]; t3[2]=lst[c];
            uint32_t T3 = (1u<<t3[0])|(1u<<t3[1])|(1u<<t3[2]);
            uint32_t cn3 = g->adj[t3[0]] & g->adj[t3[1]] & g->adj[t3[2]];
            uint32_t both = cn2 & cn3 & ~T3 & ~T2;
            if (POP(both) >= 3) return 1;   /* T1 inside both */
        }
    }
    return 0;
}
/* collect all cliques of size `want` as bitmasks (vertices ascending) */
static void collect_k(const G *g, uint32_t cand, uint32_t cur, int size, int want,
                      uint32_t *out, int *nout, int cap){
    if (size == want){
        if (*nout < cap) out[(*nout)++] = cur;
        else { fprintf(stderr, "clique overflow\n"); exit(1); }
        return;
    }
    while (cand){
        int v = CTZ(cand); cand &= cand - 1;
        collect_k(g, cand & g->adj[v], cur | (1u << v), size + 1, want, out, nout, cap);
    }
}

static int is_minimal(const G *g, int n){
    uint32_t full = (1u << n) - 1;
    for (int u = 0; u < n; u++)
        for (int v = u + 1; v < n; v++)
            if ((g->adj[u] >> v) & 1)
                if ((g->adj[u] | g->adj[v] | (1u<<u) | (1u<<v)) == full)
                    return 0;
    return 1;
}

/* ------------------------------- main ------------------------------- */
static const long EXPECTED_ALL[15] = { /* index n, d=5 totals (Table 3) */
    0,0,0,0, 7, 14, 38, 106, 402, 1817, 11132, 86053, 803299, 7623096, -1 };
static const long EXPECTED_MIN[15] = { /* minimal (Table 2, d=5) */
    0,0,0,0, 2, 3, 4, 6, 9, 14, 25, 46, 106, 242, 653 };

int main(int argc, char **argv){
    int maxn = argc > 1 ? atoi(argv[1]) : 13;
    if (maxn > 14){ fprintf(stderr, "maxn <= 14\n"); return 1; }
    clock_t t0 = clock();

    size_t lvl_start = 0, lvl_len = 1;
    { G g; memset(&g, 0, sizeof g); pool_push(&g); }   /* K_1 */
    int n = 1;

    while (n < maxn && lvl_len > 0){
        int m = n + 1;
        size_t next_start = pool_len;
        /* hash table sized for expected level */
        size_t est = (n + 1 <= 14 && EXPECTED_ALL[n+1] > 0) ?
                     (size_t)EXPECTED_ALL[n+1] * 2 + 1024 : lvl_len * 16 + 1024;
        ht_init(est * 2);
        hn_cap = est * 4 + 1024;
        hnodes = (HNode*)malloc(hn_cap * sizeof(HNode));
        if (!hnodes){ fprintf(stderr, "OOM hn\n"); return 1; }
        long ncand = 0;

        for (size_t pi = lvl_start; pi < lvl_start + lvl_len; pi++){
            G parent = pool[pi];              /* copy: pool may realloc */
            int degs[MAXV], mind = MAXV;
            for (int i = 0; i < n; i++){ degs[i] = POP(parent.adj[i]); if (degs[i] < mind) mind = degs[i]; }
            int lo = m - 7; if (lo < 0) lo = 0;
            int hi = mind + 1; if (hi > n) hi = n;
            /* 6-cliques of parent (for K7 check) */
            static uint32_t cl6[200000]; int ncl6 = 0;
            if (m >= 7)
                collect_k(&parent, (1u << n) - 1, 0, 0, 6, cl6, &ncl6, 200000);
            /* iterate subsets S of {0..n-1} of size k */
            for (int k = lo; k <= hi; k++){
                int c[MAXV];
                for (int i = 0; i < k; i++) c[i] = i;
                while (1){
                    uint32_t S = 0;
                    for (int i = 0; i < k; i++) S |= 1u << c[i];
                    /* min-degree of new vertex */
                    int okd = 1;
                    for (int u = 0; u < n && okd; u++)
                        if (degs[u] + ((S >> u) & 1) < k) okd = 0;
                    if (okd){
                        /* alpha<=2: non-neighbours of v pairwise adjacent */
                        uint32_t nonN = ((1u << n) - 1) & ~S;
                        int okA = 1;
                        uint32_t t = nonN;
                        while (t && okA){
                            int i = CTZ(t); t &= t - 1;
                            if ((nonN & ~parent.adj[i] & ~(1u << i) & (~0u << (i+1)))) okA = 0;
                        }
                        if (okA){
                            /* K7: 6-clique of parent inside S */
                            int okK = 1;
                            for (int q = 0; q < ncl6; q++)
                                if ((cl6[q] & S) == cl6[q]){ okK = 0; break; }
                            if (okK){
                                /* K333 through v: T2,T3 in S w/ complete bipartite,
                                   plus >=2 more common nbrs (any vertices) */
                                G g2 = parent;
                                g2.adj[n] = S;
                                uint32_t tt = S;
                                while (tt){ int i = CTZ(tt); tt &= tt-1; g2.adj[i] |= 1u << n; }
                                int bad = 0;
                                int ls[MAXV], lk = 0; uint32_t s2 = S;
                                while (s2){ int i = CTZ(s2); s2 &= s2-1; ls[lk++] = i; }
                                for (int a = 0; a < lk && !bad; a++)
                                for (int b = a+1; b < lk && !bad; b++)
                                for (int cc = b+1; cc < lk && !bad; cc++){
                                    uint32_t T2 = (1u<<ls[a])|(1u<<ls[b])|(1u<<ls[cc]);
                                    uint32_t cn2 = g2.adj[ls[a]] & g2.adj[ls[b]] & g2.adj[ls[cc]] & ~T2;
                                    uint32_t cn2S = cn2 & S;
                                    if (POP(cn2S) < 3) continue;
                                    int l3[MAXV], k3 = 0; uint32_t s3 = cn2S;
                                    while (s3){ int i = CTZ(s3); s3 &= s3-1; l3[k3++] = i; }
                                    for (int x = 0; x < k3 && !bad; x++)
                                    for (int y = x+1; y < k3 && !bad; y++)
                                    for (int z = y+1; z < k3 && !bad; z++){
                                        uint32_t T3 = (1u<<l3[x])|(1u<<l3[y])|(1u<<l3[z]);
                                        uint32_t cn3 = g2.adj[l3[x]] & g2.adj[l3[y]] & g2.adj[l3[z]];
                                        uint32_t both = cn2 & cn3 & ~T2 & ~T3 & ~(1u << n);
                                        if (POP(both) >= 2) bad = 1;
                                    }
                                }
                                if (!bad){
                                    ncand++;
                                    uint64_t h = invariant_hash(&g2, m);
                                    size_t slot = (size_t)(h % ht_size);
                                    HNode *nd = ht[slot]; int dup = 0;
                                    for (; nd; nd = nd->next)
                                        if (nd->h == h && isomorphic(&g2, &pool[nd->gi], m)){ dup = 1; break; }
                                    if (!dup){
                                        size_t gi = pool_push(&g2);
                                        HNode *nn2 = hnode_new();
                                        nn2->h = h; nn2->gi = gi;
                                        nn2->next = ht[slot]; ht[slot] = nn2;
                                    }
                                }
                            }
                        }
                    }
                    /* next combination */
                    int i = k - 1;
                    while (i >= 0 && c[i] == i + n - k) i--;
                    if (i < 0) break;
                    c[i]++;
                    for (int j = i + 1; j < k; j++) c[j] = c[j-1] + 1;
                }
            }
        }
        free(ht); free(hnodes); ht = NULL; hnodes = NULL;
        n = m;
        lvl_start = next_start; lvl_len = pool_len - next_start;
        /* validate + count minimal */
        long nmin = 0; long bad = 0;
        for (size_t i = lvl_start; i < pool_len; i++){
            if (is_minimal(&pool[i], n)) nmin++;
            if (n <= 11 || (i % 997) == 0){
                if (has_indep3(&pool[i], n) || has_k7(&pool[i], n) || has_k333(&pool[i], n)) bad++;
            }
        }
        long expA = (n <= 14) ? EXPECTED_ALL[n] : -1;
        long expM = (n <= 14) ? EXPECTED_MIN[n] : -1;
        printf("n=%2d: %8zu abstract (expected %ld)  minimal %ld (expected %ld)"
               "  cand %ld  bad %ld  t=%.1fs\n",
               n, lvl_len, expA, nmin, expM, ncand, bad,
               (double)(clock() - t0) / CLOCKS_PER_SEC);
        fflush(stdout);
        if (bad){ fprintf(stderr, "VALIDATION FAILURE\n"); return 1; }
        /* dump level to file for downstream use */
        if (n >= 12){
            char fn[64]; snprintf(fn, sizeof fn, "aeq_d5_n%d.bin", n);
            FILE *f = fopen(fn, "wb");
            if (f){
                uint32_t hdr[2] = { (uint32_t)n, (uint32_t)lvl_len };
                fwrite(hdr, 4, 2, f);
                for (size_t i = lvl_start; i < pool_len; i++)
                    fwrite(pool[i].adj, 4, n, f);
                fclose(f);
            }
        }
    }
    printf("done, %.1fs\n", (double)(clock() - t0) / CLOCKS_PER_SEC);
    return 0;
}
