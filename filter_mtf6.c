/* Filter maximal triangle-free graphs (multicode stream from triangleramsey)
 * down to minimal abstract almost-equidistant graphs in R^6:
 *   input:  H on n vertices, maximal triangle-free (verified from scratch);
 *   accept: G = complement(H) has no K_8 (alpha(H) <= 7) and no
 *           K_{1,3,3,3} subgraph (BPSSV Lemma 11, even d).
 * Accepted G: alpha(G) <= 2 and minimal, automatically.
 *
 * Output: one line per accepted graph: n then n adjacency masks of G.
 * Usage: filter_mtf6 <file.mc|-> [n] > out.txt
 * Expected (BPSSV Table 2, d=6): n=19: 3,971,787 of 13,734,745.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#ifdef _MSC_VER
#include <intrin.h>
static inline int POP(uint32_t x){ return (int)__popcnt(x); }
static inline int CTZ(uint32_t x){ unsigned long i; _BitScanForward(&i, x); return (int)i; }
#else
static inline int POP(uint32_t x){ return __builtin_popcount(x); }
static inline int CTZ(uint32_t x){ return __builtin_ctz(x); }
#endif

#define MAXV 30

static int is_triangle_free(const uint32_t *adj, int n){
    for (int v = 0; v < n; v++){
        uint32_t t = adj[v];
        while (t){
            int u = CTZ(t); t &= t - 1;
            if (adj[u] & adj[v] & (~0u << (u + 1))) return 0;
        }
    }
    return 1;
}
static int is_maximal_tf(const uint32_t *adj, int n){
    for (int u = 0; u < n; u++)
        for (int v = u + 1; v < n; v++)
            if (!((adj[u] >> v) & 1) && !(adj[u] & adj[v])) return 0;
    return 1;
}
static int has_clique(const uint32_t *adj, uint32_t cand, int size, int want){
    if (size == want) return 1;
    while (cand){
        int v = CTZ(cand); cand &= cand - 1;
        if (has_clique(adj, cand & adj[v], size + 1, want)) return 1;
    }
    return 0;
}
/* K_{1,3,3,3} subgraph: disjoint triples T1,T2,T3, all 27 cross edges,
   plus an apex adjacent to all nine. Enumerate complete-bipartite triple
   pairs (T2,T3); B = common neighbours of T2 u T3 (minus the triples);
   need an apex a in B with >=3 further B-vertices adjacent to a (=T1). */
static int has_k1333(const uint32_t *adj, int n){
    int t2[3];
    for (t2[0] = 0; t2[0] < n; t2[0]++)
    for (t2[1] = t2[0]+1; t2[1] < n; t2[1]++)
    for (t2[2] = t2[1]+1; t2[2] < n; t2[2]++){
        uint32_t T2 = (1u<<t2[0])|(1u<<t2[1])|(1u<<t2[2]);
        uint32_t cn2 = adj[t2[0]] & adj[t2[1]] & adj[t2[2]] & ~T2;
        if (POP(cn2) < 7) continue;      /* need T3(3) + T1(3) + apex(1) */
        int lst[MAXV], k = 0; uint32_t s = cn2;
        while (s){ int i = CTZ(s); s &= s-1; lst[k++] = i; }
        for (int a = 0; a < k; a++)
        for (int b = a+1; b < k; b++)
        for (int c = b+1; c < k; c++){
            uint32_t T3 = (1u<<lst[a])|(1u<<lst[b])|(1u<<lst[c]);
            uint32_t cn3 = adj[lst[a]] & adj[lst[b]] & adj[lst[c]];
            uint32_t B = cn2 & cn3 & ~T3 & ~T2;
            if (POP(B) < 4) continue;    /* apex + triple T1 */
            uint32_t bb = B;
            while (bb){
                int ap = CTZ(bb); bb &= bb - 1;
                if (POP(B & adj[ap] & ~(1u << ap)) >= 3) return 1;
            }
        }
    }
    return 0;
}

int main(int argc, char **argv){
    if (argc < 2){ fprintf(stderr, "usage: filter_mtf6 <file.mc|-> [n]\n"); return 1; }
    FILE *f = strcmp(argv[1], "-") ? fopen(argv[1], "rb") : stdin;
    if (!f){ perror("open"); return 1; }
    int expn = argc > 2 ? atoi(argv[2]) : -1;
    long total = 0, kept = 0, k8 = 0, k1333 = 0, badver = 0;
    for (;;){
        int c = fgetc(f);
        if (c == EOF) break;
        int n = c;
        if (n < 3 || n >= MAXV){ fprintf(stderr, "bad n=%d\n", n); return 1; }
        if (expn > 0 && n != expn){ fprintf(stderr, "unexpected n=%d\n", n); return 1; }
        uint32_t H[MAXV]; memset(H, 0, sizeof H);
        for (int x = 1; x <= n - 1; x++){
            for (;;){
                int y = fgetc(f);
                if (y == EOF){ fprintf(stderr, "EOF mid-graph\n"); return 1; }
                if (y == 0) break;
                H[x-1] |= 1u << (y-1); H[y-1] |= 1u << (x-1);
            }
        }
        total++;
        if (!is_triangle_free(H, n) || !is_maximal_tf(H, n)){ badver++; continue; }
        uint32_t G[MAXV];
        uint32_t full = (1u << n) - 1;
        for (int i = 0; i < n; i++) G[i] = full & ~H[i] & ~(1u << i);
        if (has_clique(G, full, 0, 8)){ k8++; continue; }
        if (has_k1333(G, n)){ k1333++; continue; }
        kept++;
        printf("%d", n);
        for (int i = 0; i < n; i++) printf(" %u", G[i]);
        printf("\n");
    }
    if (f != stdin) fclose(f);
    fprintf(stderr, "total %ld  kept %ld  rejected: K8 %ld, K1333 %ld, badver %ld\n",
            total, kept, k8, k1333, badver);
    if (badver){ fprintf(stderr, "VERIFICATION FAILURE\n"); return 1; }
    return 0;
}
