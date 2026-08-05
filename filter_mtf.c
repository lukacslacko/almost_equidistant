/* Filter maximal triangle-free graphs (multicode stream from triangleramsey)
 * down to minimal abstract almost-equidistant graphs in R^5:
 *   input:  H on n vertices, maximal triangle-free (verified independently);
 *   accept: G = complement(H) has no K_7 (alpha(H) <= 6) and no K_{3,3,3}
 *           subgraph.
 * Every accepted G automatically satisfies: alpha(G) <= 2 (H triangle-free)
 * and G minimal (H maximal). Verified here from scratch for every graph.
 *
 * Output: one line per accepted graph: n then the n adjacency bitmasks of G
 * (decimal). Usage: filter_mtf <file.mc> [n] > out.txt
 * Prints stats to stderr. Expected (BPSSV Table 2, d=5):
 *   n=17: 12654 of 164796; n=18: 8825 of 1337848; n=19: 340 of 13734745;
 *   n=20: 8 of 178587364; n=21: 0 of 2911304940.
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

#define MAXV 24

static int is_triangle_free(const uint32_t *adj, int n){
    for (int v = 0; v < n; v++){
        uint32_t nb = adj[v];
        uint32_t t = nb;
        while (t){
            int u = CTZ(t); t &= t - 1;
            if (adj[u] & nb & (~0u << (u + 1))) return 0;
        }
    }
    return 1;
}
static int is_maximal_tf(const uint32_t *adj, int n){
    /* every non-adjacent pair has a common neighbour */
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
static int has_k333(const uint32_t *adj, int n){
    int t2[3];
    for (t2[0] = 0; t2[0] < n; t2[0]++)
    for (t2[1] = t2[0]+1; t2[1] < n; t2[1]++)
    for (t2[2] = t2[1]+1; t2[2] < n; t2[2]++){
        uint32_t T2 = (1u<<t2[0])|(1u<<t2[1])|(1u<<t2[2]);
        uint32_t cn2 = adj[t2[0]] & adj[t2[1]] & adj[t2[2]] & ~T2;
        if (POP(cn2) < 6) continue;
        int lst[MAXV], k = 0; uint32_t s = cn2;
        while (s){ int i = CTZ(s); s &= s-1; lst[k++] = i; }
        for (int a = 0; a < k; a++)
        for (int b = a+1; b < k; b++)
        for (int c = b+1; c < k; c++){
            uint32_t T3 = (1u<<lst[a])|(1u<<lst[b])|(1u<<lst[c]);
            uint32_t cn3 = adj[lst[a]] & adj[lst[b]] & adj[lst[c]];
            if (POP(cn2 & cn3 & ~T3 & ~T2) >= 3) return 1;
        }
    }
    return 0;
}

int main(int argc, char **argv){
    if (argc < 2){ fprintf(stderr, "usage: filter_mtf <file.mc> [expected_n]\n"); return 1; }
    FILE *f = strcmp(argv[1], "-") ? fopen(argv[1], "rb") : stdin;
    if (!f){ perror("open"); return 1; }
    int expn = argc > 2 ? atoi(argv[2]) : -1;
    long total = 0, kept = 0, k7 = 0, k333 = 0, badver = 0;
    long omega_hist[8] = {0};
    for (;;){
        int c = fgetc(f);
        if (c == EOF) break;
        int n = c;
        if (n < 3 || n >= MAXV){ fprintf(stderr, "bad n=%d at graph %ld\n", n, total); return 1; }
        if (expn > 0 && n != expn){ fprintf(stderr, "unexpected n=%d\n", n); return 1; }
        uint32_t H[MAXV]; memset(H, 0, sizeof H);
        for (int x = 1; x <= n - 1; x++){
            for (;;){
                int y = fgetc(f);
                if (y == EOF){ fprintf(stderr, "EOF mid-graph\n"); return 1; }
                if (y == 0) break;
                int a = x - 1, b = y - 1;
                H[a] |= 1u << b; H[b] |= 1u << a;
            }
        }
        total++;
        /* verify H really is maximal triangle-free */
        if (!is_triangle_free(H, n) || !is_maximal_tf(H, n)){ badver++; continue; }
        /* G = complement */
        uint32_t G[MAXV];
        uint32_t full = (1u << n) - 1;
        for (int i = 0; i < n; i++) G[i] = full & ~H[i] & ~(1u << i);
        /* omega(G) <= 6, i.e. no K7 */
        if (has_clique(G, full, 0, 7)){ k7++; continue; }
        int om = 6;
        while (om > 1 && !has_clique(G, full, 0, om)) om--;
        if (has_k333(G, n)){ k333++; omega_hist[7]++; continue; }
        omega_hist[om]++;
        kept++;
        printf("%d", n);
        for (int i = 0; i < n; i++) printf(" %u", G[i]);
        printf("\n");
    }
    if (f != stdin) fclose(f);
    fprintf(stderr, "total %ld  kept %ld  rejected: K7 %ld, K333 %ld, badver %ld\n",
            total, kept, k7, k333, badver);
    fprintf(stderr, "omega histogram of kept (5=K5 max, 6=K6): ");
    for (int i = 4; i <= 6; i++) fprintf(stderr, "omega=%d: %ld  ", i, omega_hist[i]);
    fprintf(stderr, "\n");
    if (badver){ fprintf(stderr, "VERIFICATION FAILURE\n"); return 1; }
    return 0;
}
