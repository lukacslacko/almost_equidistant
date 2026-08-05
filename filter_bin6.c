/* Read a .bin level dump (hdr: n, count; then count*n uint32 masks) of
 * complement-triangle-free graphs and count how many are abstract
 * almost-equidistant in R^6 (no K_8, no K_{1,3,3,3}).
 * Usage: filter_bin6 tf_n11.bin */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

static inline int POP(uint32_t x){ return __builtin_popcount(x); }
static inline int CTZ(uint32_t x){ return __builtin_ctz(x); }
#define MAXV 15

static int has_clique(const uint32_t *adj, uint32_t cand, int size, int want){
    if (size == want) return 1;
    while (cand){
        int v = CTZ(cand); cand &= cand - 1;
        if (has_clique(adj, cand & adj[v], size + 1, want)) return 1;
    }
    return 0;
}
static int has_k1333(const uint32_t *adj, int n){
    int t2[3];
    for (t2[0] = 0; t2[0] < n; t2[0]++)
    for (t2[1] = t2[0]+1; t2[1] < n; t2[1]++)
    for (t2[2] = t2[1]+1; t2[2] < n; t2[2]++){
        uint32_t T2 = (1u<<t2[0])|(1u<<t2[1])|(1u<<t2[2]);
        uint32_t cn2 = adj[t2[0]] & adj[t2[1]] & adj[t2[2]] & ~T2;
        if (POP(cn2) < 7) continue;
        int lst[MAXV], k = 0; uint32_t s = cn2;
        while (s){ int i = CTZ(s); s &= s-1; lst[k++] = i; }
        for (int a = 0; a < k; a++)
        for (int b = a+1; b < k; b++)
        for (int c = b+1; c < k; c++){
            uint32_t T3 = (1u<<lst[a])|(1u<<lst[b])|(1u<<lst[c]);
            uint32_t cn3 = adj[lst[a]] & adj[lst[b]] & adj[lst[c]];
            uint32_t B = cn2 & cn3 & ~T3 & ~T2;
            if (POP(B) < 4) continue;
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
    FILE *f = fopen(argv[1], "rb");
    if (!f){ perror("open"); return 1; }
    uint32_t hdr[2];
    if (fread(hdr, 4, 2, f) != 2) return 1;
    int n = (int)hdr[0]; long cnt = (long)hdr[1];
    long kept = 0, k8 = 0, k13 = 0;
    for (long i = 0; i < cnt; i++){
        uint32_t adj[MAXV];
        if (fread(adj, 4, n, f) != (size_t)n){ fprintf(stderr, "short read\n"); return 1; }
        uint32_t full = (1u << n) - 1;
        if (has_clique(adj, full, 0, 8)){ k8++; continue; }
        if (has_k1333(adj, n)){ k13++; continue; }
        kept++;
    }
    printf("n=%d total=%ld kept=%ld rejected: K8=%ld K1333=%ld\n",
           n, cnt, kept, k8, k13);
    return 0;
}
