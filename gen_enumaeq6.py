"""Generate enumaeq6.c from enumaeq5.c with the d=6 rule changes."""
src = open('enumaeq5.c', encoding='utf-8').read()

src = src.replace('abstract almost-equidistant graphs in R^5',
                  'abstract almost-equidistant graphs in R^6')
src = src.replace('(2) no K_7,', '(2) no K_8,')
src = src.replace('(3) no K_{3,3,3} as a (not necessarily induced) subgraph.',
                  '(3) no K_{1,3,3,3} as a (not necessarily induced) subgraph.')
src = src.replace(''' * Expected totals (BPSSV Table 3, d=5), n = 4..13:
 *   7, 14, 38, 106, 402, 1817, 11132, 86053, 803299, 7623096
 * Expected minimal counts (BPSSV Table 2, d=5), n = 10..13: 25, 46, 106, 242''',
''' * Expected totals (BPSSV Table 3, d=6), n = 4..13:
 *   7, 14, 38, 107, 409, 1888, 12064, 103333, 1217849, 19170728
 * Expected minimal counts (BPSSV Table 2, d=6), n = 10..13: 29, 54, 130, 339''')
src = src.replace('cc -O2 -o enumaeq5 enumaeq5.c', 'cc -O2 -o enumaeq6 enumaeq6.c')
src = src.replace('enumaeq5 [maxn]', 'enumaeq6 [maxn]')

src = src.replace('''static const long EXPECTED_ALL[15] = { /* index n, d=5 totals (Table 3) */
    0,0,0,0, 7, 14, 38, 106, 402, 1817, 11132, 86053, 803299, 7623096, -1 };
static const long EXPECTED_MIN[15] = { /* minimal (Table 2, d=5) */
    0,0,0,0, 2, 3, 4, 6, 9, 14, 25, 46, 106, 242, 653 };''',
'''static const long EXPECTED_ALL[15] = { /* index n, d=6 totals (Table 3) */
    0,0,0,0, 7, 14, 38, 107, 409, 1888, 12064, 103333, 1217849, 19170728, -1 };
static const long EXPECTED_MIN[15] = { /* minimal (Table 2, d=6) */
    0,0,0,0, 2, 3, 4, 6, 10, 15, 29, 54, 130, 339, 1052 };''')

src = src.replace('int lo = m - 7; if (lo < 0) lo = 0;',
                  'int lo = m - 8; if (lo < 0) lo = 0;')

src = src.replace('''            static uint32_t cl6[200000]; int ncl6 = 0;
            if (m >= 7)
                collect_k(&parent, (1u << n) - 1, 0, 0, 6, cl6, &ncl6, 200000);''',
'''            static uint32_t cl7[400000]; int ncl7 = 0;
            if (m >= 8)
                collect_k(&parent, (1u << n) - 1, 0, 0, 7, cl7, &ncl7, 400000);''')
src = src.replace('''                            int okK = 1;
                            for (int q = 0; q < ncl6; q++)
                                if ((cl6[q] & S) == cl6[q]){ okK = 0; break; }''',
'''                            int okK = 1;
                            for (int q = 0; q < ncl7; q++)
                                if ((cl7[q] & S) == cl7[q]){ okK = 0; break; }''')

src = src.replace('''static int has_k7(const G *g, int n){
    uint32_t full = (1u << n) - 1;
    return clique_ext(g, full, 0, 7);
}''',
'''static int has_k8(const G *g, int n){
    uint32_t full = (1u << n) - 1;
    return clique_ext(g, full, 0, 8);
}''')

old_k333 = '''/* K_{1,3,3,3} subgraph: disjoint triples T1,T2,T3, all 27 cross edges */
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
}'''
new_k1333 = '''/* K_{1,3,3,3} subgraph: disjoint triples T1,T2,T3 with all 27 cross
   edges plus an apex adjacent to all nine */
static int has_k1333(const G *g, int n){
    int t2[3];
    for (t2[0] = 0; t2[0] < n; t2[0]++)
    for (t2[1] = t2[0]+1; t2[1] < n; t2[1]++)
    for (t2[2] = t2[1]+1; t2[2] < n; t2[2]++){
        uint32_t T2 = (1u<<t2[0])|(1u<<t2[1])|(1u<<t2[2]);
        uint32_t cn2 = g->adj[t2[0]] & g->adj[t2[1]] & g->adj[t2[2]] & ~T2;
        if (POP(cn2) < 7) continue;
        int lst[MAXV], k = 0; uint32_t s = cn2;
        while (s){ int i = CTZ(s); s &= s-1; lst[k++] = i; }
        for (int a = 0; a < k; a++)
        for (int b = a+1; b < k; b++)
        for (int c = b+1; c < k; c++){
            uint32_t T3 = (1u<<lst[a])|(1u<<lst[b])|(1u<<lst[c]);
            uint32_t cn3 = g->adj[lst[a]] & g->adj[lst[b]] & g->adj[lst[c]];
            uint32_t B = cn2 & cn3 & ~T3 & ~T2;
            if (POP(B) < 4) continue;
            uint32_t bb = B;
            while (bb){
                int ap = CTZ(bb); bb &= bb - 1;
                if (POP(B & g->adj[ap] & ~(1u << ap)) >= 3) return 1;
            }
        }
    }
    return 0;
}'''
# note: source still names it has_k333 with the K_{3,3,3} comment; match on
# the function as it exists in enumaeq5.c
old_k333 = old_k333.replace('K_{1,3,3,3} subgraph: disjoint triples T1,T2,T3, all 27 cross edges',
                            'K_{3,3,3} subgraph: disjoint triples T1,T2,T3, all 27 cross edges')
assert old_k333 in src, "k333 block not found"
src = src.replace(old_k333, new_k1333)

old_inc = '''                                /* K333 through v: T2,T3 in S w/ complete bipartite,
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
                                if (!bad){'''
new_inc = '''                                /* K1333 involving v: parent is
                                   K1333-free, so any occurrence in the
                                   extended graph passes through v; run
                                   the full checker on the extension */
                                G g2 = parent;
                                g2.adj[n] = S;
                                uint32_t tt = S;
                                while (tt){ int i = CTZ(tt); tt &= tt-1; g2.adj[i] |= 1u << n; }
                                int bad = has_k1333(&g2, m);
                                if (!bad){'''
assert old_inc in src, "incremental block not found"
src = src.replace(old_inc, new_inc)

src = src.replace('if (has_indep3(&pool[i], n) || has_k7(&pool[i], n) || has_k333(&pool[i], n)) bad++;',
                  'if (has_indep3(&pool[i], n) || has_k8(&pool[i], n) || has_k1333(&pool[i], n)) bad++;')
src = src.replace('snprintf(fn, sizeof fn, "aeq_d5_n%d.bin", n);',
                  'snprintf(fn, sizeof fn, "aeq_d6_n%d.bin", n);')

open('enumaeq6.c', 'w', encoding='utf-8', newline='\n').write(src)
print("enumaeq6.c written")
