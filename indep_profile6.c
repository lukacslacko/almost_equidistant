/* indep_profile6.c — INDEPENDENT reimplementation of the codex first-
 * generation exact graph filters for the d=6 level-19 corpus, written
 * from d6_theory_filters.md alone (deliberately not from profile_d6.c),
 * as the owed cross-check before their kills enter the certified ledger.
 *
 * Rules implemented (a graph is rejected by a rule if ANY seed clique of
 * the required size makes the rule's necessary condition fail):
 *   1. link bounds: common unit neighbours of a unit K_m form an almost-
 *      equidistant set in a sphere of dimension 7-m => bounds
 *      K2:16 K3:12 K4:10 K5:4 K6:2 K7:0.
 *   2. K7 facet-reflection: an outside vertex adjacent to exactly six of
 *      a K7 is the reflection -(4/3)q_i; two of the same omitted type
 *      coincide; two of different types have distance^2 16/9, so a
 *      required edge between them is impossible.
 *   3. K7 defect-support CSP (items 2-4 of the discrete layer; only
 *      vertices with |D_x| <= 2 branch).
 *   4. K7 disjoint-edge bounded cover: L-edges = required edges with
 *      disjoint seed-defect masks; Z = {c_x = 0} must be an L-vertex-
 *      cover, only |D_x| >= 3 vertices are eligible, |Z| <= 7.
 *   5. K7 tight-cover matching: if the minimum eligible cover is exactly
 *      7, some size-7 eligible cover must have defect masks with a
 *      perfect matching into the seven seed coordinates.
 *   6. K6 two-light-ray CSP (K6-only graphs): exists Z0 (eligible
 *      |D_x|>=3, inside initially-non-bipartite components of L, |Z0|<=6,
 *      masks matchable) and a 2-colouring of the non-bipartite components
 *      of L-Z0 such that each colour's masks plus Z0's masks match
 *      injectively into the six seed coordinates.
 *
 * Output: aggregate tallies to stderr (to be compared against
 * d6_profile.json schema 3) and one verdict byte per graph to a file:
 * bit0..bit5 = rules above, bit6 = clique number 7 (else 6), bit7 = any.
 * Usage: indep_profile6 <corpus> <verdicts.bin> [start [end]]
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

typedef uint32_t u32;
#define MAXV 24

static inline int POP(u32 x){ return __builtin_popcount(x); }
static inline int CTZ(u32 x){ return __builtin_ctz(x); }

static int N;                 /* vertices (19) */
static u32 A[MAXV];           /* adjacency */
static u32 FULL;

/* ---------------- clique enumeration ---------------- */
static int cliques[4096][8], ncl;   /* vertex lists of current size */
static void cliq_rec(u32 cand, int *cur, int size, int want){
    if (size == want){
        if (ncl < 4096){
            memcpy(cliques[ncl], cur, want * sizeof(int));
            ncl++;
        }
        return;
    }
    while (cand){
        int v = CTZ(cand); cand &= cand - 1;
        cur[size] = v;
        cliq_rec(cand & A[v], cur, size + 1, want);
    }
}
static int find_cliques(int want){
    int cur[8];
    ncl = 0;
    cliq_rec(FULL, cur, 0, want);
    return ncl;
}
static int has_clique_size(int want){
    /* existence only */
    int cur[8];
    ncl = 4096;              /* suppress storing */
    /* quick recursive existence */
    ncl = 0;
    cliq_rec(FULL, cur, 0, want);   /* stores up to cap; fine */
    return ncl > 0;
}

/* ---------------- rule 1: link bounds ---------------- */
static const int LINKB[8] = {0,0,16,12,10,4,2,0};  /* index by clique size */
static int rule_link(void){
    for (int m = 2; m <= 7; m++){
        if (find_cliques(m) == 0) continue;
        if (ncl >= 4096){ fprintf(stderr, "clique overflow m=%d\n", m); exit(3); }
        for (int c = 0; c < ncl; c++){
            u32 cn = FULL;
            for (int i = 0; i < m; i++) cn &= A[cliques[c][i]];
            if (POP(cn) > LINKB[m]) return 1;
        }
    }
    return 0;
}

/* ---------------- K7 machinery ---------------- */
static int nk7;
static int k7s[4096][8];
static u32 Dm[MAXV];          /* defect mask (7 bits) per outside vertex */
static int outv[MAXV], nout;

static void setup_seed7(const int *Q){
    u32 qm = 0;
    for (int i = 0; i < 7; i++) qm |= 1u << Q[i];
    nout = 0;
    for (int v = 0; v < N; v++){
        if (qm & (1u << v)) continue;
        u32 d = 0;
        for (int i = 0; i < 7; i++)
            if (!((A[Q[i]] >> v) & 1)) d |= 1u << i;
        Dm[nout] = d;
        outv[nout] = v;
        nout++;
    }
}

static int rule_reflection_seed(void){
    /* type i = outside vertex with defect exactly {i} */
    int type_of[MAXV];
    for (int x = 0; x < nout; x++)
        type_of[x] = (POP(Dm[x]) == 1) ? CTZ(Dm[x]) : -1;
    for (int x = 0; x < nout; x++){
        if (type_of[x] < 0) continue;
        for (int y = x + 1; y < nout; y++){
            if (type_of[y] < 0) continue;
            if (type_of[x] == type_of[y]) return 1;      /* coincide */
            if ((A[outv[x]] >> outv[y]) & 1) return 1;   /* dist^2 16/9 */
        }
    }
    return 0;
}

/* rule 3: defect-support CSP, branching on |D|<=2 vertices */
static int csp_bx[MAXV], csp_nb;
static u32 csp_sup[MAXV];
static int csp_ok(int k){
    if (k == csp_nb) return 1;
    int x = csp_bx[k];
    u32 D = Dm[x];
    /* iterate nonempty subsets of D */
    for (u32 S = D;; S = (S - 1) & D){
        if (S){
            csp_sup[x] = S;
            /* check against previously assigned */
            int good = 1;
            for (int j = 0; j < k && good; j++){
                int y = csp_bx[j];
                u32 T = csp_sup[y];
                if (POP(S) == 1 && POP(T) == 1){
                    if (S == T) good = 0;                    /* item 2 */
                    else if ((A[outv[x]] >> outv[y]) & 1) good = 0;
                }
            }
            if (good && POP(S) == 2){
                /* item 4: at most two vertices share a fixed 2-support */
                int tot = 1;
                for (int l = 0; l < k; l++)
                    if (csp_sup[csp_bx[l]] == S) tot++;
                if (tot > 2) good = 0;
            }
            if (good && csp_ok(k + 1)) return 1;
        }
        if (S == 0) break;
    }
    return 0;
}
static int rule_defect_csp_seed(void){
    csp_nb = 0;
    for (int x = 0; x < nout; x++)
        if (POP(Dm[x]) <= 2) csp_bx[csp_nb++] = x;
    /* item 3 (<=2 singleton supports) follows from item 2; not separate */
    return csp_ok(0) ? 0 : 1;
}

/* rules 4+5: bounded cover and tight-cover matching */
static u32 Ledge[MAXV];       /* L adjacency over outside indices */
static u32 eligible;          /* bitmask over outside indices */

static void build_L(void){
    eligible = 0;
    for (int x = 0; x < nout; x++){
        Ledge[x] = 0;
        if (POP(Dm[x]) >= 3) eligible |= 1u << x;
    }
    for (int x = 0; x < nout; x++)
        for (int y = x + 1; y < nout; y++)
            if (((A[outv[x]] >> outv[y]) & 1) && !(Dm[x] & Dm[y])){
                Ledge[x] |= 1u << y;
                Ledge[y] |= 1u << x;
            }
}
/* min eligible vertex cover via edge branching; returns >7 as 8 */
static int cover_rec(u32 covered, int used, int best){
    if (used >= best) return best;
    /* find an uncovered L-edge */
    int fx = -1, fy = -1;
    for (int x = 0; x < nout && fx < 0; x++){
        if ((covered >> x) & 1) continue;
        u32 m = Ledge[x] & ~covered;
        if (m){ fx = x; fy = CTZ(m); }
    }
    if (fx < 0) return used;
    int res = best;
    if ((eligible >> fx) & 1){
        int r = cover_rec(covered | (1u << fx), used + 1, res);
        if (r < res) res = r;
    }
    if ((eligible >> fy) & 1){
        int r = cover_rec(covered | (1u << fy), used + 1, res);
        if (r < res) res = r;
    }
    return res;
}
static int zcap = 7;      /* --zcap3: sharpened |Z|<=3 bound (theory doc,
                             n=19 twelve-outside case only) */
static int mincover(void){
    return cover_rec(0, 0, 8);
}
/* matching masks into 7 coordinates: standard augmenting paths */
static int mat_to[8];
static u32 mat_masks[16]; static int mat_n;
static int mat_seen[8];
static int aug(int x){
    u32 m = mat_masks[x];
    while (m){
        int c = CTZ(m); m &= m - 1;
        if (mat_seen[c]) continue;
        mat_seen[c] = 1;
        if (mat_to[c] < 0 || aug(mat_to[c])){
            mat_to[c] = x;
            return 1;
        }
    }
    return 0;
}
static int perfect_match(u32 *masks, int k, int ncoord){
    mat_n = k;
    for (int i = 0; i < k; i++) mat_masks[i] = masks[i];
    for (int c = 0; c < ncoord; c++) mat_to[c] = -1;
    for (int x = 0; x < k; x++){
        for (int c = 0; c < ncoord; c++) mat_seen[c] = 0;
        if (!aug(x)) return 0;
    }
    return 1;
}
/* enumerate ALL eligible size-7 covers directly (C(|eligible|,7) <= 792)
   and test the perfect-matching refinement */
static int any_cover7_matches(void){
    int el[MAXV], en = 0;
    u32 t = eligible;
    while (t){ el[en++] = CTZ(t); t &= t - 1; }
    if (en < 7) return 1;   /* mincover()==7 implies en>=7; defensive */
    int idx[8];
    for (int i = 0; i < 7; i++) idx[i] = i;
    for (;;){
        u32 ch = 0;
        for (int i = 0; i < 7; i++) ch |= 1u << el[idx[i]];
        /* is ch a vertex cover of L? */
        int isc = 1;
        for (int x = 0; x < nout && isc; x++){
            if ((ch >> x) & 1) continue;
            if (Ledge[x] & ~ch) isc = 0;
        }
        if (isc){
            u32 masks[8];
            for (int i = 0; i < 7; i++) masks[i] = Dm[el[idx[i]]];
            if (perfect_match(masks, 7, 7)) return 1;
        }
        int i = 6;
        while (i >= 0 && idx[i] == en - 7 + i) i--;
        if (i < 0) break;
        idx[i]++;
        for (int j = i + 1; j < 7; j++) idx[j] = idx[j-1] + 1;
    }
    return 0;
}
static int rule_cover_seed(int *tight_reject){
    build_L();
    int mc = mincover();
    *tight_reject = 0;
    if (mc > 7) return 1;
    if (mc == 7 && !any_cover7_matches()) *tight_reject = 1;
    if (zcap == 3 && nout == 12 && mc > 3) return 1;
    return 0;
}
/* ---------------- rule 6: K6 two-light-ray CSP ---------------- */
static int nk6out;
static int outv6[MAXV];
static u32 Dm6[MAXV];
static u32 L6[MAXV];

static void setup_seed6(const int *Q){
    u32 qm = 0;
    for (int i = 0; i < 6; i++) qm |= 1u << Q[i];
    nk6out = 0;
    for (int v = 0; v < N; v++){
        if (qm & (1u << v)) continue;
        u32 d = 0;
        for (int i = 0; i < 6; i++)
            if (!((A[Q[i]] >> v) & 1)) d |= 1u << i;
        Dm6[nk6out] = d;
        outv6[nk6out] = v;
        nk6out++;
    }
    for (int x = 0; x < nk6out; x++){
        L6[x] = 0;
    }
    for (int x = 0; x < nk6out; x++)
        for (int y = x + 1; y < nk6out; y++)
            if (((A[outv6[x]] >> outv6[y]) & 1) && !(Dm6[x] & Dm6[y])){
                L6[x] |= 1u << y;
                L6[y] |= 1u << x;
            }
}
/* connected components of L6 minus a removed set; bipartite check */
static int comp_id[MAXV], comp_bip[MAXV], ncomp;
static void components(u32 removed){
    int color[MAXV];
    for (int x = 0; x < nk6out; x++){ comp_id[x] = -1; color[x] = 0; }
    ncomp = 0;
    for (int s = 0; s < nk6out; s++){
        if (comp_id[s] >= 0 || ((removed >> s) & 1)) continue;
        /* BFS */
        int queue[MAXV], qh = 0, qt = 0;
        queue[qt++] = s; comp_id[s] = ncomp; color[s] = 0;
        int bip = 1;
        while (qh < qt){
            int x = queue[qh++];
            u32 m = L6[x] & ~removed;
            while (m){
                int y = CTZ(m); m &= m - 1;
                if (comp_id[y] < 0){
                    comp_id[y] = ncomp; color[y] = color[x] ^ 1;
                    queue[qt++] = y;
                } else if (color[y] == color[x]) bip = 0;
            }
        }
        comp_bip[ncomp] = bip;
        ncomp++;
    }
}
static int lightray_feasible_Z0(u32 Z0){
    /* Z0 masks must match injectively (assertion, kept) */
    u32 zmasks[8]; int zn = 0;
    u32 t = Z0;
    while (t){ int x = CTZ(t); t &= t - 1;
        if (zn >= 6) return 0;
        zmasks[zn++] = Dm6[x]; }
    if (zn > 6) return 0;
    if (!perfect_match(zmasks, zn, 6)) return 0;
    components(Z0);
    /* collect non-bipartite components */
    int nb[MAXV], nnb = 0;
    for (int c = 0; c < ncomp; c++)
        if (!comp_bip[c]) nb[nnb++] = c;
    if (nnb == 0) return 1;
    if (nnb > 12) return 0;   /* cannot happen: 13 outside */
    /* 2-colourings of the non-bipartite components */
    for (u32 col = 0; col < (1u << nnb); col++){
        int ok = 1;
        for (int side = 0; side < 2 && ok; side++){
            u32 masks[16]; int k = 0;
            for (int i = 0; i < zn; i++) masks[k++] = zmasks[i];
            for (int i = 0; i < nnb && ok; i++){
                if (((col >> i) & 1) != (u32)side) continue;
                for (int x = 0; x < nk6out; x++){
                    if (((Z0 >> x) & 1)) continue;
                    if (comp_id[x] == nb[i]){
                        if (k >= 16){ ok = 0; break; }
                        masks[k++] = Dm6[x];
                    }
                }
            }
            if (ok && k > 6) ok = 0;
            if (ok && !perfect_match(masks, k, 6)) ok = 0;
        }
        if (ok) return 1;
    }
    return 0;
}
static int rule_lightray_seed(void){
    /* eligible Z0 vertices: |D|>=3 and inside an initially non-bipartite
       component (doc: bipartite-component vertices never need Z0) */
    components(0);
    u32 elig = 0;
    for (int x = 0; x < nk6out; x++)
        if (POP(Dm6[x]) >= 3 && !comp_bip[comp_id[x]]) elig |= 1u << x;
    /* try Z0 = empty first, then all eligible subsets of size <= 6 */
    if (lightray_feasible_Z0(0)) return 0;
    int el[MAXV], en = 0;
    u32 t = elig;
    while (t){ el[en++] = CTZ(t); t &= t - 1; }
    for (int k = 1; k <= 6 && k <= en; k++){
        int idx[8];
        for (int i = 0; i < k; i++) idx[i] = i;
        for (;;){
            u32 Z0 = 0;
            for (int i = 0; i < k; i++) Z0 |= 1u << el[idx[i]];
            if (lightray_feasible_Z0(Z0)) return 0;
            int i = k - 1;
            while (i >= 0 && idx[i] == en - k + i) i--;
            if (i < 0) break;
            idx[i]++;
            for (int j = i + 1; j < k; j++) idx[j] = idx[j-1] + 1;
        }
    }
    return 1;   /* every choice failed: reject */
}

/* ------------------------------ main ------------------------------ */
int main(int argc, char **argv){
    if (argc < 3){
        fprintf(stderr, "usage: indep_profile6 <corpus> <verdicts.bin> [start [end]]\n");
        return 1;
    }
    FILE *f = fopen(argv[1], "r");
    if (!f){ perror(argv[1]); return 1; }
    FILE *vf = fopen(argv[2], "wb");
    if (!vf){ perror(argv[2]); return 1; }
    long start = 0, end = -1;
    {
        int ai = 3;
        if (ai < argc && !strcmp(argv[ai], "--zcap3")){ zcap = 3; ai++; }
        if (ai < argc) start = atol(argv[ai++]);
        if (ai < argc) end = atol(argv[ai++]);
    }

    long tally[8] = {0};   /* link, refl, csp, cover, tight, lightray */
    long total = 0, cl7 = 0, cl6 = 0, any7 = 0, any6 = 0, anyall = 0;
    long line = 0;
    char buf[512];
    while (fgets(buf, sizeof buf, f)){
        if (line < start){ line++; continue; }
        if (end >= 0 && line >= end) break;
        char *p = buf;
        N = (int)strtol(p, &p, 10);
        if (N < 8 || N > MAXV){ fprintf(stderr, "bad n at line %ld\n", line); return 1; }
        for (int i = 0; i < N; i++) A[i] = (u32)strtoul(p, &p, 10);
        FULL = (1u << N) - 1;
        /* validate: symmetric, loopless, alpha<=2 */
        for (int i = 0; i < N; i++){
            if ((A[i] >> i) & 1){ fprintf(stderr, "loop line %ld\n", line); return 1; }
            for (int j = 0; j < N; j++)
                if (((A[i] >> j) & 1) != ((A[j] >> i) & 1)){
                    fprintf(stderr, "asym line %ld\n", line); return 1; }
        }
        for (int i = 0; i < N; i++)
            for (int j = i + 1; j < N; j++){
                if ((A[i] >> j) & 1) continue;
                u32 nonnbr = FULL & ~A[i] & ~A[j] & ~(1u << i) & ~(1u << j);
                u32 m = nonnbr;
                int bad = 0;
                while (m){
                    int v = CTZ(m); m &= m - 1;
                    if (!((A[i] >> v) & 1) && !((A[j] >> v) & 1)){ bad = 1; break; }
                }
                if (bad){ fprintf(stderr, "alpha>2 line %ld\n", line); return 1; }
            }

        int is7 = 0;
        int nk7l = find_cliques(7);
        if (nk7l > 0) is7 = 1;
        if (ncl >= 4096){ fprintf(stderr, "K7 overflow line %ld\n", line); return 1; }
        int k7list[4096][8]; int nk7c = 0;
        if (is7){
            nk7c = ncl;
            memcpy(k7list, cliques, sizeof(int) * 8 * (size_t)ncl);
        }
        unsigned char verdict = 0;
        if (is7) verdict |= 1u << 6;

        if (rule_link()) verdict |= 1u << 0;

        if (is7){
            int rref = 0, rcsp = 0, rcov = 0, rtig = 0;
            for (int c = 0; c < nk7c && (!rref || !rcsp || !rcov || !rtig); c++){
                int Q[8];
                memcpy(Q, k7list[c], sizeof(int) * 7);
                setup_seed7(Q);
                if (!rref && rule_reflection_seed()) rref = 1;
                if (!rcsp && rule_defect_csp_seed()) rcsp = 1;
                if (!rcov || !rtig){
                    int tr = 0;
                    int cv = rule_cover_seed(&tr);
                    if (cv) rcov = 1;
                    if (tr) rtig = 1;
                }
            }
            if (rref) verdict |= 1u << 1;
            if (rcsp) verdict |= 1u << 2;
            if (rcov) verdict |= 1u << 3;
            if (rtig) verdict |= 1u << 4;
        } else {
            int nk6 = find_cliques(6);
            if (ncl >= 4096){ fprintf(stderr, "K6 overflow line %ld\n", line); return 1; }
            int rlr = 0;
            int k6list[4096][8]; int nk6c = nk6;
            memcpy(k6list, cliques, sizeof(int) * 8 * (size_t)ncl);
            for (int c = 0; c < nk6c && !rlr; c++){
                setup_seed6(k6list[c]);
                if (rule_lightray_seed()) rlr = 1;
            }
            if (rlr) verdict |= 1u << 5;
        }
        int any = (verdict & 0x3f) != 0;
        if (any) verdict |= 1u << 7;
        fputc(verdict, vf);
        total++;
        if (is7){ cl7++; if (any) any7++; } else { cl6++; if (any) any6++; }
        if (any) anyall++;
        for (int b = 0; b < 6; b++) if ((verdict >> b) & 1) tally[b]++;
        line++;
        if (total % 100000 == 0)
            fprintf(stderr, "... %ld processed (cl7 %ld cl6 %ld) rejected %ld\n",
                    total, cl7, cl6, anyall);
    }
    fclose(f); fclose(vf);
    fprintf(stderr,
        "INDEP PROFILE: total %ld  clique7 %ld  clique6 %ld\n"
        "  link %ld\n  K7 reflection %ld\n  K7 defect-CSP %ld\n"
        "  K7 bounded-cover %ld\n  K7 tight-cover-matching %ld\n"
        "  K6 two-light-ray %ld\n"
        "  union rejected %ld (K7-class %ld, K6-class %ld)\n"
        "  residue %ld (K7-class %ld, K6-class %ld)\n",
        total, cl7, cl6, tally[0], tally[1], tally[2], tally[3], tally[4],
        tally[5], anyall, any7, any6,
        total - anyall, cl7 - any7, cl6 - any6);
    return 0;
}
