/* Exact graph-only profiler for the d=6, n=19 candidate corpus.
 *
 * Input is the text format written by filter_mtf6.c:
 *     19 adj[0] ... adj[18]
 * one graph per line.  The optional second argument is the uncompressed
 * append-only kill log ("index decomposition" per line).  Candidate indices
 * are zero based, exactly as in reproduce6.py.
 *
 * Reported rejection rules are necessary conditions for realizability in R^6:
 * for every unit K_m, the number of common unit neighbours is at most
 *   m=2:16, m=3:12, m=4:10, m=5:4, m=6:2, m=7:0.
 * See d6_theory_filters.md for the proof.  This program uses only integer
 * bit operations; there is no numerical geometry.
 *
 * Build: cc -O3 -o profile_d6 profile_d6.c
 * Usage: ./profile_d6 aeq_d6_n19.txt [killed_d6_n19.log]
 */
#include <errno.h>
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define N 19
#define EXPECTED 3971787u

typedef struct {
    uint64_t total;
    uint64_t omega[9];
    uint64_t link[8];
    uint64_t cumulative;
    uint64_t k7_reflection;
    uint64_t k7_defect_csp;
    uint64_t cumulative_exact;
    uint64_t cumulative_exact_by_omega[9];
} Stats;

static inline int pop(uint32_t x) { return __builtin_popcount(x); }
static inline int ctz(uint32_t x) { return __builtin_ctz(x); }

static int has_clique_rec(const uint32_t adj[N], uint32_t cand, int need) {
    if (need == 0) return 1;
    if (pop(cand) < need) return 0;
    while (cand) {
        int v = ctz(cand);
        cand &= cand - 1;
        if (has_clique_rec(adj, cand & adj[v], need - 1)) return 1;
        if (pop(cand) < need) return 0;
    }
    return 0;
}

/* Find a K_want whose full common neighbourhood is larger than bound.
 * chosen_common is the intersection of the neighbourhoods of vertices
 * already chosen; cand imposes increasing labels and avoids duplicates. */
static int bad_link_rec(const uint32_t adj[N], uint32_t cand,
                        uint32_t chosen_common, int depth,
                        int want, int bound) {
    if (depth == want) return pop(chosen_common) > bound;
    if (pop(cand) < want - depth) return 0;
    while (cand) {
        int v = ctz(cand);
        cand &= cand - 1;
        if (bad_link_rec(adj, cand & adj[v], chosen_common & adj[v],
                         depth + 1, want, bound)) return 1;
        if (pop(cand) < want - depth) return 0;
    }
    return 0;
}

static int bad_link(const uint32_t adj[N], int m, int bound) {
    return bad_link_rec(adj, (1u << N) - 1, (1u << N) - 1,
                        0, m, bound);
}

/* A K7 is a full-dimensional regular simplex.  An outside vertex joined to
 * exactly six seed vertices is forced to be the nonduplicate reflection of
 * the omitted simplex vertex in the opposite K6 facet.  Reflections of two
 * different seed vertices are not at unit distance; two of the same type
 * coincide.  Either situation contradicts the corresponding required edge
 * or distinctness. */
static int bad_reflection_seed(const uint32_t adj[N], uint32_t seed) {
    int forced_vertex[7];
    for (int i = 0; i < 7; ++i) forced_vertex[i] = -1;
    uint32_t outside = ((1u << N) - 1) & ~seed;
    while (outside) {
        int v = ctz(outside);
        outside &= outside - 1;
        uint32_t defects = seed & ~adj[v];
        if (pop(defects) != 1) continue;
        int type_vertex = ctz(defects);
        int type = pop(seed & ((1u << type_vertex) - 1));
        if (forced_vertex[type] >= 0) return 1; /* forced collision */
        for (int j = 0; j < 7; ++j)
            if (forced_vertex[j] >= 0 &&
                    ((adj[v] >> forced_vertex[j]) & 1u))
                return 1; /* distinct reflections have non-unit distance */
        forced_vertex[type] = v;
    }
    return 0;
}

static int bad_k7_reflection_rec(const uint32_t adj[N], uint32_t cand,
                                 uint32_t seed, int depth) {
    if (depth == 7) return bad_reflection_seed(adj, seed);
    if (pop(cand) < 7 - depth) return 0;
    while (cand) {
        int v = ctz(cand);
        cand &= cand - 1;
        if (bad_k7_reflection_rec(adj, cand & adj[v], seed | (1u << v),
                                  depth + 1)) return 1;
        if (pop(cand) < 7 - depth) return 0;
    }
    return 0;
}

static int bad_k7_reflection(const uint32_t adj[N]) {
    return bad_k7_reflection_rec(adj, (1u << N) - 1, 0, 0);
}

typedef struct {
    int v[12];
    uint8_t allowed[12];
    int count;
} SmallDefects;

static int defect_dfs(const uint32_t adj[N], const SmallDefects *sd, int at,
                      uint8_t singleton_types, int singleton_count,
                      int singleton_vertices[2], uint8_t pair_count[128]) {
    if (at == sd->count) return 1;
    int v = sd->v[at];
    uint8_t a = sd->allowed[at];

    /* Actual support may be either singleton contained in the allowed pair. */
    uint8_t choices = a;
    while (choices) {
        uint8_t bit = choices & (uint8_t)(-choices);
        choices ^= bit;
        if (singleton_count >= 2 || (singleton_types & bit)) continue;
        int ok = 1;
        for (int j = 0; j < singleton_count; ++j)
            if ((adj[v] >> singleton_vertices[j]) & 1u) ok = 0;
        if (!ok) continue;
        singleton_vertices[singleton_count] = v;
        if (defect_dfs(adj, sd, at + 1, singleton_types | bit,
                       singleton_count + 1, singleton_vertices, pair_count))
            return 1;
    }

    /* Or the actual support is the full allowed pair.  Points with the same
     * two-defect support are pairwise unit and lie on a circle of radius
     * sqrt(3/5), hence there are at most two of them. */
    if (pop(a) == 2 && pair_count[a] < 2) {
        ++pair_count[a];
        if (defect_dfs(adj, sd, at + 1, singleton_types, singleton_count,
                       singleton_vertices, pair_count)) {
            --pair_count[a];
            return 1;
        }
        --pair_count[a];
    }
    return 0;
}

static int bad_defect_seed(const uint32_t adj[N], uint32_t seed) {
    int seed_pos[N];
    int k = 0;
    for (int v = 0; v < N; ++v)
        seed_pos[v] = (seed >> v) & 1u ? k++ : -1;
    SmallDefects sd = {{0}, {0}, 0};
    uint32_t outside = ((1u << N) - 1) & ~seed;
    while (outside) {
        int v = ctz(outside);
        outside &= outside - 1;
        uint32_t am = seed & ~adj[v];
        int sz = pop(am);
        if (sz == 0) return 1; /* K8, retained as a defensive check */
        if (sz > 2) continue;  /* choosing the whole support is CSP-safe */
        uint8_t a = 0;
        while (am) {
            int q = ctz(am);
            am &= am - 1;
            a |= (uint8_t)(1u << seed_pos[q]);
        }
        sd.v[sd.count] = v;
        sd.allowed[sd.count] = a;
        ++sd.count;
    }
    /* Forced singleton variables first, then pairs with fewest choices. */
    for (int i = 0; i < sd.count; ++i)
        for (int j = i + 1; j < sd.count; ++j)
            if (pop(sd.allowed[j]) < pop(sd.allowed[i])) {
                int tv = sd.v[i]; sd.v[i] = sd.v[j]; sd.v[j] = tv;
                uint8_t ta = sd.allowed[i];
                sd.allowed[i] = sd.allowed[j]; sd.allowed[j] = ta;
            }
    int singleton_vertices[2] = {-1, -1};
    uint8_t pair_count[128] = {0};
    return !defect_dfs(adj, &sd, 0, 0, 0, singleton_vertices, pair_count);
}

static int bad_k7_defect_rec(const uint32_t adj[N], uint32_t cand,
                             uint32_t seed, int depth) {
    if (depth == 7) return bad_defect_seed(adj, seed);
    if (pop(cand) < 7 - depth) return 0;
    while (cand) {
        int v = ctz(cand);
        cand &= cand - 1;
        if (bad_k7_defect_rec(adj, cand & adj[v], seed | (1u << v),
                              depth + 1)) return 1;
        if (pop(cand) < 7 - depth) return 0;
    }
    return 0;
}

static int bad_k7_defect(const uint32_t adj[N]) {
    return bad_k7_defect_rec(adj, (1u << N) - 1, 0, 0);
}

static int validate(const uint32_t adj[N]) {
    const uint32_t full = (1u << N) - 1;
    for (int i = 0; i < N; ++i) {
        if (adj[i] & ~full) return 0;
        if (adj[i] & (1u << i)) return 0;
        for (int j = i + 1; j < N; ++j)
            if (((adj[i] >> j) & 1u) != ((adj[j] >> i) & 1u)) return 0;
    }
    /* Complement triangle-free, equivalently alpha(G) <= 2. */
    for (int i = 0; i < N; ++i) {
        uint32_t non = full & ~adj[i] & ~(1u << i);
        while (non) {
            int j = ctz(non);
            non &= non - 1;
            if (non & ~adj[j] & ~(1u << j)) return 0;
        }
    }
    return 1;
}

static void profile_one(Stats *s, const uint32_t adj[N]) {
    static const int bounds[8] = {0, 0, 16, 12, 10, 4, 2, 0};
    int rejected = 0;
    ++s->total;
    int omega = has_clique_rec(adj, (1u << N) - 1, 7) ? 7 :
                (has_clique_rec(adj, (1u << N) - 1, 6) ? 6 : 5);
    ++s->omega[omega];
    for (int m = 2; m <= 7; ++m) {
        int bad = bad_link(adj, m, bounds[m]);
        if (bad) {
            ++s->link[m];
            rejected = 1;
        }
    }
    if (rejected) ++s->cumulative;
    int bad_refl = omega == 7 && bad_k7_reflection(adj);
    if (bad_refl) ++s->k7_reflection;
    int bad_csp = omega == 7 && bad_k7_defect(adj);
    if (bad_csp) ++s->k7_defect_csp;
    if (rejected || bad_refl || bad_csp) {
        ++s->cumulative_exact;
        ++s->cumulative_exact_by_omega[omega];
    }
}

static uint8_t *load_killed(const char *path, uint64_t *count) {
    FILE *f = fopen(path, "r");
    if (!f) {
        fprintf(stderr, "cannot open kill log %s: %s\n", path,
                strerror(errno));
        exit(2);
    }
    uint8_t *bits = calloc(EXPECTED, 1);
    if (!bits) { fprintf(stderr, "out of memory\n"); exit(2); }
    unsigned idx;
    int dec;
    *count = 0;
    while (fscanf(f, "%u %d", &idx, &dec) == 2) {
        if (idx >= EXPECTED) {
            fprintf(stderr, "bad kill-log index %u\n", idx);
            exit(2);
        }
        if (!bits[idx]) { bits[idx] = 1; ++*count; }
    }
    if (!feof(f)) {
        fprintf(stderr, "malformed kill log near record %" PRIu64 "\n", *count);
        exit(2);
    }
    fclose(f);
    return bits;
}

static void print_stats(const char *name, const Stats *s) {
    printf("    \"%s\": {\n", name);
    printf("      \"total\": %" PRIu64 ",\n", s->total);
    printf("      \"clique_number\": {\"5_or_less\": %" PRIu64
           ", \"6\": %" PRIu64 ", \"7\": %" PRIu64 "},\n",
           s->omega[5], s->omega[6], s->omega[7]);
    printf("      \"link_rejections\": {\"K2\": %" PRIu64
           ", \"K3\": %" PRIu64 ", \"K4\": %" PRIu64
           ", \"K5\": %" PRIu64 ", \"K6\": %" PRIu64
           ", \"K7\": %" PRIu64 "},\n",
           s->link[2], s->link[3], s->link[4], s->link[5],
           s->link[6], s->link[7]);
    printf("      \"cumulative_link_rejections\": %" PRIu64 ",\n",
           s->cumulative);
    printf("      \"K7_reflection_rejections\": %" PRIu64 ",\n",
           s->k7_reflection);
    printf("      \"K7_defect_CSP_rejections\": %" PRIu64 ",\n",
           s->k7_defect_csp);
    printf("      \"cumulative_exact_rejections\": %" PRIu64 ",\n",
           s->cumulative_exact);
    printf("      \"exact_rejections_by_clique_number\": {\"5_or_less\": %"
           PRIu64 ", \"6\": %" PRIu64 ", \"7\": %" PRIu64 "},\n",
           s->cumulative_exact_by_omega[5],
           s->cumulative_exact_by_omega[6],
           s->cumulative_exact_by_omega[7]);
    printf("      \"residue_by_clique_number\": {\"5_or_less\": %" PRIu64
           ", \"6\": %" PRIu64 ", \"7\": %" PRIu64 "}\n",
           s->omega[5] - s->cumulative_exact_by_omega[5],
           s->omega[6] - s->cumulative_exact_by_omega[6],
           s->omega[7] - s->cumulative_exact_by_omega[7]);
    printf("    }");
}

int main(int argc, char **argv) {
    if (argc < 2 || argc > 3) {
        fprintf(stderr, "usage: %s aeq_d6_n19.txt [killed_d6_n19.log]\n",
                argv[0]);
        return 2;
    }
    FILE *f = fopen(argv[1], "r");
    if (!f) { perror(argv[1]); return 2; }
    uint64_t nkilled = 0;
    uint8_t *killed = argc == 3 ? load_killed(argv[2], &nkilled) : NULL;
    Stats all = {0}, certified = {0}, deferred = {0};
    unsigned idx = 0;
    for (;;) {
        int n;
        int got = fscanf(f, "%d", &n);
        if (got == EOF) break;
        if (got != 1 || n != N) {
            fprintf(stderr, "bad candidate header at index %u\n", idx);
            return 2;
        }
        uint32_t adj[N];
        for (int i = 0; i < N; ++i)
            if (fscanf(f, "%u", &adj[i]) != 1) {
                fprintf(stderr, "truncated candidate at index %u\n", idx);
                return 2;
            }
        if (!validate(adj)) {
            fprintf(stderr, "invalid graph at index %u\n", idx);
            return 2;
        }
        profile_one(&all, adj);
        if (killed && killed[idx]) profile_one(&certified, adj);
        else if (killed) profile_one(&deferred, adj);
        ++idx;
        if ((idx % 500000u) == 0)
            fprintf(stderr, "profiled %u candidates\n", idx);
    }
    fclose(f);
    if (idx != EXPECTED) {
        fprintf(stderr, "candidate count %u, expected %u\n", idx, EXPECTED);
        return 2;
    }
    if (killed && nkilled != certified.total) {
        fprintf(stderr, "kill-log/profile mismatch: log=%" PRIu64
                " profile=%" PRIu64 "\n", nkilled, certified.total);
        return 2;
    }
    printf("{\n  \"schema\": 1,\n  \"n\": 19,\n  \"dimension\": 6,\n");
    printf("  \"candidate_count_expected\": %u,\n", EXPECTED);
    printf("  \"kill_log_unique_indices\": %" PRIu64 ",\n", nkilled);
    printf("  \"populations\": {\n");
    print_stats("all", &all);
    if (killed) {
        printf(",\n"); print_stats("certified", &certified);
        printf(",\n"); print_stats("deferred", &deferred);
    }
    printf("\n  }\n}\n");
    free(killed);
    return 0;
}
