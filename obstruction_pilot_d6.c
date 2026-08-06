/* Heuristic obstruction-containment coverage pilot for d=6.
 *
 * This program does NOT certify any of the selected n=14 patterns as a
 * geometric obstruction.  It only proves ordinary (non-induced) graph
 * containment on a deterministic, stratified sample, in order to rank the
 * possible payoff of later rigorous certification.
 *
 * Inputs:
 *   - the complete 1052-graph n=14 corpus;
 *   - the LM triage output, from which column 4 > 1e-6 selects 89 patterns;
 *   - the complete 3971787-graph n=19 corpus;
 *   - the old certified-kill log.
 *
 * The pull-boundary strata reproduce the exact link/reflection union in
 * d6_profile.json at commit 9292f29be9cd26e265ed44f7473cf36bf931e68c.
 * The measured equality of the defect-CSP and reflection decisions on this
 * corpus is used only to name the pull-boundary "exact residue" strata; the
 * population counts are checked against that committed profile.
 *
 * Containment means edges of a pattern map to required edges of a target;
 * pattern nonedges impose no condition.  Search uses only sound degree,
 * arc-consistency, and all-different pruning.  Every reported hit is checked
 * edge by edge.  A per-pair node-limit exhaustion is TIMEOUT and is counted
 * conservatively as no hit.
 *
 * Build:
 *   cc -O3 -std=c11 -Wall -Wextra -pedantic -pthread \
 *      -o obstruction_pilot_d6.bin obstruction_pilot_d6.c -lm
 *
 * Usage:
 *   ./obstruction_pilot_d6.bin N14.txt LM.out N19.txt KILL.log OUT.json \
 *       [sample_per_stratum=160] [node_limit=50000] [threads=12] \
 *       [seed=20260806]
 */

#define _POSIX_C_SOURCE 200809L
#include <errno.h>
#include <inttypes.h>
#include <math.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define NP 14
#define NT 19
#define EXPECTED_TARGETS 3971787u
#define EXPECTED_KILLED 3055474u
#define EXPECTED_N14 1052
#define EXPECTED_SELECTED 89
#define MAX_PATTERNS 128
#define NSTRATA 5
#define FULL19 ((1u << NT) - 1u)

enum {
    CERTIFIED_EXACT_REJECTED = 0,
    CERTIFIED_EXACT_RESIDUE = 1,
    DEFERRED_K7_EXACT_REJECTED = 2,
    DEFERRED_K7_EXACT_RESIDUE = 3,
    DEFERRED_K6_EXACT_RESIDUE = 4
};

static const char *const stratum_names[NSTRATA] = {
    "certified_K7_preexisting_exact_rejected",
    "certified_K7_preexisting_exact_residue",
    "deferred_K7_preexisting_exact_rejected",
    "deferred_K7_preexisting_exact_residue",
    "deferred_K6_preexisting_exact_residue"
};

/* Counts in d6_profile.json at the review boundary. */
static const uint64_t expected_strata[NSTRATA] = {
    1004938, 2050536, 12466, 728028, 175819
};

typedef struct {
    uint32_t adj[NP];
    uint8_t deg[NP];
    int source_index;
    int edges;
    double lm_residual;
} Pattern;

typedef struct {
    uint32_t adj[NT];
    uint32_t corpus_index;
    uint8_t stratum;
} Target;

typedef struct {
    Target *items;
    size_t used;
    size_t cap;
    uint64_t seen;
    uint64_t rng;
} Reservoir;

typedef struct {
    uint64_t hit[2];
    uint64_t timeout[2];
    uint64_t nodes;
} TargetResult;

static inline int pop32(uint32_t x) { return __builtin_popcount(x); }
static inline int ctz32(uint32_t x) { return __builtin_ctz(x); }

static double monotime(void) {
    struct timespec ts;
    if (clock_gettime(CLOCK_MONOTONIC, &ts) != 0) {
        perror("clock_gettime");
        exit(2);
    }
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

static void *xcalloc(size_t n, size_t z) {
    void *p = calloc(n, z);
    if (!p) {
        fprintf(stderr, "out of memory allocating %zu bytes\n", n * z);
        exit(2);
    }
    return p;
}

static uint64_t splitmix64(uint64_t *state) {
    uint64_t z = (*state += UINT64_C(0x9e3779b97f4a7c15));
    z = (z ^ (z >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    z = (z ^ (z >> 27)) * UINT64_C(0x94d049bb133111eb);
    return z ^ (z >> 31);
}

/* Lemire's unbiased reduction of a 64-bit random integer to [0,bound). */
static uint64_t rand_bounded(uint64_t *state, uint64_t bound) {
    uint64_t x, low, threshold;
    __uint128_t product;
    if (bound == 0) return 0;
    x = splitmix64(state);
    product = (__uint128_t)x * (__uint128_t)bound;
    low = (uint64_t)product;
    if (low < bound) {
        threshold = (uint64_t)(-bound) % bound;
        while (low < threshold) {
            x = splitmix64(state);
            product = (__uint128_t)x * (__uint128_t)bound;
            low = (uint64_t)product;
        }
    }
    return (uint64_t)(product >> 64);
}

static int validate_graph(const uint32_t *adj, int n) {
    uint32_t full = (1u << n) - 1u;
    for (int i = 0; i < n; ++i) {
        if (adj[i] & ~full) return 0;
        if (adj[i] & (1u << i)) return 0;
        for (int j = i + 1; j < n; ++j)
            if (((adj[i] >> j) & 1u) != ((adj[j] >> i) & 1u))
                return 0;
    }
    return 1;
}

static int load_n14(const char *path, Pattern all[EXPECTED_N14]) {
    FILE *f = fopen(path, "r");
    int count = 0;
    if (!f) {
        fprintf(stderr, "cannot open %s: %s\n", path, strerror(errno));
        exit(2);
    }
    for (;;) {
        int n;
        int got = fscanf(f, "%d", &n);
        if (got == EOF) break;
        if (got != 1 || n != NP || count >= EXPECTED_N14) {
            fprintf(stderr, "bad n=14 corpus header at graph %d\n", count);
            exit(2);
        }
        memset(&all[count], 0, sizeof(all[count]));
        for (int i = 0; i < NP; ++i)
            if (fscanf(f, "%u", &all[count].adj[i]) != 1) {
                fprintf(stderr, "truncated n=14 graph %d\n", count);
                exit(2);
            }
        if (!validate_graph(all[count].adj, NP)) {
            fprintf(stderr, "invalid n=14 graph %d\n", count);
            exit(2);
        }
        all[count].source_index = count;
        for (int i = 0; i < NP; ++i) {
            all[count].deg[i] = (uint8_t)pop32(all[count].adj[i]);
            all[count].edges += all[count].deg[i];
        }
        all[count].edges /= 2;
        ++count;
    }
    fclose(f);
    if (count != EXPECTED_N14) {
        fprintf(stderr, "n=14 count %d, expected %d\n", count,
                EXPECTED_N14);
        exit(2);
    }
    return count;
}

static int select_lm_patterns(const char *path,
                              Pattern all[EXPECTED_N14],
                              Pattern selected[MAX_PATTERNS]) {
    FILE *f = fopen(path, "r");
    char line[1024];
    uint8_t seen[EXPECTED_N14] = {0};
    int count = 0, lines = 0;
    if (!f) {
        fprintf(stderr, "cannot open %s: %s\n", path, strerror(errno));
        exit(2);
    }
    while (fgets(line, sizeof(line), f)) {
        int idx, n;
        double residual_any, residual_distinct;
        if (sscanf(line, "%d %d %lf %lf", &idx, &n,
                   &residual_any, &residual_distinct) != 4 ||
                idx < 0 || idx >= EXPECTED_N14 || n != NP || seen[idx]) {
            fprintf(stderr, "malformed/duplicate LM row %d\n", lines);
            exit(2);
        }
        (void)residual_any;
        seen[idx] = 1;
        all[idx].lm_residual = residual_distinct;
        if (residual_distinct > 1e-6) {
            if (count >= MAX_PATTERNS) {
                fprintf(stderr, "too many selected patterns\n");
                exit(2);
            }
            selected[count++] = all[idx];
        }
        ++lines;
    }
    fclose(f);
    if (lines != EXPECTED_N14 || count != EXPECTED_SELECTED) {
        fprintf(stderr, "LM rows/selected = %d/%d, expected %d/%d\n",
                lines, count, EXPECTED_N14, EXPECTED_SELECTED);
        exit(2);
    }
    return count;
}

static uint8_t *load_killed(const char *path, uint64_t *unique_count) {
    FILE *f = fopen(path, "r");
    uint8_t *bits = xcalloc(EXPECTED_TARGETS, 1);
    char line[256];
    uint64_t line_number = 0;
    *unique_count = 0;
    if (!f) {
        fprintf(stderr, "cannot open %s: %s\n", path, strerror(errno));
        exit(2);
    }
    while (fgets(line, sizeof(line), f)) {
        unsigned idx;
        int decomposition;
        char extra;
        int fields;
        ++line_number;
        if (!strchr(line, '\n') && !feof(f)) {
            fprintf(stderr, "kill-log line %" PRIu64 " is too long\n",
                    line_number);
            exit(2);
        }
        /* The third conversion catches any non-whitespace trailing field,
         * preventing a dangling/misaligned token from becoming the next
         * record.  A final valid line need not end in a newline. */
        fields = sscanf(line, " %u %d %c", &idx, &decomposition, &extra);
        if (fields != 2) {
            fprintf(stderr, "malformed kill-log line %" PRIu64 "\n",
                    line_number);
            exit(2);
        }
        (void)decomposition;
        if (idx >= EXPECTED_TARGETS) {
            fprintf(stderr, "kill-log index out of range: %u\n", idx);
            exit(2);
        }
        if (!bits[idx]) {
            bits[idx] = 1;
            ++*unique_count;
        }
    }
    if (ferror(f)) {
        fprintf(stderr, "error reading kill log\n");
        exit(2);
    }
    fclose(f);
    if (*unique_count != EXPECTED_KILLED) {
        fprintf(stderr, "unique killed count %" PRIu64 ", expected %u\n",
                *unique_count, EXPECTED_KILLED);
        exit(2);
    }
    return bits;
}

/* Pull-boundary exact-stratum classification. */
static int bad_link_rec(const uint32_t adj[NT], uint32_t cand,
                        uint32_t common, int depth, int want, int bound) {
    if (depth == want) return pop32(common) > bound;
    if (pop32(cand) < want - depth) return 0;
    while (cand) {
        int v = ctz32(cand);
        cand &= cand - 1;
        if (bad_link_rec(adj, cand & adj[v], common & adj[v], depth + 1,
                         want, bound))
            return 1;
        if (pop32(cand) < want - depth) return 0;
    }
    return 0;
}

static int bad_k5_link(const uint32_t adj[NT]) {
    return bad_link_rec(adj, FULL19, FULL19, 0, 5, 4);
}

static int bad_reflection_seed(const uint32_t adj[NT], uint32_t seed) {
    int forced[7];
    uint32_t outside = FULL19 & ~seed;
    for (int i = 0; i < 7; ++i) forced[i] = -1;
    while (outside) {
        uint32_t defects;
        int v = ctz32(outside), type_vertex, type;
        outside &= outside - 1;
        defects = seed & ~adj[v];
        if (pop32(defects) != 1) continue;
        type_vertex = ctz32(defects);
        type = pop32(seed & ((1u << type_vertex) - 1u));
        if (forced[type] >= 0) return 1;
        for (int j = 0; j < 7; ++j)
            if (forced[j] >= 0 && ((adj[v] >> forced[j]) & 1u))
                return 1;
        forced[type] = v;
    }
    return 0;
}

typedef struct {
    int has_k7;
    int bad_reflection;
} K7Scan;

static void scan_k7_rec(const uint32_t adj[NT], uint32_t cand,
                        uint32_t seed, int depth, K7Scan *out) {
    if (out->bad_reflection) return;
    if (depth == 7) {
        out->has_k7 = 1;
        if (bad_reflection_seed(adj, seed)) out->bad_reflection = 1;
        return;
    }
    if (pop32(cand) < 7 - depth) return;
    while (cand && !out->bad_reflection) {
        int v = ctz32(cand);
        cand &= cand - 1;
        scan_k7_rec(adj, cand & adj[v], seed | (1u << v), depth + 1,
                    out);
        if (pop32(cand) < 7 - depth) return;
    }
}

static K7Scan scan_k7(const uint32_t adj[NT]) {
    K7Scan out = {0, 0};
    scan_k7_rec(adj, FULL19, 0, 0, &out);
    return out;
}

static void reservoir_offer(Reservoir *r, const Target *t) {
    uint64_t j;
    ++r->seen;
    if (r->used < r->cap) {
        r->items[r->used++] = *t;
        return;
    }
    j = rand_bounded(&r->rng, r->seen);
    if (j < r->cap) r->items[j] = *t;
}

static int compare_target(const void *a, const void *b) {
    const Target *x = a, *y = b;
    if (x->stratum != y->stratum)
        return (int)x->stratum - (int)y->stratum;
    if (x->corpus_index < y->corpus_index) return -1;
    if (x->corpus_index > y->corpus_index) return 1;
    return 0;
}

static Target *sample_targets(const char *corpus_path, const uint8_t *killed,
                              size_t per_stratum, uint64_t seed,
                              size_t *sample_count, uint64_t populations[NSTRATA],
                              double *scan_seconds) {
    Reservoir reservoirs[NSTRATA];
    FILE *f = fopen(corpus_path, "r");
    uint32_t idx = 0;
    double t0 = monotime();
    if (!f) {
        fprintf(stderr, "cannot open %s: %s\n", corpus_path,
                strerror(errno));
        exit(2);
    }
    memset(populations, 0, NSTRATA * sizeof(populations[0]));
    for (int s = 0; s < NSTRATA; ++s) {
        reservoirs[s].items = xcalloc(per_stratum, sizeof(Target));
        reservoirs[s].used = 0;
        reservoirs[s].cap = per_stratum;
        reservoirs[s].seen = 0;
        reservoirs[s].rng = seed ^
            (UINT64_C(0xd1b54a32d192ed03) * (uint64_t)(s + 1));
    }
    for (;;) {
        Target target;
        K7Scan k7;
        int n, got, bad, stratum;
        got = fscanf(f, "%d", &n);
        if (got == EOF) break;
        if (got != 1 || n != NT || idx >= EXPECTED_TARGETS) {
            fprintf(stderr, "bad n=19 header at index %u\n", idx);
            exit(2);
        }
        memset(&target, 0, sizeof(target));
        target.corpus_index = idx;
        for (int i = 0; i < NT; ++i)
            if (fscanf(f, "%u", &target.adj[i]) != 1) {
                fprintf(stderr, "truncated n=19 graph %u\n", idx);
                exit(2);
            }
        if (!validate_graph(target.adj, NT)) {
            fprintf(stderr, "invalid n=19 graph %u\n", idx);
            exit(2);
        }
        k7 = scan_k7(target.adj);
        bad = k7.bad_reflection || bad_k5_link(target.adj);
        if (killed[idx]) {
            if (!k7.has_k7) {
                fprintf(stderr, "unexpected certified K6-only graph %u\n",
                        idx);
                exit(2);
            }
            stratum = bad ? CERTIFIED_EXACT_REJECTED :
                            CERTIFIED_EXACT_RESIDUE;
        } else if (k7.has_k7) {
            stratum = bad ? DEFERRED_K7_EXACT_REJECTED :
                            DEFERRED_K7_EXACT_RESIDUE;
        } else {
            if (bad) {
                fprintf(stderr, "unexpected exact-rejected K6 graph %u\n",
                        idx);
                exit(2);
            }
            stratum = DEFERRED_K6_EXACT_RESIDUE;
        }
        target.stratum = (uint8_t)stratum;
        ++populations[stratum];
        reservoir_offer(&reservoirs[stratum], &target);
        ++idx;
        if ((idx % 500000u) == 0)
            fprintf(stderr, "sample scan: %u/%u targets\n", idx,
                    EXPECTED_TARGETS);
    }
    fclose(f);
    if (idx != EXPECTED_TARGETS) {
        fprintf(stderr, "n=19 count %u, expected %u\n", idx,
                EXPECTED_TARGETS);
        exit(2);
    }
    for (int s = 0; s < NSTRATA; ++s) {
        if (populations[s] != expected_strata[s] ||
                reservoirs[s].seen != populations[s]) {
            fprintf(stderr, "stratum %s count %" PRIu64
                    ", expected %" PRIu64 "\n", stratum_names[s],
                    populations[s], expected_strata[s]);
            exit(2);
        }
    }
    *sample_count = 0;
    for (int s = 0; s < NSTRATA; ++s) *sample_count += reservoirs[s].used;
    Target *sample = xcalloc(*sample_count, sizeof(Target));
    size_t at = 0;
    for (int s = 0; s < NSTRATA; ++s) {
        memcpy(sample + at, reservoirs[s].items,
               reservoirs[s].used * sizeof(Target));
        at += reservoirs[s].used;
        free(reservoirs[s].items);
    }
    qsort(sample, *sample_count, sizeof(Target), compare_target);
    *scan_seconds = monotime() - t0;
    return sample;
}

/* ---------- Conservative non-induced subgraph monomorphism ---------- */

typedef struct {
    const Pattern *pattern;
    const Target *target;
    uint32_t base_domain[NP];
    int8_t map[NP];
    uint32_t used;
    uint64_t nodes;
    uint64_t limit;
} MatchState;

/* Kuhn matching for a small family of target-vertex domains. */
static int augment_domain(int item, const uint32_t *sets,
                          int owner[NT], uint32_t *seen) {
    uint32_t candidates = sets[item] & ~*seen;
    while (candidates) {
        int t = ctz32(candidates);
        candidates &= candidates - 1;
        *seen |= 1u << t;
        if (owner[t] < 0 || augment_domain(owner[t], sets, owner, seen)) {
            owner[t] = item;
            return 1;
        }
    }
    return 0;
}

static int domains_have_sdr(const uint32_t *sets, int count) {
    int order[NP], owner[NT];
    for (int i = 0; i < count; ++i) order[i] = i;
    for (int i = 0; i < count; ++i)
        for (int j = i + 1; j < count; ++j)
            if (pop32(sets[order[j]]) < pop32(sets[order[i]])) {
                int tmp = order[i]; order[i] = order[j]; order[j] = tmp;
            }
    for (int t = 0; t < NT; ++t) owner[t] = -1;
    for (int oi = 0; oi < count; ++oi) {
        uint32_t seen = 0;
        if (!augment_domain(order[oi], sets, owner, &seen)) return 0;
    }
    return 1;
}

/* Necessary local condition for p -> t: all neighbours of p must map to
 * distinct neighbours of t with adequate total target degree. */
static int root_neighbour_sdr(const MatchState *m, int p, int t,
                              const uint8_t target_deg[NT]) {
    uint32_t sets[NP];
    int count = 0;
    uint32_t neighbours = m->pattern->adj[p];
    while (neighbours) {
        int q = ctz32(neighbours);
        uint32_t choices = 0, target_neighbours = m->target->adj[t];
        neighbours &= neighbours - 1;
        while (target_neighbours) {
            int u = ctz32(target_neighbours);
            target_neighbours &= target_neighbours - 1;
            if (target_deg[u] >= m->pattern->deg[q]) choices |= 1u << u;
        }
        if (!choices) return 0;
        sets[count++] = choices;
    }
    return domains_have_sdr(sets, count);
}

static int initialize_domains(MatchState *m) {
    uint8_t target_deg[NT];
    for (int t = 0; t < NT; ++t)
        target_deg[t] = (uint8_t)pop32(m->target->adj[t]);
    for (int p = 0; p < NP; ++p) {
        uint32_t domain = 0;
        for (int t = 0; t < NT; ++t) {
            if (target_deg[t] < m->pattern->deg[p]) continue;
            if (root_neighbour_sdr(m, p, t, target_deg)) domain |= 1u << t;
        }
        m->base_domain[p] = domain;
        if (!domain) return 0;
    }
    return 1;
}

static int build_domains(const MatchState *m, uint32_t domain[NP],
                         int *unmapped_count) {
    int changed;
    *unmapped_count = 0;
    for (int p = 0; p < NP; ++p) {
        if (m->map[p] >= 0) {
            domain[p] = 1u << m->map[p];
            continue;
        }
        uint32_t d = m->base_domain[p] & ~m->used;
        for (int q = 0; q < NP; ++q)
            if (m->map[q] >= 0 && ((m->pattern->adj[p] >> q) & 1u))
                d &= m->target->adj[(int)m->map[q]];
        if (!d) return 0;
        domain[p] = d;
        ++*unmapped_count;
    }
    /* Edge arc consistency among unassigned pattern vertices. */
    do {
        changed = 0;
        for (int p = 0; p < NP; ++p) {
            uint32_t old, keep = 0;
            if (m->map[p] >= 0) continue;
            old = domain[p];
            while (old) {
                int t = ctz32(old);
                int supported = 1;
                uint32_t qs = m->pattern->adj[p];
                old &= old - 1;
                while (qs) {
                    int q = ctz32(qs);
                    qs &= qs - 1;
                    if (m->map[q] < 0 &&
                            !(domain[q] & m->target->adj[t])) {
                        supported = 0;
                        break;
                    }
                }
                if (supported) keep |= 1u << t;
            }
            if (!keep) return 0;
            if (keep != domain[p]) {
                domain[p] = keep;
                changed = 1;
            }
        }
    } while (changed);

    uint32_t sets[NP];
    int count = 0;
    for (int p = 0; p < NP; ++p)
        if (m->map[p] < 0) sets[count++] = domain[p];
    return domains_have_sdr(sets, count);
}

static int branch_neighbour_sdr(const MatchState *m, int p, int t,
                                const uint32_t domain[NP]) {
    uint32_t sets[NP];
    int count = 0;
    uint32_t qs = m->pattern->adj[p];
    while (qs) {
        int q = ctz32(qs);
        qs &= qs - 1;
        if (m->map[q] >= 0) continue;
        sets[count] = domain[q] & m->target->adj[t] & ~(1u << t);
        if (!sets[count]) return 0;
        ++count;
    }
    return domains_have_sdr(sets, count);
}

static int verify_mapping(const MatchState *m) {
    uint32_t used = 0;
    for (int p = 0; p < NP; ++p) {
        int t = m->map[p];
        if (t < 0 || t >= NT || (used & (1u << t))) return 0;
        used |= 1u << t;
    }
    for (int p = 0; p < NP; ++p) {
        uint32_t qs = m->pattern->adj[p];
        while (qs) {
            int q = ctz32(qs);
            qs &= qs - 1;
            if (!((m->target->adj[(int)m->map[p]] >> m->map[q]) & 1u))
                return 0;
        }
    }
    return 1;
}

enum MatchResult { MATCH_NO = 0, MATCH_YES = 1, MATCH_TIMEOUT = 2 };

static enum MatchResult match_rec(MatchState *m) {
    uint32_t domain[NP];
    int unmapped, chosen = -1, chosen_size = NT + 1, chosen_frontier = -1;
    if (++m->nodes > m->limit) return MATCH_TIMEOUT;
    if (!build_domains(m, domain, &unmapped)) return MATCH_NO;
    if (unmapped == 0) {
        if (!verify_mapping(m)) {
            fprintf(stderr, "internal error: invalid reported embedding\n");
            exit(3);
        }
        return MATCH_YES;
    }
    for (int p = 0; p < NP; ++p) {
        int size, frontier = 0;
        if (m->map[p] >= 0) continue;
        size = pop32(domain[p]);
        for (int q = 0; q < NP; ++q)
            if (m->map[q] < 0 && ((m->pattern->adj[p] >> q) & 1u))
                ++frontier;
        if (size < chosen_size ||
                (size == chosen_size && frontier > chosen_frontier) ||
                (size == chosen_size && frontier == chosen_frontier &&
                 m->pattern->deg[p] > m->pattern->deg[chosen])) {
            chosen = p;
            chosen_size = size;
            chosen_frontier = frontier;
        }
    }
    uint32_t candidates = domain[chosen];
    while (candidates) {
        int best_t = -1, best_score = -1;
        uint32_t tmp = candidates;
        while (tmp) {
            int t = ctz32(tmp), score = 0;
            tmp &= tmp - 1;
            for (int q = 0; q < NP; ++q)
                if (m->map[q] < 0 && q != chosen &&
                        ((m->pattern->adj[chosen] >> q) & 1u))
                    score += pop32(domain[q] & m->target->adj[t]);
            if (score > best_score) {
                best_t = t;
                best_score = score;
            }
        }
        candidates &= ~(1u << best_t);
        if (!branch_neighbour_sdr(m, chosen, best_t, domain)) continue;
        m->map[chosen] = (int8_t)best_t;
        m->used |= 1u << best_t;
        enum MatchResult result = match_rec(m);
        if (result == MATCH_YES || result == MATCH_TIMEOUT) return result;
        m->used &= ~(1u << best_t);
        m->map[chosen] = -1;
    }
    return MATCH_NO;
}

static enum MatchResult embeds(const Pattern *pattern, const Target *target,
                               uint64_t limit, uint64_t *nodes) {
    MatchState m;
    memset(&m, 0, sizeof(m));
    m.pattern = pattern;
    m.target = target;
    m.limit = limit;
    for (int p = 0; p < NP; ++p) m.map[p] = -1;
    if (!initialize_domains(&m)) {
        *nodes = 0;
        return MATCH_NO;
    }
    enum MatchResult result = match_rec(&m);
    *nodes = m.nodes;
    return result;
}

typedef struct {
    const Pattern *patterns;
    int pattern_count;
    const Target *targets;
    TargetResult *results;
    size_t target_count;
    volatile size_t next_target;
    volatile size_t finished;
    uint64_t node_limit;
    pthread_mutex_t print_mutex;
    double start_time;
} Work;

static void *worker_main(void *arg) {
    Work *w = arg;
    for (;;) {
        size_t ti = __sync_fetch_and_add(&w->next_target, 1);
        if (ti >= w->target_count) break;
        TargetResult *tr = &w->results[ti];
        for (int pi = 0; pi < w->pattern_count; ++pi) {
            uint64_t nodes = 0;
            enum MatchResult result = embeds(&w->patterns[pi],
                                             &w->targets[ti],
                                             w->node_limit, &nodes);
            tr->nodes += nodes;
            if (result == MATCH_YES)
                tr->hit[pi >> 6] |= UINT64_C(1) << (pi & 63);
            else if (result == MATCH_TIMEOUT)
                tr->timeout[pi >> 6] |= UINT64_C(1) << (pi & 63);
        }
        size_t done = __sync_add_and_fetch(&w->finished, 1);
        if (done % 50 == 0 || done == w->target_count) {
            pthread_mutex_lock(&w->print_mutex);
            fprintf(stderr, "containment: %zu/%zu targets (%.1f s)\n",
                    done, w->target_count, monotime() - w->start_time);
            pthread_mutex_unlock(&w->print_mutex);
        }
    }
    return NULL;
}

static int bit_is_set(const uint64_t bits[2], int p) {
    return (int)((bits[p >> 6] >> (p & 63)) & 1u);
}

static int target_has_any_hit(const TargetResult *r) {
    return r->hit[0] != 0 || r->hit[1] != 0;
}

static void json_string(FILE *out, const char *s) {
    fputc('"', out);
    while (*s) {
        unsigned char c = (unsigned char)*s++;
        if (c == '"' || c == '\\') fprintf(out, "\\%c", c);
        else if (c == '\n') fputs("\\n", out);
        else if (c < 0x20) fprintf(out, "\\u%04x", c);
        else fputc(c, out);
    }
    fputc('"', out);
}

static void print_count_by_stratum(FILE *out, const uint64_t count[NSTRATA]) {
    fputc('{', out);
    for (int s = 0; s < NSTRATA; ++s) {
        if (s) fputs(", ", out);
        json_string(out, stratum_names[s]);
        fprintf(out, ": %" PRIu64, count[s]);
    }
    fputc('}', out);
}

static void wilson95(uint64_t hits, uint64_t total,
                     double *low, double *high) {
    const double z = 1.959963984540054;
    if (!total) { *low = 0.0; *high = 1.0; return; }
    double n = (double)total, p = (double)hits / n;
    double denominator = 1.0 + z * z / n;
    double center = (p + z * z / (2.0 * n)) / denominator;
    double half = z * sqrt(p * (1.0 - p) / n +
                           z * z / (4.0 * n * n)) / denominator;
    *low = center - half;
    *high = center + half;
}

static void print_greedy(FILE *out, const Pattern *patterns, int pattern_count,
                         const Target *targets, const TargetResult *results,
                         size_t target_count, int stratum_filter) {
    uint8_t *covered = xcalloc(target_count, 1);
    uint8_t chosen[MAX_PATTERNS] = {0};
    uint64_t cumulative = 0;
    int first = 1;
    fputc('[', out);
    for (int step = 0; step < pattern_count; ++step) {
        int best = -1;
        uint64_t best_marginal = 0;
        for (int p = 0; p < pattern_count; ++p) {
            uint64_t marginal = 0;
            if (chosen[p]) continue;
            for (size_t i = 0; i < target_count; ++i)
                if (!covered[i] &&
                        (stratum_filter < 0 ||
                         targets[i].stratum == stratum_filter) &&
                        bit_is_set(results[i].hit, p))
                    ++marginal;
            if (marginal > best_marginal ||
                    (marginal == best_marginal && marginal > 0 &&
                     (best < 0 || patterns[p].source_index <
                                  patterns[best].source_index))) {
                best = p;
                best_marginal = marginal;
            }
        }
        if (best < 0 || best_marginal == 0) break;
        chosen[best] = 1;
        for (size_t i = 0; i < target_count; ++i)
            if (!covered[i] &&
                    (stratum_filter < 0 ||
                     targets[i].stratum == stratum_filter) &&
                    bit_is_set(results[i].hit, best))
                covered[i] = 1;
        cumulative += best_marginal;
        if (!first) fputs(",", out);
        fprintf(out, "\n        {\"step\": %d, \"pattern_index\": %d, "
                "\"marginal_hits\": %" PRIu64 ", "
                "\"cumulative_hits\": %" PRIu64 "}",
                step + 1, patterns[best].source_index, best_marginal,
                cumulative);
        first = 0;
    }
    if (!first) fputc('\n', out);
    fputs("      ]", out);
    free(covered);
}

static void write_json(const char *path, int argc, char **argv,
                       const char *n14_path, const char *lm_path,
                       const char *n19_path, const char *kill_path,
                       const Pattern *patterns, int pattern_count,
                       const Target *targets, const TargetResult *results,
                       size_t target_count, const uint64_t populations[NSTRATA],
                       size_t per_stratum, uint64_t node_limit, int threads,
                       uint64_t seed, double scan_seconds, double match_seconds,
                       double cpu_seconds) {
    FILE *out = fopen(path, "w");
    uint64_t sample_by[NSTRATA] = {0}, union_by[NSTRATA] = {0};
    uint64_t total_nodes = 0, total_hits = 0, total_timeouts = 0;
    double weighted_old_residue, weighted_low, weighted_high;
    if (!out) {
        fprintf(stderr, "cannot create %s: %s\n", path, strerror(errno));
        exit(2);
    }
    for (size_t i = 0; i < target_count; ++i) {
        int s = targets[i].stratum;
        ++sample_by[s];
        if (target_has_any_hit(&results[i])) ++union_by[s];
        total_nodes += results[i].nodes;
        total_hits += (uint64_t)__builtin_popcountll(results[i].hit[0]) +
                      (uint64_t)__builtin_popcountll(results[i].hit[1]);
        total_timeouts +=
            (uint64_t)__builtin_popcountll(results[i].timeout[0]) +
            (uint64_t)__builtin_popcountll(results[i].timeout[1]);
    }
    fprintf(out, "{\n  \"schema\": 1,\n"
            "  \"claim_status\": \"HEURISTIC COVERAGE ONLY; the 89 LM-selected patterns are not certified geometric obstructions\",\n"
            "  \"containment_semantics\": \"non-induced: pattern edges map to target edges; pattern nonedges are unconstrained\",\n"
            "  \"timeout_policy\": \"every target-pattern timeout is recorded and counted conservatively as no hit\",\n"
            "  \"sampling_design\": \"five equal-size deterministic reservoir strata; aggregate sample fractions and greedy orders are unweighted across deliberately unequal populations\",\n"
            "  \"profile_boundary_commit\": \"9292f29be9cd26e265ed44f7473cf36bf931e68c\",\n");
    fputs("  \"command\": [", out);
    for (int i = 0; i < argc; ++i) {
        if (i) fputs(", ", out);
        json_string(out, argv[i]);
    }
    fputs("],\n  \"inputs\": {\n", out);
    fputs("    \"n14_corpus\": ", out); json_string(out, n14_path);
    fputs(",\n    \"lm_triage\": ", out); json_string(out, lm_path);
    fputs(",\n    \"n19_corpus\": ", out); json_string(out, n19_path);
    fputs(",\n    \"kill_log\": ", out); json_string(out, kill_path);
    fputs("\n  },\n", out);
    fprintf(out, "  \"parameters\": {\"sample_per_stratum\": %zu, "
            "\"node_limit_per_target_pattern\": %" PRIu64 ", "
            "\"threads\": %d, \"seed\": %" PRIu64 "},\n",
            per_stratum, node_limit, threads, seed);
    fprintf(out, "  \"input_counts\": {\"n14\": %d, "
            "\"lm_selected\": %d, \"n19\": %u, "
            "\"old_certified\": %u},\n", EXPECTED_N14, pattern_count,
            EXPECTED_TARGETS, EXPECTED_KILLED);
    fputs("  \"pull_boundary_population_counts\": ", out);
    print_count_by_stratum(out, populations);
    fputs(",\n  \"sample_counts\": ", out);
    print_count_by_stratum(out, sample_by);
    fputs(",\n  \"runtime_seconds\": {", out);
    fprintf(out, "\"classification_and_sampling_wall\": %.6f, "
            "\"containment_wall\": %.6f, \"total_wall\": %.6f, "
            "\"process_cpu\": %.6f},\n", scan_seconds, match_seconds,
            scan_seconds + match_seconds, cpu_seconds);
    fprintf(out, "  \"search_totals\": {\"target_pattern_calls\": %zu, "
            "\"hits\": %" PRIu64 ", \"timeouts\": %" PRIu64 ", "
            "\"exhaustive_no_hits\": %" PRIu64 ", "
            "\"search_nodes\": %" PRIu64 "},\n",
            target_count * (size_t)pattern_count, total_hits, total_timeouts,
            (uint64_t)(target_count * (size_t)pattern_count) - total_hits -
                total_timeouts, total_nodes);
    uint64_t union_total = 0;
    for (int s = 0; s < NSTRATA; ++s) union_total += union_by[s];
    fprintf(out, "  \"unweighted_stratified_union_coverage\": {"
            "\"interpretation\": \"descriptive only; not a population estimate\", "
            "\"hits\": %" PRIu64 ", "
            "\"sample\": %zu, \"fraction\": %.9f, \"by_stratum\": ",
            union_total, target_count,
            target_count ? (double)union_total / (double)target_count : 0.0);
    print_count_by_stratum(out, union_by);
    fputs("},\n", out);
    {
        const int strata[2] = {DEFERRED_K7_EXACT_RESIDUE,
                               DEFERRED_K6_EXACT_RESIDUE};
        uint64_t population_total = populations[strata[0]] +
                                    populations[strata[1]];
        weighted_old_residue = 0.0;
        weighted_low = 0.0;
        weighted_high = 0.0;
        for (int j = 0; j < 2; ++j) {
            int s = strata[j];
            double low, high;
            double weight = (double)populations[s] /
                            (double)population_total;
            wilson95(union_by[s], sample_by[s], &low, &high);
            weighted_old_residue += weight *
                ((double)union_by[s] / (double)sample_by[s]);
            weighted_low += weight * low;
            weighted_high += weight * high;
        }
        fprintf(out, "  \"population_weighted_preexisting_deferred_residue_estimate\": {"
                "\"population\": %" PRIu64 ", \"point_estimate\": %.12f, "
                "\"approximate_95pct_low\": %.12f, "
                "\"approximate_95pct_high\": %.12f, "
                "\"uncertainty_method\": \"population-weighted stratumwise Wilson score intervals; descriptive, not simultaneous\"},\n",
                population_total, weighted_old_residue, weighted_low,
                weighted_high);
    }
    fputs("  \"individual_patterns\": [\n", out);
    for (int p = 0; p < pattern_count; ++p) {
        uint64_t hit_by[NSTRATA] = {0}, timeout_by[NSTRATA] = {0};
        uint64_t hits = 0, timeouts = 0;
        for (size_t i = 0; i < target_count; ++i) {
            int s = targets[i].stratum;
            if (bit_is_set(results[i].hit, p)) { ++hits; ++hit_by[s]; }
            if (bit_is_set(results[i].timeout, p)) {
                ++timeouts; ++timeout_by[s];
            }
        }
        fprintf(out, "    {\"ordinal\": %d, \"pattern_index\": %d, "
                "\"edges\": %d, \"lm_distinct_residual\": %.17g, "
                "\"hits\": %" PRIu64 ", \"timeouts\": %" PRIu64
                ", \"hits_by_stratum\": ", p,
                patterns[p].source_index, patterns[p].edges,
                patterns[p].lm_residual, hits, timeouts);
        print_count_by_stratum(out, hit_by);
        fputs(", \"timeouts_by_stratum\": ", out);
        print_count_by_stratum(out, timeout_by);
        fprintf(out, "}%s\n", p + 1 == pattern_count ? "" : ",");
    }
    fputs("  ],\n  \"greedy_marginal_coverage_unweighted_stratified\": {\n"
          "    \"interpretation\": \"ranked on equal-size strata; not population weighted\",\n"
          "    \"all_strata\": ", out);
    print_greedy(out, patterns, pattern_count, targets, results, target_count,
                 -1);
    for (int s = 0; s < NSTRATA; ++s) {
        fputs(",\n    ", out); json_string(out, stratum_names[s]);
        fputs(": ", out);
        print_greedy(out, patterns, pattern_count, targets, results,
                     target_count, s);
    }
    fputs("\n  },\n  \"sample_audit\": [\n", out);
    for (size_t i = 0; i < target_count; ++i) {
        fprintf(out, "    {\"corpus_index\": %u, \"stratum\": ",
                targets[i].corpus_index);
        json_string(out, stratum_names[targets[i].stratum]);
        fprintf(out, ", \"hit_mask_hex_high_low\": "
                "\"%016" PRIx64 "%016" PRIx64 "\", "
                "\"timeout_mask_hex_high_low\": "
                "\"%016" PRIx64 "%016" PRIx64 "\", "
                "\"search_nodes\": %" PRIu64 "}%s\n",
                results[i].hit[1], results[i].hit[0],
                results[i].timeout[1], results[i].timeout[0],
                results[i].nodes, i + 1 == target_count ? "" : ",");
    }
    fputs("  ]\n}\n", out);
    if (fclose(out) != 0) {
        fprintf(stderr, "error closing %s\n", path);
        exit(2);
    }
}

static uint64_t parse_u64(const char *s, const char *what) {
    char *end = NULL;
    errno = 0;
    unsigned long long x;
    if (!s[0] || s[0] == '-') {
        fprintf(stderr, "invalid %s: %s\n", what, s);
        exit(2);
    }
    x = strtoull(s, &end, 10);
    if (errno || !end || *end) {
        fprintf(stderr, "invalid %s: %s\n", what, s);
        exit(2);
    }
    return (uint64_t)x;
}

int main(int argc, char **argv) {
    Pattern all[EXPECTED_N14], patterns[MAX_PATTERNS];
    uint8_t *killed;
    uint64_t unique_killed, populations[NSTRATA], seed = 20260806;
    uint64_t node_limit = 50000;
    size_t per_stratum = 160, target_count;
    int threads = 12, pattern_count;
    Target *targets;
    TargetResult *results;
    Work work;
    pthread_t *worker_threads;
    double scan_seconds, match_seconds, match_start;
    clock_t cpu_start;

    if (argc < 6 || argc > 10) {
        fprintf(stderr, "usage: %s N14.txt LM.out N19.txt KILL.log OUT.json "
                "[sample_per_stratum] [node_limit] [threads] [seed]\n",
                argv[0]);
        return 2;
    }
    if (argc >= 7) {
        uint64_t parsed = parse_u64(argv[6], "sample size");
        if (parsed > EXPECTED_TARGETS || parsed > SIZE_MAX) {
            fprintf(stderr, "sample size is too large: %" PRIu64 "\n",
                    parsed);
            return 2;
        }
        per_stratum = (size_t)parsed;
    }
    if (argc >= 8) node_limit = parse_u64(argv[7], "node limit");
    if (argc >= 9) {
        uint64_t parsed = parse_u64(argv[8], "thread count");
        if (parsed > 128) {
            fprintf(stderr, "thread count is too large: %" PRIu64 "\n",
                    parsed);
            return 2;
        }
        threads = (int)parsed;
    }
    if (argc >= 10) seed = parse_u64(argv[9], "seed");
    if (per_stratum == 0 || node_limit == 0 || threads < 1 || threads > 128) {
        fprintf(stderr, "sample, node limit, and thread count must be positive"
                " (threads <= 128)\n");
        return 2;
    }
    cpu_start = clock();
    load_n14(argv[1], all);
    pattern_count = select_lm_patterns(argv[2], all, patterns);
    killed = load_killed(argv[4], &unique_killed);
    (void)unique_killed;
    fprintf(stderr, "loaded %d LM-selected patterns; scanning corpus\n",
            pattern_count);
    targets = sample_targets(argv[3], killed, per_stratum, seed,
                             &target_count, populations, &scan_seconds);
    free(killed);
    fprintf(stderr, "sampled %zu targets in %.1f s; matching with %d threads\n",
            target_count, scan_seconds, threads);
    results = xcalloc(target_count, sizeof(TargetResult));
    memset(&work, 0, sizeof(work));
    work.patterns = patterns;
    work.pattern_count = pattern_count;
    work.targets = targets;
    work.results = results;
    work.target_count = target_count;
    work.node_limit = node_limit;
    work.start_time = monotime();
    if (pthread_mutex_init(&work.print_mutex, NULL) != 0) {
        fprintf(stderr, "pthread_mutex_init failed\n");
        return 2;
    }
    worker_threads = xcalloc((size_t)threads, sizeof(pthread_t));
    match_start = monotime();
    int created = 0;
    for (int i = 0; i < threads; ++i) {
        int err = pthread_create(&worker_threads[i], NULL, worker_main, &work);
        if (err != 0) {
            fprintf(stderr, "pthread_create failed at worker %d: %s\n", i,
                    strerror(err));
            for (int j = 0; j < created; ++j)
                (void)pthread_join(worker_threads[j], NULL);
            pthread_mutex_destroy(&work.print_mutex);
            free(worker_threads);
            free(results);
            free(targets);
            return 2;
        }
        ++created;
    }
    for (int i = 0; i < created; ++i) {
        int err = pthread_join(worker_threads[i], NULL);
        if (err != 0) {
            fprintf(stderr, "pthread_join failed at worker %d: %s\n", i,
                    strerror(err));
            pthread_mutex_destroy(&work.print_mutex);
            free(worker_threads);
            free(results);
            free(targets);
            return 2;
        }
    }
    match_seconds = monotime() - match_start;
    pthread_mutex_destroy(&work.print_mutex);
    free(worker_threads);
    write_json(argv[5], argc, argv, argv[1], argv[2], argv[3], argv[4],
               patterns, pattern_count, targets, results, target_count,
               populations, per_stratum, node_limit, threads, seed,
               scan_seconds, match_seconds,
               (double)(clock() - cpu_start) / (double)CLOCKS_PER_SEC);
    fprintf(stderr, "wrote %s (scan %.1f s, match %.1f s)\n", argv[5],
            scan_seconds, match_seconds);
    free(results);
    free(targets);
    return 0;
}
