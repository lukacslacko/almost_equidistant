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
 * Build: cc -O3 -Wall -Wextra -Werror -pthread -o profile_d6 profile_d6.c
 * Usage: ./profile_d6 aeq_d6_n19.txt
 *        [killed_d6_n19.log|- [limit [workers [decisions.tsv]]]]
 */
#include <errno.h>
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <unistd.h>

#ifndef N
#define N 19
#endif
#ifndef EXPECTED
#define EXPECTED 3971787u
#endif
#if N < 7 || N > 19
#error "profile_d6 requires 7 <= N <= 19"
#endif
#define NWIT 5
#define BATCH_GRAPHS 65536u

typedef struct {
    uint32_t index;
    uint32_t seed;
    uint32_t auxiliary;
} Witness;

typedef struct {
    uint64_t total;
    uint64_t omega[9];
    uint64_t link[8];
    uint64_t cumulative;
    uint64_t k7_reflection;
    uint64_t k7_defect_csp;
    uint64_t k7_clique_hall;
    uint64_t k7_bounded_cover;
    uint64_t k7_tight_cover_matching;
    uint64_t k6_clique_hall;
    uint64_t k6_two_light_ray;
    uint64_t previous_exact;
    uint64_t previous_exact_by_omega[9];
    uint64_t cumulative_exact;
    uint64_t cumulative_exact_by_omega[9];
    Witness k7_hall_witness[NWIT];
    Witness k7_cover_witness[NWIT];
    Witness k7_tight_cover_witness[NWIT];
    Witness k6_hall_witness[NWIT];
    Witness k6_two_light_ray_witness[NWIT];
} Stats;

typedef struct {
    int omega;
    uint8_t link_mask;
    int bad_reflection;
    int bad_defect_csp;
    int bad_k7_hall;
    int bad_k7_cover;
    int bad_k7_tight_cover;
    int bad_k6_hall;
    int bad_k6_two_light_ray;
    uint32_t k7_hall_seed, k7_hall_aux;
    uint32_t k7_cover_seed, k7_cover_aux;
    uint32_t k7_tight_cover_seed, k7_tight_cover_aux;
    uint32_t k6_hall_seed, k6_hall_aux;
    uint32_t k6_two_light_ray_seed, k6_two_light_ray_aux;
} ProfileResult;

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

/* Convert graph non-neighbours in a simplex seed to coordinate masks. */
static void seed_defects(const uint32_t adj[N], uint32_t seed, int k,
                         uint8_t defects[N], uint32_t *outside_out) {
    int pos[N];
    int next = 0;
    for (int v = 0; v < N; ++v)
        pos[v] = ((seed >> v) & 1u) ? next++ : -1;
    if (next != k) { fprintf(stderr, "internal seed-size error\n"); exit(2); }
    uint32_t outside = ((1u << N) - 1) & ~seed;
    uint32_t todo = outside;
    memset(defects, 0, N * sizeof(*defects));
    while (todo) {
        int v = ctz(todo);
        todo &= todo - 1;
        uint32_t d = seed & ~adj[v];
        while (d) {
            int q = ctz(d);
            d &= d - 1;
            defects[v] |= (uint8_t)(1u << pos[q]);
        }
    }
    *outside_out = outside;
}

/* Hall failure, phrased equivalently as a clique whose union of allowed
 * coordinates has size smaller than the clique (plus virtual-coordinate
 * slack=1 for a K6 seed).  Enumerating the at most 2^13 outside cliques is
 * exact and covers every Hall subset. */
static int hall_clique_rec(const uint32_t adj[N], const uint8_t defects[N],
                           uint32_t cand, uint32_t chosen,
                           int depth, uint8_t coord_union, int slack,
                           uint32_t *aux) {
    while (cand) {
        int v = ctz(cand);
        cand &= cand - 1;
        uint8_t u = coord_union | defects[v];
        uint32_t ch = chosen | (1u << v);
        if (depth + 1 > pop(u) + slack) {
            /* Low 19 bits: failing outside clique; high bits: coordinate
             * union.  K7 uses at most 7 high bits, K6 at most 6. */
            *aux = ch | ((uint32_t)u << N);
            return 1;
        }
        if (hall_clique_rec(adj, defects, cand & adj[v], ch,
                            depth + 1, u, slack, aux)) return 1;
    }
    return 0;
}

static int bad_hall_seed(const uint32_t adj[N], uint32_t seed, int k,
                         int virtual_slack, uint32_t *aux) {
    uint8_t defects[N];
    uint32_t outside;
    seed_defects(adj, seed, k, defects, &outside);
    return hall_clique_rec(adj, defects, outside, 0, 0, 0,
                           virtual_slack, aux);
}

/* Is there an eligible vertex cover of L using at most budget vertices? */
static int cover_dfs(const uint32_t ladj[N], uint32_t outside,
                     uint32_t eligible, uint32_t chosen, int budget) {
    int u = -1, v = -1;
    uint32_t left = outside & ~chosen;
    uint32_t t = left;
    while (t && u < 0) {
        int x = ctz(t);
        t &= t - 1;
        uint32_t nb = ladj[x] & left;
        if (nb) { u = x; v = ctz(nb); }
    }
    if (u < 0) return 1;
    if (budget == 0) return 0;
    if ((eligible >> u) & 1u)
        if (cover_dfs(ladj, outside, eligible, chosen | (1u << u),
                      budget - 1)) return 1;
    if ((eligible >> v) & 1u)
        if (cover_dfs(ladj, outside, eligible, chosen | (1u << v),
                      budget - 1)) return 1;
    return 0;
}

static int is_cover(const uint32_t ladj[N], uint32_t outside,
                    uint32_t chosen) {
    uint32_t left = outside & ~chosen;
    uint32_t t = left;
    while (t) {
        int v = ctz(t);
        t &= t - 1;
        if (ladj[v] & left) return 0;
    }
    return 1;
}

/* A nonsingular 7-by-7 matrix supported in the chosen allowed defect masks
 * requires a perfect matching.  Dynamic programming over the 128 coordinate
 * subsets is exact and deliberately ignores coefficient values. */
static int defect_perfect_matching(const uint8_t defects[N],
                                   uint32_t chosen) {
    uint8_t reachable[128] = {0};
    reachable[0] = 1;
    while (chosen) {
        int v = ctz(chosen);
        chosen &= chosen - 1;
        uint8_t next[128] = {0};
        for (int used = 0; used < 128; ++used) {
            if (!reachable[used]) continue;
            uint8_t available = defects[v] & (uint8_t)~used;
            while (available) {
                int q = __builtin_ctz((unsigned)available);
                available &= (uint8_t)(available - 1);
                next[used | (1 << q)] = 1;
            }
        }
        memcpy(reachable, next, sizeof(reachable));
    }
    return reachable[127] != 0;
}

static int tight_cover_with_matching(const uint32_t ladj[N],
                                     uint32_t outside, uint32_t eligible,
                                     const uint8_t defects[N]) {
    uint32_t chosen = eligible;
    for (;;) {
        if (pop(chosen) == 7 && is_cover(ladj, outside, chosen) &&
                defect_perfect_matching(defects, chosen)) return 1;
        if (chosen == 0) break;
        chosen = (chosen - 1) & eligible;
    }
    return 0;
}

typedef struct { int bound, tight; uint32_t bound_aux, tight_aux; } CoverFlags;

static CoverFlags cover_seed_flags(const uint32_t adj[N], uint32_t seed) {
    CoverFlags flags = {0};
    uint8_t defects[N];
    uint32_t outside;
    seed_defects(adj, seed, 7, defects, &outside);
    uint32_t eligible = 0;
    uint32_t ladj[N] = {0};
    uint32_t t = outside;
    while (t) {
        int u = ctz(t);
        t &= t - 1;
        if (pop(defects[u]) >= 3) eligible |= 1u << u;
        uint32_t vv = t & adj[u];
        while (vv) {
            int v = ctz(vv);
            vv &= vv - 1;
            if ((defects[u] & defects[v]) == 0) {
                ladj[u] |= 1u << v;
                ladj[v] |= 1u << u;
            }
        }
    }
    if (!cover_dfs(ladj, outside, eligible, 0, 7)) {
        flags.bound = 1;
        flags.bound_aux = eligible;
        return flags;
    }
    if (!cover_dfs(ladj, outside, eligible, 0, 6) &&
            !tight_cover_with_matching(ladj, outside, eligible, defects)) {
        flags.tight = 1;
        flags.tight_aux = eligible;
    }
    return flags;
}

typedef struct {
    int hall, cover, tight_cover;
    uint32_t hall_seed, hall_aux;
    uint32_t cover_seed, cover_aux;
    uint32_t tight_cover_seed, tight_cover_aux;
} K7NewFlags;

static void scan_k7_new_rec(const uint32_t adj[N], uint32_t cand,
                            uint32_t seed, int depth, K7NewFlags *f) {
    if (f->hall && f->cover && f->tight_cover) return;
    if (depth == 7) {
        uint32_t aux = 0;
        if (!f->hall && bad_hall_seed(adj, seed, 7, 0, &aux)) {
            f->hall = 1; f->hall_seed = seed; f->hall_aux = aux;
        }
        if (!f->cover || !f->tight_cover) {
            CoverFlags cover = cover_seed_flags(adj, seed);
            if (!f->cover && cover.bound) {
                f->cover = 1;
                f->cover_seed = seed;
                f->cover_aux = cover.bound_aux;
            }
            if (!f->tight_cover && cover.tight) {
                f->tight_cover = 1;
                f->tight_cover_seed = seed;
                f->tight_cover_aux = cover.tight_aux;
            }
        }
        return;
    }
    if (pop(cand) < 7 - depth) return;
    while (cand) {
        int v = ctz(cand);
        cand &= cand - 1;
        scan_k7_new_rec(adj, cand & adj[v], seed | (1u << v),
                        depth + 1, f);
        if (f->hall && f->cover && f->tight_cover) return;
        if (pop(cand) < 7 - depth) return;
    }
}

static K7NewFlags scan_k7_new(const uint32_t adj[N]) {
    K7NewFlags f = {0};
    scan_k7_new_rec(adj, (1u << N) - 1, 0, 0, &f);
    return f;
}

/* Find a matching from the selected outside vertices into the six allowed
 * seed coordinates.  The tiny augmenting-path implementation is exact and
 * is called repeatedly by the K6 two-light-ray CSP. */
static int defect_matching6_augment(int v, const uint8_t defects[N],
                                    int owner[6], uint8_t *seen) {
    uint8_t choices = defects[v] & (uint8_t)~*seen;
    while (choices) {
        int q = __builtin_ctz((unsigned)choices);
        choices &= (uint8_t)(choices - 1);
        *seen |= (uint8_t)(1u << q);
        if (owner[q] < 0 ||
                defect_matching6_augment(owner[q], defects, owner, seen)) {
            owner[q] = v;
            return 1;
        }
    }
    return 0;
}

static int defect_matching6(const uint8_t defects[N], uint32_t chosen) {
    if (pop(chosen) > 6) return 0;
    int order[N], count = 0;
    while (chosen) {
        int v = ctz(chosen);
        chosen &= chosen - 1;
        int at = count;
        while (at > 0 &&
                pop(defects[order[at - 1]]) > pop(defects[v])) {
            order[at] = order[at - 1];
            --at;
        }
        order[at] = v;
        ++count;
    }
    int owner[6];
    for (int q = 0; q < 6; ++q) owner[q] = -1;
    for (int i = 0; i < count; ++i) {
        uint8_t seen = 0;
        if (!defect_matching6_augment(order[i], defects, owner, &seen))
            return 0;
    }
    return 1;
}

/* Return the connected non-bipartite components of L after deleting a set
 * of vertices.  All masks use the original graph labels. */
static int k6_nonbipartite_components(const uint32_t ladj[N],
                                      uint32_t remaining,
                                      uint32_t components[N]) {
    int count = 0;
    while (remaining) {
        uint32_t root = remaining & (uint32_t)(-(int32_t)remaining);
        uint32_t colored = root, color_one = 0, frontier = root;
        int nonbipartite = 0;
        while (frontier) {
            uint32_t bit = frontier & (uint32_t)(-(int32_t)frontier);
            frontier ^= bit;
            int v = ctz(bit);
            uint32_t neighbors = ladj[v] & remaining;
            uint32_t same = (color_one & bit) ? color_one :
                            (colored & ~color_one);
            if (neighbors & same) nonbipartite = 1;
            uint32_t fresh = neighbors & ~colored;
            if (!(color_one & bit)) color_one |= fresh;
            colored |= fresh;
            frontier |= fresh;
        }
        remaining &= ~colored;
        if (nonbipartite) components[count++] = colored;
    }
    return count;
}

/* Test one proposed zero-factor set Z0.  Every non-bipartite component of
 * L-Z0 must choose one of the two Lorentz-lightlike directions.  In either
 * direction its vertices together with Z0 support pairwise orthogonal
 * vectors in R^6, hence require a six-coordinate matching. */
static int k6_two_light_ray_z0(const uint32_t ladj[N], uint32_t outside,
                               const uint8_t defects[N], uint32_t z0) {
    if (!defect_matching6(defects, z0)) return 0;
    uint32_t components[N];
    int ncomp = k6_nonbipartite_components(ladj, outside & ~z0,
                                           components);
    if (ncomp == 0) return 1;
    int capacity = 6 - pop(z0);
    int odd_vertices = 0;
    for (int i = 0; i < ncomp; ++i) {
        int size = pop(components[i]);
        if (size > capacity) return 0;
        odd_vertices += size;
    }
    if (odd_vertices > 2 * capacity) return 0;

    /* Swapping the two light rays is a symmetry, so fix component zero in
     * the first bin and enumerate only the remaining component colors. */
    uint32_t colorings = 1u << (ncomp - 1);
    for (uint32_t coloring = 0; coloring < colorings; ++coloring) {
        uint32_t bins[2] = {z0 | components[0], z0};
        for (int i = 1; i < ncomp; ++i)
            bins[(coloring >> (i - 1)) & 1u] |= components[i];
        if (pop(bins[0]) <= 6 && pop(bins[1]) <= 6 &&
                defect_matching6(defects, bins[0]) &&
                defect_matching6(defects, bins[1])) return 1;
    }
    return 0;
}

static int k6_choose_z0(const uint32_t ladj[N], uint32_t outside,
                        const uint8_t defects[N], uint32_t available,
                        uint32_t chosen, int need) {
    if (need == 0)
        return k6_two_light_ray_z0(ladj, outside, defects, chosen);
    while (pop(available) >= need) {
        uint32_t bit = available & (uint32_t)(-(int32_t)available);
        available ^= bit;
        if (k6_choose_z0(ladj, outside, defects, available,
                         chosen | bit, need - 1)) return 1;
    }
    return 0;
}

static int bad_k6_two_light_ray_seed(const uint32_t adj[N], uint32_t seed,
                                     uint32_t *aux) {
    uint8_t defects[N];
    uint32_t outside;
    seed_defects(adj, seed, 6, defects, &outside);
    uint32_t ladj[N] = {0};
    uint32_t eligible = 0;
    uint32_t todo = outside;
    while (todo) {
        int u = ctz(todo);
        todo &= todo - 1;
        if (pop(defects[u]) >= 3) eligible |= 1u << u;
        uint32_t vv = todo & adj[u];
        while (vv) {
            int v = ctz(vv);
            vv &= vv - 1;
            if ((defects[u] & defects[v]) == 0) {
                ladj[u] |= 1u << v;
                ladj[v] |= 1u << u;
            }
        }
    }

    /* Zeros in an initially bipartite component are never needed: deleting
     * them leaves that component bipartite and omitting them only relaxes the
     * two support matchings. */
    uint32_t initial_components[N], initial_odd_union = 0;
    int initial_count = k6_nonbipartite_components(
        ladj, outside, initial_components);
    for (int i = 0; i < initial_count; ++i)
        initial_odd_union |= initial_components[i];
    eligible &= initial_odd_union;
    *aux = eligible;

    int maximum = pop(eligible);
    if (maximum > 6) maximum = 6;
    for (int size = 0; size <= maximum; ++size)
        if (k6_choose_z0(ladj, outside, defects, eligible, 0, size))
            return 0;
    return 1;
}

typedef struct {
    int hall, two_light_ray;
    uint32_t hall_seed, hall_aux;
    uint32_t two_light_ray_seed, two_light_ray_aux;
} K6NewFlags;

static void scan_k6_new_rec(const uint32_t adj[N], uint32_t cand,
                            uint32_t seed, int depth, K6NewFlags *f) {
    if (f->hall && f->two_light_ray) return;
    if (depth == 6) {
        uint32_t aux = 0;
        if (!f->hall && bad_hall_seed(adj, seed, 6, 1, &aux)) {
            f->hall = 1; f->hall_seed = seed; f->hall_aux = aux;
        }
        if (!f->two_light_ray &&
                bad_k6_two_light_ray_seed(adj, seed, &aux)) {
            f->two_light_ray = 1;
            f->two_light_ray_seed = seed;
            f->two_light_ray_aux = aux;
        }
        return;
    }
    if (pop(cand) < 6 - depth) return;
    while (cand) {
        int v = ctz(cand);
        cand &= cand - 1;
        scan_k6_new_rec(adj, cand & adj[v], seed | (1u << v),
                        depth + 1, f);
        if (f->hall && f->two_light_ray) return;
        if (pop(cand) < 6 - depth) return;
    }
}

static K6NewFlags scan_k6_new(const uint32_t adj[N]) {
    K6NewFlags f = {0};
    scan_k6_new_rec(adj, (1u << N) - 1, 0, 0, &f);
    return f;
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

static ProfileResult profile_graph(const uint32_t adj[N]) {
    static const int bounds[8] = {0, 0, 16, 12, 10, 4, 2, 0};
    ProfileResult r = {0};
    r.omega = has_clique_rec(adj, (1u << N) - 1, 7) ? 7 :
              (has_clique_rec(adj, (1u << N) - 1, 6) ? 6 : 5);
    for (int m = 2; m <= 7; ++m) {
        if (bad_link(adj, m, bounds[m])) r.link_mask |= (uint8_t)(1u << m);
    }
    if (r.omega == 7) {
        r.bad_reflection = bad_k7_reflection(adj);
        r.bad_defect_csp = bad_k7_defect(adj);
        K7NewFlags f = scan_k7_new(adj);
        r.bad_k7_hall = f.hall;
        r.bad_k7_cover = f.cover;
        r.bad_k7_tight_cover = f.tight_cover;
        r.k7_hall_seed = f.hall_seed; r.k7_hall_aux = f.hall_aux;
        r.k7_cover_seed = f.cover_seed; r.k7_cover_aux = f.cover_aux;
        r.k7_tight_cover_seed = f.tight_cover_seed;
        r.k7_tight_cover_aux = f.tight_cover_aux;
    } else if (r.omega == 6) {
        K6NewFlags f = scan_k6_new(adj);
        r.bad_k6_hall = f.hall;
        r.bad_k6_two_light_ray = f.two_light_ray;
        r.k6_hall_seed = f.hall_seed; r.k6_hall_aux = f.hall_aux;
        r.k6_two_light_ray_seed = f.two_light_ray_seed;
        r.k6_two_light_ray_aux = f.two_light_ray_aux;
    }
    return r;
}

static void add_witness(Witness out[NWIT], uint64_t previous_count,
                        unsigned idx, uint32_t seed, uint32_t auxiliary) {
    if (previous_count < NWIT) {
        out[previous_count].index = idx;
        out[previous_count].seed = seed;
        out[previous_count].auxiliary = auxiliary;
    }
}

static void stats_add(Stats *s, const ProfileResult *r, unsigned idx) {
    ++s->total;
    ++s->omega[r->omega];
    int link_any = r->link_mask != 0;
    for (int m = 2; m <= 7; ++m)
        if ((r->link_mask >> m) & 1u) ++s->link[m];
    if (link_any) ++s->cumulative;
    if (r->bad_reflection) ++s->k7_reflection;
    if (r->bad_defect_csp) ++s->k7_defect_csp;
    if (r->bad_k7_hall) {
        add_witness(s->k7_hall_witness, s->k7_clique_hall, idx,
                    r->k7_hall_seed, r->k7_hall_aux);
        ++s->k7_clique_hall;
    }
    if (r->bad_k7_cover) {
        add_witness(s->k7_cover_witness, s->k7_bounded_cover, idx,
                    r->k7_cover_seed, r->k7_cover_aux);
        ++s->k7_bounded_cover;
    }
    if (r->bad_k7_tight_cover) {
        add_witness(s->k7_tight_cover_witness,
                    s->k7_tight_cover_matching, idx,
                    r->k7_tight_cover_seed, r->k7_tight_cover_aux);
        ++s->k7_tight_cover_matching;
    }
    if (r->bad_k6_hall) {
        add_witness(s->k6_hall_witness, s->k6_clique_hall, idx,
                    r->k6_hall_seed, r->k6_hall_aux);
        ++s->k6_clique_hall;
    }
    if (r->bad_k6_two_light_ray) {
        add_witness(s->k6_two_light_ray_witness, s->k6_two_light_ray, idx,
                    r->k6_two_light_ray_seed, r->k6_two_light_ray_aux);
        ++s->k6_two_light_ray;
    }
    int previous = link_any || r->bad_reflection || r->bad_defect_csp;
    if (previous) {
        ++s->previous_exact;
        ++s->previous_exact_by_omega[r->omega];
    }
    int exact = previous || r->bad_k7_hall || r->bad_k7_cover ||
                r->bad_k7_tight_cover || r->bad_k6_hall ||
                r->bad_k6_two_light_ray;
    if (exact) {
        ++s->cumulative_exact;
        ++s->cumulative_exact_by_omega[r->omega];
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
    char line[256];
    uint64_t line_number = 0;
    *count = 0;
    while (fgets(line, sizeof(line), f)) {
        unsigned idx;
        int dec;
        char extra;
        ++line_number;
        if (!strchr(line, '\n') && !feof(f)) {
            fprintf(stderr, "kill-log line %" PRIu64 " is too long\n",
                    line_number);
            exit(2);
        }
        if (sscanf(line, " %u %d %c", &idx, &dec, &extra) != 2) {
            fprintf(stderr, "malformed kill-log line %" PRIu64 "\n",
                    line_number);
            exit(2);
        }
        if (idx >= EXPECTED) {
            fprintf(stderr, "bad kill-log index %u\n", idx);
            exit(2);
        }
        (void)dec;
        if (!bits[idx]) { bits[idx] = 1; ++*count; }
    }
    if (ferror(f)) { fprintf(stderr, "error reading kill log\n"); exit(2); }
    fclose(f);
    return bits;
}

static void print_witness_array(const Witness w[NWIT], uint64_t count) {
    putchar('[');
    uint64_t n = count < NWIT ? count : NWIT;
    for (uint64_t i = 0; i < n; ++i) {
        if (i) printf(", ");
        printf("{\"index\":%u,\"seed\":%u,\"auxiliary\":%u}",
               w[i].index, w[i].seed, w[i].auxiliary);
    }
    putchar(']');
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
    printf("      \"K7_clique_Hall_rejections\": %" PRIu64 ",\n",
           s->k7_clique_hall);
    printf("      \"K7_disjoint_edge_bounded_cover_rejections\": %" PRIu64
           ",\n", s->k7_bounded_cover);
    printf("      \"K7_tight_cover_matching_rejections\": %" PRIu64
           ",\n", s->k7_tight_cover_matching);
    printf("      \"K6_clique_Hall_rejections\": %" PRIu64 ",\n",
           s->k6_clique_hall);
    printf("      \"K6_odd_component_two_light_ray_rejections\": %" PRIu64
           ",\n", s->k6_two_light_ray);
    printf("      \"cumulative_previous_exact_rejections\": %" PRIu64
           ",\n", s->previous_exact);
    printf("      \"cumulative_exact_rejections\": %" PRIu64 ",\n",
           s->cumulative_exact);
    printf("      \"exact_rejections_by_clique_number\": {\"5_or_less\": %"
           PRIu64 ", \"6\": %" PRIu64 ", \"7\": %" PRIu64 "},\n",
           s->cumulative_exact_by_omega[5],
           s->cumulative_exact_by_omega[6],
           s->cumulative_exact_by_omega[7]);
    printf("      \"residue_by_clique_number\": {\"5_or_less\": %" PRIu64
           ", \"6\": %" PRIu64 ", \"7\": %" PRIu64 "},\n",
           s->omega[5] - s->cumulative_exact_by_omega[5],
           s->omega[6] - s->cumulative_exact_by_omega[6],
           s->omega[7] - s->cumulative_exact_by_omega[7]);
    printf("      \"new_rule_witnesses\": {\n");
    printf("        \"K7_clique_Hall\": ");
    print_witness_array(s->k7_hall_witness, s->k7_clique_hall);
    printf(",\n        \"K7_disjoint_edge_bounded_cover\": ");
    print_witness_array(s->k7_cover_witness, s->k7_bounded_cover);
    printf(",\n        \"K7_tight_cover_matching\": ");
    print_witness_array(s->k7_tight_cover_witness,
                        s->k7_tight_cover_matching);
    printf(",\n        \"K6_clique_Hall\": ");
    print_witness_array(s->k6_hall_witness, s->k6_clique_hall);
    printf(",\n        \"K6_odd_component_two_light_ray\": ");
    print_witness_array(s->k6_two_light_ray_witness,
                        s->k6_two_light_ray);
    printf("\n      }\n");
    printf("    }");
}

static int witness_cmp(const void *a, const void *b) {
    const Witness *x = a, *y = b;
    return x->index < y->index ? -1 : x->index > y->index;
}

static void merge_witnesses(Witness dst[NWIT], uint64_t dst_count,
                            const Witness src[NWIT], uint64_t src_count) {
    Witness tmp[2 * NWIT];
    int n = 0;
    int nd = dst_count < NWIT ? (int)dst_count : NWIT;
    int ns = src_count < NWIT ? (int)src_count : NWIT;
    for (int i = 0; i < nd; ++i) tmp[n++] = dst[i];
    for (int i = 0; i < ns; ++i) tmp[n++] = src[i];
    qsort(tmp, (size_t)n, sizeof(*tmp), witness_cmp);
    int keep = n < NWIT ? n : NWIT;
    for (int i = 0; i < keep; ++i) dst[i] = tmp[i];
}

static void stats_merge(Stats *dst, const Stats *src) {
    merge_witnesses(dst->k7_hall_witness, dst->k7_clique_hall,
                    src->k7_hall_witness, src->k7_clique_hall);
    merge_witnesses(dst->k7_cover_witness, dst->k7_bounded_cover,
                    src->k7_cover_witness, src->k7_bounded_cover);
    merge_witnesses(dst->k7_tight_cover_witness,
                    dst->k7_tight_cover_matching,
                    src->k7_tight_cover_witness,
                    src->k7_tight_cover_matching);
    merge_witnesses(dst->k6_hall_witness, dst->k6_clique_hall,
                    src->k6_hall_witness, src->k6_clique_hall);
    merge_witnesses(dst->k6_two_light_ray_witness,
                    dst->k6_two_light_ray,
                    src->k6_two_light_ray_witness,
                    src->k6_two_light_ray);
    dst->total += src->total;
    for (int i = 0; i < 9; ++i) {
        dst->omega[i] += src->omega[i];
        dst->previous_exact_by_omega[i] += src->previous_exact_by_omega[i];
        dst->cumulative_exact_by_omega[i] += src->cumulative_exact_by_omega[i];
    }
    for (int i = 0; i < 8; ++i) dst->link[i] += src->link[i];
    dst->cumulative += src->cumulative;
    dst->k7_reflection += src->k7_reflection;
    dst->k7_defect_csp += src->k7_defect_csp;
    dst->k7_clique_hall += src->k7_clique_hall;
    dst->k7_bounded_cover += src->k7_bounded_cover;
    dst->k7_tight_cover_matching += src->k7_tight_cover_matching;
    dst->k6_clique_hall += src->k6_clique_hall;
    dst->k6_two_light_ray += src->k6_two_light_ray;
    dst->previous_exact += src->previous_exact;
    dst->cumulative_exact += src->cumulative_exact;
}

typedef struct {
    uint32_t (*graphs)[N];
    const uint8_t *killed;
    ProfileResult *results;
    unsigned base, start, end;
    Stats all, certified, deferred;
} Worker;

static void *profile_worker(void *arg) {
    Worker *w = arg;
    for (unsigned local = w->start; local < w->end; ++local) {
        unsigned idx = w->base + local;
        ProfileResult r = profile_graph(w->graphs[local]);
        if (w->results) w->results[local] = r;
        stats_add(&w->all, &r, idx);
        if (w->killed && w->killed[idx]) stats_add(&w->certified, &r, idx);
        else if (w->killed) stats_add(&w->deferred, &r, idx);
    }
    return NULL;
}

static int write_decision_header(FILE *f) {
    return fprintf(f, "index\tomega\tlink_mask\treflection\tdefect_csp\t"
                   "K7_Hall\tK7_cover\tK7_tight_cover\tK6_Hall\t"
                   "K6_two_light_ray\t"
                   "K7_Hall_seed\tK7_Hall_aux\tK7_cover_seed\tK7_cover_aux\t"
                   "K7_tight_cover_seed\tK7_tight_cover_aux\t"
                   "K6_Hall_seed\tK6_Hall_aux\t"
                   "K6_two_light_ray_seed\tK6_two_light_ray_aux\n") >= 0;
}

static int write_decision(FILE *f, unsigned idx, const ProfileResult *r) {
    return fprintf(f, "%u\t%d\t%u\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t"
                   "%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\n",
                   idx, r->omega, r->link_mask, r->bad_reflection,
                   r->bad_defect_csp, r->bad_k7_hall, r->bad_k7_cover,
                   r->bad_k7_tight_cover, r->bad_k6_hall,
                   r->bad_k6_two_light_ray,
                   r->k7_hall_seed, r->k7_hall_aux,
                   r->k7_cover_seed, r->k7_cover_aux,
                   r->k7_tight_cover_seed, r->k7_tight_cover_aux,
                   r->k6_hall_seed, r->k6_hall_aux,
                   r->k6_two_light_ray_seed,
                   r->k6_two_light_ray_aux) >= 0;
}

static int publish_decisions(FILE *temporary, const char *path) {
    if (fflush(temporary) != 0 || fseek(temporary, 0, SEEK_SET) != 0) {
        fprintf(stderr, "cannot rewind temporary decisions stream: %s\n",
                strerror(errno));
        return 0;
    }
    FILE *output = fopen(path, "w");
    if (!output) {
        perror(path);
        return 0;
    }
    char buffer[1 << 16];
    int ok = 1;
    size_t got;
    while ((got = fread(buffer, 1, sizeof(buffer), temporary)) != 0) {
        if (fwrite(buffer, 1, got, output) != got) {
            fprintf(stderr, "error writing decisions file %s\n", path);
            ok = 0;
            break;
        }
    }
    if (ferror(temporary)) {
        fprintf(stderr, "error reading temporary decisions stream\n");
        ok = 0;
    }
    if (fclose(output) != 0) {
        fprintf(stderr, "error closing decisions file %s\n", path);
        ok = 0;
    }
    return ok;
}

int main(int argc, char **argv) {
    if (argc < 2 || argc > 6) {
        fprintf(stderr, "usage: %s aeq_d6_n19.txt "
                "[killed_d6_n19.log|- [limit [workers [decisions.tsv]]]]\n",
                argv[0]);
        return 2;
    }
    FILE *f = fopen(argv[1], "r");
    if (!f) { perror(argv[1]); return 2; }
    uint64_t nkilled = 0;
    uint8_t *killed = argc >= 3 && strcmp(argv[2], "-") != 0 ?
                      load_killed(argv[2], &nkilled) : NULL;
    unsigned limit = argc >= 4 ? (unsigned)strtoul(argv[3], NULL, 10) : 0;
    long online = sysconf(_SC_NPROCESSORS_ONLN);
    int nworkers = argc >= 5 ? atoi(argv[4]) : (online > 0 ? (int)online : 1);
    if (nworkers < 1) nworkers = 1;
    if (nworkers > 64) nworkers = 64;
    unsigned goal = limit ? limit : EXPECTED;
    unsigned batch_capacity = goal < BATCH_GRAPHS ? goal : BATCH_GRAPHS;
    uint32_t (*graphs)[N] = malloc((size_t)batch_capacity * sizeof(*graphs));
    if (!graphs) { fprintf(stderr, "cannot allocate graph corpus\n"); return 2; }
    Worker *workers = calloc((size_t)nworkers, sizeof(*workers));
    pthread_t *threads = malloc((size_t)nworkers * sizeof(*threads));
    FILE *decision_tmp = argc >= 6 ? tmpfile() : NULL;
    ProfileResult *results = argc >= 6 ?
                             calloc(batch_capacity, sizeof(*results)) : NULL;
    if (!workers || !threads || (argc >= 6 && !results)) {
        fprintf(stderr, "worker allocation failed\n"); return 2;
    }
    if (argc >= 6 && (!decision_tmp || !write_decision_header(decision_tmp))) {
        fprintf(stderr, "cannot create temporary decisions stream\n");
        return 2;
    }
    Stats all = {0}, certified = {0}, deferred = {0};
    unsigned idx = 0;
    int reached_eof = 0, failed = 0;
    while (idx < goal && !reached_eof) {
        unsigned count = 0;
        unsigned wanted = goal - idx;
        if (wanted > batch_capacity) wanted = batch_capacity;
        while (count < wanted) {
            int n;
            int got = fscanf(f, "%d", &n);
            if (got == EOF) {
                reached_eof = 1;
                break;
            }
            if (got != 1 || n != N) {
                fprintf(stderr, "bad candidate header at index %u\n",
                        idx + count);
                failed = 1;
                break;
            }
            for (int i = 0; i < N; ++i) {
                if (fscanf(f, "%u", &graphs[count][i]) != 1) {
                    fprintf(stderr, "truncated candidate at index %u\n",
                            idx + count);
                    failed = 1;
                    break;
                }
            }
            if (failed) break;
            if (!validate(graphs[count])) {
                fprintf(stderr, "invalid graph at index %u\n", idx + count);
                failed = 1;
                break;
            }
            ++count;
        }
        if (failed || count == 0) break;

        int active_workers = nworkers;
        if (active_workers > (int)count) active_workers = (int)count;
        memset(workers, 0, (size_t)nworkers * sizeof(*workers));
        int created = 0;
        for (int t = 0; t < active_workers; ++t) {
            workers[t].graphs = graphs;
            workers[t].killed = killed;
            workers[t].results = results;
            workers[t].base = idx;
            workers[t].start = (unsigned)((uint64_t)count * (uint64_t)t /
                                          (uint64_t)active_workers);
            workers[t].end = (unsigned)((uint64_t)count * (uint64_t)(t + 1) /
                                        (uint64_t)active_workers);
            int err = pthread_create(&threads[t], NULL, profile_worker,
                                     &workers[t]);
            if (err) {
                fprintf(stderr, "pthread_create failed at worker %d: %s\n",
                        t, strerror(err));
                failed = 1;
                break;
            }
            ++created;
        }
        int join_failed = 0;
        for (int t = 0; t < created; ++t) {
            int err = pthread_join(threads[t], NULL);
            if (err) {
                fprintf(stderr, "pthread_join failed at worker %d: %s\n",
                        t, strerror(err));
                join_failed = 1;
            }
        }
        if (failed || join_failed) {
            failed = 1;
            break;
        }
        for (int t = 0; t < created; ++t) {
            stats_merge(&all, &workers[t].all);
            stats_merge(&certified, &workers[t].certified);
            stats_merge(&deferred, &workers[t].deferred);
        }
        if (decision_tmp) {
            for (unsigned i = 0; i < count; ++i) {
                if (!write_decision(decision_tmp, idx + i, &results[i])) {
                    fprintf(stderr, "error writing temporary decisions stream\n");
                    failed = 1;
                    break;
                }
            }
        }
        idx += count;
        if (failed) break;
    }
    if (!failed && !limit) {
        if (idx != EXPECTED) {
            fprintf(stderr, "candidate count %u, expected %u\n", idx,
                    EXPECTED);
            failed = 1;
        } else {
            int extra;
            int got = fscanf(f, " %d", &extra);
            if (got != EOF) {
                fprintf(stderr, "trailing data after candidate %u\n", idx - 1);
                failed = 1;
            } else if (ferror(f)) {
                fprintf(stderr, "error checking corpus end\n");
                failed = 1;
            }
        }
    }
    fclose(f);
    if (failed) {
        if (decision_tmp) fclose(decision_tmp);
        free(results); free(threads); free(workers); free(graphs);
        free(killed);
        return 2;
    }
    if (nworkers > (int)idx) nworkers = (int)idx;
    fprintf(stderr, "profiling %u candidates with %d workers\n", idx, nworkers);
    if (!limit && killed && nkilled != certified.total) {
        fprintf(stderr, "kill-log/profile mismatch: log=%" PRIu64
                " profile=%" PRIu64 "\n", nkilled, certified.total);
        return 2;
    }
    if (decision_tmp && !publish_decisions(decision_tmp, argv[5])) {
        fclose(decision_tmp);
        free(results); free(threads); free(workers); free(graphs);
        free(killed);
        return 2;
    }
    printf("{\n  \"schema\": 3,\n  \"n\": %d,\n  \"dimension\": 6,\n", N);
    printf("  \"candidate_count_expected\": %u,\n", EXPECTED);
    printf("  \"candidate_count_profiled\": %u,\n", idx);
    printf("  \"worker_count\": %d,\n", nworkers);
    printf("  \"kill_log_unique_indices\": %" PRIu64 ",\n", nkilled);
    printf("  \"populations\": {\n");
    print_stats("all", &all);
    if (killed) {
        printf(",\n"); print_stats("certified", &certified);
        printf(",\n"); print_stats("deferred", &deferred);
    }
    printf("\n  }\n}\n");
    if (decision_tmp) fclose(decision_tmp);
    free(results); free(threads); free(workers); free(graphs);
    free(killed);
    return 0;
}
