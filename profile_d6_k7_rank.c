/* Exact standalone evaluator for the post-cover K7 filters at d=6, n=19.
 *
 * Input:
 *     corpus_index 19 adj[0] ... adj[18]
 *
 * Build:
 *   cc -O3 -Wall -Wextra -Werror -pthread -o profile_d6_k7_rank \
 *      profile_d6_k7_rank.c
 * Usage:
 *   ./profile_d6_k7_rank selected.txt [workers [decisions.tsv]]
 *
 * On macOS production runs revoke inherited background priority and set
 * worker QoS.  Sandboxed tests may explicitly set
 * D6_ALLOW_BACKGROUND_TEST_ONLY=1; the evaluator prints a warning then.
 */

#include <errno.h>
#include <inttypes.h>
#include <pthread.h>
#include <stdatomic.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#ifdef __APPLE__
#include <pthread/qos.h>
#include <sys/resource.h>
#endif

#define N 19
#define OUT 12
#define COORD 7
#define LOCAL_FULL ((uint16_t)((1u << OUT) - 1u))
#define GLOBAL_FULL ((uint32_t)((1u << N) - 1u))
#define SUPPORT_CACHE_SIZE (1u << 16)
#define WORK_CHUNK 1

typedef struct {
    uint32_t index;
    uint32_t adj[N];
} Record;

typedef struct {
    uint32_t index;
    uint32_t first_cap_seed;
    uint32_t first_support_seed;
    uint32_t first_clique_seed;
    uint32_t first_degree_seed;
    uint32_t first_mask_seed;
    uint32_t first_basis_seed;
    uint32_t first_k_seed;
    uint32_t first_b_seed;
    uint32_t first_joint_seed;
    uint32_t seeds;
    uint32_t covers;
    uint32_t cover_size[4];
    uint32_t support_cover_failures;
    uint32_t clique_cover_failures;
    uint32_t degree_cover_failures;
    uint32_t mask_cover_failures;
    uint32_t basis_cover_failures;
    uint32_t k_cover_failures;
    uint32_t b_cover_failures;
    uint32_t base_joint_cover_failures;
    uint32_t joint_cover_failures;
    uint8_t applicable;
    uint8_t cap3_rejected;
    uint8_t support_rejected;
    uint8_t clique_rejected;
    uint8_t degree_rejected;
    uint8_t mask_rejected;
    uint8_t basis_rejected;
    uint8_t k_rejected;
    uint8_t b_rejected;
    uint8_t base_joint_rejected;
    uint8_t joint_rejected;
    uint8_t internal_error;
} Result;

typedef struct {
    uint32_t key;
    uint8_t value;
} SupportCacheEntry;

typedef struct {
    SupportCacheEntry *support_cache;
    uint64_t support_nodes;
    uint64_t zero_forcing_initial_sets;
} EvalContext;

typedef struct {
    const Record *records;
    Result *results;
    size_t count;
    _Atomic size_t next;
} Scheduler;

typedef struct {
    Scheduler *scheduler;
    uint64_t support_nodes;
    uint64_t zero_forcing_initial_sets;
    int qos_error;
    int failed;
} Worker;

static inline int pop16(uint16_t value) {
    return __builtin_popcount((unsigned)value);
}

static inline int pop32(uint32_t value) {
    return __builtin_popcount(value);
}

static inline int ctz16(uint16_t value) {
    return __builtin_ctz((unsigned)value);
}

static inline int ctz32(uint32_t value) {
    return __builtin_ctz(value);
}

static int validate_graph(const uint32_t adj[N]) {
    for (int u = 0; u < N; ++u) {
        if (adj[u] & ~GLOBAL_FULL) return 0;
        if (adj[u] & (1u << u)) return 0;
        for (int v = u + 1; v < N; ++v) {
            if (((adj[u] >> v) & 1u) != ((adj[v] >> u) & 1u)) return 0;
        }
    }
    /* Complement triangle-free, equivalently alpha(G)<=2. */
    for (int u = 0; u < N; ++u) {
        uint32_t non = GLOBAL_FULL & ~adj[u] & ~(1u << u);
        uint32_t todo = non;
        while (todo) {
            int v = ctz32(todo);
            todo &= todo - 1;
            if (todo & ~adj[v] & ~(1u << v)) return 0;
        }
    }
    return 1;
}

static int augment_allowed(int vertex, const uint8_t defects[OUT],
                           int owner[COORD], uint8_t *seen) {
    uint8_t choices = defects[vertex] & (uint8_t)~*seen;
    while (choices) {
        uint8_t bit = choices & (uint8_t)(-choices);
        choices ^= bit;
        int coordinate = __builtin_ctz((unsigned)bit);
        *seen |= bit;
        if (owner[coordinate] < 0 ||
                augment_allowed(owner[coordinate], defects, owner, seen)) {
            owner[coordinate] = vertex;
            return 1;
        }
    }
    return 0;
}

static int allowed_term_rank(const uint8_t defects[OUT], uint16_t selected) {
    int owner[COORD];
    for (int coordinate = 0; coordinate < COORD; ++coordinate)
        owner[coordinate] = -1;
    int rank = 0;
    while (selected) {
        int vertex = ctz16(selected);
        selected &= (uint16_t)(selected - 1);
        uint8_t seen = 0;
        if (augment_allowed(vertex, defects, owner, &seen)) ++rank;
    }
    return rank;
}

static int augment_support(int column, const uint8_t supports[3],
                           int owner[COORD], uint8_t *seen) {
    uint8_t choices = supports[column] & (uint8_t)~*seen;
    while (choices) {
        uint8_t bit = choices & (uint8_t)(-choices);
        choices ^= bit;
        int coordinate = __builtin_ctz((unsigned)bit);
        *seen |= bit;
        if (owner[coordinate] < 0 ||
                augment_support(owner[coordinate], supports, owner, seen)) {
            owner[coordinate] = column;
            return 1;
        }
    }
    return 0;
}

static int supports_match(const uint8_t supports[3], int count) {
    int owner[COORD];
    for (int coordinate = 0; coordinate < COORD; ++coordinate)
        owner[coordinate] = -1;
    for (int column = 0; column < count; ++column) {
        uint8_t seen = 0;
        if (!augment_support(column, supports, owner, &seen)) return 0;
    }
    return 1;
}

static int support_search(const uint8_t allowed[3], int count, int at,
                          uint8_t supports[3], EvalContext *context) {
    ++context->support_nodes;
    if (at == count) return supports_match(supports, count);
    for (uint8_t support = allowed[at]; support;
         support = (uint8_t)((support - 1u) & allowed[at])) {
        if (pop16(support) < 3) continue;
        int compatible = 1;
        for (int previous = 0; previous < at; ++previous) {
            if (pop16((uint16_t)(support & supports[previous])) == 1) {
                compatible = 0;
                break;
            }
        }
        if (!compatible) continue;
        supports[at] = support;
        if (support_search(allowed, count, at + 1, supports, context))
            return 1;
    }
    return 0;
}

static int support_possible(const uint8_t defects[OUT], uint16_t zmask,
                            EvalContext *context) {
    int count = pop16(zmask);
    if (count == 0) return 1;
    uint8_t allowed[3] = {0, 0, 0};
    int at = 0;
    while (zmask) {
        int vertex = ctz16(zmask);
        zmask &= (uint16_t)(zmask - 1);
        allowed[at++] = defects[vertex];
    }
    for (int i = 1; i < count; ++i) {
        uint8_t value = allowed[i];
        int j = i;
        while (j > 0 && allowed[j - 1] > value) {
            allowed[j] = allowed[j - 1];
            --j;
        }
        allowed[j] = value;
    }
    uint32_t packed = (uint32_t)count;
    for (int i = 0; i < count; ++i)
        packed |= (uint32_t)allowed[i] << (2 + 7 * i);
    uint32_t slot = (packed * UINT32_C(2654435761)) >> 16;
    SupportCacheEntry *entry = &context->support_cache[slot];
    if (entry->key == packed) return entry->value == 1;
    uint8_t supports[3] = {0, 0, 0};
    int answer = support_search(allowed, count, 0, supports, context);
    entry->key = packed;
    entry->value = (uint8_t)(answer ? 1 : 2);
    return answer;
}

static uint16_t forcing_closure(const uint16_t adj[OUT], uint16_t active,
                                uint16_t black) {
    for (;;) {
        int changed = 0;
        uint16_t sources = black;
        while (sources) {
            int vertex = ctz16(sources);
            sources &= (uint16_t)(sources - 1);
            uint16_t white = adj[vertex] & active & (uint16_t)~black;
            if (white && !(white & (uint16_t)(white - 1))) {
                black |= white;
                changed = 1;
                break;
            }
        }
        if (!changed) return black;
    }
}

static int forcing_combination(const uint16_t adj[OUT], uint16_t active,
                               uint16_t available, uint16_t chosen, int need,
                               EvalContext *context) {
    if (need == 0) {
        ++context->zero_forcing_initial_sets;
        return forcing_closure(adj, active, chosen) == active;
    }
    while (pop16(available) >= need) {
        uint16_t bit = available & (uint16_t)(-available);
        available ^= bit;
        if (forcing_combination(adj, active, available,
                                (uint16_t)(chosen | bit), need - 1,
                                context)) return 1;
    }
    return 0;
}

static int zero_forcing_number(const uint16_t adj[OUT], uint16_t active,
                               EvalContext *context) {
    int count = pop16(active);
    for (int size = 0; size <= count; ++size) {
        if (forcing_combination(adj, active, active, 0, size, context))
            return size;
    }
    return count;
}

static uint16_t component_from(const uint16_t adj[OUT], uint16_t active,
                               uint16_t root) {
    uint16_t component = 0;
    uint16_t frontier = root;
    while (frontier) {
        uint16_t bit = frontier & (uint16_t)(-frontier);
        frontier ^= bit;
        if (component & bit) continue;
        component |= bit;
        int vertex = ctz16(bit);
        frontier |= adj[vertex] & active & (uint16_t)~component;
    }
    return component;
}

static int component_nullity_cap(const uint16_t fadj[OUT], uint16_t active,
                                 EvalContext *context) {
    uint16_t remaining = active;
    int components = 0;
    int maximum_component_zf = 0;
    while (remaining) {
        uint16_t root = remaining & (uint16_t)(-remaining);
        uint16_t component = component_from(fadj, remaining, root);
        remaining &= (uint16_t)~component;
        if (pop16(component) <= 1) continue;
        ++components;
        int zf = zero_forcing_number(fadj, component, context);
        if (zf > maximum_component_zf) maximum_component_zf = zf;
    }
    if (components == 0) return 0;
    int indefinite = components - 1 + maximum_component_zf;
    return components > indefinite ? components : indefinite;
}

static int is_vertex_cover(const uint16_t ladj[OUT], uint16_t zmask) {
    uint16_t remaining = LOCAL_FULL & (uint16_t)~zmask;
    uint16_t todo = remaining;
    while (todo) {
        int vertex = ctz16(todo);
        todo &= (uint16_t)(todo - 1);
        if (ladj[vertex] & remaining) return 0;
    }
    return 1;
}

static int has_clique_local(const uint16_t adj[OUT], uint16_t candidates,
                            int need) {
    if (need == 0) return 1;
    if (pop16(candidates) < need) return 0;
    while (candidates) {
        int vertex = ctz16(candidates);
        candidates &= (uint16_t)(candidates - 1);
        if (has_clique_local(adj, candidates & adj[vertex], need - 1))
            return 1;
        if (pop16(candidates) < need) return 0;
    }
    return 0;
}

/* Equality case for a size-U_K clique C.  The fixed 0/1 cross Gram
 * columns are the G-neighbourhood masks in C.  Positivity in the exact
 * inverse of J+diag(r_c^2) forces these mask conditions. */
static int equality_mask_clique_compatible(const uint16_t adj[OUT],
                                           uint16_t active,
                                           uint16_t clique) {
    uint16_t remaining = active & (uint16_t)~clique;
    uint16_t todo = remaining;
    while (todo) {
        int vertex = ctz16(todo);
        todo &= (uint16_t)(todo - 1);
        uint16_t mask = adj[vertex] & clique;
        if (mask == 0 || mask == clique) return 0;
    }
    todo = remaining;
    while (todo) {
        int first = ctz16(todo);
        todo &= (uint16_t)(todo - 1);
        uint16_t first_mask = adj[first] & clique;
        uint16_t later = todo;
        while (later) {
            int second = ctz16(later);
            later &= (uint16_t)(later - 1);
            uint16_t second_mask = adj[second] & clique;
            uint16_t intersection = first_mask & second_mask;
            if (first_mask == second_mask || intersection == 0) return 0;
            if (adj[first] & (uint16_t)(1u << second)) {
                if ((uint16_t)(first_mask | second_mask) == clique)
                    return 0;
            } else if (intersection == first_mask ||
                       intersection == second_mask) {
                return 0;
            }
        }
    }
    return 1;
}

static int equality_mask_cliques_compatible(const uint16_t adj[OUT],
                                            uint16_t active,
                                            uint16_t candidates,
                                            uint16_t clique, int need) {
    if (need == 0)
        return equality_mask_clique_compatible(adj, active, clique);
    if (pop16(candidates) < need) return 1;
    while (candidates) {
        uint16_t bit = candidates & (uint16_t)(-candidates);
        int vertex = ctz16(bit);
        candidates ^= bit;
        if (!equality_mask_cliques_compatible(
                adj, active, candidates & adj[vertex],
                (uint16_t)(clique | bit), need - 1))
            return 0;
        if (pop16(candidates) < need) return 1;
    }
    return 1;
}

typedef struct {
    int64_t numerator;
    int64_t denominator;
} Rational;

static unsigned __int128 abs128(__int128 value) {
    if (value >= 0) return (unsigned __int128)value;
    return (unsigned __int128)(-(value + 1)) + 1;
}

static unsigned __int128 gcd128(unsigned __int128 first,
                                unsigned __int128 second) {
    while (second) {
        unsigned __int128 remainder = first % second;
        first = second;
        second = remainder;
    }
    return first;
}

static Rational rational_make(__int128 numerator, __int128 denominator,
                              int *ok) {
    Rational zero = {0, 1};
    if (!*ok) return zero;
    if (denominator == 0) {
        *ok = 0;
        return zero;
    }
    if (numerator == 0) return zero;
    if (denominator < 0) {
        numerator = -numerator;
        denominator = -denominator;
    }
    unsigned __int128 divisor = gcd128(abs128(numerator),
                                       (unsigned __int128)denominator);
    numerator /= (__int128)divisor;
    denominator /= (__int128)divisor;
    if (numerator > INT64_MAX || numerator < INT64_MIN ||
        denominator > INT64_MAX) {
        *ok = 0;
        return zero;
    }
    Rational answer = {(int64_t)numerator, (int64_t)denominator};
    return answer;
}

static Rational rational_add(Rational first, Rational second, int *ok) {
    return rational_make(
        (__int128)first.numerator * second.denominator +
            (__int128)second.numerator * first.denominator,
        (__int128)first.denominator * second.denominator, ok);
}

static Rational rational_subtract(Rational first, Rational second,
                                   int *ok) {
    return rational_make(
        (__int128)first.numerator * second.denominator -
            (__int128)second.numerator * first.denominator,
        (__int128)first.denominator * second.denominator, ok);
}

static Rational rational_multiply(Rational first, Rational second,
                                   int *ok) {
    return rational_make(
        (__int128)first.numerator * second.numerator,
        (__int128)first.denominator * second.denominator, ok);
}

static Rational rational_divide(Rational first, Rational second, int *ok) {
    return rational_make(
        (__int128)first.numerator * second.denominator,
        (__int128)first.denominator * second.numerator, ok);
}

static Rational rational_negate(Rational value, int *ok) {
    return rational_make(-(__int128)value.numerator, value.denominator, ok);
}

static int rational_equal(Rational first, Rational second) {
    return first.numerator == second.numerator &&
           first.denominator == second.denominator;
}

/* Return 1 for compatible, 0 for an exact incompatibility, and -1 only
 * if the checked rational arithmetic exceeds its deliberately wide bounds. */
static int basis_kernel_clique_compatible(const uint16_t adj[OUT],
                                          uint16_t active,
                                          uint16_t clique) {
    int clique_vertices[COORD];
    int clique_size = 0;
    uint16_t todo = clique;
    while (todo) {
        clique_vertices[clique_size++] = ctz16(todo);
        todo &= (uint16_t)(todo - 1);
    }
    int remainder_vertices[OUT];
    int remainder_size = 0;
    todo = active & (uint16_t)~clique;
    while (todo) {
        remainder_vertices[remainder_size++] = ctz16(todo);
        todo &= (uint16_t)(todo - 1);
    }
    if (remainder_size == 0) return 1;

    Rational matrix[COORD][OUT];
    Rational zero = {0, 1};
    Rational one = {1, 1};
    for (int row = 0; row < COORD; ++row)
        for (int column = 0; column < OUT; ++column)
            matrix[row][column] = zero;
    for (int row = 0; row < clique_size; ++row) {
        for (int column = 0; column < remainder_size; ++column) {
            if (adj[clique_vertices[row]] &
                (uint16_t)(1u << remainder_vertices[column]))
                matrix[row][column] = one;
        }
    }

    int ok = 1;
    int pivots[COORD];
    int rank = 0;
    for (int column = 0; column < remainder_size && rank < clique_size;
         ++column) {
        int selected = rank;
        while (selected < clique_size &&
               matrix[selected][column].numerator == 0)
            ++selected;
        if (selected == clique_size) continue;
        if (selected != rank) {
            for (int other = 0; other < remainder_size; ++other) {
                Rational swap = matrix[rank][other];
                matrix[rank][other] = matrix[selected][other];
                matrix[selected][other] = swap;
            }
        }
        Rational pivot = matrix[rank][column];
        for (int other = 0; other < remainder_size; ++other)
            matrix[rank][other] = rational_divide(
                matrix[rank][other], pivot, &ok);
        for (int row = 0; row < clique_size; ++row) {
            if (row == rank || matrix[row][column].numerator == 0) continue;
            Rational multiplier = matrix[row][column];
            for (int other = 0; other < remainder_size; ++other) {
                Rational product = rational_multiply(
                    multiplier, matrix[rank][other], &ok);
                matrix[row][other] = rational_subtract(
                    matrix[row][other], product, &ok);
            }
        }
        if (!ok) return -1;
        pivots[rank++] = column;
    }
    int nullity = remainder_size - rank;
    if (nullity == 0) return 1;

    Rational kernel[OUT][OUT];
    for (int row = 0; row < OUT; ++row)
        for (int column = 0; column < OUT; ++column)
            kernel[row][column] = zero;
    int free_at = 0;
    for (int column = 0; column < remainder_size; ++column) {
        int is_pivot = 0;
        for (int row = 0; row < rank; ++row) {
            if (pivots[row] == column) {
                is_pivot = 1;
                break;
            }
        }
        if (is_pivot) continue;
        kernel[column][free_at] = one;
        for (int row = 0; row < rank; ++row) {
            Rational value = rational_negate(matrix[row][column], &ok);
            kernel[pivots[row]][free_at] = value;
        }
        ++free_at;
    }
    if (free_at != nullity) return -1;

    for (int row = 0; row < remainder_size; ++row) {
        Rational action[OUT];
        for (int column = 0; column < nullity; ++column) {
            action[column] = zero;
            for (int other = 0; other < remainder_size; ++other) {
                if (adj[remainder_vertices[row]] &
                    (uint16_t)(1u << remainder_vertices[other])) {
                    action[column] = rational_add(
                        action[column], kernel[other][column], &ok);
                }
            }
        }
        if (!ok) return -1;
        int nonzero = 0;
        while (nonzero < nullity &&
               kernel[row][nonzero].numerator == 0)
            ++nonzero;
        if (nonzero == nullity) {
            for (int column = 0; column < nullity; ++column)
                if (action[column].numerator != 0) return 0;
            continue;
        }
        Rational negative_action = rational_negate(action[nonzero], &ok);
        Rational diagonal = rational_divide(
            negative_action, kernel[row][nonzero], &ok);
        if (!ok) return -1;
        if (diagonal.numerator <= diagonal.denominator) return 0;
        for (int column = 0; column < nullity; ++column) {
            Rational expected = rational_multiply(
                diagonal, kernel[row][column], &ok);
            Rational target = rational_negate(action[column], &ok);
            if (!ok) return -1;
            if (!rational_equal(expected, target)) return 0;
        }
    }
    return 1;
}

static int basis_kernel_cliques_compatible(const uint16_t adj[OUT],
                                           uint16_t active,
                                           uint16_t candidates,
                                           uint16_t clique, int need) {
    if (need == 0)
        return basis_kernel_clique_compatible(adj, active, clique);
    if (pop16(candidates) < need) return 1;
    while (candidates) {
        uint16_t bit = candidates & (uint16_t)(-candidates);
        int vertex = ctz16(bit);
        candidates ^= bit;
        int answer = basis_kernel_cliques_compatible(
            adj, active, candidates & adj[vertex],
            (uint16_t)(clique | bit), need - 1);
        if (answer <= 0) return answer;
        if (pop16(candidates) < need) return 1;
    }
    return 1;
}

typedef struct {
    uint32_t covers;
    uint32_t cover_size[4];
    uint32_t support_failures;
    uint32_t clique_failures;
    uint32_t degree_failures;
    uint32_t mask_failures;
    uint32_t basis_failures;
    uint32_t k_failures;
    uint32_t b_failures;
    uint32_t base_joint_failures;
    uint32_t joint_failures;
    int support_pass;
    int clique_pass;
    int degree_pass;
    int mask_pass;
    int basis_pass;
    int k_pass;
    int b_pass;
    int base_joint_pass;
    int joint_pass;
    int internal_error;
} SeedResult;

static SeedResult evaluate_seed(const uint32_t adj[N], uint32_t seed,
                                EvalContext *context) {
    SeedResult result = {0};
    int seed_vertices[COORD];
    int seed_count = 0;
    for (int vertex = 0; vertex < N; ++vertex) {
        if (seed & (1u << vertex)) seed_vertices[seed_count++] = vertex;
    }
    if (seed_count != COORD) {
        result.internal_error = 1;
        return result;
    }
    int outside_vertices[OUT];
    int outside_count = 0;
    for (int vertex = 0; vertex < N; ++vertex) {
        if (!(seed & (1u << vertex)))
            outside_vertices[outside_count++] = vertex;
    }
    if (outside_count != OUT) {
        result.internal_error = 1;
        return result;
    }

    uint8_t defects[OUT] = {0};
    uint16_t gadj[OUT] = {0};
    uint16_t ladj[OUT] = {0};
    uint16_t eligible = 0;
    for (int i = 0; i < OUT; ++i) {
        for (int coordinate = 0; coordinate < COORD; ++coordinate) {
            if (!(adj[outside_vertices[i]] &
                  (1u << seed_vertices[coordinate]))) {
                defects[i] |= (uint8_t)(1u << coordinate);
            }
        }
        if (pop16(defects[i]) >= 3)
            eligible |= (uint16_t)(1u << i);
    }
    for (int i = 0; i < OUT; ++i) {
        for (int j = i + 1; j < OUT; ++j) {
            if (adj[outside_vertices[i]] & (1u << outside_vertices[j])) {
                gadj[i] |= (uint16_t)(1u << j);
                gadj[j] |= (uint16_t)(1u << i);
                if (!(defects[i] & defects[j])) {
                    ladj[i] |= (uint16_t)(1u << j);
                    ladj[j] |= (uint16_t)(1u << i);
                }
            }
        }
    }
    int nu_t = allowed_term_rank(defects, LOCAL_FULL);

    for (uint16_t zmask = 0; zmask <= LOCAL_FULL; ++zmask) {
        int zsize = pop16(zmask);
        if (zsize > 3 || (zmask & (uint16_t)~eligible)) continue;
        if (!is_vertex_cover(ladj, zmask)) continue;
        ++result.covers;
        ++result.cover_size[zsize];
        uint16_t active = LOCAL_FULL & (uint16_t)~zmask;

        /* Cover plus alpha<=2 makes G[N] exactly the mask-intersection
         * graph.  Keep this as a production invariant check. */
        uint16_t pairs = active;
        while (pairs) {
            int i = ctz16(pairs);
            pairs &= (uint16_t)(pairs - 1);
            uint16_t later = pairs;
            while (later) {
                int j = ctz16(later);
                later &= (uint16_t)(later - 1);
                int intersects = (defects[i] & defects[j]) != 0;
                int edge = (gadj[i] & (1u << j)) != 0;
                if (intersects != edge) result.internal_error = 1;
            }
        }

        int support_failed = !support_possible(defects, zmask, context);
        int nu_n = allowed_term_rank(defects, active);
        int uk = COORD - zsize;
        if (nu_n < uk) uk = nu_n;
        int combined = nu_t - zsize;
        if (combined < 0) combined = 0;
        if (combined < uk) uk = combined;
        int nsize = pop16(active);

        /* Every required clique principal K block is positive definite. */
        int clique_failed = uk < nsize &&
                            has_clique_local(gadj, active, uk + 1);

        uint16_t fadj[OUT] = {0};
        int maximum_f_degree = 0;
        for (int i = 0; i < OUT; ++i) {
            if (active & (1u << i)) {
                fadj[i] = active & (uint16_t)~gadj[i] &
                           (uint16_t)~(1u << i);
                int degree = pop16(fadj[i]);
                if (degree > maximum_f_degree) maximum_f_degree = degree;
            }
        }
        /* F-neighbors form a required clique in v_i^perp. */
        int degree_failed = nsize > 0 && maximum_f_degree > uk - 1;

        int mask_failed = uk > 0 && uk <= nsize &&
            !equality_mask_cliques_compatible(gadj, active, active, 0, uk);
        int basis_state = uk > 0 && uk <= nsize ?
            basis_kernel_cliques_compatible(gadj, active, active, 0, uk) : 1;
        if (basis_state < 0) result.internal_error = 1;
        int basis_failed = basis_state == 0;

        int k_zf = zero_forcing_number(gadj, active, context);
        int k_failed = nsize - k_zf > uk;

        int b_upper = uk + 1;
        if (b_upper > nsize) b_upper = nsize;
        if (b_upper > 8 - zsize) b_upper = 8 - zsize;
        int nullity_cap = component_nullity_cap(fadj, active, context);
        int b_failed = nsize - nullity_cap > b_upper;
        int base_joint_failed = support_failed || k_failed || b_failed;
        int joint_failed = base_joint_failed || clique_failed ||
                           degree_failed || mask_failed || basis_failed;

        if (support_failed) ++result.support_failures;
        else result.support_pass = 1;
        if (clique_failed) ++result.clique_failures;
        else result.clique_pass = 1;
        if (degree_failed) ++result.degree_failures;
        else result.degree_pass = 1;
        if (mask_failed) ++result.mask_failures;
        else result.mask_pass = 1;
        if (basis_failed) ++result.basis_failures;
        else result.basis_pass = 1;
        if (k_failed) ++result.k_failures;
        else result.k_pass = 1;
        if (b_failed) ++result.b_failures;
        else result.b_pass = 1;
        if (base_joint_failed) ++result.base_joint_failures;
        else result.base_joint_pass = 1;
        if (joint_failed) ++result.joint_failures;
        else result.joint_pass = 1;
    }
    return result;
}

static void scan_k7_seeds(const uint32_t adj[N], uint32_t candidates,
                          uint32_t seed, int depth, Result *graph,
                          EvalContext *context) {
    if (depth == COORD) {
        graph->applicable = 1;
        ++graph->seeds;
        SeedResult current = evaluate_seed(adj, seed, context);
        graph->covers += current.covers;
        for (int size = 0; size <= 3; ++size)
            graph->cover_size[size] += current.cover_size[size];
        graph->support_cover_failures += current.support_failures;
        graph->clique_cover_failures += current.clique_failures;
        graph->degree_cover_failures += current.degree_failures;
        graph->mask_cover_failures += current.mask_failures;
        graph->basis_cover_failures += current.basis_failures;
        graph->k_cover_failures += current.k_failures;
        graph->b_cover_failures += current.b_failures;
        graph->base_joint_cover_failures += current.base_joint_failures;
        graph->joint_cover_failures += current.joint_failures;
        graph->internal_error |= (uint8_t)current.internal_error;
        if (current.covers == 0) {
            if (!graph->cap3_rejected) graph->first_cap_seed = seed;
            graph->cap3_rejected = 1;
        } else {
            if (!current.support_pass) {
                if (!graph->support_rejected)
                    graph->first_support_seed = seed;
                graph->support_rejected = 1;
            }
            if (!current.clique_pass) {
                if (!graph->clique_rejected)
                    graph->first_clique_seed = seed;
                graph->clique_rejected = 1;
            }
            if (!current.degree_pass) {
                if (!graph->degree_rejected)
                    graph->first_degree_seed = seed;
                graph->degree_rejected = 1;
            }
            if (!current.mask_pass) {
                if (!graph->mask_rejected)
                    graph->first_mask_seed = seed;
                graph->mask_rejected = 1;
            }
            if (!current.basis_pass) {
                if (!graph->basis_rejected)
                    graph->first_basis_seed = seed;
                graph->basis_rejected = 1;
            }
            if (!current.k_pass) {
                if (!graph->k_rejected) graph->first_k_seed = seed;
                graph->k_rejected = 1;
            }
            if (!current.b_pass) {
                if (!graph->b_rejected) graph->first_b_seed = seed;
                graph->b_rejected = 1;
            }
        }
        if (current.covers == 0 || !current.base_joint_pass)
            graph->base_joint_rejected = 1;
        if (current.covers == 0 || !current.joint_pass) {
            if (!graph->joint_rejected) graph->first_joint_seed = seed;
            graph->joint_rejected = 1;
        }
        return;
    }
    if (pop32(candidates) < COORD - depth) return;
    while (candidates) {
        int vertex = ctz32(candidates);
        candidates &= candidates - 1;
        scan_k7_seeds(adj, candidates & adj[vertex],
                      seed | (1u << vertex), depth + 1, graph, context);
        if (pop32(candidates) < COORD - depth) return;
    }
}

static Result evaluate_graph(const Record *record, EvalContext *context) {
    Result result = {0};
    result.index = record->index;
    scan_k7_seeds(record->adj, GLOBAL_FULL, 0, 0, &result, context);
    return result;
}

static void *worker_main(void *opaque) {
    Worker *worker = opaque;
#ifdef __APPLE__
    worker->qos_error = pthread_set_qos_class_self_np(
        QOS_CLASS_USER_INITIATED, 0);
    if (worker->qos_error) {
        worker->failed = 1;
        return NULL;
    }
#endif
    EvalContext context = {0};
    context.support_cache = calloc(SUPPORT_CACHE_SIZE,
                                   sizeof(*context.support_cache));
    if (!context.support_cache) {
        worker->failed = 1;
        return NULL;
    }
    for (;;) {
        size_t start = atomic_fetch_add_explicit(
            &worker->scheduler->next, WORK_CHUNK, memory_order_relaxed);
        if (start >= worker->scheduler->count) break;
        size_t end = start + WORK_CHUNK;
        if (end > worker->scheduler->count) end = worker->scheduler->count;
        for (size_t i = start; i < end; ++i) {
            worker->scheduler->results[i] = evaluate_graph(
                &worker->scheduler->records[i], &context);
        }
    }
    worker->support_nodes = context.support_nodes;
    worker->zero_forcing_initial_sets = context.zero_forcing_initial_sets;
    free(context.support_cache);
    return NULL;
}

static double monotonic_seconds(void) {
    struct timespec now;
    if (clock_gettime(CLOCK_MONOTONIC, &now) != 0) return 0.0;
    return (double)now.tv_sec + (double)now.tv_nsec * 1e-9;
}

static int write_decisions(const char *path, const Result *results,
                           size_t count) {
    FILE *stream = fopen(path, "w");
    if (!stream) {
        perror(path);
        return 0;
    }
    fprintf(stream,
            "index\tseeds\tcovers\tcover0\tcover1\tcover2\tcover3\t"
            "support_cover_failures\tclique_cover_failures\t"
            "degree_cover_failures\tmask_cover_failures\t"
            "basis_kernel_cover_failures\t"
            "K_cover_failures\tB_cover_failures\t"
            "base_joint_cover_failures\tjoint_cover_failures\tcap3_rejected\t"
            "support_rejected\t"
            "clique_rejected\tdegree_rejected\tmask_rejected\t"
            "basis_kernel_rejected\t"
            "K_rejected\tB_rejected\t"
            "base_joint_rejected\tjoint_rejected\tfirst_cap_seed\t"
            "first_support_seed\t"
            "first_clique_seed\tfirst_degree_seed\tfirst_mask_seed\t"
            "first_basis_kernel_seed\t"
            "first_K_seed\t"
            "first_B_seed\tfirst_joint_seed\tinternal_error\n");
    for (size_t i = 0; i < count; ++i) {
        const Result *r = &results[i];
        fprintf(stream,
                "%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t"
                "%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t"
                "%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\n",
                r->index, r->seeds, r->covers,
                r->cover_size[0], r->cover_size[1], r->cover_size[2],
                r->cover_size[3], r->support_cover_failures,
                r->clique_cover_failures, r->degree_cover_failures,
                r->mask_cover_failures,
                r->basis_cover_failures,
                r->k_cover_failures, r->b_cover_failures,
                r->base_joint_cover_failures, r->joint_cover_failures,
                r->cap3_rejected,
                r->support_rejected, r->clique_rejected,
                r->degree_rejected, r->mask_rejected,
                r->basis_rejected,
                r->k_rejected, r->b_rejected,
                r->base_joint_rejected, r->joint_rejected, r->first_cap_seed,
                r->first_support_seed, r->first_clique_seed,
                r->first_degree_seed, r->first_mask_seed,
                r->first_basis_seed, r->first_k_seed,
                r->first_b_seed, r->first_joint_seed, r->internal_error);
    }
    if (ferror(stream) || fclose(stream) != 0) {
        fprintf(stderr, "error writing decisions file %s\n", path);
        return 0;
    }
    return 1;
}

int main(int argc, char **argv) {
    if (argc < 2 || argc > 4) {
        fprintf(stderr, "usage: %s selected.txt [workers [decisions.tsv]]\n",
                argv[0]);
        return 2;
    }
    FILE *input = fopen(argv[1], "r");
    if (!input) {
        perror(argv[1]);
        return 2;
    }
    size_t capacity = 1024;
    size_t count = 0;
    Record *records = malloc(capacity * sizeof(*records));
    if (!records) {
        fprintf(stderr, "record allocation failed\n");
        return 2;
    }
    for (;;) {
        uint32_t index;
        int n;
        int got = fscanf(input, "%u %d", &index, &n);
        if (got == EOF) break;
        if (got != 2 || n != N) {
            fprintf(stderr, "bad input header at record %zu\n", count);
            return 2;
        }
        if (count == capacity) {
            size_t next = capacity * 2;
            Record *grown = realloc(records, next * sizeof(*records));
            if (!grown) {
                fprintf(stderr, "record allocation failed\n");
                return 2;
            }
            records = grown;
            capacity = next;
        }
        records[count].index = index;
        for (int vertex = 0; vertex < N; ++vertex) {
            if (fscanf(input, "%u", &records[count].adj[vertex]) != 1) {
                fprintf(stderr, "truncated input at record %zu\n", count);
                return 2;
            }
        }
        if (!validate_graph(records[count].adj)) {
            fprintf(stderr, "invalid graph at corpus index %u\n", index);
            return 2;
        }
        ++count;
    }
    if (ferror(input)) {
        fprintf(stderr, "error reading %s\n", argv[1]);
        return 2;
    }
    fclose(input);
    if (count == 0) {
        fprintf(stderr, "input contains no graphs\n");
        return 2;
    }

#ifdef __APPLE__
    if (setpriority(PRIO_DARWIN_PROCESS, getpid(), 0) != 0) {
        const char *test_override = getenv("D6_ALLOW_BACKGROUND_TEST_ONLY");
        if (!test_override || strcmp(test_override, "1") != 0) {
            perror("could not revoke inherited Darwin background priority");
            return 2;
        }
        fprintf(stderr,
                "WARNING: test-only override retained Darwin background "
                "priority\n");
    }
#endif

    long online = sysconf(_SC_NPROCESSORS_ONLN);
    int workers_requested = argc >= 3 ? atoi(argv[2]) :
                            (online > 0 ? (int)online : 1);
    if (workers_requested < 1) workers_requested = 1;
    if (workers_requested > 64) workers_requested = 64;
    if ((size_t)workers_requested > count) workers_requested = (int)count;
    Result *results = calloc(count, sizeof(*results));
    Worker *workers = calloc((size_t)workers_requested, sizeof(*workers));
    pthread_t *threads = malloc((size_t)workers_requested * sizeof(*threads));
    if (!results || !workers || !threads) {
        fprintf(stderr, "worker allocation failed\n");
        return 2;
    }
    Scheduler scheduler = {
        .records = records,
        .results = results,
        .count = count,
        .next = 0,
    };

    double started = monotonic_seconds();
    int created = 0;
    for (int worker = 0; worker < workers_requested; ++worker) {
        workers[worker].scheduler = &scheduler;
        int error = pthread_create(&threads[worker], NULL, worker_main,
                                   &workers[worker]);
        if (error) {
            fprintf(stderr, "pthread_create failed: %s\n", strerror(error));
            break;
        }
        ++created;
    }
    int failed = created != workers_requested;
    for (int worker = 0; worker < created; ++worker) {
        int error = pthread_join(threads[worker], NULL);
        if (error) {
            fprintf(stderr, "pthread_join failed: %s\n", strerror(error));
            failed = 1;
        }
        if (workers[worker].failed) failed = 1;
        if (workers[worker].qos_error) {
            fprintf(stderr,
                    "worker %d could not set user-initiated QoS: %s\n",
                    worker, strerror(workers[worker].qos_error));
        }
    }
    double elapsed = monotonic_seconds() - started;
    if (failed) return 2;

    uint64_t seeds = 0, covers = 0, support_failures = 0;
    uint64_t clique_failures = 0, degree_failures = 0;
    uint64_t mask_failures = 0;
    uint64_t basis_failures = 0;
    uint64_t k_failures = 0, b_failures = 0, base_joint_failures = 0;
    uint64_t joint_failures = 0;
    uint64_t cover_size[4] = {0, 0, 0, 0};
    uint64_t applicable = 0, cap_rejected = 0, support_rejected = 0;
    uint64_t clique_rejected = 0, degree_rejected = 0;
    uint64_t mask_rejected = 0;
    uint64_t basis_rejected = 0;
    uint64_t k_rejected = 0, b_rejected = 0, base_joint_rejected = 0;
    uint64_t joint_rejected = 0;
    uint64_t internal_errors = 0, support_nodes = 0, zf_initials = 0;
    for (size_t i = 0; i < count; ++i) {
        const Result *r = &results[i];
        seeds += r->seeds;
        covers += r->covers;
        for (int size = 0; size <= 3; ++size)
            cover_size[size] += r->cover_size[size];
        support_failures += r->support_cover_failures;
        clique_failures += r->clique_cover_failures;
        degree_failures += r->degree_cover_failures;
        mask_failures += r->mask_cover_failures;
        basis_failures += r->basis_cover_failures;
        k_failures += r->k_cover_failures;
        b_failures += r->b_cover_failures;
        base_joint_failures += r->base_joint_cover_failures;
        joint_failures += r->joint_cover_failures;
        applicable += r->applicable;
        cap_rejected += r->cap3_rejected;
        support_rejected += r->support_rejected;
        clique_rejected += r->clique_rejected;
        degree_rejected += r->degree_rejected;
        mask_rejected += r->mask_rejected;
        basis_rejected += r->basis_rejected;
        k_rejected += r->k_rejected;
        b_rejected += r->b_rejected;
        base_joint_rejected += r->base_joint_rejected;
        joint_rejected += r->joint_rejected;
        internal_errors += r->internal_error;
    }
    for (int worker = 0; worker < workers_requested; ++worker) {
        support_nodes += workers[worker].support_nodes;
        zf_initials += workers[worker].zero_forcing_initial_sets;
    }
    if (internal_errors) {
        fprintf(stderr, "internal exact-pattern checks failed on %" PRIu64
                " graphs\n", internal_errors);
        return 2;
    }
    if (argc >= 4 && !write_decisions(argv[3], results, count)) return 2;

    printf("{\n");
    printf("  \"schema\": 1,\n");
    printf("  \"filter\": \"K7_cap3_support_rank_equality\",\n");
    printf("  \"graphs\": %zu,\n", count);
    printf("  \"workers\": %d,\n", workers_requested);
    printf("  \"wall_seconds\": %.9f,\n", elapsed);
    printf("  \"applicable_K7_graphs\": %" PRIu64 ",\n", applicable);
    printf("  \"seeds_checked\": %" PRIu64 ",\n", seeds);
    printf("  \"eligible_cap3_covers_checked\": %" PRIu64 ",\n", covers);
    printf("  \"cover_size_histogram\": {\"0\": %" PRIu64
           ", \"1\": %" PRIu64 ", \"2\": %" PRIu64
           ", \"3\": %" PRIu64 "},\n",
           cover_size[0], cover_size[1], cover_size[2], cover_size[3]);
    printf("  \"cover_failures\": {\"support\": %" PRIu64
           ", \"clique_K\": %" PRIu64 ", \"degree_F\": %" PRIu64
           ", \"equality_masks\": %" PRIu64
           ", \"basis_kernel\": %" PRIu64
           ", \"subspace_K\": %" PRIu64 ", \"component_B\": %" PRIu64
           ", \"base_joint\": %" PRIu64 ", \"enhanced_joint\": %" PRIu64
           "},\n",
           support_failures, clique_failures, degree_failures, mask_failures,
           basis_failures, k_failures, b_failures, base_joint_failures,
           joint_failures);
    printf("  \"graph_rejections\": {\"cap3\": %" PRIu64
           ", \"support_individual_post_cap3\": %" PRIu64
           ", \"clique_K_individual_post_cap3\": %" PRIu64
           ", \"degree_F_individual_post_cap3\": %" PRIu64
           ", \"equality_masks_individual_post_cap3\": %" PRIu64
           ", \"basis_kernel_individual_post_cap3\": %" PRIu64
           ", \"subspace_K_individual_post_cap3\": %" PRIu64
           ", \"component_B_individual_post_cap3\": %" PRIu64
           ", \"base_joint_existential\": %" PRIu64
           ", \"enhanced_joint_existential\": %" PRIu64 "},\n",
           cap_rejected, support_rejected, clique_rejected, degree_rejected,
           mask_rejected, basis_rejected, k_rejected, b_rejected,
           base_joint_rejected, joint_rejected);
    printf("  \"survivors\": %zu,\n", count - (size_t)joint_rejected);
    printf("  \"support_CSP_nodes\": %" PRIu64 ",\n", support_nodes);
    printf("  \"zero_forcing_initial_sets_checked\": %" PRIu64 "\n",
           zf_initials);
    printf("}\n");

    free(threads);
    free(workers);
    free(results);
    free(records);
    return 0;
}
