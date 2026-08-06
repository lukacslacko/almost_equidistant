/* Exact non-induced subgraph monomorphism kernel for the 13-vertex
 * pattern-954 obstruction core in 19-vertex targets.
 *
 * Pattern edges must map to target edges. Pattern nonedges impose no
 * condition. The search uses only degree bounds, edge arc consistency,
 * all-different matching, and exhaustive branching. A node-limit exhaustion
 * is TIMEOUT, never NO_HIT. Every HIT returns a full mapping for independent
 * edge-by-edge verification.
 *
 * Build on POSIX:
 *   cc -O3 -std=c11 -Wall -Wextra -pedantic -shared \
 *      -o d6_n14_pattern_954_containment.dylib \
 *      d6_n14_pattern_954_containment.c
 */

#include <stdint.h>
#include <string.h>

#ifdef _WIN32
#define EXPORT __declspec(dllexport)
#else
#define EXPORT __attribute__((visibility("default")))
#endif

#define NP 13
#define NT 19

enum MatchResult {
    MATCH_NO_HIT = 0,
    MATCH_HIT = 1,
    MATCH_TIMEOUT = 2,
    MATCH_INVALID = 3
};

typedef struct {
    const uint32_t *pattern;
    const uint32_t *target;
    uint8_t pattern_degree[NP];
    uint8_t target_degree[NT];
    uint32_t base_domain[NP];
    int8_t map[NP];
    uint32_t used;
    uint64_t nodes;
    uint64_t limit;
} MatchState;

static int pop32(uint32_t value) {
    return __builtin_popcount(value);
}

static int ctz32(uint32_t value) {
    return __builtin_ctz(value);
}

static int validate_graph(const uint32_t *adjacency, int order) {
    uint32_t full = ((uint32_t)1 << order) - 1u;
    for (int vertex = 0; vertex < order; ++vertex) {
        if (adjacency[vertex] & ~full) return 0;
        if (adjacency[vertex] & ((uint32_t)1 << vertex)) return 0;
        for (int other = 0; other < vertex; ++other) {
            int left = !!(adjacency[vertex] & ((uint32_t)1 << other));
            int right = !!(adjacency[other] & ((uint32_t)1 << vertex));
            if (left != right) return 0;
        }
    }
    return 1;
}

/* Kuhn matching for a small family of target-vertex domains. */
static int augment_domain(int item, const uint32_t *sets,
                          int owner[NT], uint32_t *seen) {
    uint32_t candidates = sets[item] & ~*seen;
    while (candidates) {
        int target = ctz32(candidates);
        candidates &= candidates - 1;
        *seen |= (uint32_t)1 << target;
        if (owner[target] < 0 ||
                augment_domain(owner[target], sets, owner, seen)) {
            owner[target] = item;
            return 1;
        }
    }
    return 0;
}

static int domains_have_sdr(const uint32_t *sets, int count) {
    int order[NP], owner[NT];
    for (int item = 0; item < count; ++item) order[item] = item;
    for (int first = 0; first < count; ++first) {
        for (int second = first + 1; second < count; ++second) {
            if (pop32(sets[order[second]]) < pop32(sets[order[first]])) {
                int temporary = order[first];
                order[first] = order[second];
                order[second] = temporary;
            }
        }
    }
    for (int target = 0; target < NT; ++target) owner[target] = -1;
    for (int position = 0; position < count; ++position) {
        uint32_t seen = 0;
        if (!augment_domain(order[position], sets, owner, &seen)) return 0;
    }
    return 1;
}

/* Necessary local condition for pattern -> target: distinct target
 * neighbours with sufficient total degrees must exist. */
static int root_neighbour_sdr(const MatchState *state,
                              int pattern_vertex, int target_vertex) {
    uint32_t sets[NP];
    int count = 0;
    uint32_t neighbours = state->pattern[pattern_vertex];
    while (neighbours) {
        int pattern_neighbour = ctz32(neighbours);
        uint32_t choices = 0;
        uint32_t target_neighbours = state->target[target_vertex];
        neighbours &= neighbours - 1;
        while (target_neighbours) {
            int candidate = ctz32(target_neighbours);
            target_neighbours &= target_neighbours - 1;
            if (state->target_degree[candidate] >=
                    state->pattern_degree[pattern_neighbour]) {
                choices |= (uint32_t)1 << candidate;
            }
        }
        if (!choices) return 0;
        sets[count++] = choices;
    }
    return domains_have_sdr(sets, count);
}

static int initialize_domains(MatchState *state) {
    for (int vertex = 0; vertex < NP; ++vertex) {
        state->pattern_degree[vertex] =
            (uint8_t)pop32(state->pattern[vertex]);
    }
    for (int vertex = 0; vertex < NT; ++vertex) {
        state->target_degree[vertex] =
            (uint8_t)pop32(state->target[vertex]);
    }
    for (int pattern_vertex = 0; pattern_vertex < NP; ++pattern_vertex) {
        uint32_t domain = 0;
        for (int target_vertex = 0; target_vertex < NT; ++target_vertex) {
            if (state->target_degree[target_vertex] <
                    state->pattern_degree[pattern_vertex]) continue;
            if (root_neighbour_sdr(state, pattern_vertex, target_vertex)) {
                domain |= (uint32_t)1 << target_vertex;
            }
        }
        state->base_domain[pattern_vertex] = domain;
        if (!domain) return 0;
    }
    return 1;
}

static int build_domains(const MatchState *state, uint32_t domain[NP],
                         int *unmapped_count) {
    int changed;
    *unmapped_count = 0;
    for (int pattern_vertex = 0; pattern_vertex < NP; ++pattern_vertex) {
        if (state->map[pattern_vertex] >= 0) {
            domain[pattern_vertex] =
                (uint32_t)1 << state->map[pattern_vertex];
            continue;
        }
        uint32_t current =
            state->base_domain[pattern_vertex] & ~state->used;
        for (int neighbour = 0; neighbour < NP; ++neighbour) {
            if (state->map[neighbour] >= 0 &&
                    (state->pattern[pattern_vertex] &
                     ((uint32_t)1 << neighbour))) {
                current &= state->target[(int)state->map[neighbour]];
            }
        }
        if (!current) return 0;
        domain[pattern_vertex] = current;
        ++*unmapped_count;
    }

    /* Edge arc consistency among unassigned pattern vertices. */
    do {
        changed = 0;
        for (int pattern_vertex = 0; pattern_vertex < NP;
                ++pattern_vertex) {
            uint32_t candidates, keep = 0;
            if (state->map[pattern_vertex] >= 0) continue;
            candidates = domain[pattern_vertex];
            while (candidates) {
                int target_vertex = ctz32(candidates);
                int supported = 1;
                uint32_t neighbours = state->pattern[pattern_vertex];
                candidates &= candidates - 1;
                while (neighbours) {
                    int pattern_neighbour = ctz32(neighbours);
                    neighbours &= neighbours - 1;
                    if (state->map[pattern_neighbour] < 0 &&
                            !(domain[pattern_neighbour] &
                              state->target[target_vertex])) {
                        supported = 0;
                        break;
                    }
                }
                if (supported) keep |= (uint32_t)1 << target_vertex;
            }
            if (!keep) return 0;
            if (keep != domain[pattern_vertex]) {
                domain[pattern_vertex] = keep;
                changed = 1;
            }
        }
    } while (changed);

    uint32_t sets[NP];
    int count = 0;
    for (int vertex = 0; vertex < NP; ++vertex) {
        if (state->map[vertex] < 0) sets[count++] = domain[vertex];
    }
    return domains_have_sdr(sets, count);
}

static int branch_neighbour_sdr(const MatchState *state,
                                int pattern_vertex, int target_vertex,
                                const uint32_t domain[NP]) {
    uint32_t sets[NP];
    int count = 0;
    uint32_t neighbours = state->pattern[pattern_vertex];
    while (neighbours) {
        int pattern_neighbour = ctz32(neighbours);
        neighbours &= neighbours - 1;
        if (state->map[pattern_neighbour] >= 0) continue;
        sets[count] = domain[pattern_neighbour] &
            state->target[target_vertex] &
            ~((uint32_t)1 << target_vertex);
        if (!sets[count]) return 0;
        ++count;
    }
    return domains_have_sdr(sets, count);
}

static int verify_mapping(const MatchState *state) {
    uint32_t used = 0;
    for (int pattern_vertex = 0; pattern_vertex < NP; ++pattern_vertex) {
        int target_vertex = state->map[pattern_vertex];
        if (target_vertex < 0 || target_vertex >= NT ||
                (used & ((uint32_t)1 << target_vertex))) return 0;
        used |= (uint32_t)1 << target_vertex;
    }
    for (int pattern_vertex = 0; pattern_vertex < NP; ++pattern_vertex) {
        uint32_t neighbours = state->pattern[pattern_vertex];
        while (neighbours) {
            int pattern_neighbour = ctz32(neighbours);
            neighbours &= neighbours - 1;
            if (!(state->target[(int)state->map[pattern_vertex]] &
                    ((uint32_t)1 << state->map[pattern_neighbour]))) {
                return 0;
            }
        }
    }
    return 1;
}

static enum MatchResult match_recursive(MatchState *state) {
    uint32_t domain[NP];
    int unmapped, chosen = -1, chosen_size = NT + 1;
    int chosen_frontier = -1;
    if (++state->nodes > state->limit) return MATCH_TIMEOUT;
    if (!build_domains(state, domain, &unmapped)) return MATCH_NO_HIT;
    if (unmapped == 0) {
        return verify_mapping(state) ? MATCH_HIT : MATCH_INVALID;
    }

    for (int pattern_vertex = 0; pattern_vertex < NP; ++pattern_vertex) {
        int size, frontier = 0;
        if (state->map[pattern_vertex] >= 0) continue;
        size = pop32(domain[pattern_vertex]);
        for (int neighbour = 0; neighbour < NP; ++neighbour) {
            if (state->map[neighbour] < 0 &&
                    (state->pattern[pattern_vertex] &
                     ((uint32_t)1 << neighbour))) ++frontier;
        }
        if (size < chosen_size ||
                (size == chosen_size && frontier > chosen_frontier) ||
                (size == chosen_size && frontier == chosen_frontier &&
                 (chosen < 0 || state->pattern_degree[pattern_vertex] >
                                  state->pattern_degree[chosen]))) {
            chosen = pattern_vertex;
            chosen_size = size;
            chosen_frontier = frontier;
        }
    }

    uint32_t candidates = domain[chosen];
    while (candidates) {
        int best_target = -1, best_score = -1;
        uint32_t scan = candidates;
        while (scan) {
            int target = ctz32(scan), score = 0;
            scan &= scan - 1;
            for (int neighbour = 0; neighbour < NP; ++neighbour) {
                if (state->map[neighbour] < 0 && neighbour != chosen &&
                        (state->pattern[chosen] &
                         ((uint32_t)1 << neighbour))) {
                    score += pop32(domain[neighbour] & state->target[target]);
                }
            }
            if (score > best_score) {
                best_target = target;
                best_score = score;
            }
        }
        candidates &= ~((uint32_t)1 << best_target);
        if (!branch_neighbour_sdr(
                state, chosen, best_target, domain)) continue;
        state->map[chosen] = (int8_t)best_target;
        state->used |= (uint32_t)1 << best_target;
        enum MatchResult result = match_recursive(state);
        if (result != MATCH_NO_HIT) return result;
        state->used &= ~((uint32_t)1 << best_target);
        state->map[chosen] = -1;
    }
    return MATCH_NO_HIT;
}

EXPORT int pattern954_match(const uint32_t pattern[NP],
                            const uint32_t target[NT],
                            uint64_t node_limit,
                            uint64_t *nodes_out,
                            int8_t mapping_out[NP]) {
    MatchState state;
    if (!pattern || !target || !nodes_out || !mapping_out || !node_limit) {
        return MATCH_INVALID;
    }
    if (!validate_graph(pattern, NP) || !validate_graph(target, NT)) {
        return MATCH_INVALID;
    }
    memset(&state, 0, sizeof(state));
    state.pattern = pattern;
    state.target = target;
    state.limit = node_limit;
    for (int vertex = 0; vertex < NP; ++vertex) {
        state.map[vertex] = -1;
        mapping_out[vertex] = -1;
    }
    if (!initialize_domains(&state)) {
        *nodes_out = 0;
        return MATCH_NO_HIT;
    }
    enum MatchResult result = match_recursive(&state);
    *nodes_out = state.nodes;
    if (result == MATCH_HIT) {
        memcpy(mapping_out, state.map, NP * sizeof(mapping_out[0]));
    }
    return result;
}
