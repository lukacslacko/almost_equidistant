/* Heuristic-only R^6 realization/rigidity screen for the 18-deletion corpus.
 *
 * This deliberately reuses the campaign's existing LM implementation without
 * changing it.  The input format extends lm5.c's line format by appending
 *
 *     seed_flag [18*6 seed coordinates]
 *
 * after the 18 adjacency masks.  A seed is evaluated separately from the
 * random-start screen so that known positive controls cannot be mistaken for
 * discoveries.  Nothing printed by this program is a non-realizability
 * certificate; rank is floating-point complete-pivot elimination at an LM
 * endpoint and is meaningful only as numerical triage near zero residual.
 */

#define D 6
#define main inherited_lm5_main
#include "lm5.c"
#undef main

#include <errno.h>

static double minimum_distance(const Graph *g, const double *p){
    double best = 1e300;
    for (int i = 0; i < g->n; i++)
        for (int j = i + 1; j < g->n; j++){
            double squared = 0;
            for (int axis = 0; axis < D; axis++){
                double difference = p[i * D + axis] - p[j * D + axis];
                squared += difference * difference;
            }
            if (squared < best) best = squared;
        }
    return sqrt(best);
}

/* Numerical rank of the rigidity matrix using complete-pivot elimination.
 * The relative threshold is recorded by the Python run manifest. */
static int rigidity_rank(const Graph *g, const double *p, double relative_tolerance){
    static double matrix[MAXE][NVMAX];
    int rows = g->m, columns = g->n * D;
    memset(matrix, 0, sizeof(matrix));
    double scale = 0;
    for (int edge = 0; edge < rows; edge++){
        int first = g->ei[edge], second = g->ej[edge];
        for (int axis = 0; axis < D; axis++){
            double difference = p[first * D + axis] - p[second * D + axis];
            matrix[edge][first * D + axis] = difference;
            matrix[edge][second * D + axis] = -difference;
            if (fabs(difference) > scale) scale = fabs(difference);
        }
    }
    if (scale == 0) return 0;
    double threshold = relative_tolerance * scale;
    int rank = 0;
    while (rank < rows && rank < columns){
        int pivot_row = -1, pivot_column = -1;
        double pivot_size = 0;
        for (int row = rank; row < rows; row++)
            for (int column = rank; column < columns; column++){
                double candidate = fabs(matrix[row][column]);
                if (candidate > pivot_size){
                    pivot_size = candidate;
                    pivot_row = row;
                    pivot_column = column;
                }
            }
        if (pivot_size <= threshold) break;
        if (pivot_row != rank)
            for (int column = 0; column < columns; column++){
                double temporary = matrix[rank][column];
                matrix[rank][column] = matrix[pivot_row][column];
                matrix[pivot_row][column] = temporary;
            }
        if (pivot_column != rank)
            for (int row = 0; row < rows; row++){
                double temporary = matrix[row][rank];
                matrix[row][rank] = matrix[row][pivot_column];
                matrix[row][pivot_column] = temporary;
            }
        double pivot = matrix[rank][rank];
        for (int row = rank + 1; row < rows; row++){
            if (matrix[row][rank] == 0) continue;
            double factor = matrix[row][rank] / pivot;
            matrix[row][rank] = 0;
            for (int column = rank + 1; column < columns; column++)
                matrix[row][column] -= factor * matrix[rank][column];
        }
        rank++;
    }
    return rank;
}

static int parse_long(const char *text, long *answer){
    char *end = NULL;
    errno = 0;
    long value = strtol(text, &end, 10);
    if (errno || !end || *end) return 0;
    *answer = value;
    return 1;
}

int main(int argc, char **argv){
    if (argc != 6){
        fprintf(stderr,
                "usage: d6_18_lm_screen <graphfile> <start> <end> <restarts> <seedbase>\n");
        return 1;
    }
    long start, end;
    if (!parse_long(argv[2], &start) || !parse_long(argv[3], &end) ||
        start < 0 || end < start){
        fprintf(stderr, "invalid half-open index range\n");
        return 1;
    }
    int restarts = atoi(argv[4]);
    uint64_t seedbase = strtoull(argv[5], NULL, 10);
    if (restarts < 0){
        fprintf(stderr, "restarts must be nonnegative\n");
        return 1;
    }
    FILE *input = fopen(argv[1], "r");
    if (!input){ perror("open"); return 1; }
    char line[32768];
    long index = -1, emitted = 0;
    while (fgets(line, sizeof(line), input)){
        index++;
        if (index < start) continue;
        if (index >= end) break;
        Graph graph;
        graph.m = 0;
        char *token = strtok(line, " \t\r\n");
        if (!token){ fprintf(stderr, "empty graph line %ld\n", index); return 2; }
        graph.n = atoi(token);
        if (graph.n <= 0 || graph.n > MAXN){
            fprintf(stderr, "invalid n on graph line %ld\n", index); return 2;
        }
        uint32_t adjacency[MAXN];
        for (int vertex = 0; vertex < graph.n; vertex++){
            token = strtok(NULL, " \t\r\n");
            if (!token){ fprintf(stderr, "short adjacency line %ld\n", index); return 2; }
            adjacency[vertex] = (uint32_t)strtoul(token, NULL, 10);
        }
        for (int first = 0; first < graph.n; first++)
            for (int second = first + 1; second < graph.n; second++)
                if ((adjacency[first] >> second) & 1){
                    if (graph.m >= MAXE){ fprintf(stderr, "too many edges\n"); return 2; }
                    graph.ei[graph.m] = first;
                    graph.ej[graph.m] = second;
                    graph.m++;
                }
        token = strtok(NULL, " \t\r\n");
        if (!token){ fprintf(stderr, "missing seed flag line %ld\n", index); return 2; }
        int has_seed = atoi(token);
        static double seed_point[NVMAX];
        if (has_seed){
            for (int coordinate = 0; coordinate < graph.n * D; coordinate++){
                token = strtok(NULL, " \t\r\n");
                if (!token){ fprintf(stderr, "short seed line %ld\n", index); return 2; }
                seed_point[coordinate] = strtod(token, NULL);
            }
        }
        if (strtok(NULL, " \t\r\n")){
            fprintf(stderr, "unexpected trailing token line %ld\n", index); return 2;
        }

        double seed_residual = -1, seed_minimum_distance = -1;
        int seed_rank = -1;
        if (has_seed){
            seed_residual = residual(&graph, seed_point, NULL);
            seed_minimum_distance = minimum_distance(&graph, seed_point);
            seed_rank = rigidity_rank(&graph, seed_point, 1e-10);
        }

        double best = 1e300, best_distinct = 1e300, best_minimum_distance = 0;
        int near_solutions = 0, distinct_near_solutions = 0;
        static double point[NVMAX], best_distinct_point[NVMAX];
        const double sigmas[3] = {0.35, 0.5, 0.75};
        for (int restart = 0; restart < restarts; restart++){
            rng_s = seedbase * 1000003ULL + (uint64_t)index * 7919ULL +
                    (uint64_t)restart + 1;
            if (restart % 2 == 0){
                place_init(adjacency, graph.n, point);
            } else {
                double sigma = sigmas[(restart / 2) % 3];
                for (int coordinate = 0; coordinate < graph.n * D; coordinate++)
                    point[coordinate] = sigma * gauss_rand();
            }
            double value = lm_run(&graph, point, 400);
            if (value < best) best = value;
            if (value < 1e-20) near_solutions++;
            double separation = minimum_distance(&graph, point);
            if (separation > 1e-3 && value < best_distinct){
                best_distinct = value;
                best_minimum_distance = separation;
                memcpy(best_distinct_point, point,
                       sizeof(double) * graph.n * D);
            }
            if (separation > 1e-3 && value < 1e-20)
                distinct_near_solutions++;
        }
        int endpoint_rank = best_distinct < 1e299
            ? rigidity_rank(&graph, best_distinct_point, 1e-10) : -1;
        printf("%ld\t%d\t%d\t%d\t%.17g\t%.17g\t%d\t%.17g\t%.17g\t%.17g\t%d\t%d\t%d\n",
               index, graph.n, graph.m, has_seed,
               seed_residual, seed_minimum_distance, seed_rank,
               best, best_distinct, best_minimum_distance, endpoint_rank,
               near_solutions, distinct_near_solutions);
        fflush(stdout);
        emitted++;
        if (best_distinct < 1e-20 && best_minimum_distance > 1e-3){
            /* A machine-readable numerical witness for reconstruction.  It is
             * still only binary64 output, never an exact certificate. */
            fprintf(stderr, "HEURISTIC-WITNESS\t%ld\t%.17g\t%.17g\t%d",
                    index, best_distinct, best_minimum_distance, endpoint_rank);
            for (int coordinate = 0; coordinate < graph.n * D; coordinate++)
                fprintf(stderr, "\t%.17g", best_distinct_point[coordinate]);
            fputc('\n', stderr);
            fflush(stderr);
        }
    }
    fclose(input);
    if (emitted != end - start){
        fprintf(stderr, "input ended after %ld of %ld requested rows\n",
                emitted, end - start);
        return 3;
    }
    return 0;
}
