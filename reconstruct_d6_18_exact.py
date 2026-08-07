#!/usr/bin/env python3
"""Exactify the two nonstandard 18-point frameworks found numerically.

The numerical sweep supplies only discovery coordinates.  This producer rounds
their squared distances to small rational numbers, then emits a self-contained
exact Euclidean-distance-matrix certificate.  It also proves that neither
exact configuration can be extended by a nineteenth point while preserving
almost-equidistance: the possible non-unit neighbours of a new point form a
unit clique, so it suffices to solve the equal-sphere equations for complements
of all maximal unit cliques.

The companion verifier recomputes every exact claim from the rational distance
matrices and does not import this producer.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import os
import tempfile
from collections import Counter
from fractions import Fraction
from pathlib import Path
from typing import Sequence


ROOT = Path(__file__).resolve().parent
DEFAULT_NUMERICAL_REPORT = ROOT / "d6_18_numerical_screen_report.json"
DEFAULT_CORPUS = ROOT / "d6_residue_18_deletions.json"
DEFAULT_OUTPUT = ROOT / "d6_18_exact_reconstruction.json"
TARGETS = {
    10261: [0, 1, 2, 3, 4, 6],
    10887: [0, 1, 2, 3, 4, 5],
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path, value: object) -> None:
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=path.name + ".tmp.", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
            json.dump(value, stream, indent=2, sort_keys=True)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary_name, path)
    except BaseException:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise


def encode(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else (
        f"{value.numerator}/{value.denominator}"
    )


def encode_matrix(matrix: Sequence[Sequence[Fraction]]) -> list[list[str]]:
    return [[encode(value) for value in row] for row in matrix]


def determinant(matrix: Sequence[Sequence[Fraction]]) -> Fraction:
    data = [list(row) for row in matrix]
    answer = Fraction(1)
    for column in range(len(data)):
        pivot = next(
            (row for row in range(column, len(data)) if data[row][column]), None
        )
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            data[column], data[pivot] = data[pivot], data[column]
            answer = -answer
        value = data[column][column]
        answer *= value
        for row in range(column + 1, len(data)):
            factor = data[row][column] / value
            for other in range(column + 1, len(data)):
                data[row][other] -= factor * data[column][other]
    return answer


def inverse(matrix: Sequence[Sequence[Fraction]]) -> list[list[Fraction]]:
    size = len(matrix)
    data = [
        list(matrix[row])
        + [Fraction(int(row == column)) for column in range(size)]
        for row in range(size)
    ]
    for column in range(size):
        pivot = next(
            row for row in range(column, size) if data[row][column]
        )
        data[column], data[pivot] = data[pivot], data[column]
        value = data[column][column]
        data[column] = [entry / value for entry in data[column]]
        for row in range(size):
            if row == column or not data[row][column]:
                continue
            factor = data[row][column]
            data[row] = [
                left - factor * right
                for left, right in zip(data[row], data[column])
            ]
    return [row[size:] for row in data]


def matvec(
    matrix: Sequence[Sequence[Fraction]], vector: Sequence[Fraction]
) -> list[Fraction]:
    return [sum((value * coordinate for value, coordinate in zip(row, vector)),
                Fraction(0)) for row in matrix]


def centered_gram(distances: Sequence[Sequence[Fraction]]) -> list[list[Fraction]]:
    size = len(distances)
    row_means = [sum(row, Fraction(0)) / size for row in distances]
    overall = sum(row_means, Fraction(0)) / size
    return [
        [
            -(distances[first][second] - row_means[first]
              - row_means[second] + overall) / 2
            for second in range(size)
        ]
        for first in range(size)
    ]


def standard_base_vectors() -> list[tuple[int, ...]]:
    return [
        tuple(1 if word & (1 << axis) else -1 for axis in range(5))
        for word in range(32)
        if word.bit_count() % 2 == 1
    ]


def switching_distance_matrix(
    replacements: dict[int, int],
) -> list[list[Fraction]]:
    """Exact distances after q -> (a_sign-q)/2.

    A vector is represented by rational coefficients of the 16 standard base
    vectors plus a rational multiple of a_+, whose squared norm is 3/8 and
    which is orthogonal to the base five-flat.
    """
    base = standard_base_vectors()
    base_gram = [
        [
            Fraction(sum(left * right for left, right in zip(first, second)), 8)
            for second in base
        ]
        for first in base
    ]
    vectors: list[tuple[dict[int, Fraction], Fraction]] = []
    for vertex in range(16):
        if vertex in replacements:
            vectors.append((
                {vertex: Fraction(-1, 2)},
                Fraction(replacements[vertex], 2),
            ))
        else:
            vectors.append(({vertex: Fraction(1)}, Fraction(0)))
    vectors.extend([
        ({}, Fraction(1)),
        ({}, Fraction(-1)),
    ])

    def inner_product(
        first: tuple[dict[int, Fraction], Fraction],
        second: tuple[dict[int, Fraction], Fraction],
    ) -> Fraction:
        return (
            sum(
                coefficient_first * coefficient_second * base_gram[i][j]
                for i, coefficient_first in first[0].items()
                for j, coefficient_second in second[0].items()
            )
            + first[1] * second[1] * Fraction(3, 8)
        )

    norms = [inner_product(vector, vector) for vector in vectors]
    distances = [[Fraction(0) for _ in range(18)] for _ in range(18)]
    for first in range(18):
        for second in range(first):
            value = (
                norms[first] + norms[second]
                - 2 * inner_product(vectors[first], vectors[second])
            )
            distances[first][second] = distances[second][first] = value
    return distances


def weighted_isomorphism(
    source: Sequence[Sequence[Fraction]],
    target: Sequence[Sequence[Fraction]],
) -> list[int]:
    signatures_source = [Counter(row) for row in source]
    signatures_target = [Counter(row) for row in target]
    mapping: dict[int, int] = {}
    used: set[int] = set()

    def visit() -> bool:
        if len(mapping) == 18:
            return True
        best_source = -1
        best_options: list[int] | None = None
        for source_vertex in range(18):
            if source_vertex in mapping:
                continue
            options = [
                target_vertex
                for target_vertex in range(18)
                if target_vertex not in used
                and signatures_source[source_vertex]
                    == signatures_target[target_vertex]
                and all(
                    source[source_vertex][other_source]
                    == target[target_vertex][other_target]
                    for other_source, other_target in mapping.items()
                )
            ]
            if not options:
                return False
            if best_options is None or len(options) < len(best_options):
                best_source, best_options = source_vertex, options
        assert best_options is not None
        for target_vertex in best_options:
            mapping[best_source] = target_vertex
            used.add(target_vertex)
            if visit():
                return True
            used.remove(target_vertex)
            del mapping[best_source]
        return False

    if not visit():
        raise AssertionError("rational EDM does not match switching construction")
    return [mapping[vertex] for vertex in range(18)]


def maximal_cliques(rows: Sequence[int]) -> list[tuple[int, ...]]:
    size = len(rows)
    answer = []

    def search(current: int, possible: int, excluded: int) -> None:
        if not possible and not excluded:
            answer.append(tuple(
                vertex for vertex in range(size) if current & (1 << vertex)
            ))
            return
        union = possible | excluded
        pivot = max(
            (vertex for vertex in range(size) if union & (1 << vertex)),
            key=lambda vertex: (possible & rows[vertex]).bit_count(),
            default=-1,
        )
        candidates = possible & ~(rows[pivot] if pivot >= 0 else 0)
        while candidates:
            bit = candidates & -candidates
            vertex = bit.bit_length() - 1
            search(
                current | bit,
                possible & rows[vertex],
                excluded & rows[vertex],
            )
            possible ^= bit
            excluded |= bit
            candidates ^= bit

    search(0, (1 << size) - 1, 0)
    return sorted(answer)


def solve_linear(
    matrix: Sequence[Sequence[Fraction]], right: Sequence[Fraction]
) -> dict:
    variables = len(matrix[0])
    data = [list(row) + [value] for row, value in zip(matrix, right)]
    pivots = []
    rank = 0
    for column in range(variables):
        pivot = next(
            (row for row in range(rank, len(data)) if data[row][column]), None
        )
        if pivot is None:
            continue
        data[rank], data[pivot] = data[pivot], data[rank]
        value = data[rank][column]
        data[rank] = [entry / value for entry in data[rank]]
        for row in range(len(data)):
            if row == rank or not data[row][column]:
                continue
            factor = data[row][column]
            data[row] = [
                left - factor * right_entry
                for left, right_entry in zip(data[row], data[rank])
            ]
        pivots.append(column)
        rank += 1
    inconsistent = any(
        all(not entry for entry in row[:-1]) and row[-1]
        for row in data
    )
    result = {
        "coefficient_rank": rank,
        "augmented_rank": rank + int(inconsistent),
        "status": "INCONSISTENT" if inconsistent else "CONSISTENT",
    }
    if not inconsistent and rank == variables:
        solution = [Fraction(0)] * variables
        for row, column in enumerate(pivots):
            solution[column] = data[row][-1]
        result["unique_solution"] = solution
    elif not inconsistent:
        result["status"] = "CONSISTENT_UNDERDETERMINED"
    return result


def rational_distance_matrix(points: Sequence[Sequence[float]]) -> tuple[
    list[list[Fraction]], float
]:
    distances = [[Fraction(0) for _ in range(18)] for _ in range(18)]
    maximum_error = 0.0
    for first in range(18):
        for second in range(first):
            observed = sum(
                (points[first][axis] - points[second][axis]) ** 2
                for axis in range(6)
            )
            exact = Fraction(observed).limit_denominator(64)
            maximum_error = max(maximum_error, abs(observed - float(exact)))
            if abs(observed - float(exact)) > 1e-10:
                raise AssertionError("distance failed small-rational reconstruction")
            distances[first][second] = distances[second][first] = exact
    return distances, maximum_error


def reconstruct(index: int, witness: dict, corpus_record: dict) -> dict:
    points = witness["coordinates_vertex_major_binary64"]
    distances, reconstruction_error = rational_distance_matrix(points)
    gram = centered_gram(distances)
    basis = TARGETS[index]
    principal = [[gram[first][second] for second in basis] for first in basis]
    principal_inverse = inverse(principal)
    leading_minors = [
        determinant([row[:size] for row in principal[:size]])
        for size in range(1, 7)
    ]
    if not all(value > 0 for value in leading_minors):
        raise AssertionError("principal Gram block is not positive definite")
    # Exact column-space factorization B = B[:,J] B[J,J]^-1 B[J,:].
    for first in range(18):
        left = [gram[first][vertex] for vertex in basis]
        coefficients = matvec(principal_inverse, left)
        for second in range(18):
            reconstructed = sum(
                (coefficients[position] * gram[basis[position]][second]
                 for position in range(6)),
                Fraction(0),
            )
            if reconstructed != gram[first][second]:
                raise AssertionError("rank-six Gram factorization failed")

    unit_rows = [
        sum(
            1 << second
            for second in range(18)
            if second != first and distances[first][second] == 1
        )
        for first in range(18)
    ]
    required_edges = sum(mask.bit_count() for mask in corpus_record["adjacency"]) // 2
    if any(
        distances[first][second] != 1
        for first in range(18)
        for second in range(first)
        if corpus_record["adjacency"][first] & (1 << second)
    ):
        raise AssertionError("a required edge did not reconstruct to unit length")
    bad_triples = [
        triple
        for triple in itertools.combinations(range(18), 3)
        if not any(
            distances[first][second] == 1
            for first, second in itertools.combinations(triple, 2)
        )
    ]
    if bad_triples:
        raise AssertionError("reconstruction is not almost equidistant")
    if any(distances[first][second] <= 0
           for first in range(18) for second in range(first)):
        raise AssertionError("reconstruction contains a collision")

    cliques = maximal_cliques(unit_rows)
    extension_checks = []
    for clique in cliques:
        unit_neighbours = [
            vertex for vertex in range(18) if vertex not in clique
        ]
        matrix = [
            [-2 * gram[vertex][basis_vertex] for basis_vertex in basis]
            + [Fraction(1)]
            for vertex in unit_neighbours
        ]
        right = [Fraction(1) - gram[vertex][vertex]
                 for vertex in unit_neighbours]
        solution = solve_linear(matrix, right)
        record = {
            "maximal_nonunit_clique": list(clique),
            "forced_unit_neighbours": unit_neighbours,
            "linear_status": solution["status"],
            "coefficient_rank": solution["coefficient_rank"],
            "augmented_rank": solution["augmented_rank"],
        }
        if solution["status"] == "CONSISTENT":
            values = solution["unique_solution"]
            coefficients, asserted_norm = values[:6], values[6]
            actual_norm = sum(
                (coefficients[first] * principal[first][second]
                 * coefficients[second]
                 for first in range(6) for second in range(6)),
                Fraction(0),
            )
            record.update({
                "basis_coefficients": [encode(value) for value in coefficients],
                "linear_norm_variable": encode(asserted_norm),
                "quadratic_norm": encode(actual_norm),
                "norm_discrepancy": encode(asserted_norm - actual_norm),
                "extension_exists": asserted_norm == actual_norm,
            })
        elif solution["status"] == "CONSISTENT_UNDERDETERMINED":
            raise AssertionError("unexpected underdetermined extension system")
        extension_checks.append(record)
    if any(record.get("extension_exists") for record in extension_checks):
        raise AssertionError("configuration unexpectedly extends")

    replacements = {0: 1} if index == 10887 else {0: 1, 7: -1}
    switched_distances = switching_distance_matrix(replacements)
    switching_permutation = weighted_isomorphism(distances, switched_distances)
    if any(
        distances[first][second]
        != switched_distances[switching_permutation[first]][
            switching_permutation[second]
        ]
        for first in range(18) for second in range(18)
    ):
        raise AssertionError("switching permutation failed exact distance check")

    distance_histogram = Counter(
        distances[first][second]
        for first in range(18) for second in range(first)
    )
    clique_histogram = Counter(map(len, cliques))
    consistent = [
        record for record in extension_checks
        if record["linear_status"] == "CONSISTENT"
    ]
    return {
        "name": f"nonstandard_exact_18_from_deletion_class_{index}",
        "deletion_class_index": index,
        "class_id": corpus_record["class_id"],
        "discovery_numerical_maximum_squared_distance_rounding_error": (
            reconstruction_error
        ),
        "exact_squared_distance_matrix": encode_matrix(distances),
        "exact_centered_gram_matrix": encode_matrix(gram),
        "gram_certificate": {
            "centroid": "origin",
            "principal_basis_indices": basis,
            "principal_matrix": encode_matrix(principal),
            "principal_determinant": encode(determinant(principal)),
            "leading_principal_minors": [encode(value) for value in leading_minors],
            "factorization": "B = B[:,J] inverse(B[J,J]) B[J,:]",
            "positive_semidefinite": True,
            "rank": 6,
        },
        "exact_geometry": {
            "points": 18,
            "affine_dimension": 6,
            "all_points_distinct": True,
            "almost_equidistant": True,
            "unit_graph_adjacency": unit_rows,
            "unit_edges": sum(row.bit_count() for row in unit_rows) // 2,
            "required_edges": required_edges,
            "extra_unit_edges_beyond_required_support": (
                sum(row.bit_count() for row in unit_rows) // 2 - required_edges
            ),
            "squared_distance_histogram": {
                encode(value): count
                for value, count in sorted(distance_histogram.items())
            },
        },
        "switching_construction": {
            "standard_base": (
                "q_w=(sign vector w)/(2 sqrt(2)) in R5 x {0}, for the "
                "16 odd-parity sign vectors; ||q_w||^2=5/8"
            ),
            "standard_apices": (
                "a_+=(0,...,0,sqrt(3/8)), a_-=-a_+; ||a_+||^2=3/8"
            ),
            "replacement_rule": "replace q by x_q=(a_sign-q)/2",
            "replacements_in_switching_labels": [
                {
                    "base_label": vertex,
                    "apex_sign": "+" if sign > 0 else "-",
                }
                for vertex, sign in sorted(replacements.items())
            ],
            "canonical_vertex_to_switching_label": switching_permutation,
            "exact_distance_matrix_match": True,
        },
        "global_nonextension_certificate": {
            "argument": (
                "For any added point x, its non-unit neighbours in the exact "
                "18-set must form a unit clique. Extend that clique to a "
                "maximal unit clique C. Then x is unit from every vertex "
                "outside C. All maximal C are enumerated below; their exact "
                "linear sphere equations are inconsistent, or their unique "
                "linear solution violates s=y^T B[J,J] y."
            ),
            "maximal_cliques_checked": len(cliques),
            "maximal_clique_size_histogram": {
                str(size): count for size, count in sorted(clique_histogram.items())
            },
            "linear_inconsistent": sum(
                record["linear_status"] == "INCONSISTENT"
                for record in extension_checks
            ),
            "linear_consistent_but_norm_failed": len(consistent),
            "extension_exists": False,
            "checks": extension_checks,
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--numerical-report", type=Path, default=DEFAULT_NUMERICAL_REPORT
    )
    parser.add_argument("--corpus", type=Path, default=DEFAULT_CORPUS)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    numerical = json.loads(args.numerical_report.read_text(encoding="utf-8"))
    corpus = json.loads(args.corpus.read_text(encoding="utf-8"))
    witnesses = {
        record["index"]: record
        for record in numerical["summary"][
            "random_distinct_near_solution_witnesses"
        ]
    }
    if set(witnesses) != set(TARGETS):
        raise RuntimeError("unexpected numerical witness set")
    configurations = [
        reconstruct(index, witnesses[index], corpus["unique_deletions"][index])
        for index in sorted(TARGETS)
    ]
    unit_edge_counts = [
        record["exact_geometry"]["unit_edges"] for record in configurations
    ]
    if len(set(unit_edge_counts + [112])) != 3:
        raise AssertionError("unit-edge counts did not separate the configurations")
    # The distance-2 graph on the 16 base points is the Clebsch graph.  These
    # exact parameters are the only combinatorial facts used in the switching
    # family proof below.
    base = standard_base_vectors()
    clebsch = [
        sum(
            1 << second
            for second in range(16)
            if first != second
            and sum(
                (base[first][axis] - base[second][axis]) ** 2
                for axis in range(5)
            ) == 16
        )
        for first in range(16)
    ]
    clebsch_degrees = sorted(set(row.bit_count() for row in clebsch))
    adjacent_common = sorted(set(
        (clebsch[first] & clebsch[second]).bit_count()
        for first in range(16) for second in range(first)
        if clebsch[first] & (1 << second)
    ))
    nonadjacent_common = sorted(set(
        (clebsch[first] & clebsch[second]).bit_count()
        for first in range(16) for second in range(first)
        if not clebsch[first] & (1 << second)
    ))
    if (clebsch_degrees, adjacent_common, nonadjacent_common) != ([5], [0], [2]):
        raise AssertionError("unexpected Clebsch parameters")
    report = {
        "schema": "d6-18-exact-nonstandard-reconstruction-v1",
        "status": "PASS",
        "source": {
            "numerical_report": str(args.numerical_report.resolve()),
            "numerical_report_sha256": sha256_file(args.numerical_report),
            "deletion_corpus": str(args.corpus.resolve()),
            "deletion_corpus_sha256": sha256_file(args.corpus),
            "producer": str(Path(__file__).resolve()),
            "producer_sha256": sha256_file(Path(__file__).resolve()),
        },
        "exact_conclusions": [
            "The two displayed rational squared-distance matrices are Euclidean, positive semidefinite of affine dimension 6, and have 18 distinct points.",
            "Each displayed configuration is almost equidistant.",
            "Neither displayed configuration admits any nineteenth distinct point preserving almost-equidistance.",
            "Their unit-edge counts are 110 and 111, versus 112 for the standard odd-halfcube-plus-apices configuration, so all three are pairwise nonisometric.",
        ],
        "nonclaims": [
            "This does not classify every 18-point almost-equidistant configuration in R6.",
            "This does not reject any other 18-deletion support or any 19-point candidate graph.",
            "The numerical coordinates are discovery provenance only; the exact claims depend solely on the rational matrices and exact verifier.",
        ],
        "pairwise_nonisometry_witness": {
            "invariant": "number of unit-distance pairs",
            "nonstandard_counts": unit_edge_counts,
            "standard18_count": 112,
        },
        "switching_family_theorem": {
            "operation": (
                "Starting from the standard odd-halfcube-plus-apices set, "
                "replace selected base vertices q by (a_sign-q)/2."
            ),
            "exact_distance_identities": [
                "For an unchanged base point p, the new squared distance is 1 if ||p-q||^2=1 and 1/2 if ||p-q||^2=2.",
                "The replacement has squared distance 1/4 from its assigned apex and 1 from the opposite apex.",
                "Two replacements with the same sign have squared distance ||q-r||^2/4; with opposite signs they have squared distance 3/8+||q-r||^2/4. None of these values is 1.",
            ],
            "clebsch_distance2_graph": {
                "vertices": 16,
                "degree_values": clebsch_degrees,
                "common_neighbours_for_adjacent_pair": adjacent_common,
                "common_neighbours_for_nonadjacent_pair": nonadjacent_common,
                "strongly_regular_parameters": [16, 5, 0, 2],
            },
            "classification_within_this_switching_operation": {
                "maximum_replacements": 2,
                "reason_maximum_two": (
                    "Every pair of replacements is non-unit, so three "
                    "replacements themselves form a bad triple."
                ),
                "zero_replacements": "the standard 112-unit-edge type",
                "one_replacement": (
                    "always valid; the Clebsch neighbourhood of q is "
                    "independent because lambda=0; gives the 111-edge type"
                ),
                "two_replacements": (
                    "valid exactly for opposite apex signs and a Clebsch "
                    "edge qr (standard squared distance 2); same signs fail "
                    "with their common assigned apex, while a Clebsch "
                    "nonedge has two common Clebsch neighbours and creates "
                    "a bad triple; gives the 110-edge type"
                ),
                "orbit_note": (
                    "Even sign flips act transitively on odd sign vectors; "
                    "after composing coordinate permutations with the even "
                    "sign flips that return a chosen vertex, its stabilizer "
                    "acts as S5 and is transitive on distance-2 neighbours; "
                    "swapping the apices "
                    "swaps signs. Thus there is one type for each of 0,1,2 "
                    "valid replacements within this operation."
                ),
                "number_of_types": 3,
            },
            "scope": (
                "This classifies only the displayed switching family, not all "
                "18-point almost-equidistant configurations in R6."
            ),
        },
        "configurations": configurations,
    }
    atomic_json(args.output, report)
    print(json.dumps({
        "status": report["status"],
        "output": str(args.output),
        "sha256": sha256_file(args.output),
        "configurations": [
            {
                "index": record["deletion_class_index"],
                "unit_edges": record["exact_geometry"]["unit_edges"],
                "maximal_cliques_checked": record[
                    "global_nonextension_certificate"
                ]["maximal_cliques_checked"],
                "linear_consistent_but_norm_failed": record[
                    "global_nonextension_certificate"
                ]["linear_consistent_but_norm_failed"],
            }
            for record in configurations
        ],
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
