#!/usr/bin/env python3
"""Independent exact checker for the two nonstandard d=6 18-sets."""

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
DEFAULT_CERTIFICATE = ROOT / "d6_18_exact_reconstruction.json"
DEFAULT_CORPUS = ROOT / "d6_residue_18_deletions.json"
DEFAULT_OUTPUT = ROOT / "d6_18_exact_reconstruction_verification.json"
PRODUCER = ROOT / "reconstruct_d6_18_exact.py"


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


def decode(value: str) -> Fraction:
    return Fraction(value)


def decode_matrix(matrix: Sequence[Sequence[str]]) -> list[list[Fraction]]:
    return [[decode(value) for value in row] for row in matrix]


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
            (row for row in range(column, size) if data[row][column]), None
        )
        if pivot is None:
            raise AssertionError("singular claimed principal block")
        data[column], data[pivot] = data[pivot], data[column]
        value = data[column][column]
        data[column] = [entry / value for entry in data[column]]
        for row in range(size):
            if row == column:
                continue
            factor = data[row][column]
            if factor:
                data[row] = [
                    left - factor * right
                    for left, right in zip(data[row], data[column])
                ]
    return [row[size:] for row in data]


def centered_gram(distances: Sequence[Sequence[Fraction]]) -> list[list[Fraction]]:
    size = len(distances)
    means = [sum(row, Fraction(0)) / size for row in distances]
    overall = sum(means, Fraction(0)) / size
    return [
        [
            -(distances[first][second] - means[first] - means[second] + overall)
            / 2
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


def switching_distances(replacements: dict[int, int]) -> list[list[Fraction]]:
    base = standard_base_vectors()
    base_gram = [[
        Fraction(sum(a * b for a, b in zip(first, second)), 8)
        for second in base
    ] for first in base]
    vectors = []
    for vertex in range(16):
        vectors.append(
            ({vertex: Fraction(-1, 2)}, Fraction(replacements[vertex], 2))
            if vertex in replacements
            else ({vertex: Fraction(1)}, Fraction(0))
        )
    vectors += [({}, Fraction(1)), ({}, Fraction(-1))]

    def dot(first: tuple, second: tuple) -> Fraction:
        return (
            sum(
                left * right * base_gram[i][j]
                for i, left in first[0].items()
                for j, right in second[0].items()
            )
            + first[1] * second[1] * Fraction(3, 8)
        )

    norms = [dot(vector, vector) for vector in vectors]
    answer = [[Fraction(0) for _ in range(18)] for _ in range(18)]
    for first in range(18):
        for second in range(first):
            answer[first][second] = answer[second][first] = (
                norms[first] + norms[second] - 2 * dot(vectors[first], vectors[second])
            )
    return answer


def maximal_cliques(rows: Sequence[int]) -> list[tuple[int, ...]]:
    # Independent lexicographic branch-and-bound enumeration.  This is not the
    # pivoting Bron--Kerbosch implementation used by the producer.
    size = len(rows)
    maximal = []

    def visit(chosen: tuple[int, ...], candidates: tuple[int, ...]) -> None:
        extended = False
        for position, vertex in enumerate(candidates):
            remaining = tuple(
                other for other in candidates[position + 1 :]
                if rows[vertex] & (1 << other)
            )
            visit(chosen + (vertex,), remaining)
            extended = True
        if not extended and chosen:
            mask = sum(1 << vertex for vertex in chosen)
            if not any(
                (rows[vertex] & mask) == mask
                for vertex in range(size) if not mask & (1 << vertex)
            ):
                maximal.append(chosen)

    visit(tuple(), tuple(range(size)))
    return sorted(set(maximal))


def solve_linear(
    matrix: Sequence[Sequence[Fraction]], right: Sequence[Fraction]
) -> tuple[str, int, int, list[Fraction] | None]:
    variables = len(matrix[0])
    data = [list(row) + [value] for row, value in zip(matrix, right)]
    pivot_columns = []
    rank = 0
    for column in range(variables):
        pivot = next(
            (row for row in range(rank, len(data)) if data[row][column]), None
        )
        if pivot is None:
            continue
        data[rank], data[pivot] = data[pivot], data[rank]
        value = data[rank][column]
        for entry in range(column, variables + 1):
            data[rank][entry] /= value
        for row in range(rank + 1, len(data)):
            factor = data[row][column]
            if factor:
                for entry in range(column, variables + 1):
                    data[row][entry] -= factor * data[rank][entry]
        pivot_columns.append(column)
        rank += 1
    inconsistent = any(
        all(not row[column] for column in range(variables)) and row[-1]
        for row in data
    )
    if inconsistent:
        return "INCONSISTENT", rank, rank + 1, None
    if rank < variables:
        return "CONSISTENT_UNDERDETERMINED", rank, rank, None
    # Back substitution, unlike the producer's all-row RREF.
    solution = [Fraction(0)] * variables
    for row in range(rank - 1, -1, -1):
        column = pivot_columns[row]
        solution[column] = data[row][-1] - sum(
            data[row][other] * solution[other]
            for other in range(column + 1, variables)
        )
    return "CONSISTENT", rank, rank, solution


def check_configuration(record: dict, corpus_record: dict) -> dict:
    distances = decode_matrix(record["exact_squared_distance_matrix"])
    if len(distances) != 18 or any(len(row) != 18 for row in distances):
        raise AssertionError("distance matrix is not 18 by 18")
    if any(distances[first][second] != distances[second][first]
           for first in range(18) for second in range(18)):
        raise AssertionError("distance matrix is not symmetric")
    if any(distances[first][first] for first in range(18)):
        raise AssertionError("distance diagonal is nonzero")
    if any(distances[first][second] <= 0
           for first in range(18) for second in range(first)):
        raise AssertionError("points are not exactly distinct")

    gram = centered_gram(distances)
    if gram != decode_matrix(record["exact_centered_gram_matrix"]):
        raise AssertionError("stored centered Gram matrix is wrong")
    if any(sum(row, Fraction(0)) for row in gram):
        raise AssertionError("centered Gram matrix does not annihilate one")
    for first in range(18):
        for second in range(18):
            recovered = gram[first][first] + gram[second][second] - 2 * gram[first][second]
            if recovered != distances[first][second]:
                raise AssertionError("Gram matrix does not reproduce distances")

    basis = record["gram_certificate"]["principal_basis_indices"]
    if len(basis) != 6 or len(set(basis)) != 6:
        raise AssertionError("invalid principal basis")
    principal = [[gram[first][second] for second in basis] for first in basis]
    minors = [
        determinant([row[:size] for row in principal[:size]])
        for size in range(1, 7)
    ]
    if not all(value > 0 for value in minors):
        raise AssertionError("principal block is not positive definite")
    principal_inverse = inverse(principal)
    for first in range(18):
        for second in range(18):
            factored = sum(
                gram[first][basis[left]] * principal_inverse[left][right]
                * gram[basis[right]][second]
                for left in range(6) for right in range(6)
            )
            if factored != gram[first][second]:
                raise AssertionError("Gram rank-six factorization failed")

    unit_rows = [
        sum(1 << second for second in range(18)
            if first != second and distances[first][second] == 1)
        for first in range(18)
    ]
    if unit_rows != record["exact_geometry"]["unit_graph_adjacency"]:
        raise AssertionError("stored unit graph is wrong")
    for first in range(18):
        for second in range(first):
            if (corpus_record["adjacency"][first] & (1 << second)
                    and distances[first][second] != 1):
                raise AssertionError("required support edge is not unit")
    bad_triples = sum(
        not any(distances[first][second] == 1
                for first, second in itertools.combinations(triple, 2))
        for triple in itertools.combinations(range(18), 3)
    )
    if bad_triples:
        raise AssertionError("exact configuration has a bad triple")

    construction = record["switching_construction"]
    replacements = {
        item["base_label"]: 1 if item["apex_sign"] == "+" else -1
        for item in construction["replacements_in_switching_labels"]
    }
    switched = switching_distances(replacements)
    permutation = construction["canonical_vertex_to_switching_label"]
    if sorted(permutation) != list(range(18)):
        raise AssertionError("switching-label map is not a permutation")
    if any(
        distances[first][second]
        != switched[permutation[first]][permutation[second]]
        for first in range(18) for second in range(18)
    ):
        raise AssertionError("exact EDM does not match the switching formula")

    cliques = maximal_cliques(unit_rows)
    stored_checks = {
        tuple(item["maximal_nonunit_clique"]): item
        for item in record["global_nonextension_certificate"]["checks"]
    }
    if set(stored_checks) != set(cliques):
        raise AssertionError("maximal-clique list is incomplete or has extras")
    consistent = 0
    discrepancy_histogram = Counter()
    for clique in cliques:
        neighbours = [vertex for vertex in range(18) if vertex not in clique]
        matrix = [
            [-2 * gram[vertex][basis_vertex] for basis_vertex in basis]
            + [Fraction(1)]
            for vertex in neighbours
        ]
        right = [Fraction(1) - gram[vertex][vertex] for vertex in neighbours]
        status, rank, augmented_rank, solution = solve_linear(matrix, right)
        stored = stored_checks[clique]
        if (status, rank, augmented_rank) != (
            stored["linear_status"], stored["coefficient_rank"],
            stored["augmented_rank"],
        ):
            raise AssertionError("stored sphere-system rank/status is wrong")
        if status == "CONSISTENT":
            consistent += 1
            coefficients, asserted_norm = solution[:6], solution[6]
            actual_norm = sum(
                coefficients[first] * principal[first][second]
                * coefficients[second]
                for first in range(6) for second in range(6)
            )
            discrepancy = asserted_norm - actual_norm
            if not discrepancy:
                raise AssertionError("an exact extension exists")
            discrepancy_histogram[discrepancy] += 1
            if Fraction(stored["norm_discrepancy"]) != discrepancy:
                raise AssertionError("stored norm discrepancy is wrong")
        elif status != "INCONSISTENT":
            raise AssertionError("underdetermined case needs a polynomial check")

    unit_edges = sum(row.bit_count() for row in unit_rows) // 2
    if unit_edges != record["exact_geometry"]["unit_edges"]:
        raise AssertionError("stored unit-edge count is wrong")
    if len(cliques) != record["global_nonextension_certificate"][
        "maximal_cliques_checked"
    ]:
        raise AssertionError("stored maximal-clique count is wrong")
    if consistent != record["global_nonextension_certificate"][
        "linear_consistent_but_norm_failed"
    ]:
        raise AssertionError("stored consistent-system count is wrong")
    return {
        "deletion_class_index": record["deletion_class_index"],
        "unit_edges": unit_edges,
        "required_edges": sum(
            mask.bit_count() for mask in corpus_record["adjacency"]
        ) // 2,
        "affine_dimension": 6,
        "distinct_points": 18,
        "bad_triples": bad_triples,
        "maximal_cliques_recomputed": len(cliques),
        "maximal_clique_size_histogram": {
            str(size): count
            for size, count in sorted(Counter(map(len, cliques)).items())
        },
        "linear_consistent_norm_failures": consistent,
        "norm_discrepancy_histogram": {
            str(value): count
            for value, count in sorted(discrepancy_histogram.items())
        },
        "globally_nonextendable": True,
        "switching_replacements": len(replacements),
        "switching_formula_recomputed": True,
    }


def verify(certificate_path: Path, corpus_path: Path) -> dict:
    certificate = json.loads(certificate_path.read_text(encoding="utf-8"))
    corpus = json.loads(corpus_path.read_text(encoding="utf-8"))
    if certificate["schema"] != "d6-18-exact-nonstandard-reconstruction-v1":
        raise AssertionError("unexpected certificate schema")
    if sha256_file(PRODUCER) != certificate["source"]["producer_sha256"]:
        raise AssertionError("producer source hash mismatch")
    if sha256_file(corpus_path) != certificate["source"]["deletion_corpus_sha256"]:
        raise AssertionError("deletion corpus hash mismatch")
    checked = []
    for record in certificate["configurations"]:
        index = record["deletion_class_index"]
        checked.append(check_configuration(
            record, corpus["unique_deletions"][index]
        ))
    counts = [item["unit_edges"] for item in checked]
    if len(set(counts + [112])) != 3:
        raise AssertionError("unit-edge invariant does not prove nonisometry")
    base = standard_base_vectors()
    clebsch = [
        sum(
            1 << second for second in range(16)
            if first != second
            and sum(
                (base[first][axis] - base[second][axis]) ** 2
                for axis in range(5)
            ) == 16
        )
        for first in range(16)
    ]
    degrees = sorted(set(row.bit_count() for row in clebsch))
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
    if (degrees, adjacent_common, nonadjacent_common) != ([5], [0], [2]):
        raise AssertionError("switching-family Clebsch parameters failed")
    stored_clebsch = certificate["switching_family_theorem"][
        "clebsch_distance2_graph"
    ]
    if stored_clebsch["strongly_regular_parameters"] != [16, 5, 0, 2]:
        raise AssertionError("stored switching-family parameters are wrong")
    return {
        "schema": "d6-18-exact-nonstandard-reconstruction-verification-v1",
        "status": "PASS",
        "certificate_path": str(certificate_path.resolve()),
        "certificate_sha256": sha256_file(certificate_path),
        "producer_sha256": sha256_file(PRODUCER),
        "arithmetic": "fractions.Fraction exact rational arithmetic",
        "configurations_checked": checked,
        "exact_conclusions": {
            "new_exact_18_point_configurations": 2,
            "each_is_almost_equidistant_in_R6": True,
            "each_is_globally_nonextendable_to_19": True,
            "pairwise_nonisometric_with_each_other_and_standard18": True,
            "switching_family_formula_verified": True,
            "switching_family_clebsch_parameters": [16, 5, 0, 2],
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--certificate", type=Path, default=DEFAULT_CERTIFICATE)
    parser.add_argument("--corpus", type=Path, default=DEFAULT_CORPUS)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    try:
        result = verify(args.certificate, args.corpus)
    except Exception as error:
        result = {
            "schema": "d6-18-exact-nonstandard-reconstruction-verification-v1",
            "status": "FAIL",
            "error": f"{type(error).__name__}: {error}",
            "certificate_path": str(args.certificate.resolve()),
            "certificate_sha256": sha256_file(args.certificate),
        }
    atomic_json(args.output, result)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
