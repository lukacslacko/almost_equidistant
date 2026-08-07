#!/usr/bin/env python3
"""Exact geometry and nonextension certificate for the standard d=6 18-set.

The package concerns one explicit realization only: the 16 odd halfcube
vertices in a five-flat plus the two orthogonal apices.  It proves exact
rigidity of the full 112-edge unit framework and proves that no distinct point
can be appended to these coordinates while preserving almost-equidistance.
It does not classify all 16- or 18-point almost-equidistant realizations.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import os
import time
from pathlib import Path
from typing import Sequence


ROOT = Path(__file__).resolve().parent


def stable_hash(value: object) -> str:
    raw = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")
    return hashlib.sha256(raw).hexdigest()


def sign_vectors() -> list[tuple[int, ...]]:
    return [
        tuple(1 if word & (1 << coordinate) else -1 for coordinate in range(5))
        for word in range(32)
        if word.bit_count() % 2 == 1
    ]


def squared_difference(first: Sequence[int], second: Sequence[int]) -> int:
    return sum((a - b) ** 2 for a, b in zip(first, second))


def add_edge(rows: list[int], first: int, second: int) -> None:
    rows[first] |= 1 << second
    rows[second] |= 1 << first


def full_unit_graph() -> tuple[int, ...]:
    base = sign_vectors()
    rows = [0] * 18
    for first in range(16):
        for second in range(first):
            # Coordinates are scaled by 2 sqrt(2); unit squared distance is 8.
            if squared_difference(base[first], base[second]) == 8:
                add_edge(rows, first, second)
    for apex in (16, 17):
        for vertex in range(16):
            add_edge(rows, apex, vertex)
    return tuple(rows)


def edge_count(rows: Sequence[int]) -> int:
    return sum(row.bit_count() for row in rows) // 2


def base_edge_count(rows: Sequence[int]) -> int:
    mask = (1 << 16) - 1
    return sum((rows[vertex] & mask).bit_count() for vertex in range(16)) // 2


def maximal_cliques(rows: Sequence[int]) -> list[tuple[int, ...]]:
    answer = []
    n = len(rows)

    def bron_kerbosch(r: int, p: int, x: int) -> None:
        if not p and not x:
            answer.append(tuple(vertex for vertex in range(n) if r & (1 << vertex)))
            return
        union = p | x
        pivot = -1
        if union:
            pivot = max(
                (vertex for vertex in range(n) if union & (1 << vertex)),
                key=lambda vertex: (p & rows[vertex]).bit_count(),
            )
        candidates = p & ~(rows[pivot] if pivot >= 0 else 0)
        while candidates:
            bit = candidates & -candidates
            vertex = bit.bit_length() - 1
            bron_kerbosch(r | bit, p & rows[vertex], x & rows[vertex])
            p ^= bit
            x |= bit
            candidates ^= bit

    bron_kerbosch(0, (1 << n) - 1, 0)
    return sorted(answer)


def modular_rank(matrix: Sequence[Sequence[int]], prime: int) -> int:
    if not matrix:
        return 0
    data = [[entry % prime for entry in row] for row in matrix]
    rows = len(data)
    columns = len(data[0])
    rank = 0
    for column in range(columns):
        pivot = next(
            (row for row in range(rank, rows) if data[row][column]), None
        )
        if pivot is None:
            continue
        data[rank], data[pivot] = data[pivot], data[rank]
        inverse = pow(data[rank][column], prime - 2, prime)
        data[rank] = [(entry * inverse) % prime for entry in data[rank]]
        for row in range(rows):
            if row == rank or not data[row][column]:
                continue
            factor = data[row][column]
            data[row] = [
                (left - factor * right) % prime
                for left, right in zip(data[row], data[rank])
            ]
        rank += 1
        if rank == rows:
            break
    return rank


def rigidity_matrix_mod13(rows: Sequence[int]) -> list[list[int]]:
    # Scaled coordinates lie in Q(sqrt(3)); modulo 13 choose sqrt(3)=4.
    base = [list(vector) + [0] for vector in sign_vectors()]
    coordinates = base + [[0, 0, 0, 0, 0, 4], [0, 0, 0, 0, 0, -4]]
    matrix = []
    for first in range(18):
        for second in range(first):
            if not rows[first] & (1 << second):
                continue
            row = [0] * (18 * 6)
            difference = [
                coordinates[first][axis] - coordinates[second][axis]
                for axis in range(6)
            ]
            row[first * 6 : first * 6 + 6] = difference
            row[second * 6 : second * 6 + 6] = [-value for value in difference]
            matrix.append(row)
    return matrix


def base_rigidity_matrix() -> list[list[int]]:
    base = sign_vectors()
    rows = full_unit_graph()
    matrix = []
    for first in range(16):
        for second in range(first):
            if not rows[first] & (1 << second):
                continue
            row = [0] * (16 * 5)
            difference = [
                base[first][axis] - base[second][axis] for axis in range(5)
            ]
            row[first * 5 : first * 5 + 5] = difference
            row[second * 5 : second * 5 + 5] = [-value for value in difference]
            matrix.append(row)
    return matrix


def affine_rank_mod_prime(points: Sequence[Sequence[int]], prime: int) -> int:
    if len(points) <= 1:
        return 0
    origin = points[0]
    differences = [
        [(point[axis] - origin[axis]) % prime for axis in range(len(origin))]
        for point in points[1:]
    ]
    return modular_rank(differences, prime)


def affine_spanning_audit() -> dict:
    base = sign_vectors()
    prime = 1_000_000_007
    checked = 0
    for subset in itertools.combinations(range(16), 11):
        checked += 1
        if affine_rank_mod_prime([base[index] for index in subset], prime) != 5:
            raise AssertionError(f"11-subset failed to span: {subset}")
    ten_point_hyperplane = [
        index for index, point in enumerate(base) if sum(point) == 1
    ]
    if len(ten_point_hyperplane) != 10:
        raise AssertionError("the explicit ten-point hyperplane has wrong size")
    if affine_rank_mod_prime(
        [base[index] for index in ten_point_hyperplane], prime
    ) != 4:
        raise AssertionError("the ten-point witness is not a hyperplane section")
    return {
        "prime": prime,
        "all_11_subsets_checked": checked,
        "all_11_subsets_affinely_span_R5": True,
        "maximum_hyperplane_intersection": 10,
        "ten_point_hyperplane_equation_scaled": "s1+s2+s3+s4+s5=1",
        "ten_point_hyperplane_vertices": ten_point_hyperplane,
    }


def build_report() -> dict:
    base = sign_vectors()
    rows = full_unit_graph()
    base_distance_histogram = {}
    for first in range(16):
        for second in range(first):
            distance = squared_difference(base[first], base[second])
            base_distance_histogram[str(distance)] = (
                base_distance_histogram.get(str(distance), 0) + 1
            )
    cliques = maximal_cliques(rows)
    base_cliques = maximal_cliques(rows[:16])
    base_clique_number = max(map(len, base_cliques))
    if base_clique_number != 5:
        raise AssertionError("unexpected halfcube clique number")
    full_rank = modular_rank(rigidity_matrix_mod13(rows), 13)
    base_rank = modular_rank(base_rigidity_matrix(), 1_000_000_007)
    if full_rank != 87 or base_rank != 65:
        raise AssertionError("unexpected rigidity rank")
    affine = affine_spanning_audit()

    proof_steps = [
        "For a new point x, its non-unit neighbours in the 16-point base form a unit clique: otherwise x and two such vertices form a bad triple.",
        "The base unit graph has clique number 5, so x is unit from at least 11 base vertices.",
        "Every 11 base vertices affinely span the five-flat. Since all base norms are 5/8, subtracting their unit-sphere equations forces the horizontal projection of x to be zero.",
        "One remaining sphere equation gives x_6^2=3/8, so x equals one of the two existing apices and violates distinctness.",
    ]
    return {
        "schema": "d6-standard18-exact-geometry-v1",
        "status": "PASS",
        "coordinates": {
            "scale": "listed coordinates equal (2*sqrt(2)) times actual coordinates",
            "base": [list(point) + [0] for point in base],
            "apices": [
                [0, 0, 0, 0, 0, "+sqrt(3)"],
                [0, 0, 0, 0, 0, "-sqrt(3)"],
            ],
            "actual_base_norm_squared": "5/8",
            "actual_apex_height_squared": "3/8",
        },
        "unit_graph": {
            "adjacency": list(rows),
            "adjacency_sha256": stable_hash(list(rows)),
            "edges": edge_count(rows),
            "base_edges": base_edge_count(rows),
            "base_scaled_squared_distance_histogram": base_distance_histogram,
            "maximal_clique_size_histogram": {
                str(size): sum(len(clique) == size for clique in cliques)
                for size in sorted(set(map(len, cliques)))
            },
            "K6_seeds": sum(len(clique) == 6 for clique in cliques),
            "K7_seeds": sum(len(clique) == 7 for clique in cliques),
            "base_clique_number": base_clique_number,
        },
        "rigidity": {
            "full_framework": {
                "variables": 108,
                "euclidean_motion_dimension": 21,
                "rank_upper_bound": 87,
                "rank_lower_bound_mod_13": full_rank,
                "sqrt3_mod_13": 4,
                "exact_rank": 87,
                "self_stress_dimension": edge_count(rows) - 87,
            },
            "base_framework_R5": {
                "variables": 80,
                "euclidean_motion_dimension": 15,
                "exact_rank": 65,
                "self_stress_dimension": base_edge_count(rows) - 65,
            },
        },
        "affine_spanning": affine,
        "nonextension": {
            "result": "no distinct nineteenth point extends these coordinates",
            "proof_steps": proof_steps,
            "minimum_forced_unit_base_neighbours": 11,
            "forced_horizontal_projection": [0, 0, 0, 0, 0],
            "forced_height_squared": "3/8",
            "terminal": "collision with an existing apex",
        },
        "semantics": {
            "scope": "the displayed standard coordinates and their full actual unit graph",
            "nonedges": "unconstrained and may also be unit in abstract supports",
            "classification_claim": False,
        },
        "nonclaims": [
            "This does not prove uniqueness of 16-point almost-equidistant sets in R5.",
            "This does not prove uniqueness of 18-point almost-equidistant sets in R6.",
            "Local rigidity does not prove global uniqueness of an abstract support embedding.",
            "Failure to extend this component alone does not reject a 19-vertex parent support.",
        ],
    }


def atomic_json(path: Path, value: object) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path, default=ROOT / "d6_standard18_geometry_report.json"
    )
    args = parser.parse_args()
    report = build_report()
    atomic_json(args.output, report)
    print(
        json.dumps(
            {
                "output": str(args.output),
                "status": "PASS",
                "full_rigidity_rank": report["rigidity"]["full_framework"][
                    "exact_rank"
                ],
                "base_clique_number": report["unit_graph"]["base_clique_number"],
                "all_11_subsets_checked": report["affine_spanning"][
                    "all_11_subsets_checked"
                ],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
