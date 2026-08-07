#!/usr/bin/env python3
"""Independent exact checker for the standard 18-point geometry report."""

from __future__ import annotations

import argparse
import ast
import hashlib
import itertools
import json
import os
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parent
PRODUCER = "d6_standard18_geometry.py"
EXPECTED_PRODUCER_SHA256 = (
    "6e4a1af59a6a6063cade592cd3897ce3e0a0adbd2a46b28d4652a93c05303536"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_hash(value: object) -> str:
    raw = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")
    return hashlib.sha256(raw).hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def base_points() -> list[tuple[int, ...]]:
    result = []
    for word in range(32):
        if word.bit_count() & 1:
            result.append(
                tuple(2 * int(bool(word & (1 << coordinate))) - 1 for coordinate in range(5))
            )
    return result


def graph() -> tuple[int, ...]:
    points = base_points()
    rows = [0] * 18
    for first in range(16):
        for second in range(first):
            if sum(
                (points[first][axis] - points[second][axis]) ** 2
                for axis in range(5)
            ) == 8:
                rows[first] |= 1 << second
                rows[second] |= 1 << first
    for pole in (16, 17):
        for vertex in range(16):
            rows[pole] |= 1 << vertex
            rows[vertex] |= 1 << pole
    return tuple(rows)


def is_clique(rows: tuple[int, ...], vertices: tuple[int, ...]) -> bool:
    return all(
        rows[first] & (1 << second)
        for position, first in enumerate(vertices)
        for second in vertices[:position]
    )


def row_rank_mod(matrix: list[list[int]], prime: int) -> int:
    rows = [[entry % prime for entry in row] for row in matrix]
    if not rows:
        return 0
    pivot_row = 0
    for column in range(len(rows[0])):
        pivot = None
        for candidate in range(pivot_row, len(rows)):
            if rows[candidate][column]:
                pivot = candidate
                break
        if pivot is None:
            continue
        rows[pivot_row], rows[pivot] = rows[pivot], rows[pivot_row]
        scale = pow(rows[pivot_row][column], -1, prime)
        for index in range(column, len(rows[0])):
            rows[pivot_row][index] = rows[pivot_row][index] * scale % prime
        for candidate in range(pivot_row + 1, len(rows)):
            factor = rows[candidate][column]
            if factor:
                for index in range(column, len(rows[0])):
                    rows[candidate][index] = (
                        rows[candidate][index]
                        - factor * rows[pivot_row][index]
                    ) % prime
        pivot_row += 1
        if pivot_row == len(rows):
            break
    return pivot_row


def rigidity_rank(rows: tuple[int, ...]) -> int:
    points = [list(point) + [0] for point in base_points()]
    points += [[0, 0, 0, 0, 0, 4], [0, 0, 0, 0, 0, -4]]
    matrix = []
    for first in range(18):
        for second in range(first):
            if not rows[first] & (1 << second):
                continue
            equation = [0] * 108
            for axis in range(6):
                difference = points[first][axis] - points[second][axis]
                equation[6 * first + axis] = difference
                equation[6 * second + axis] = -difference
            matrix.append(equation)
    return row_rank_mod(matrix, 13)


def base_rigidity_rank(rows: tuple[int, ...]) -> int:
    points = base_points()
    matrix = []
    for first in range(16):
        for second in range(first):
            if not rows[first] & (1 << second):
                continue
            equation = [0] * 80
            for axis in range(5):
                difference = points[first][axis] - points[second][axis]
                equation[5 * first + axis] = difference
                equation[5 * second + axis] = -difference
            matrix.append(equation)
    return row_rank_mod(matrix, 1_000_000_007)


def affine_rank(points: list[tuple[int, ...]]) -> int:
    first = points[0]
    return row_rank_mod(
        [
            [point[axis] - first[axis] for axis in range(5)]
            for point in points[1:]
        ],
        1_000_000_007,
    )


def verify(path: Path, expected_report_sha256: str | None) -> dict:
    checks = {}
    report_hash = sha256(path)
    if expected_report_sha256 is not None:
        require(report_hash == expected_report_sha256, "report hash")
    require(sha256(ROOT / PRODUCER) == EXPECTED_PRODUCER_SHA256, "producer hash")
    report = json.loads(path.read_text(encoding="utf-8"))
    require(report.get("schema") == "d6-standard18-exact-geometry-v1", "schema")

    points = base_points()
    rows = graph()
    require(report["unit_graph"]["adjacency"] == list(rows), "adjacency")
    require(report["unit_graph"]["adjacency_sha256"] == stable_hash(list(rows)), "adjacency hash")
    base_edges = sum(
        bool(rows[first] & (1 << second))
        for first in range(16)
        for second in range(first)
    )
    all_edges = sum(row.bit_count() for row in rows) // 2
    require((base_edges, all_edges) == (80, 112), "edge counts")
    distance_histogram = {8: 0, 16: 0}
    for first in range(16):
        for second in range(first):
            distance = sum(
                (points[first][axis] - points[second][axis]) ** 2
                for axis in range(5)
            )
            require(distance in distance_histogram, "unexpected base distance")
            distance_histogram[distance] += 1
    require(distance_histogram == {8: 80, 16: 40}, "distance histogram")
    checks["coordinates_and_unit_graph"] = True

    base_k5 = sum(
        is_clique(rows, vertices)
        for vertices in itertools.combinations(range(16), 5)
    )
    base_k6 = any(
        is_clique(rows, vertices)
        for vertices in itertools.combinations(range(16), 6)
    )
    full_k6 = sum(
        is_clique(rows, vertices)
        for vertices in itertools.combinations(range(18), 6)
    )
    full_k7 = any(
        is_clique(rows, vertices)
        for vertices in itertools.combinations(range(18), 7)
    )
    require(base_k5 > 0 and not base_k6, "base clique number")
    require(full_k6 == 32 and not full_k7, "full clique counts")
    require(report["unit_graph"]["base_clique_number"] == 5, "reported clique")
    checks["clique_bounds"] = True

    require(rigidity_rank(rows) == 87, "full rigidity rank")
    require(base_rigidity_rank(rows) == 65, "base rigidity rank")
    require(report["rigidity"]["full_framework"]["exact_rank"] == 87, "reported full rank")
    require(report["rigidity"]["base_framework_R5"]["exact_rank"] == 65, "reported base rank")
    checks["exact_rigidity"] = True

    subsets_checked = 0
    for subset in itertools.combinations(range(16), 11):
        subsets_checked += 1
        require(affine_rank([points[index] for index in subset]) == 5, "11-subset span")
    section = [point for point in points if sum(point) == 1]
    require(len(section) == 10 and affine_rank(section) == 4, "ten-point section")
    require(report["affine_spanning"]["all_11_subsets_checked"] == subsets_checked, "subset count")
    checks["affine_spanning"] = True

    nonextension = report.get("nonextension")
    require(isinstance(nonextension, dict), "nonextension section")
    require(nonextension.get("minimum_forced_unit_base_neighbours") == 11, "forced neighbours")
    require(nonextension.get("forced_height_squared") == "3/8", "forced height")
    require(nonextension.get("terminal") == "collision with an existing apex", "terminal")
    semantics = report.get("semantics")
    require(semantics.get("classification_claim") is False, "classification overclaim")
    require(len(report.get("nonclaims", [])) >= 4, "missing nonclaims")
    checks["nonextension_scope"] = True

    syntax = ast.parse((ROOT / PRODUCER).read_text(encoding="utf-8"))
    imports = []
    for node in ast.walk(syntax):
        if isinstance(node, ast.Import):
            imports.extend(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imports.append(node.module)
    require(not any(name.startswith("d6_") for name in imports), "producer imported campaign kernel")
    checks["producer_import_independence"] = True

    return {
        "schema": "d6-standard18-exact-geometry-verification-v1",
        "status": "PASS",
        "report": {"path": path.name, "sha256": report_hash},
        "producer_sha256": EXPECTED_PRODUCER_SHA256,
        "checks": checks,
        "exact_results": {
            "full_rigidity_rank": 87,
            "base_rigidity_rank": 65,
            "base_clique_number": 5,
            "eleven_subsets_checked": subsets_checked,
            "standard_coordinates_nonextendable": True,
        },
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
        "--report", type=Path, default=ROOT / "d6_standard18_geometry_report.json"
    )
    parser.add_argument("--expected-report-sha256")
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "d6_standard18_geometry_verification.json",
    )
    args = parser.parse_args()
    result = verify(args.report.resolve(), args.expected_report_sha256)
    atomic_json(args.output, result)
    print(json.dumps({"status": "PASS", **result["exact_results"], "output": str(args.output)}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
