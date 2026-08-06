#!/usr/bin/env python3
"""Independent proof and containment checker for K7 Mobius-cycle cores.

This checker imports neither the production containment program nor the K7
search stack.  It rebuilds both required-edge patterns, verifies the exact
order-six Mobius calculation in Q(sqrt(7)), and independently exhausts all
K7 seeds, coordinate cycles, and role assignments in every pinned target.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from typing import Iterator, Sequence


ROOT = Path(__file__).resolve().parent
TARGETS = ROOT / ".runs/d6_n14_pattern_954_targets_12839.tsv"
SELECTION = ROOT / "d6_k7_rankone_tetrad_full_selection.json"
RUNNER = ROOT / "d6_k7_mobius_cycle_containment.py"
PROOF = ROOT / "d6_k7_mobius_cycle_obstructions.md"
TEST = ROOT / "test_d6_k7_mobius_cycle_containment.py"
REPORT = ROOT / "d6_k7_mobius_cycle_containment_report.json"
TARGETS_SHA256 = (
    "9355854172c324f9d93cc4085020974fb9622a010b2e8fe19f5a954d17a7277d"
)
SELECTION_SHA256 = (
    "86e4f19a203bca111a47924ae05461ada5b21b6341be59d38b1428bebe1f9479"
)
INDICES_SHA256 = (
    "320a68221ef597a1252d0c16cba6b30b4bef2bb6c13b73ce529e860ed431d497"
)
TARGET_COUNT = 12_839
TARGET_ORDER = 19


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_hash(value: object) -> str:
    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


@dataclass(frozen=True)
class QRoot7:
    rational: Fraction = Fraction(0)
    root: Fraction = Fraction(0)

    def __add__(self, other: "QRoot7") -> "QRoot7":
        return QRoot7(self.rational + other.rational, self.root + other.root)

    def __sub__(self, other: "QRoot7") -> "QRoot7":
        return QRoot7(self.rational - other.rational, self.root - other.root)

    def __neg__(self) -> "QRoot7":
        return QRoot7(-self.rational, -self.root)

    def __mul__(self, other: "QRoot7") -> "QRoot7":
        return QRoot7(
            self.rational * other.rational + 7 * self.root * other.root,
            self.rational * other.root + self.root * other.rational,
        )


ZERO = QRoot7()
ONE = QRoot7(Fraction(1))
ROOT7 = QRoot7(Fraction(0), Fraction(1))
Matrix = tuple[tuple[QRoot7, QRoot7], tuple[QRoot7, QRoot7]]


def matrix_multiply(left: Matrix, right: Matrix) -> Matrix:
    return tuple(
        tuple(
            left[row][0] * right[0][column]
            + left[row][1] * right[1][column]
            for column in range(2)
        )
        for row in range(2)
    )  # type: ignore[return-value]


def matrix_power(matrix: Matrix, exponent: int) -> Matrix:
    answer: Matrix = ((ONE, ZERO), (ZERO, ONE))
    factor = matrix
    while exponent:
        if exponent & 1:
            answer = matrix_multiply(answer, factor)
        factor = matrix_multiply(factor, factor)
        exponent >>= 1
    return answer


def fixed_point_discriminant(matrix: Matrix) -> QRoot7:
    a, b = matrix[0]
    c, d = matrix[1]
    difference = d - a
    four = QRoot7(Fraction(4))
    return difference * difference + four * b * c


def verify_mobius_algebra() -> dict:
    matrix: Matrix = ((ONE, -ROOT7), (ROOT7, QRoot7(Fraction(-4))))
    expected = [-3, -27, -108, -243, -243]
    observed = []
    for exponent, value in enumerate(expected, 1):
        discriminant = fixed_point_discriminant(
            matrix_power(matrix, exponent)
        )
        if discriminant != QRoot7(Fraction(value)):
            raise ValueError(f"wrong fixed-point discriminant for T^{exponent}")
        observed.append(value)
    sixth = matrix_power(matrix, 6)
    scalar: Matrix = (
        (QRoot7(Fraction(-27)), ZERO),
        (ZERO, QRoot7(Fraction(-27))),
    )
    if sixth != scalar:
        raise ValueError("Mobius matrix does not have projective order six")
    return {
        "fixed_point_discriminants_T1_through_T5": observed,
        "M6": [[[-27, 1], [0, 1]], [[0, 1], [-27, 1]]],
        "forbidden_cycle_lengths_checked": [3, 4],
        "retained_negative_control_cycle_length": 6,
    }


def add_edge(adjacency: list[int], first: int, second: int) -> None:
    adjacency[first] |= 1 << second
    adjacency[second] |= 1 << first


def expected_pattern(length: int) -> dict:
    order = 7 + length
    adjacency = [0] * order
    for first, second in combinations(range(7), 2):
        add_edge(adjacency, first, second)
    types = [(position, (position + 1) % length)
             for position in range(length)]
    for position, defect_type in enumerate(types):
        role = 7 + position
        for seed in range(7):
            if seed not in defect_type:
                add_edge(adjacency, role, seed)
        add_edge(adjacency, role, 7 + ((position + 1) % length))
    return {
        "schema": 1,
        "kind": f"d6_k7_two_defect_{length}_cycle_obstruction",
        "order": order,
        "edges": sum(row.bit_count() for row in adjacency) // 2,
        "adjacency": adjacency,
        "adjacency_sha256": stable_hash(adjacency),
        "seed_vertices": list(range(7)),
        "role_vertices": list(range(7, order)),
        "role_defect_upper_bounds": [list(pair) for pair in types],
        "required_role_cycle": [
            [7 + position, 7 + ((position + 1) % length)]
            for position in range(length)
        ],
        "nonedges_used_as_distance_constraints": False,
    }


def validate_graph(adjacency: Sequence[int], order: int) -> tuple[int, ...]:
    if len(adjacency) != order:
        raise ValueError("wrong graph order")
    full = (1 << order) - 1
    answer = tuple(adjacency)
    for vertex, row in enumerate(answer):
        if type(row) is not int or row & ~full or row & (1 << vertex):
            raise ValueError("invalid graph row")
        for other in range(vertex):
            if bool(row & (1 << other)) != bool(
                answer[other] & (1 << vertex)
            ):
                raise ValueError("asymmetric graph")
    return answer


def target_rows() -> list[dict]:
    if sha256(TARGETS) != TARGETS_SHA256:
        raise ValueError("target hash mismatch")
    if sha256(SELECTION) != SELECTION_SHA256:
        raise ValueError("selection hash mismatch")
    selection = json.loads(SELECTION.read_text(encoding="utf-8"))
    indices = selection["selected_indices"]
    if len(indices) != TARGET_COUNT or stable_hash(indices) != INDICES_SHA256:
        raise ValueError("selection indices mismatch")
    rows = []
    with TARGETS.open("r", encoding="ascii") as stream:
        for line_number, line in enumerate(stream, 1):
            fields = [int(field) for field in line.split()]
            if len(fields) != TARGET_ORDER + 1:
                raise ValueError(f"bad target row {line_number}")
            rows.append({
                "index": fields[0],
                "adjacency": validate_graph(fields[1:], TARGET_ORDER),
            })
    if [row["index"] for row in rows] != indices:
        raise ValueError("target rows do not follow the selection")
    return rows


def vertices(mask: int) -> list[int]:
    output = []
    while mask:
        bit = mask & -mask
        mask ^= bit
        output.append(bit.bit_length() - 1)
    return output


def seven_cliques(adjacency: Sequence[int]) -> Iterator[int]:
    """Independent iterative clique enumerator."""

    stack = [((1 << len(adjacency)) - 1, 7, 0)]
    while stack:
        candidates, needed, chosen = stack.pop()
        if needed == 0:
            yield chosen
            continue
        choices = []
        remaining = candidates
        while remaining.bit_count() >= needed:
            bit = remaining & -remaining
            remaining ^= bit
            choices.append((remaining & adjacency[bit.bit_length() - 1], bit))
        for following, bit in reversed(choices):
            stack.append((following, needed - 1, chosen | bit))


def cycles(length: int) -> Iterator[tuple[int, ...]]:
    if length == 3:
        for triple in combinations(range(7), 3):
            yield triple
    elif length == 4:
        for a, b, c, d in combinations(range(7), 4):
            yield (a, b, c, d)
            yield (a, b, d, c)
            yield (a, c, b, d)
    else:
        raise ValueError("unexpected cycle length")


def domains_for(
    adjacency: Sequence[int], seed: Sequence[int], cycle: Sequence[int]
) -> list[int]:
    full = (1 << len(adjacency)) - 1
    seed_mask = sum(1 << vertex for vertex in seed)
    outside = full ^ seed_mask
    domains = []
    for position, coordinate in enumerate(cycle):
        next_coordinate = cycle[(position + 1) % len(cycle)]
        required = seed_mask & ~(1 << seed[coordinate])
        required &= ~(1 << seed[next_coordinate])
        domain = 0
        for vertex in vertices(outside):
            if adjacency[vertex] & required == required:
                domain |= 1 << vertex
        domains.append(domain)
    return domains


def triangle_roles(adjacency: Sequence[int], domains: Sequence[int]) -> list[int] | None:
    for first in vertices(domains[0]):
        seconds = domains[1] & adjacency[first] & ~(1 << first)
        for second in vertices(seconds):
            thirds = domains[2] & adjacency[first] & adjacency[second]
            thirds &= ~(1 << first) & ~(1 << second)
            if thirds:
                return [first, second, vertices(thirds)[0]]
    return None


def quadrilateral_roles(
    adjacency: Sequence[int], domains: Sequence[int]
) -> list[int] | None:
    for first in vertices(domains[0]):
        seconds = domains[1] & adjacency[first] & ~(1 << first)
        for second in vertices(seconds):
            thirds = domains[2] & adjacency[second]
            thirds &= ~(1 << first) & ~(1 << second)
            for third in vertices(thirds):
                fourths = domains[3] & adjacency[third] & adjacency[first]
                fourths &= ~(1 << first) & ~(1 << second) & ~(1 << third)
                if fourths:
                    return [first, second, third, vertices(fourths)[0]]
    return None


def verify_mapping(
    pattern: Sequence[int], target: Sequence[int], mapping: Sequence[int]
) -> None:
    if len(mapping) != len(pattern) or len(set(mapping)) != len(mapping):
        raise ValueError("reported mapping is not injective")
    for first, row in enumerate(pattern):
        for second in vertices(row):
            if not target[mapping[first]] & (1 << mapping[second]):
                raise ValueError("reported mapping loses a required edge")


def independent_embedding(
    adjacency: Sequence[int], length: int
) -> tuple[list[int] | None, int, int]:
    seeds_tested = coordinate_cycles = 0
    for seed_mask in seven_cliques(adjacency):
        seeds_tested += 1
        seed = vertices(seed_mask)
        for coordinate_cycle in cycles(length):
            coordinate_cycles += 1
            domains = domains_for(adjacency, seed, coordinate_cycle)
            if any(not domain for domain in domains):
                continue
            roles = (
                triangle_roles(adjacency, domains)
                if length == 3
                else quadrilateral_roles(adjacency, domains)
            )
            if roles is None:
                continue
            remainder = [coordinate for coordinate in range(7)
                         if coordinate not in coordinate_cycle]
            mapping = [seed[value] for value in coordinate_cycle]
            mapping += [seed[value] for value in remainder]
            mapping += roles
            verify_mapping(expected_pattern(length)["adjacency"], adjacency, mapping)
            return mapping, seeds_tested, coordinate_cycles
    return None, seeds_tested, coordinate_cycles


def verify(report_path: Path, report_sha256: str) -> dict:
    if sha256(report_path) != report_sha256:
        raise ValueError("containment report hash mismatch")
    report = json.loads(report_path.read_text(encoding="utf-8"))
    if (
        report.get("schema") != 1
        or report.get("kind") != "d6_k7_mobius_cycle_containment_screen"
        or report.get("status") != "COMPLETE"
        or report.get("required_edge_only") is not True
        or report.get("targets") != TARGET_COUNT
        or report.get("target_indices_sha256") != INDICES_SHA256
    ):
        raise ValueError("containment report header mismatch")
    if report.get("controls", {}).get("status") != "PASS":
        raise ValueError("production controls did not pass")
    source_hashes = report.get("provenance", {}).get("source_hashes", {})
    expected_sources = {
        RUNNER.name: sha256(RUNNER),
        PROOF.name: sha256(PROOF),
        Path(__file__).name: sha256(Path(__file__)),
        TEST.name: sha256(TEST),
        SELECTION.name: SELECTION_SHA256,
        TARGETS.relative_to(ROOT).as_posix(): TARGETS_SHA256,
    }
    if source_hashes != expected_sources:
        raise ValueError("containment report source provenance mismatch")
    algebra = verify_mobius_algebra()
    for length in (3, 4):
        if report["patterns"][str(length)]["pattern"] != expected_pattern(length):
            raise ValueError(f"reported length-{length} pattern mismatch")

    rows = target_rows()
    independently_found = {3: [], 4: []}
    seeds = {3: 0, 4: 0}
    coordinate_counts = {3: 0, 4: 0}
    for ordinal, row in enumerate(rows):
        for length in (3, 4):
            mapping, local_seeds, local_cycles = independent_embedding(
                row["adjacency"], length
            )
            seeds[length] += local_seeds
            coordinate_counts[length] += local_cycles
            if mapping is not None:
                independently_found[length].append({
                    "ordinal": ordinal,
                    "index": row["index"],
                    "mapping": mapping,
                })
    for length in (3, 4):
        production = report["patterns"][str(length)]
        production_hits = production["hits"]
        for item in production_hits:
            row = rows[item["ordinal"]]
            if row["index"] != item["index"]:
                raise ValueError("production HIT ordinal/index mismatch")
            verify_mapping(
                expected_pattern(length)["adjacency"],
                row["adjacency"],
                item["mapping"],
            )
        if [item["index"] for item in production_hits] != [
            item["index"] for item in independently_found[length]
        ]:
            raise ValueError(f"independent length-{length} hit set mismatch")
        if production["hit_count"] != len(independently_found[length]):
            raise ValueError(f"length-{length} hit count mismatch")
        expected_cycles = seeds[length] * (35 if length == 3 else 105)
        if coordinate_counts[length] != expected_cycles:
            raise ValueError("independent coordinate-cycle completeness mismatch")
        if production["stats"]["k7_seeds"] != seeds[length]:
            raise ValueError("production/independent K7 seed count mismatch")
        if production["stats"]["coordinate_cycles"] != expected_cycles:
            raise ValueError("production coordinate-cycle count mismatch")
    combined = sorted({
        item["index"]
        for length in (3, 4)
        for item in independently_found[length]
    }, key={row["index"]: position for position, row in enumerate(rows)}.get)
    residue = [row["index"] for row in rows if row["index"] not in set(combined)]
    if (
        report["combined_hit_indices"] != combined
        or report["combined_hit_indices_sha256"] != stable_hash(combined)
        or report["combined_residue_indices"] != residue
        or report["combined_residue_indices_sha256"] != stable_hash(residue)
    ):
        raise ValueError("combined containment complement mismatch")
    return {
        "schema": 1,
        "kind": "d6_k7_mobius_cycle_containment_independent_verification",
        "status": "PASS",
        "report": {"path": report_path.name, "sha256": report_sha256},
        "required_edge_only": True,
        "mobius_algebra": algebra,
        "targets": len(rows),
        "K7_seeds_reenumerated": seeds[3],
        "triangle_coordinate_cycles_reenumerated": coordinate_counts[3],
        "quadrilateral_coordinate_cycles_reenumerated": coordinate_counts[4],
        "triangle_hits": len(independently_found[3]),
        "quadrilateral_hits": len(independently_found[4]),
        "combined_hits": len(combined),
        "residue": len(residue),
        "residue_indices_sha256": stable_hash(residue),
        "verifier_source_sha256": sha256(Path(__file__)),
    }


def atomic_json(path: Path, value: object) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    temporary.replace(path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, default=REPORT)
    parser.add_argument("--report-sha256", required=True)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "d6_k7_mobius_cycle_containment_verification.json",
    )
    args = parser.parse_args()
    result = verify(args.report, args.report_sha256)
    atomic_json(args.output, result)
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
