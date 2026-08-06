#!/usr/bin/env python3
"""Independent exact checker for the pattern-954 algebraic obstruction.

This checker does not import the extractor or any K7 theory implementation.
It reconstructs the 13-vertex core from the hash-pinned n=14 corpus, checks
the required-edge hypotheses of the simplex argument, verifies the arithmetic
in Q(sqrt(7)), and independently reconstructs the 12,839 target selection and
its transient containment input.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from typing import Sequence


ROOT = Path(__file__).resolve().parent
EXPECTED_INPUT_SHA256 = (
    "e302b818747a08e07fbd137211c05886a7e290cd4b561c1e977d4e0830b953a0"
)
EXPECTED_EXTRACTOR_SHA256 = (
    "972d522d1efd3f38cb7f2d95cdb3e0c2ffd4008268071ab0b0dab18343d1d320"
)
EXPECTED_N14_SHA256 = (
    "0e3d74c081b272731848d655da0c68cfba09435e39a2fd2107c32c2bd3b378f0"
)
EXPECTED_RANK_SHA256 = (
    "7a0a350142930e4e13830ea44c4fd217f25aeb07bbd6286e51e316f72712c4f9"
)
EXPECTED_PARENT_SHA256 = (
    "814c80651746ce05048f72f0cfa7e49a921554abaf167bd5c8cf4063a34d35c0"
)
EXPECTED_DUAL_REPORT_SHA256 = (
    "c2de7e06bbe4f3d163a4f747d42721f7f60da985b867ce74669664431738e345"
)
EXPECTED_DUAL_DECISIONS_SHA256 = (
    "724a97928422fc8705946ed64e3ce9afd37579d3229e5b3ebb7c4f1144c002ec"
)
EXPECTED_TARGET_SHA256 = (
    "9355854172c324f9d93cc4085020974fb9622a010b2e8fe19f5a954d17a7277d"
)
EXPECTED_CURRENT_INDICES_SHA256 = (
    "320a68221ef597a1252d0c16cba6b30b4bef2bb6c13b73ce529e860ed431d497"
)
EXPECTED_CORE_ADJACENCY = (
    2842, 6777, 3496, 8039, 4931, 7374, 8122,
    5476, 7901, 7515, 7148, 5999, 4090,
)
EXPECTED_ORIGINAL_VERTICES = (0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13)
EXPECTED_SEED = (3, 6, 8, 9, 10, 11, 12)
EXPECTED_ROLES = {"A": 0, "B": 2, "C": 4, "D": 7, "E": 1, "F": 5}
EXPECTED_DEFECT_BOUNDS = {
    "A": {1, 4, 6},
    "C": {0, 4, 5},
    "E": {2, 4},
    "B": {1, 3, 6},
    "D": {0, 3, 5},
    "F": {2, 3},
}
EXPECTED_TRIANGLES = (("A", "C", "E"), ("B", "D", "F"))
EXPECTED_BRIDGE = ("E", "F")


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


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def validate_graph(adjacency: Sequence[int]) -> None:
    count = len(adjacency)
    full = (1 << count) - 1
    for vertex, row in enumerate(adjacency):
        require(isinstance(row, int) and not isinstance(row, bool), "noninteger row")
        require(not (row & ~full), "out-of-range adjacency bit")
        require(not (row & (1 << vertex)), "adjacency loop")
        for other in range(vertex):
            require(
                bool(row & (1 << other))
                == bool(adjacency[other] & (1 << vertex)),
                "asymmetric adjacency",
            )


def is_clique(adjacency: Sequence[int], vertices: Sequence[int]) -> bool:
    return len(vertices) == len(set(vertices)) and all(
        adjacency[first] & (1 << second)
        for first, second in combinations(vertices, 2)
    )


def clique_masks(adjacency: Sequence[int], size: int) -> list[int]:
    answer: list[int] = []

    def visit(candidates: int, need: int, chosen: int) -> None:
        if need == 0:
            answer.append(chosen)
            return
        while candidates.bit_count() >= need:
            bit = candidates & -candidates
            candidates ^= bit
            vertex = bit.bit_length() - 1
            visit(candidates & adjacency[vertex], need - 1, chosen | bit)

    visit((1 << len(adjacency)) - 1, size, 0)
    return answer


def source_adjacency(corpus: Path) -> tuple[int, ...]:
    require(sha256(corpus) == EXPECTED_N14_SHA256, "n=14 corpus hash mismatch")
    rows = corpus.read_text(encoding="ascii").splitlines()
    require(len(rows) == 1_052, "n=14 corpus population mismatch")
    fields = tuple(map(int, rows[954].split()))
    require(fields[0] == 14 and len(fields) == 15, "malformed pattern row")
    adjacency = fields[1:]
    validate_graph(adjacency)
    require(len(clique_masks(adjacency, 7)) == 1, "source K7 is not unique")
    return adjacency


def reconstruct_core(source: Sequence[int]) -> tuple[int, ...]:
    position = {vertex: index for index, vertex in enumerate(EXPECTED_ORIGINAL_VERTICES)}
    removed = {frozenset((0, 3)), frozenset((5, 8))}
    adjacency = [0] * 13
    for first, second in combinations(EXPECTED_ORIGINAL_VERTICES, 2):
        if not (source[first] & (1 << second)):
            continue
        if frozenset((first, second)) in removed:
            continue
        left, right = position[first], position[second]
        adjacency[left] |= 1 << right
        adjacency[right] |= 1 << left
    validate_graph(adjacency)
    return tuple(adjacency)


def defect_coordinates(
    adjacency: Sequence[int], seed: Sequence[int], vertex: int
) -> set[int]:
    return {
        coordinate
        for coordinate, seed_vertex in enumerate(seed)
        if not (adjacency[vertex] & (1 << seed_vertex))
    }


def verify_obstruction_hypotheses(
    adjacency: Sequence[int], seed: Sequence[int], roles: dict[str, int]
) -> dict:
    """Check only required-edge hypotheses; extra edges are permitted."""

    validate_graph(adjacency)
    require(len(seed) == 7 and is_clique(adjacency, seed), "seed is not a K7")
    require(set(roles) == set(EXPECTED_ROLES), "role set changed")
    labels = list(seed) + list(roles.values())
    require(len(labels) == len(set(labels)), "seed/role labels are not distinct")
    require(all(0 <= label < len(adjacency) for label in labels), "bad role label")
    actual_defects = {
        role: defect_coordinates(adjacency, seed, vertex)
        for role, vertex in roles.items()
    }
    for role, actual in actual_defects.items():
        require(
            actual <= EXPECTED_DEFECT_BOUNDS[role],
            f"role {role} has an unallowed defect coordinate",
        )
    for triangle in EXPECTED_TRIANGLES:
        require(
            is_clique(adjacency, [roles[role] for role in triangle]),
            f"required role triangle {triangle} is absent",
        )
    first, second = (roles[role] for role in EXPECTED_BRIDGE)
    require(adjacency[first] & (1 << second), "required E-F bridge is absent")

    # The coordinate-overlap upper bounds used by every dot-product equation.
    require(
        EXPECTED_DEFECT_BOUNDS["A"] & EXPECTED_DEFECT_BOUNDS["C"]
        == EXPECTED_DEFECT_BOUNDS["A"] & EXPECTED_DEFECT_BOUNDS["E"]
        == EXPECTED_DEFECT_BOUNDS["C"] & EXPECTED_DEFECT_BOUNDS["E"]
        == {4},
        "first triangle does not have singleton overlap 4",
    )
    require(
        EXPECTED_DEFECT_BOUNDS["B"] & EXPECTED_DEFECT_BOUNDS["D"]
        == EXPECTED_DEFECT_BOUNDS["B"] & EXPECTED_DEFECT_BOUNDS["F"]
        == EXPECTED_DEFECT_BOUNDS["D"] & EXPECTED_DEFECT_BOUNDS["F"]
        == {3},
        "second triangle does not have singleton overlap 3",
    )
    require(
        EXPECTED_DEFECT_BOUNDS["E"] & EXPECTED_DEFECT_BOUNDS["F"] == {2},
        "bridge does not have singleton overlap 2",
    )
    return {role: sorted(value) for role, value in actual_defects.items()}


# Elements a + b*sqrt(7), represented as rational pairs (a,b).
Q7 = tuple[Fraction, Fraction]


def q7(a: int | Fraction = 0, b: int | Fraction = 0) -> Q7:
    return Fraction(a), Fraction(b)


def q7_add(left: Q7, right: Q7) -> Q7:
    return left[0] + right[0], left[1] + right[1]


def q7_sub(left: Q7, right: Q7) -> Q7:
    return left[0] - right[0], left[1] - right[1]


def q7_mul(left: Q7, right: Q7) -> Q7:
    return (
        left[0] * right[0] + 7 * left[1] * right[1],
        left[0] * right[1] + left[1] * right[0],
    )


def q7_scale(value: Q7, scalar: Fraction) -> Q7:
    return value[0] * scalar, value[1] * scalar


def verify_algebra() -> list[dict]:
    """Check the two-support formula and all four terminal sign cases."""

    root = q7(0, 1)
    one = q7(1, 0)
    records = []
    for sigma in (-1, 1):
        # From e^2=(sqrt(7)-e-sigma)^2.  The other factor would require
        # sqrt(7)=sigma, whose square would say 7=1.
        e = q7_scale(q7_sub(root, q7(sigma)), Fraction(1, 2))
        right = q7_sub(q7_sub(root, e), q7(sigma))
        require(q7_mul(e, e) == q7_mul(right, right), "E diagonal formula failed")
        require(7 != sigma * sigma, "spurious diagonal factor was not excluded")
        for tau in (-1, 1):
            f = q7_scale(q7_sub(root, q7(tau)), Fraction(1, 2))
            right_f = q7_sub(q7_sub(root, f), q7(tau))
            require(
                q7_mul(f, f) == q7_mul(right_f, right_f),
                "F diagonal formula failed",
            )
            product = q7_mul(e, f)
            require(product != one, "a terminal sign case unexpectedly survives")
            left_terminal = 3 + sigma * tau
            rational_coefficient = sigma + tau
            # Equality of a rational and sqrt(7) times an integer is checked
            # directly in the quadratic basis (1,sqrt(7)).
            terminal_difference = q7(left_terminal, -rational_coefficient)
            require(terminal_difference != q7(), "terminal equation vanished")
            records.append({
                "sigma": sigma,
                "tau": tau,
                "E_coordinate_2": [e[0].numerator, e[0].denominator,
                                     e[1].numerator, e[1].denominator],
                "F_coordinate_2": [f[0].numerator, f[0].denominator,
                                     f[1].numerator, f[1].denominator],
                "bridge_product_minus_one": [
                    (product[0] - 1).numerator,
                    (product[0] - 1).denominator,
                    product[1].numerator,
                    product[1].denominator,
                ],
            })
    require(len(records) == 4, "sign-case enumeration is incomplete")
    return records


def reconstruct_targets(payload: dict) -> tuple[list[int], dict[int, tuple[int, ...]]]:
    targets = payload["current_K7_targets"]
    parent_path = ROOT / targets["parent_selection"]
    report_path = ROOT / targets["full_dual_report"]
    decisions_path = ROOT / targets["full_dual_decisions"]
    rank_path = ROOT / targets["rank_survivors"]
    target_path = ROOT / targets["transient_tsv"]
    require(sha256(parent_path) == EXPECTED_PARENT_SHA256, "parent hash mismatch")
    require(sha256(report_path) == EXPECTED_DUAL_REPORT_SHA256, "dual report hash mismatch")
    require(
        sha256(decisions_path) == EXPECTED_DUAL_DECISIONS_SHA256,
        "dual decisions hash mismatch",
    )
    require(sha256(rank_path) == EXPECTED_RANK_SHA256, "rank input hash mismatch")
    require(sha256(target_path) == EXPECTED_TARGET_SHA256, "target TSV hash mismatch")

    parent = json.loads(parent_path.read_text(encoding="utf-8"))["selected_indices"]
    with gzip.open(decisions_path, "rt", encoding="ascii", newline="") as stream:
        decisions = list(csv.DictReader(stream, delimiter="\t"))
    require([int(row["index"]) for row in decisions] == parent, "decision order mismatch")
    selected = [int(row["index"]) for row in decisions if row["status"] == "SURVIVOR"]
    require(len(selected) == 12_839, "wrong current selection size")
    require(stable_hash(selected) == EXPECTED_CURRENT_INDICES_SHA256, "index hash mismatch")

    rank_payload = json.loads(rank_path.read_text(encoding="utf-8"))
    rank_rows = {
        int(row["index"]): tuple(map(int, row["adjacency"]))
        for row in rank_payload["graphs"]
    }
    observed: list[int] = []
    target_rows: dict[int, tuple[int, ...]] = {}
    for line_number, line in enumerate(target_path.read_text(encoding="ascii").splitlines(), 1):
        fields = tuple(map(int, line.split()))
        require(len(fields) == 20, f"malformed target row {line_number}")
        index, adjacency = fields[0], fields[1:]
        require(index not in target_rows, "target TSV repeats an index")
        validate_graph(adjacency)
        require(rank_rows.get(index) == adjacency, "target adjacency differs from rank input")
        observed.append(index)
        target_rows[index] = adjacency
    require(observed == selected, "target TSV order/coverage mismatch")
    return selected, target_rows


def verify(input_path: Path, corpus_path: Path) -> dict:
    require(sha256(input_path) == EXPECTED_INPUT_SHA256, "input JSON hash mismatch")
    require(
        sha256(ROOT / "d6_n14_pattern_954_extract.py")
        == EXPECTED_EXTRACTOR_SHA256,
        "extractor hash mismatch",
    )
    payload = json.loads(input_path.read_text(encoding="utf-8"))
    require(payload.get("schema") == 1, "input schema mismatch")
    require(
        payload.get("kind") == "d6_n14_pattern_954_algebraic_core",
        "input kind mismatch",
    )
    source = source_adjacency(corpus_path)
    require(payload["source_pattern"]["adjacency"] == list(source), "source copy mismatch")
    core = reconstruct_core(source)
    require(core == EXPECTED_CORE_ADJACENCY, "independent core reconstruction changed")
    recorded = payload["obstruction_core"]
    require(recorded["adjacency"] == list(core), "recorded core adjacency mismatch")
    require(recorded["original_vertex_labels"] == list(EXPECTED_ORIGINAL_VERTICES),
            "core original-label map mismatch")
    require(recorded["edges"] == 54, "core edge count mismatch")
    require(recorded["nonedges_used_as_distance_constraints"] is False,
            "certificate prescribes a nonedge distance")
    require(recorded["seed"] == list(EXPECTED_SEED), "core seed mismatch")
    require(recorded["roles"] == EXPECTED_ROLES, "core roles mismatch")
    actual_defects = verify_obstruction_hypotheses(core, EXPECTED_SEED, EXPECTED_ROLES)
    require(
        actual_defects
        == {name: sorted(value) for name, value in EXPECTED_DEFECT_BOUNDS.items()},
        "core defect masks are not the exact advertised bounds",
    )
    signs = verify_algebra()
    selected, _ = reconstruct_targets(payload)
    return {
        "status": "PASS",
        "source_pattern_index": 954,
        "core_order": 13,
        "core_edges": 54,
        "required_edge_only": True,
        "sign_cases_checked": len(signs),
        "current_target_population": len(selected),
        "current_target_indices_sha256": stable_hash(selected),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=ROOT / "d6_n14_pattern_954_input.json",
    )
    parser.add_argument(
        "--corpus", type=Path, default=ROOT / "aeq_d6_n14.txt"
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    report = verify(args.input, args.corpus)
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
