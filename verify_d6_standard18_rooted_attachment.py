#!/usr/bin/env python3
"""Independent checker for the rooted standard-18 attachment audit.

This checker imports neither the producer nor the deletion-corpus builder.  It
reconstructs the coordinate arithmetic, rigidity matrices, colored occurrence
isomorphisms, and every attachment contradiction from pinned inputs.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import time
from fractions import Fraction
from pathlib import Path
from typing import Sequence


ROOT = Path(__file__).resolve().parent
PRODUCER = "d6_standard18_rooted_attachment.py"
EXPECTED_PRODUCER_SHA256 = (
    "ec06eaf6a09a754597d9a66b35925da35c28314dafe200d186f65bc04287cd71"
)
DELETION = "d6_residue_18_deletions.json"
DELETION_VERIFICATION = "d6_residue_18_deletions_verification.json"
V4 = "d6_current_residue_manifest_v4.json"
V4_VERIFICATION = "d6_current_residue_manifest_v4_verification.json"
EXPECTED_HASHES = {
    DELETION: "9d08f9ec579434c9601ad3908bd2ae51f10e09e9d06a7f38fbc25794a634b732",
    DELETION_VERIFICATION: "50ab5d88441b2869e60a05f3e1485db0b25008c3a816960a20a3bb21c1c4defd",
    V4: "6ab4bfffc524f5b43409d59888fb596de8130315bda3381e484bb4ce891e03e4",
    V4_VERIFICATION: "765eceb6782d131dae8c77c9b94de78735e23a070da051c8fffbd984790e1a41",
}


Q3 = tuple[Fraction, Fraction]
ZERO: Q3 = (Fraction(0), Fraction(0))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def object_hash(value: object) -> str:
    raw = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")
    return hashlib.sha256(raw).hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def decode_q3(value: object) -> Q3:
    require(
        isinstance(value, list)
        and len(value) == 4
        and all(type(entry) is int for entry in value),
        "invalid Q(sqrt(3)) encoding",
    )
    require(value[1] > 0 and value[3] > 0, "invalid Q(sqrt(3)) denominator")
    return Fraction(value[0], value[1]), Fraction(value[2], value[3])


def add(left: Q3, right: Q3) -> Q3:
    return left[0] + right[0], left[1] + right[1]


def subtract(left: Q3, right: Q3) -> Q3:
    return left[0] - right[0], left[1] - right[1]


def multiply(left: Q3, right: Q3) -> Q3:
    return (
        left[0] * right[0] + 3 * left[1] * right[1],
        left[0] * right[1] + left[1] * right[0],
    )


def dot(left: Sequence[Q3], right: Sequence[Q3]) -> Q3:
    answer = ZERO
    for first, second in zip(left, right):
        answer = add(answer, multiply(first, second))
    return answer


def distance_squared(left: Sequence[Q3], right: Sequence[Q3]) -> Q3:
    difference = [subtract(first, second) for first, second in zip(left, right)]
    return dot(difference, difference)


def field_mod(value: Q3, prime: int = 13, root: int = 4) -> int:
    rational = value[0].numerator * pow(value[0].denominator, -1, prime)
    radical = value[1].numerator * pow(value[1].denominator, -1, prime)
    return (rational + root * radical) % prime


def coordinates() -> list[list[Q3]]:
    base = []
    for word in range(32):
        if not word.bit_count() & 1:
            continue
        point = [
            (Fraction(1 if word & (1 << axis) else -1), Fraction(0))
            for axis in range(5)
        ]
        base.append(point + [ZERO])
    return base + [
        [ZERO] * 5 + [(Fraction(0), Fraction(1))],
        [ZERO] * 5 + [(Fraction(0), Fraction(-1))],
    ]


def validate_adjacency(value: object, n: int, label: str) -> tuple[int, ...]:
    require(
        isinstance(value, list)
        and len(value) == n
        and all(type(row) is int for row in value),
        f"{label}: invalid adjacency",
    )
    rows = tuple(value)
    for vertex, row in enumerate(rows):
        require(0 <= row < 1 << n, f"{label}: row range")
        require(not row & (1 << vertex), f"{label}: loop")
        for other in range(n):
            require(
                bool(row & (1 << other)) == bool(rows[other] & (1 << vertex)),
                f"{label}: asymmetry",
            )
    return rows


def remove_vertex(rows: Sequence[int], deleted: int) -> tuple[int, ...]:
    lower = (1 << deleted) - 1
    return tuple(
        (row & lower) | ((row >> (deleted + 1)) << deleted)
        for vertex, row in enumerate(rows)
        if vertex != deleted
    )


def colored_isomorphism_exists(
    source: Sequence[int],
    target: Sequence[int],
    source_marked: set[int],
    target_marked: set[int],
) -> bool:
    """Check an isomorphism that maps the distinguished subset exactly."""

    n = len(source)
    if n != len(target) or len(source_marked) != len(target_marked):
        return False
    source_degree = [row.bit_count() for row in source]
    target_degree = [row.bit_count() for row in target]
    source_key = sorted((source_degree[v], v in source_marked) for v in range(n))
    target_key = sorted((target_degree[v], v in target_marked) for v in range(n))
    if source_key != target_key:
        return False
    mapping = [-1] * n
    used = 0

    def recurse(done: int) -> bool:
        nonlocal used
        if done == n:
            return True
        vertex = max(
            (candidate for candidate in range(n) if mapping[candidate] < 0),
            key=lambda candidate: sum(
                mapping[other] >= 0 and bool(source[candidate] & (1 << other))
                for other in range(n)
            ),
        )
        for image in range(n):
            if used & (1 << image):
                continue
            if source_degree[vertex] != target_degree[image]:
                continue
            if (vertex in source_marked) != (image in target_marked):
                continue
            if any(
                bool(source[vertex] & (1 << other))
                != bool(target[image] & (1 << mapping[other]))
                for other in range(n)
                if mapping[other] >= 0
            ):
                continue
            mapping[vertex] = image
            used |= 1 << image
            if recurse(done + 1):
                return True
            used ^= 1 << image
            mapping[vertex] = -1
        return False

    return recurse(0)


def matrix_rank_mod(matrix: Sequence[Sequence[int]], prime: int = 13) -> int:
    rows = [[entry % prime for entry in row] for row in matrix]
    if not rows:
        return 0
    rank = 0
    for column in range(len(rows[0])):
        pivot = None
        for candidate in range(rank, len(rows)):
            if rows[candidate][column]:
                pivot = candidate
                break
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inverse = pow(rows[rank][column], -1, prime)
        for index in range(column, len(rows[rank])):
            rows[rank][index] = rows[rank][index] * inverse % prime
        for candidate in range(rank + 1, len(rows)):
            factor = rows[candidate][column]
            if not factor:
                continue
            for index in range(column, len(rows[candidate])):
                rows[candidate][index] = (
                    rows[candidate][index] - factor * rows[rank][index]
                ) % prime
        rank += 1
        if rank == len(rows):
            break
    return rank


def affine_rank(points: Sequence[Sequence[Q3]]) -> int:
    origin = points[0]
    matrix = [
        [field_mod(subtract(value, base)) for value, base in zip(point, origin)]
        for point in points[1:]
    ]
    return matrix_rank_mod(matrix)


def rigidity_rank(rows: Sequence[int], embedding: Sequence[int]) -> int:
    points = coordinates()
    matrix = []
    for first in range(18):
        for second in range(first):
            if not rows[first] & (1 << second):
                continue
            equation = [0] * 108
            for axis in range(6):
                value = field_mod(
                    subtract(
                        points[embedding[first]][axis],
                        points[embedding[second]][axis],
                    )
                )
                equation[6 * first + axis] = value
                equation[6 * second + axis] = -value
            matrix.append(equation)
    return matrix_rank_mod(matrix)


def verify_flex_certificate(
    certificate: object,
    rows: Sequence[int],
    embedding: Sequence[int],
) -> None:
    require(isinstance(certificate, dict), "missing rank-86 flex certificate")
    seed = certificate.get("seed_K6")
    require(
        isinstance(seed, list)
        and len(seed) == 6
        and len(set(seed)) == 6
        and all(type(vertex) is int and 0 <= vertex < 18 for vertex in seed),
        "invalid flex seed",
    )
    require(
        all(rows[first] & (1 << second) for i, first in enumerate(seed) for second in seed[:i]),
        "flex seed is not K6",
    )
    points = coordinates()
    require(
        affine_rank([points[embedding[vertex]] for vertex in seed]) == 5,
        "flex seed does not affinely span its five-flat",
    )
    encoded = certificate.get("velocity_Qsqrt3")
    require(
        isinstance(encoded, list)
        and len(encoded) == 18
        and all(isinstance(row, list) and len(row) == 6 for row in encoded),
        "invalid flex velocity shape",
    )
    velocity = [[decode_q3(value) for value in row] for row in encoded]
    require(any(value != ZERO for row in velocity for value in row), "zero flex")
    require(
        all(velocity[vertex][axis] == ZERO for vertex in seed for axis in range(6)),
        "flex does not vanish on seed",
    )
    for first in range(18):
        for second in range(first):
            if not rows[first] & (1 << second):
                continue
            edge = [
                subtract(points[embedding[first]][axis], points[embedding[second]][axis])
                for axis in range(6)
            ]
            speed = [
                subtract(velocity[first][axis], velocity[second][axis])
                for axis in range(6)
            ]
            require(dot(edge, speed) == ZERO, "invalid infinitesimal flex")


def verify(report_path: Path, output_path: Path | None) -> dict:
    require(sha256(ROOT / PRODUCER) == EXPECTED_PRODUCER_SHA256, "producer hash")
    for name, wanted in EXPECTED_HASHES.items():
        require(sha256(ROOT / name) == wanted, f"input hash {name}")
    deletion_verification = json.loads((ROOT / DELETION_VERIFICATION).read_text())
    v4_verification = json.loads((ROOT / V4_VERIFICATION).read_text())
    require(deletion_verification.get("status") == "PASS", "deletion verification")
    require(v4_verification.get("status") == "PASS", "v4 verification")
    deletion = json.loads((ROOT / DELETION).read_text())
    v4 = json.loads((ROOT / V4).read_text())
    report = json.loads(report_path.read_text())
    require(report.get("schema") == "d6-standard18-rooted-attachment-audit-v1", "schema")
    require(report.get("status") == "PASS", "status")

    inputs = report.get("inputs")
    require(isinstance(inputs, dict), "inputs")
    require(
        inputs.get("deletion_manifest")
        == {"path": DELETION, "sha256": EXPECTED_HASHES[DELETION]},
        "deletion source pin",
    )
    require(
        inputs.get("deletion_verification")
        == {
            "path": DELETION_VERIFICATION,
            "sha256": EXPECTED_HASHES[DELETION_VERIFICATION],
            "status": "PASS",
        },
        "deletion verification pin",
    )
    require(
        inputs.get("v4_manifest") == {"path": V4, "sha256": EXPECTED_HASHES[V4]},
        "v4 source pin",
    )
    require(
        inputs.get("v4_verification")
        == {
            "path": V4_VERIFICATION,
            "sha256": EXPECTED_HASHES[V4_VERIFICATION],
            "status": "PASS",
        },
        "v4 verification pin",
    )

    semantics = report.get("semantics")
    require(isinstance(semantics, dict), "semantics")
    require(semantics.get("conditional_on_standard_embedding") is True, "conditional scope")
    require(semantics.get("parent_rejections") == 0, "parent overclaim")
    require(semantics.get("classification_claim_16_point_R5") is False, "R5 overclaim")
    require(semantics.get("classification_claim_18_point_R6") is False, "R6 overclaim")
    require(
        semantics.get("global_support_embedding_uniqueness_claim") is False,
        "embedding uniqueness overclaim",
    )
    require(semantics.get("nonedges") == "unconstrained and may also be unit", "nonedge semantics")

    parent_rows = {
        record["index"]: validate_adjacency(record["adjacency"], 19, "parent")
        for class_name in ("K7", "K6_only")
        for record in v4["classes"][class_name]["graphs"]
    }
    rank_report = report.get("rigidity", {}).get("records")
    case_report = report.get("rooted_cases")
    require(isinstance(rank_report, list), "rigidity records")
    require(isinstance(case_report, list), "rooted cases")
    rank_by_class = {record.get("class_id"): record for record in rank_report}
    require(len(rank_by_class) == len(rank_report), "duplicate rigidity record")
    case_by_key = {
        (
            record.get("class_id"),
            record.get("parent_index"),
            record.get("deleted_vertex"),
            record.get("witness_ordinal"),
        ): record
        for record in case_report
    }
    require(len(case_by_key) == len(case_report), "duplicate rooted case")

    expected_case_keys = set()
    compatible_classes = 0
    compatible_occurrences = 0
    rank_counts = {86: 0, 87: 0}
    parent_ranks: dict[int, set[int]] = {}
    points = coordinates()
    for manifest_ordinal, source_record in enumerate(deletion["unique_deletions"]):
        if not source_record.get("standard18_compatible"):
            continue
        compatible_classes += 1
        class_id = source_record["class_id"]
        rows = validate_adjacency(source_record["adjacency"], 18, class_id)
        witnesses = source_record["standard18_pole_pair_witnesses"]
        require(bool(witnesses), "compatible class without witness")
        embedding = witnesses[0]["embedding_permutation"]
        require(sorted(embedding) == list(range(18)), "embedding permutation")
        rank = rigidity_rank(rows, embedding)
        require(rank in (86, 87), "unexpected rigidity rank")
        rank_counts[rank] += 1
        stored_rank = rank_by_class.get(class_id)
        require(isinstance(stored_rank, dict), "missing rigidity class")
        require(stored_rank.get("manifest_ordinal") == manifest_ordinal, "rigidity ordinal")
        require(stored_rank.get("rank_lower_bound_mod_13") == rank, "modular rank")
        require(stored_rank.get("exact_rank") == rank, "exact rank")
        require(stored_rank.get("infinitesimally_rigid") is (rank == 87), "rigidity flag")
        if rank == 86:
            verify_flex_certificate(
                stored_rank.get("rank_86_upper_bound_certificate"), rows, embedding
            )
        else:
            require(
                "rank_86_upper_bound_certificate" not in stored_rank,
                "spurious rank-86 certificate",
            )

        for occurrence in source_record["occurrences"]:
            compatible_occurrences += 1
            parent = occurrence["parent_index"]
            parent_ranks.setdefault(parent, set()).add(rank)
            deletion_rows = remove_vertex(parent_rows[parent], occurrence["deleted_vertex"])
            source_attachment = {
                vertex
                for vertex in range(18)
                if occurrence["attachment_mask"] & (1 << vertex)
            }
            for witness_ordinal, witness in enumerate(witnesses):
                key = (
                    class_id,
                    parent,
                    occurrence["deleted_vertex"],
                    witness_ordinal,
                )
                expected_case_keys.add(key)
                stored = case_by_key.get(key)
                require(isinstance(stored, dict), "missing rooted witness case")
                require(stored.get("manifest_ordinal") == manifest_ordinal, "case ordinal")
                require(
                    stored.get("attachment_degree") == occurrence["attachment_degree"],
                    "attachment degree",
                )
                require(stored.get("pole_pair") == witness["pole_pair"], "pole pair")
                mapped = stored.get("mapped_attachment")
                require(
                    isinstance(mapped, list)
                    and len(mapped) == occurrence["attachment_degree"]
                    and len(set(mapped)) == len(mapped)
                    and all(type(vertex) is int and 0 <= vertex < 18 for vertex in mapped),
                    "mapped attachment",
                )
                canonical_to_standard = witness["embedding_permutation"]
                inverse = [0] * 18
                for canonical, standard in enumerate(canonical_to_standard):
                    inverse[standard] = canonical
                canonical_marked = {inverse[standard] for standard in mapped}
                require(
                    colored_isomorphism_exists(
                        deletion_rows, rows, source_attachment, canonical_marked
                    ),
                    "attachment transport lacks a colored isomorphism",
                )
                base = sorted(vertex for vertex in mapped if vertex < 16)
                poles = [vertex for vertex in mapped if vertex >= 16]
                require(len(base) in (11, 12) and len(poles) == 1, "attachment split")
                require(stored.get("mapped_base_neighbours") == base, "stored base neighbours")
                require(stored.get("mapped_pole_neighbour") == poles[0], "stored pole")
                spanning = stored.get("spanning_base_subset")
                require(
                    isinstance(spanning, list)
                    and len(spanning) == 11
                    and set(spanning).issubset(base),
                    "spanning subset",
                )
                require(affine_rank([points[vertex] for vertex in spanning]) == 5, "base span")
                centre_encoded = stored.get("forced_equidistant_centre_Qsqrt3")
                require(
                    isinstance(centre_encoded, list) and len(centre_encoded) == 6,
                    "forced centre",
                )
                centre = [decode_q3(value) for value in centre_encoded]
                distances = {distance_squared(centre, points[vertex]) for vertex in mapped}
                require(distances == {(Fraction(16, 3), Fraction(0))}, "forced radius")
                require(
                    decode_q3(stored.get("forced_common_scaled_squared_distance_Qsqrt3"))
                    == (Fraction(16, 3), Fraction(0)),
                    "stored forced radius",
                )
                require(
                    decode_q3(stored.get("required_unit_scaled_squared_distance_Qsqrt3"))
                    == (Fraction(8), Fraction(0)),
                    "stored unit radius",
                )
                require(stored.get("contradiction") == "16/3 != 8", "contradiction")

    require(set(rank_by_class) == {
        record["class_id"]
        for record in deletion["unique_deletions"]
        if record.get("standard18_compatible")
    }, "rigidity class coverage")
    require(set(case_by_key) == expected_case_keys, "rooted case coverage")
    require(compatible_classes == 14, "compatible class count")
    require(compatible_occurrences == 39, "compatible occurrence count")
    require(len(expected_case_keys) == 42, "rooted case count")
    require(rank_counts == {86: 5, 87: 9}, "rank split")
    parents_with_87 = sorted(parent for parent, ranks in parent_ranks.items() if 87 in ranks)
    parents_only_86 = sorted(parent for parent, ranks in parent_ranks.items() if ranks == {86})
    require(len(parent_ranks) == 16, "compatible parent count")
    require(len(parents_with_87) == 14, "rank-87 parent coverage")
    require(parents_only_86 == [3_950_119, 3_950_509], "rank-86-only parents")

    summary = report.get("summary")
    require(isinstance(summary, dict), "summary")
    expected_scalars = {
        "compatible_classes": 14,
        "compatible_occurrences": 39,
        "rooted_witness_representatives": 42,
        "compatible_parents": 16,
        "rank_87_classes": 9,
        "rank_86_classes": 5,
    }
    for field, wanted in expected_scalars.items():
        require(summary.get(field) == wanted, f"summary {field}")
    require(summary.get("parents_with_rank_87_deletion") == parents_with_87, "summary rank87 parents")
    require(summary.get("parents_only_rank_86_deletions") == parents_only_86, "summary rank86 parents")
    require(summary.get("rooted_cases_sha256") == object_hash(case_report), "rooted case hash")
    require(summary.get("rigidity_records_sha256") == object_hash(rank_report), "rigidity hash")

    checks = {
        "pinned_inputs": True,
        "conditional_nonrejection_semantics": True,
        "colored_attachment_transport": True,
        "exact_Qsqrt3_attachment_contradictions": True,
        "exact_rigidity_rank_split": True,
        "complete_42_case_coverage": True,
    }
    result = {
        "schema": "d6-standard18-rooted-attachment-verification-v1",
        "status": "PASS",
        "report": {"path": report_path.name, "sha256": sha256(report_path)},
        "producer_sha256": EXPECTED_PRODUCER_SHA256,
        "checks": checks,
        "counts": {
            "compatible_classes": 14,
            "compatible_occurrences": 39,
            "rooted_witness_representatives": 42,
            "compatible_parents": 16,
            "rank_87_classes": 9,
            "rank_86_classes": 5,
        },
        "exact_result": (
            "Every recorded rooted attachment fails at every represented standard-coordinate embedding orbit; no parent rejection follows without a global embedding classification."
        ),
    }
    if output_path is not None:
        temporary = output_path.with_name(
            f".{output_path.name}.{os.getpid()}.{time.time_ns()}.tmp"
        )
        with temporary.open("w", encoding="utf-8") as stream:
            json.dump(result, stream, indent=2, sort_keys=True)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, output_path)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--report",
        type=Path,
        default=ROOT / "d6_standard18_rooted_attachment_report.json",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "d6_standard18_rooted_attachment_verification.json",
    )
    args = parser.parse_args()
    print(json.dumps(verify(args.report, args.output), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
