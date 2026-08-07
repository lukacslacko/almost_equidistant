#!/usr/bin/env python3
"""Exact rooted grow-back audit at the standard 18-point coordinates.

The 18-deletion manifest records 39 occurrences whose required unit edges
embed in the known odd-halfcube-plus-two-apices configuration.  This program
transports the deleted vertex's attachment through each recorded pole-pair
witness and proves, in Q(sqrt(3)), that the attachment cannot be restored at
those coordinates.  It also certifies the rigidity-matrix rank of the 14
compatible deletion supports at their standard embeddings.

The conclusion is deliberately conditional.  A required-edge support may
have another, noncongruent realization.  Consequently this audit rejects no
19-vertex parent and makes no uniqueness claim for 16- or 18-point sets.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import os
import time
from fractions import Fraction
from pathlib import Path
from typing import Iterable, Sequence


ROOT = Path(__file__).resolve().parent
DELETION_NAME = "d6_residue_18_deletions.json"
DELETION_VERIFICATION_NAME = "d6_residue_18_deletions_verification.json"
V4_NAME = "d6_current_residue_manifest_v4.json"
V4_VERIFICATION_NAME = "d6_current_residue_manifest_v4_verification.json"

EXPECTED_DELETION_SHA256 = (
    "9d08f9ec579434c9601ad3908bd2ae51f10e09e9d06a7f38fbc25794a634b732"
)
EXPECTED_DELETION_VERIFICATION_SHA256 = (
    "50ab5d88441b2869e60a05f3e1485db0b25008c3a816960a20a3bb21c1c4defd"
)
EXPECTED_V4_SHA256 = (
    "6ab4bfffc524f5b43409d59888fb596de8130315bda3381e484bb4ce891e03e4"
)
EXPECTED_V4_VERIFICATION_SHA256 = (
    "765eceb6782d131dae8c77c9b94de78735e23a070da051c8fffbd984790e1a41"
)

EXPECTED_COUNTS = {
    "compatible_classes": 14,
    "compatible_occurrences": 39,
    "rooted_witness_representatives": 42,
    "compatible_parents": 16,
    "rank_87_classes": 9,
    "rank_86_classes": 5,
    "parents_with_rank_87_deletion": 14,
    "parents_only_rank_86_deletions": 2,
}

RANK_86_CLASS_IDS = {
    "u18-19c7cfbd368fafcff3fb5598e1cfa875bf71ca7f0f297be10e846c92c4bd133a",
    "u18-761add09c775ad5ba966a0eef97848e076faacc647744c4079b7bee48a65d92f",
    "u18-3f9adb1f1bbf8e73d3d434a0f562c4a4e623cccbfa0b7cae1c16f548d4cc1711",
    "u18-d58d9ae1ccd1aa28d3a5c25ede20a9135c73854300eb13d64ce74eb6f04854ee",
    "u18-b6a9b81054d4363d0db724af62235050065141d531e1093b6cecda4cc9d93267",
}


# For each rank-86 framework, a non-rigid infinitesimal motion that vanishes
# on a unit K6.  Values are pairs (rational part, sqrt(3) coefficient).
# The report stores and verifies the resulting full 18-by-6 velocity array.
FLEX_DATA: dict[str, dict] = {
    "u18-19c7cfbd368fafcff3fb5598e1cfa875bf71ca7f0f297be10e846c92c4bd133a": {
        "seed": (0, 2, 3, 8, 10, 11),
        "entries": {
            1: {5: (0, Fraction(2, 3))},
            5: {5: (0, Fraction(2, 3))},
            6: {5: (0, Fraction(2, 3))},
            7: {5: (0, Fraction(2, 3))},
            9: {1: (1, 0), 2: (1, 0)},
        },
    },
    "u18-761add09c775ad5ba966a0eef97848e076faacc647744c4079b7bee48a65d92f": {
        "seed": (0, 4, 5, 6, 8, 16),
        "entries": {
            1: {5: (2, 0)},
            2: {5: (2, 0)},
            3: {5: (2, 0)},
            7: {5: (2, 0)},
            9: {4: (0, 1), 5: (1, 0)},
        },
    },
    "u18-3f9adb1f1bbf8e73d3d434a0f562c4a4e623cccbfa0b7cae1c16f548d4cc1711": {
        "seed": (0, 8, 9, 10, 11, 12),
        "entries": {
            4: {5: (2, 0)},
            5: {5: (2, 0)},
            6: {5: (2, 0)},
            7: {5: (2, 0)},
            17: {4: (0, -1), 5: (1, 0)},
        },
    },
    "u18-d58d9ae1ccd1aa28d3a5c25ede20a9135c73854300eb13d64ce74eb6f04854ee": {
        "seed": (0, 4, 5, 6, 7, 17),
        "entries": {
            8: {
                0: (0, Fraction(1, 4)),
                1: (0, Fraction(-1, 4)),
                2: (0, Fraction(-1, 4)),
                3: (0, Fraction(-1, 4)),
                4: (0, Fraction(1, 4)),
                5: (Fraction(-1, 4), 0),
            },
            9: {5: (1, 0)},
        },
    },
    "u18-b6a9b81054d4363d0db724af62235050065141d531e1093b6cecda4cc9d93267": {
        "seed": (0, 1, 8, 9, 10, 13),
        "entries": {
            3: {5: (4, 0)},
            4: {5: (4, 0)},
            5: {5: (4, 0)},
            6: {5: (4, 0)},
            7: {5: (4, 0)},
            17: {
                0: (0, -1),
                1: (0, -1),
                2: (0, 1),
                3: (0, 1),
                4: (0, -1),
                5: (1, 0),
            },
        },
    },
}


Q3 = tuple[Fraction, Fraction]
ZERO: Q3 = (Fraction(0), Fraction(0))
ONE: Q3 = (Fraction(1), Fraction(0))


def q3(rational: int | Fraction = 0, radical: int | Fraction = 0) -> Q3:
    return Fraction(rational), Fraction(radical)


def qadd(left: Q3, right: Q3) -> Q3:
    return left[0] + right[0], left[1] + right[1]


def qneg(value: Q3) -> Q3:
    return -value[0], -value[1]


def qsub(left: Q3, right: Q3) -> Q3:
    return qadd(left, qneg(right))


def qmul(left: Q3, right: Q3) -> Q3:
    return (
        left[0] * right[0] + 3 * left[1] * right[1],
        left[0] * right[1] + left[1] * right[0],
    )


def qencode(value: Q3) -> list[int]:
    return [
        value[0].numerator,
        value[0].denominator,
        value[1].numerator,
        value[1].denominator,
    ]


def qmod(value: Q3, prime: int, root: int) -> int:
    rational = value[0].numerator * pow(value[0].denominator, -1, prime)
    radical = value[1].numerator * pow(value[1].denominator, -1, prime)
    return (rational + root * radical) % prime


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


def validate_rows(rows: object, vertices: int, label: str) -> tuple[int, ...]:
    if (
        not isinstance(rows, list)
        or len(rows) != vertices
        or any(type(row) is not int for row in rows)
    ):
        raise ValueError(f"{label}: invalid adjacency container")
    result = tuple(rows)
    bound = 1 << vertices
    for vertex, row in enumerate(result):
        if row < 0 or row >= bound or row & (1 << vertex):
            raise ValueError(f"{label}: invalid row {vertex}")
        for other in range(vertices):
            if bool(row & (1 << other)) != bool(result[other] & (1 << vertex)):
                raise ValueError(f"{label}: asymmetric at {vertex},{other}")
    return result


def delete_vertex(rows: Sequence[int], deleted: int) -> tuple[int, ...]:
    low_mask = (1 << deleted) - 1
    answer = []
    for vertex, row in enumerate(rows):
        if vertex == deleted:
            continue
        answer.append((row & low_mask) | ((row >> (deleted + 1)) << deleted))
    return tuple(answer)


def graph_isomorphism(source: Sequence[int], target: Sequence[int]) -> list[int] | None:
    if len(source) != len(target):
        return None
    n = len(source)
    source_degrees = [row.bit_count() for row in source]
    target_degrees = [row.bit_count() for row in target]
    if sorted(source_degrees) != sorted(target_degrees):
        return None
    mapping = [-1] * n
    used = 0

    def search(done: int) -> bool:
        nonlocal used
        if done == n:
            return True
        best = max(
            (vertex for vertex in range(n) if mapping[vertex] < 0),
            key=lambda vertex: sum(
                mapping[other] >= 0 and bool(source[vertex] & (1 << other))
                for other in range(n)
            ),
        )
        for image in range(n):
            if used & (1 << image) or target_degrees[image] != source_degrees[best]:
                continue
            if any(
                bool(source[best] & (1 << other))
                != bool(target[image] & (1 << mapping[other]))
                for other in range(n)
                if mapping[other] >= 0
            ):
                continue
            mapping[best] = image
            used |= 1 << image
            if search(done + 1):
                return True
            used ^= 1 << image
            mapping[best] = -1
        return False

    return mapping if search(0) else None


def sign_vectors() -> list[tuple[int, ...]]:
    return [
        tuple(1 if word & (1 << axis) else -1 for axis in range(5))
        for word in range(32)
        if word.bit_count() & 1
    ]


def standard_coordinates() -> list[list[Q3]]:
    points = [[q3(value) for value in vector] + [ZERO] for vector in sign_vectors()]
    points.append([ZERO] * 5 + [q3(0, 1)])
    points.append([ZERO] * 5 + [q3(0, -1)])
    return points


def qdot(left: Sequence[Q3], right: Sequence[Q3]) -> Q3:
    answer = ZERO
    for first, second in zip(left, right):
        answer = qadd(answer, qmul(first, second))
    return answer


def squared_distance(left: Sequence[Q3], right: Sequence[Q3]) -> Q3:
    difference = [qsub(first, second) for first, second in zip(left, right)]
    return qdot(difference, difference)


def row_rank_mod(matrix: Sequence[Sequence[int]], prime: int) -> int:
    rows = [[entry % prime for entry in row] for row in matrix]
    if not rows:
        return 0
    rank = 0
    for column in range(len(rows[0])):
        pivot = next(
            (candidate for candidate in range(rank, len(rows)) if rows[candidate][column]),
            None,
        )
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inverse = pow(rows[rank][column], -1, prime)
        rows[rank] = [(entry * inverse) % prime for entry in rows[rank]]
        for candidate in range(len(rows)):
            if candidate == rank or not rows[candidate][column]:
                continue
            factor = rows[candidate][column]
            rows[candidate] = [
                (left - factor * right) % prime
                for left, right in zip(rows[candidate], rows[rank])
            ]
        rank += 1
        if rank == len(rows):
            break
    return rank


def affine_rank_mod(points: Sequence[Sequence[Q3]], prime: int = 13, root: int = 4) -> int:
    origin = points[0]
    matrix = [
        [qmod(qsub(value, base), prime, root) for value, base in zip(point, origin)]
        for point in points[1:]
    ]
    return row_rank_mod(matrix, prime)


def rigidity_rank_mod13(rows: Sequence[int], embedding: Sequence[int]) -> int:
    points = standard_coordinates()
    matrix: list[list[int]] = []
    for first in range(18):
        for second in range(first):
            if not rows[first] & (1 << second):
                continue
            equation = [0] * 108
            for axis in range(6):
                difference = qsub(
                    points[embedding[first]][axis], points[embedding[second]][axis]
                )
                value = qmod(difference, 13, 4)
                equation[6 * first + axis] = value
                equation[6 * second + axis] = -value
            matrix.append(equation)
    return row_rank_mod(matrix, 13)


def is_clique(rows: Sequence[int], vertices: Sequence[int]) -> bool:
    return all(
        rows[first] & (1 << second)
        for position, first in enumerate(vertices)
        for second in vertices[:position]
    )


def flex_array(class_id: str) -> list[list[Q3]]:
    velocities = [[ZERO for _ in range(6)] for _ in range(18)]
    for vertex, entries in FLEX_DATA[class_id]["entries"].items():
        for axis, value in entries.items():
            velocities[vertex][axis] = q3(*value)
    return velocities


def verify_flex(
    rows: Sequence[int], embedding: Sequence[int], class_id: str
) -> tuple[tuple[int, ...], list[list[Q3]]]:
    points = standard_coordinates()
    velocities = flex_array(class_id)
    seed = tuple(FLEX_DATA[class_id]["seed"])
    if not is_clique(rows, seed):
        raise AssertionError(f"{class_id}: flex seed is not a K6")
    if affine_rank_mod([points[embedding[vertex]] for vertex in seed]) != 5:
        raise AssertionError(f"{class_id}: flex seed does not span a five-flat")
    if all(value == ZERO for velocity in velocities for value in velocity):
        raise AssertionError(f"{class_id}: zero flex")
    if any(velocities[vertex][axis] != ZERO for vertex in seed for axis in range(6)):
        raise AssertionError(f"{class_id}: flex does not vanish on its K6")
    for first in range(18):
        for second in range(first):
            if not rows[first] & (1 << second):
                continue
            edge = [
                qsub(points[embedding[first]][axis], points[embedding[second]][axis])
                for axis in range(6)
            ]
            velocity = [
                qsub(velocities[first][axis], velocities[second][axis])
                for axis in range(6)
            ]
            if qdot(edge, velocity) != ZERO:
                raise AssertionError(f"{class_id}: flex violates edge {first},{second}")
    return seed, velocities


def load_inputs(
    deletion_path: Path, v4_path: Path
) -> tuple[dict, dict, dict, dict]:
    deletion_verification_path = deletion_path.with_name(DELETION_VERIFICATION_NAME)
    v4_verification_path = v4_path.with_name(V4_VERIFICATION_NAME)
    expected = (
        (deletion_path, EXPECTED_DELETION_SHA256),
        (deletion_verification_path, EXPECTED_DELETION_VERIFICATION_SHA256),
        (v4_path, EXPECTED_V4_SHA256),
        (v4_verification_path, EXPECTED_V4_VERIFICATION_SHA256),
    )
    for path, wanted in expected:
        if sha256(path) != wanted:
            raise ValueError(f"pinned input hash mismatch: {path.name}")
    deletion_verification = json.loads(deletion_verification_path.read_text())
    v4_verification = json.loads(v4_verification_path.read_text())
    if deletion_verification.get("status") != "PASS":
        raise ValueError("deletion manifest verification is not PASS")
    if v4_verification.get("status") != "PASS":
        raise ValueError("v4 verification is not PASS")
    deletion = json.loads(deletion_path.read_text())
    v4 = json.loads(v4_path.read_text())
    if deletion.get("schema") != "d6-residue-18-deletion-manifest-v1":
        raise ValueError("unexpected deletion manifest schema")
    if v4.get("schema") != "d6-current-exact-residue-v4":
        raise ValueError("unexpected v4 schema")
    return deletion, deletion_verification, v4, v4_verification


def build_report(
    deletion_path: Path = ROOT / DELETION_NAME,
    v4_path: Path = ROOT / V4_NAME,
) -> dict:
    deletion, deletion_verification, v4, v4_verification = load_inputs(
        deletion_path, v4_path
    )
    parents = {
        record["index"]: validate_rows(
            record["adjacency"], 19, f"parent {record['index']}"
        )
        for class_name in ("K7", "K6_only")
        for record in v4["classes"][class_name]["graphs"]
    }
    points = standard_coordinates()
    rooted_cases = []
    rank_records = []
    parent_ranks: dict[int, set[int]] = {}
    compatible_classes = 0
    compatible_occurrences = 0

    for ordinal, record in enumerate(deletion["unique_deletions"]):
        if not record.get("standard18_compatible"):
            continue
        compatible_classes += 1
        class_id = record["class_id"]
        rows = validate_rows(record["adjacency"], 18, class_id)
        witnesses = record["standard18_pole_pair_witnesses"]
        if not witnesses:
            raise AssertionError(f"{class_id}: missing standard witness")
        embedding = witnesses[0]["embedding_permutation"]
        if sorted(embedding) != list(range(18)):
            raise AssertionError(f"{class_id}: invalid embedding permutation")
        rank_mod13 = rigidity_rank_mod13(rows, embedding)
        expected_rank = 86 if class_id in RANK_86_CLASS_IDS else 87
        if rank_mod13 != expected_rank:
            raise AssertionError(
                f"{class_id}: rank drift {rank_mod13} != {expected_rank}"
            )
        rank_record = {
            "manifest_ordinal": ordinal,
            "class_id": class_id,
            "edges": record["edges"],
            "witness_count": len(witnesses),
            "rank_lower_bound_mod_13": rank_mod13,
            "sqrt3_mod_13": 4,
            "rank_upper_bound_from_rigid_motions": 87,
            "exact_rank": expected_rank,
            "infinitesimally_rigid": expected_rank == 87,
        }
        if expected_rank == 86:
            seed, velocities = verify_flex(rows, embedding, class_id)
            rank_record["rank_86_upper_bound_certificate"] = {
                "kind": "nontrivial infinitesimal flex vanishing on an affinely spanning unit K6",
                "seed_K6": list(seed),
                "velocity_Qsqrt3": [
                    [qencode(value) for value in velocity] for velocity in velocities
                ],
            }
        rank_records.append(rank_record)

        for occurrence in record["occurrences"]:
            compatible_occurrences += 1
            parent_index = occurrence["parent_index"]
            parent_ranks.setdefault(parent_index, set()).add(expected_rank)
            source_rows = delete_vertex(
                parents[parent_index], occurrence["deleted_vertex"]
            )
            source_to_canonical = graph_isomorphism(source_rows, rows)
            if source_to_canonical is None:
                raise AssertionError(f"{class_id}: occurrence isomorphism failed")
            attachment = [
                vertex
                for vertex in range(18)
                if occurrence["attachment_mask"] & (1 << vertex)
            ]
            for witness_ordinal, witness in enumerate(witnesses):
                canonical_to_standard = witness["embedding_permutation"]
                mapped = sorted(
                    canonical_to_standard[source_to_canonical[vertex]]
                    for vertex in attachment
                )
                base = [vertex for vertex in mapped if vertex < 16]
                poles = [vertex for vertex in mapped if vertex >= 16]
                if len(base) not in (11, 12) or len(poles) != 1:
                    raise AssertionError(
                        f"{class_id}: unexpected standard attachment split"
                    )
                spanning = base[:11]
                if affine_rank_mod([points[vertex] for vertex in spanning]) != 5:
                    raise AssertionError(f"{class_id}: base neighbours do not span")
                pole = poles[0]
                # Pole 16 is at +sqrt(3), pole 17 at -sqrt(3).  Equality
                # with the base spheres forces the opposite signed height.
                centre_sign = -1 if pole == 16 else 1
                centre = [ZERO] * 5 + [q3(0, Fraction(centre_sign, 3))]
                distances = {squared_distance(centre, points[vertex]) for vertex in mapped}
                if distances != {q3(Fraction(16, 3))}:
                    raise AssertionError(f"{class_id}: rooted radius drift")
                rooted_cases.append(
                    {
                        "case_id": stable_hash(
                            {
                                "class_id": class_id,
                                "parent_index": parent_index,
                                "deleted_vertex": occurrence["deleted_vertex"],
                                "witness_ordinal": witness_ordinal,
                            }
                        ),
                        "class_id": class_id,
                        "manifest_ordinal": ordinal,
                        "parent_class": occurrence["parent_class"],
                        "parent_index": parent_index,
                        "deleted_vertex": occurrence["deleted_vertex"],
                        "attachment_degree": occurrence["attachment_degree"],
                        "witness_ordinal": witness_ordinal,
                        "pole_pair": witness["pole_pair"],
                        "mapped_attachment": mapped,
                        "mapped_base_neighbours": base,
                        "mapped_pole_neighbour": pole,
                        "spanning_base_subset": spanning,
                        "spanning_base_affine_rank": 5,
                        "forced_equidistant_centre_Qsqrt3": [
                            qencode(value) for value in centre
                        ],
                        "forced_common_scaled_squared_distance_Qsqrt3": qencode(
                            q3(Fraction(16, 3))
                        ),
                        "required_unit_scaled_squared_distance_Qsqrt3": qencode(q3(8)),
                        "contradiction": "16/3 != 8",
                    }
                )

    parents_with_87 = sorted(
        parent for parent, ranks in parent_ranks.items() if 87 in ranks
    )
    parents_only_86 = sorted(
        parent for parent, ranks in parent_ranks.items() if ranks == {86}
    )
    observed = {
        "compatible_classes": compatible_classes,
        "compatible_occurrences": compatible_occurrences,
        "rooted_witness_representatives": len(rooted_cases),
        "compatible_parents": len(parent_ranks),
        "rank_87_classes": sum(record["exact_rank"] == 87 for record in rank_records),
        "rank_86_classes": sum(record["exact_rank"] == 86 for record in rank_records),
        "parents_with_rank_87_deletion": len(parents_with_87),
        "parents_only_rank_86_deletions": len(parents_only_86),
    }
    if observed != EXPECTED_COUNTS:
        raise AssertionError(f"count drift: {observed!r}")
    if parents_only_86 != [3_950_119, 3_950_509]:
        raise AssertionError("rank-86-only parent list drift")

    rank_records.sort(key=lambda item: item["manifest_ordinal"])
    rooted_cases.sort(
        key=lambda item: (
            item["manifest_ordinal"],
            item["parent_index"],
            item["deleted_vertex"],
            item["witness_ordinal"],
        )
    )
    return {
        "schema": "d6-standard18-rooted-attachment-audit-v1",
        "status": "PASS",
        "inputs": {
            "deletion_manifest": {
                "path": deletion_path.name,
                "sha256": EXPECTED_DELETION_SHA256,
            },
            "deletion_verification": {
                "path": DELETION_VERIFICATION_NAME,
                "sha256": EXPECTED_DELETION_VERIFICATION_SHA256,
                "status": deletion_verification["status"],
            },
            "v4_manifest": {"path": v4_path.name, "sha256": EXPECTED_V4_SHA256},
            "v4_verification": {
                "path": V4_VERIFICATION_NAME,
                "sha256": EXPECTED_V4_VERIFICATION_SHA256,
                "status": v4_verification["status"],
            },
        },
        "coordinate_model": {
            "scale": "coordinates are 2*sqrt(2) times actual coordinates",
            "base": "the 16 odd sign vectors in the first five coordinates",
            "poles": ["(0,0,0,0,0,+sqrt(3))", "(0,0,0,0,0,-sqrt(3))"],
            "unit_scaled_squared_distance": 8,
            "field": "Q(sqrt(3)); no floating-point arithmetic",
        },
        "summary": {
            **observed,
            "rooted_cases_sha256": stable_hash(rooted_cases),
            "rigidity_records_sha256": stable_hash(rank_records),
            "parents_with_rank_87_deletion": parents_with_87,
            "parents_only_rank_86_deletions": parents_only_86,
        },
        "rigidity": {
            "records": rank_records,
            "rank_87_argument": (
                "rank 87 modulo 13 (sqrt(3)=4) is an exact lower bound; "
                "the 21-dimensional rigid-motion kernel gives the matching upper bound"
            ),
            "rank_86_argument": (
                "rank 86 modulo 13 is an exact lower bound; the archived nonzero "
                "flex vanishes on an affinely spanning K6, so it is independent of "
                "all rigid motions and gives the matching upper bound"
            ),
        },
        "rooted_cases": rooted_cases,
        "exact_conclusion": (
            "For every archived occurrence and standard pole-pair embedding orbit, "
            "the deleted vertex's required unit neighbours cannot be restored at "
            "the standard coordinates."
        ),
        "semantics": {
            "conditional_on_standard_embedding": True,
            "base_automorphism_invariance": (
                "The proof uses only one pole neighbour, at least 11 base neighbours, "
                "and the fact that every 11 base vertices span; it therefore covers "
                "all base automorphisms and the pole swap represented by each witness orbit."
            ),
            "parent_rejections": 0,
            "classification_claim_16_point_R5": False,
            "classification_claim_18_point_R6": False,
            "global_support_embedding_uniqueness_claim": False,
            "nonedges": "unconstrained and may also be unit",
        },
        "nonclaims": [
            "A standard-compatible support may have a noncongruent realization.",
            "Infinitesimal rigidity is local and does not prove global uniqueness.",
            "The five rank-86 standard frameworks have one extra infinitesimal flex; this alone neither proves nor disproves a finite flex.",
            "No 19-vertex graph is rejected by this conditional audit alone.",
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
    parser.add_argument("--deletions", type=Path, default=ROOT / DELETION_NAME)
    parser.add_argument("--v4", type=Path, default=ROOT / V4_NAME)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "d6_standard18_rooted_attachment_report.json",
    )
    args = parser.parse_args()
    report = build_report(args.deletions, args.v4)
    atomic_json(args.output, report)
    print(
        json.dumps(
            {
                "status": report["status"],
                "output": str(args.output),
                "summary": report["summary"],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
