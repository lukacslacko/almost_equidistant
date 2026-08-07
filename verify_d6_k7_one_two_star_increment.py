#!/usr/bin/env python3
"""Independent checker for the K7 one/two-free star increment.

This checker imports neither the production star kernel nor its corpus
builder.  It independently rebuilds every currently passing propagated
family on the frozen 19 K7 graphs, reconstructs all correlated sign
assignments, and reclassifies every affine-line/two-free-hyperbola system in
exact Q(sqrt(7)) arithmetic.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import subprocess
import sys
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from datetime import UTC, datetime
from fractions import Fraction
from pathlib import Path
from typing import Sequence

import d6_k7_rank_reference as reference
import verify_d6_k7_one_free_conjunction as current_independent
from verify_profile_d6 import lower_bound_18_graph


ROOT = Path(__file__).resolve().parent
MANIFEST = ROOT / "d6_current_residue_manifest_v6.json"
MANIFEST_VERIFICATION = ROOT / "d6_current_residue_manifest_v6_verification.json"
CURRENT = ROOT / "d6_k7_one_free_conjunction_report.json"
CURRENT_VERIFICATION = ROOT / "d6_k7_one_free_conjunction_verification.json"
EXPECTED_UPSTREAM = {
    MANIFEST.name: (
        "c513bcc40c037af23340dce3cc91d593a101dd1d073f93feab82abd556efadc9"
    ),
    MANIFEST_VERIFICATION.name: (
        "3d11cbd20cb59b42a8c8fd5252fb3ef3dea649576982ea8017622fad3279a7bc"
    ),
    CURRENT.name: (
        "181f63015d373fb69f37f4f390cdd8ad50bc32b5a2d97a36fd0bac876d0877a2"
    ),
    CURRENT_VERIFICATION.name: (
        "ebaacfe663728eff427c9ba2280fc4261d26dfe10fad21e1a00982b9d31ac062"
    ),
}
EXPECTED_INPUT_SHA256 = (
    "af25b842ca9eaa6e9721f51d0400be2436d930b72cc33c72bc366f0ac5e41290"
)
EXPECTED_REJECTIONS = [2592657, 3785980, 3888410]
EXPECTED_BUILDER_SHA256 = (
    "3dc6003176903d457e4c2fd697eec9d1a4e036b052afd1359cb8c800b182a41a"
)
EXPECTED_KERNEL_SHA256 = (
    "d4cd813b66a0c57dd106bb9987d1417366152f4cd9162cfd9b0486e3fd0d0b5b"
)
EXPECTED_CURRENT_CHECKER_SOURCES = {
    "verify_d6_k7_double_pin_conjunction.py": (
        "41026998b065ae0600bd0ba2209f886bec465a2b9a59153114ccb34b4107f49a"
    ),
    "verify_d6_k7_full_pin_increment.py": (
        "4f1599c8ac6d7fe83db604b70bc9e551102805d7dd54ceb580eb69d5bbca4464"
    ),
    "verify_d6_k7_one_free_conjunction.py": (
        "105454f7ef2fbafcddffb4e7a83d32c344ce31ad442ad22894b6eac90d6b0848"
    ),
    "verify_d6_k7_sparse_value_full.py": (
        "6c0012e66d0b89e280f3cebed7a61f24cc1c6e128e726663a7b246e7da0edbbc"
    ),
    "verify_d6_k7_support_full.py": (
        "c6b9dce47618bd9a4247eb5a9f093682265447729e1868c4dc0677a24f51c628"
    ),
}
SOURCE_FILES = (
    "build_d6_k7_one_two_star_increment.py",
    "d6_k7_one_two_star.py",
    "d6_k7_one_free_edge.py",
    "d6_k7_correlated_one_free_edge.py",
    "d6_k7_two_defect_double_pin.py",
    "d6_k7_full_pin_odd_cycle.py",
    "d6_k7_support_propagation.py",
    "d6_k7_small_support_value.py",
    "d6_k7_rank_reference.py",
    "verify_profile_d6.py",
)

Quadratic = tuple[Fraction, Fraction]
ComponentKey = tuple[int, tuple[int, ...]]
ZERO: Quadratic = (Fraction(0), Fraction(0))
ONE: Quadratic = (Fraction(1), Fraction(0))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_hash(value: object) -> str:
    return hashlib.sha256(
        json.dumps(
            value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
        ).encode("ascii")
    ).hexdigest()


def atomic_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def q(value: int | Fraction) -> Quadratic:
    return Fraction(value), Fraction(0)


def qadd(first: Quadratic, second: Quadratic) -> Quadratic:
    return first[0] + second[0], first[1] + second[1]


def qneg(value: Quadratic) -> Quadratic:
    return -value[0], -value[1]


def qsub(first: Quadratic, second: Quadratic) -> Quadratic:
    return qadd(first, qneg(second))


def qmul(first: Quadratic, second: Quadratic) -> Quadratic:
    a, b = first
    c, d = second
    return a * c + 7 * b * d, a * d + b * c


def qscale(value: Quadratic, scalar: int | Fraction) -> Quadratic:
    factor = Fraction(scalar)
    return value[0] * factor, value[1] * factor


def qinv(value: Quadratic) -> Quadratic:
    a, b = value
    denominator = a * a - 7 * b * b
    if denominator == 0:
        raise ZeroDivisionError("zero Q(sqrt(7)) denominator")
    return a / denominator, -b / denominator


def qdiv(first: Quadratic, second: Quadratic) -> Quadratic:
    return qmul(first, qinv(second))


def qsign(value: Quadratic) -> int:
    a, b = value
    if b == 0:
        return (a > 0) - (a < 0)
    if a == 0:
        return (b > 0) - (b < 0)
    if (a > 0) == (b > 0):
        return 1 if a > 0 else -1
    comparison = a * a - 7 * b * b
    if comparison == 0:
        raise AssertionError("sqrt(7) would be rational")
    if a > 0:
        return 1 if comparison > 0 else -1
    return -1 if comparison > 0 else 1


def qjson(value: Quadratic) -> list[list[int]]:
    return [
        [value[0].numerator, value[0].denominator],
        [value[1].numerator, value[1].denominator],
    ]


def component_key(entry: current_independent.PinComponent) -> ComponentKey:
    return entry[0], entry[1]


def diagonal_two_free(
    mask_size: int, pinned_sum: int
) -> tuple[Quadratic, Quadratic]:
    if not 2 <= mask_size <= 7:
        raise ValueError("two-free mask size")
    if pinned_sum not in range(-(mask_size - 2), mask_size - 1, 2):
        raise ValueError("pinned sum")
    return (
        (Fraction(-pinned_sum), Fraction(1)),
        (
            Fraction(mask_size + 4 + pinned_sum * pinned_sum, 2),
            Fraction(-pinned_sum),
        ),
    )


def pin_keys(
    entries: Sequence[current_independent.PinComponent],
) -> dict[int, ComponentKey]:
    result = {entry[0]: component_key(entry) for entry in entries}
    if len(result) != len(entries):
        raise ValueError("duplicate pin coordinate")
    return result


def assigned_signs(
    keys: Sequence[ComponentKey], assignment: int
) -> dict[ComponentKey, int]:
    return {
        key: 1 if assignment & (1 << position) else -1
        for position, key in enumerate(keys)
    }


def one_free_values(
    mask: int,
    free_mask: int,
    keys_by_coordinate: dict[int, ComponentKey],
    signs: dict[ComponentKey, int],
) -> dict[int, Quadratic]:
    if free_mask.bit_count() != 1:
        raise ValueError("not one-free")
    values = {
        coordinate: q(signs[key])
        for coordinate, key in keys_by_coordinate.items()
    }
    coordinate = free_mask.bit_length() - 1
    pinned_sum = sum(signs[key] for key in keys_by_coordinate.values())
    values[coordinate] = current_independent.one_free_value(
        mask.bit_count(), pinned_sum
    )
    return values


def edge_line(
    center: int,
    neighbour: int,
    masks: Sequence[int],
    free_masks: Sequence[int],
    keys_by_vertex: dict[int, dict[int, ComponentKey]],
    signs: dict[ComponentKey, int],
) -> tuple[Quadratic, Quadratic, Quadratic]:
    values = one_free_values(
        masks[neighbour],
        free_masks[neighbour],
        keys_by_vertex[neighbour],
        signs,
    )
    coordinates = tuple(current_independent.bits(free_masks[center]))
    coefficients = {coordinate: ZERO for coordinate in coordinates}
    constant = ZERO
    for coordinate in current_independent.bits(masks[center] & masks[neighbour]):
        value = values[coordinate]
        if free_masks[center] & (1 << coordinate):
            coefficients[coordinate] = qadd(coefficients[coordinate], value)
        else:
            constant = qadd(
                constant,
                qscale(value, signs[keys_by_vertex[center][coordinate]]),
            )
    return (
        coefficients[coordinates[0]],
        coefficients[coordinates[1]],
        qsub(ONE, constant),
    )


def one_line_decision(
    left: Quadratic,
    right: Quadratic,
    target: Quadratic,
    center: Quadratic,
    radius: Quadratic,
) -> dict:
    if left == ZERO and right == ZERO:
        return {
            "feasible": target == ZERO,
            "reason": (
                "rank0_free_hyperbola"
                if target == ZERO
                else "zero_line_inconsistent"
            ),
            "rank": 0,
            "discriminant": None,
            "solution": None,
        }
    if left == ZERO or right == ZERO:
        coefficient = right if left == ZERO else left
        fixed = qdiv(target, coefficient)
        feasible = fixed != center or radius == ZERO
        return {
            "feasible": feasible,
            "reason": (
                "rank1_fixed_coordinate"
                if feasible
                else "rank1_axis_asymptote"
            ),
            "rank": 1,
            "discriminant": None,
            "solution": None,
        }
    shifted = qsub(target, qmul(qadd(left, right), center))
    discriminant = qsub(
        qmul(shifted, shifted), qscale(qmul(qmul(left, right), radius), 4)
    )
    sign = qsign(discriminant)
    return {
        "feasible": sign >= 0,
        "reason": (
            "rank1_negative_discriminant"
            if sign < 0
            else "rank1_tangent"
            if sign == 0
            else "rank1_secant"
        ),
        "rank": 1,
        "discriminant": qjson(discriminant),
        "solution": None,
    }


def line_system_decision(
    lines: Sequence[tuple[Quadratic, Quadratic, Quadratic]],
    center: Quadratic,
    radius: Quadratic,
) -> dict:
    nontrivial = []
    for left, right, target in lines:
        if left == ZERO and right == ZERO:
            if target != ZERO:
                return {
                    "feasible": False,
                    "reason": "zero_line_inconsistent",
                    "rank": 0,
                    "discriminant": None,
                    "solution": None,
                }
            continue
        nontrivial.append((left, right, target))
    if not nontrivial:
        return {
            "feasible": True,
            "reason": "rank0_free_hyperbola",
            "rank": 0,
            "discriminant": None,
            "solution": None,
        }

    first_left, first_right, first_target = nontrivial[0]
    solution = None
    for left, right, target in nontrivial[1:]:
        determinant = qsub(qmul(first_left, right), qmul(first_right, left))
        if determinant != ZERO and solution is None:
            solution = (
                qdiv(
                    qsub(qmul(first_target, right), qmul(first_right, target)),
                    determinant,
                ),
                qdiv(
                    qsub(qmul(first_left, target), qmul(first_target, left)),
                    determinant,
                ),
            )
            continue
        if determinant == ZERO and (
            qsub(qmul(first_left, target), qmul(first_target, left)) != ZERO
            or qsub(qmul(first_right, target), qmul(first_target, right)) != ZERO
        ):
            return {
                "feasible": False,
                "reason": (
                    "rank2_inconsistent" if solution is not None
                    else "parallel_inconsistent"
                ),
                "rank": 2 if solution is not None else 1,
                "discriminant": None,
                "solution": (
                    None if solution is None else [qjson(value) for value in solution]
                ),
            }
    if solution is None:
        return one_line_decision(
            first_left, first_right, first_target, center, radius
        )
    x_value, y_value = solution
    serialized = [qjson(value) for value in solution]
    if any(
        qadd(qmul(left, x_value), qmul(right, y_value)) != target
        for left, right, target in nontrivial
    ):
        return {
            "feasible": False,
            "reason": "rank2_inconsistent",
            "rank": 2,
            "discriminant": None,
            "solution": serialized,
        }
    feasible = qmul(qsub(x_value, center), qsub(y_value, center)) == radius
    return {
        "feasible": feasible,
        "reason": "rank2_on_diagonal" if feasible else "rank2_diagonal_mismatch",
        "rank": 2,
        "discriminant": None,
        "solution": serialized,
    }


def independent_star_locator(
    graph: Sequence[int], masks: Sequence[int]
) -> dict | None:
    required, supports = current_independent.validate_branch(graph, masks)
    pinned, data = current_independent.pinning_data(required, supports)
    free_masks = tuple(
        mask & ~pinned_mask for mask, pinned_mask in zip(supports, pinned)
    )
    for center, free_mask in enumerate(free_masks):
        if free_mask.bit_count() != 2:
            continue
        neighbours = tuple(
            neighbour
            for neighbour in current_independent.bits(required[center])
            if free_masks[neighbour].bit_count() == 1
        )
        if len(neighbours) < 2:
            continue
        relevant = (center, *neighbours)
        components = {
            vertex: tuple(
                data[(vertex, coordinate)]
                for coordinate in current_independent.bits(
                    supports[vertex] ^ free_masks[vertex]
                )
            )
            for vertex in relevant
        }
        keys_by_vertex = {
            vertex: pin_keys(components[vertex]) for vertex in relevant
        }
        keys = tuple(
            sorted(
                {
                    key
                    for vertex in relevant
                    for key in keys_by_vertex[vertex].values()
                }
            )
        )
        failures = []
        for assignment in range(1 << len(keys)):
            signs = assigned_signs(keys, assignment)
            lines = tuple(
                edge_line(
                    center,
                    neighbour,
                    supports,
                    free_masks,
                    keys_by_vertex,
                    signs,
                )
                for neighbour in neighbours
            )
            pinned_sum = sum(
                signs[key] for key in keys_by_vertex[center].values()
            )
            diagonal_center, radius = diagonal_two_free(
                supports[center].bit_count(), pinned_sum
            )
            decision = line_system_decision(lines, diagonal_center, radius)
            if decision["feasible"]:
                break
            failures.append(
                {
                    "assignment": assignment,
                    "reason": decision["reason"],
                    "rank": decision["rank"],
                    "discriminant": decision["discriminant"],
                    "solution": decision["solution"],
                }
            )
        else:
            return {
                "center": center,
                "center_mask": supports[center],
                "center_free_mask": free_mask,
                "one_free_neighbours": list(neighbours),
                "pin_components": [
                    {
                        "vertex": vertex,
                        "entries": current_independent.pins_json(
                            components[vertex]
                        ),
                    }
                    for vertex in relevant
                ],
                "sign_variables": len(keys),
                "assignments_checked": 1 << len(keys),
                "failures": failures,
            }
    return None


def evaluate_graph(payload: tuple[dict, dict]) -> dict:
    graph, archived_record = payload
    index = int(graph["index"])
    adjacency = tuple(map(int, graph["adjacency"]))
    reference.validate_graph(adjacency)
    counts: Counter[str] = Counter()
    certificates = []
    seed_records = []

    for archived_seed in archived_record["seeds"]:
        seed_mask = int(archived_seed["seed_mask"])
        outside, defects, _ladj, _eligible = (
            current_independent.support_independent.independent_seed_instance(
                adjacency, seed_mask
            )
        )
        outside = tuple(map(int, outside))
        defects = tuple(map(int, defects))
        seed = list(current_independent.support_independent.bit_positions(seed_mask))
        if seed != archived_seed["seed"]:
            raise ValueError("independent seed reconstruction mismatch")
        cover_records = []
        for archived_cover in archived_seed["current_covers"]:
            zmask = int(archived_cover["zmask"])
            zvertices = tuple(
                current_independent.support_independent.bit_positions(zmask)
            )
            nvertices = tuple(
                vertex
                for vertex in range(len(outside))
                if not zmask & (1 << vertex)
            )
            graph_n = current_independent.support_independent.independent_induced_graph(
                adjacency, tuple(outside[vertex] for vertex in nvertices)
            )
            z_allowed = tuple(defects[vertex] for vertex in zvertices)
            n_allowed = tuple(defects[vertex] for vertex in nvertices)
            current_passing = 0
            star_passing = 0
            first_passing_witness = None
            cover_certificates = []
            for fixed in current_independent.support_independent.independent_labeled_support_families(
                z_allowed
            ):
                counts["labeled_support_families"] += 1
                masks, _deletions = (
                    current_independent.sparse_independent.propagated_masks(
                        n_allowed, fixed
                    )
                )
                failure = current_independent.sparse_independent.simple_propagation_failure(
                    graph_n, masks
                )
                if failure is not None:
                    counts[f"propagation_failure:{failure}"] += 1
                    continue
                sparse = current_independent.sparse_independent.independent_small_support_check(
                    graph_n, masks
                )
                if not sparse.feasible:
                    counts["sparse_value_failure"] += 1
                    continue
                if current_independent.base_independent.independent_double_pin(
                    graph_n, masks
                ) is not None:
                    counts["double_pin_failure"] += 1
                    continue
                if current_independent.full_pin_independent.independent_full_pin(
                    graph_n, masks
                ) is not None:
                    counts["full_pin_failure"] += 1
                    continue
                if current_independent.independent_singleton_locator(
                    graph_n, masks
                ) is not None:
                    counts["singleton_one_free_failure"] += 1
                    continue
                if current_independent.independent_correlated_locator(
                    graph_n, masks
                ) is not None:
                    counts["correlated_one_free_failure"] += 1
                    continue
                current_passing += 1
                counts["current_passing_families"] += 1
                certificate = independent_star_locator(graph_n, masks)
                if certificate is None:
                    star_passing += 1
                    counts["star_passing_families"] += 1
                    if first_passing_witness is None:
                        first_passing_witness = {
                            "z_supports": list(fixed),
                            "propagated_masks": list(masks),
                        }
                    continue
                counts["star_failed_families"] += 1
                witness = {
                    "seed": seed,
                    "seed_mask": seed_mask,
                    "zmask": zmask,
                    "z_supports": list(fixed),
                    "propagated_masks": list(masks),
                    "certificate": certificate,
                }
                certificates.append(witness)
                cover_certificates.append(witness)

            expected_current = int(
                archived_cover["family_counts"].get(
                    "combined_passing_families", 0
                )
            )
            if current_passing != expected_current:
                raise ValueError("independent current family replay mismatch")
            cover_records.append(
                {
                    "zmask": zmask,
                    "current_passing_families": current_passing,
                    "star_failed_families": len(cover_certificates),
                    "star_passing_families": star_passing,
                    "status": "PASSING" if star_passing else "INFEASIBLE",
                    "first_passing_witness": first_passing_witness,
                }
            )
        if not any(cover["current_passing_families"] for cover in cover_records):
            raise ValueError("input graph already failed at a current seed")
        seed_records.append(
            {
                "seed": seed,
                "seed_mask": seed_mask,
                "status": (
                    "PASSING"
                    if any(cover["status"] == "PASSING" for cover in cover_records)
                    else "INFEASIBLE"
                ),
                "current_covers": cover_records,
            }
        )
    rejecting = next(
        (record for record in seed_records if record["status"] == "INFEASIBLE"),
        None,
    )
    return {
        "index": index,
        "decision": "REJECTED" if rejecting is not None else "SURVIVOR",
        "first_rejecting_seed_mask": (
            0 if rejecting is None else int(rejecting["seed_mask"])
        ),
        "counts": dict(sorted(counts.items())),
        "certificates": certificates,
        "seeds": seed_records,
    }


def evaluate_graph_from_independent_current(payload: tuple) -> dict:
    """Rebuild all seed/cover/current-family quantifiers, then test stars."""

    graph, current_payload, archived_current_record = payload
    rebuilt_current = current_independent.verify_graph(current_payload)
    if rebuilt_current != archived_current_record:
        raise ValueError("independent full current-layer quantifier replay differs")
    return {
        "star_record": evaluate_graph((graph, rebuilt_current)),
        "current_layer_counts": rebuilt_current["counts"],
    }


def validate_upstream_and_report(report: dict) -> tuple[list[dict], list[int], dict]:
    paths = {
        MANIFEST.name: MANIFEST,
        MANIFEST_VERIFICATION.name: MANIFEST_VERIFICATION,
        CURRENT.name: CURRENT,
        CURRENT_VERIFICATION.name: CURRENT_VERIFICATION,
    }
    for name, expected in EXPECTED_UPSTREAM.items():
        if sha256(paths[name]) != expected:
            raise ValueError(f"upstream hash mismatch: {name}")
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    manifest_check = json.loads(
        MANIFEST_VERIFICATION.read_text(encoding="utf-8")
    )
    current = json.loads(CURRENT.read_text(encoding="utf-8"))
    current_check = json.loads(
        CURRENT_VERIFICATION.read_text(encoding="utf-8")
    )
    if (
        manifest.get("schema") != "d6-current-certified-residue-v6"
        or manifest.get("status") != "COMPLETE_MIXED_CERTIFICATE_UNION"
        or manifest_check.get("status") != "PASS"
        or manifest_check.get("manifest", {}).get("sha256")
        != EXPECTED_UPSTREAM[MANIFEST.name]
        or current.get("kind") != "d6_k7_one_free_edge_seed_conjunction"
        or current.get("status") != "COMPLETE"
        or current_check.get("status") != "PASS"
        or current_check.get("report", {}).get("sha256")
        != EXPECTED_UPSTREAM[CURRENT.name]
        or current_check.get("source_sha256")
        != EXPECTED_CURRENT_CHECKER_SOURCES
    ):
        raise ValueError("upstream semantics or verification differ")
    for name, expected in EXPECTED_CURRENT_CHECKER_SOURCES.items():
        if sha256(ROOT / name) != expected:
            raise ValueError(f"independent current-layer source differs: {name}")
    current_roots = current_check.get("root_sha256", {})
    if not isinstance(current_roots, dict) or not current_roots:
        raise ValueError("current-layer verification omits root hashes")
    for name, expected in current_roots.items():
        if Path(name).name != name or sha256(ROOT / name) != expected:
            raise ValueError(f"independent current-layer root differs: {name}")
    graphs = manifest["classes"]["K7"]["graphs"]
    indices = [int(graph["index"]) for graph in graphs]
    if len(indices) != 19 or stable_hash(indices) != EXPECTED_INPUT_SHA256:
        raise ValueError("frozen 19-list differs")
    expected_semantics = {
        "candidate_nonedges_optional": True,
        "allowed_unpinned_coordinates_may_be_zero": True,
        "only_required_edges_enter_star_equations": True,
        "propagated_masks_are_support_supersets": True,
        "normalization_factor_t_nonzero_on_cover_complement_N": True,
        "positive_sqrt7_embedding_checked_exactly": True,
        "floating_point_enters_rejection": False,
        "graph_rejected_if_any_k7_seed_is_infeasible": True,
        "cap500000_campaign_used": False,
    }
    if (
        report.get("schema") != 1
        or report.get("kind")
        != "d6_k7_one_free_neighbour_two_free_center_increment"
        or report.get("status") != "COMPLETE"
        or report.get("semantics") != expected_semantics
        or report.get("upstream_artifact_sha256")
        != dict(sorted(EXPECTED_UPSTREAM.items()))
        or report.get("input", {}).get("ordered_indices") != indices
        or report.get("input", {}).get("ordered_indices_sha256")
        != stable_hash(indices)
        or report.get("input", {}).get("embedded_graphs_sha256")
        != stable_hash(graphs)
    ):
        raise ValueError("report schema, semantics, roots, or input differs")
    return graphs, indices, current


def validate_source_boundary(report: dict) -> dict:
    sources = report.get("source_sha256")
    if not isinstance(sources, dict) or set(sources) != set(SOURCE_FILES):
        raise ValueError("source boundary differs")
    if (
        sources["build_d6_k7_one_two_star_increment.py"]
        != EXPECTED_BUILDER_SHA256
        or sources["d6_k7_one_two_star.py"] != EXPECTED_KERNEL_SHA256
    ):
        raise ValueError("frozen builder/kernel hash differs")
    for name in SOURCE_FILES:
        if sha256(ROOT / name) != sources[name]:
            raise ValueError(f"working source differs: {name}")
    git = report.get("execution", {}).get("git", {})
    if not git.get("available") or git.get("branch") != "codex/dimension6":
        raise ValueError("launch git provenance unavailable or wrong branch")
    commit = git.get("commit")
    for name in SOURCE_FILES:
        blob = subprocess.run(
            ["git", "show", f"{commit}:{name}"],
            cwd=ROOT,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        ).stdout
        if hashlib.sha256(blob).hexdigest() != sources[name]:
            raise ValueError(f"committed source differs: {name}")
    porcelain = "\n".join(git.get("porcelain_lines", ()))
    if (
        hashlib.sha256(porcelain.encode("utf-8")).hexdigest()
        != git.get("porcelain_sha256")
        or git.get("dirty") != bool(git.get("porcelain_lines"))
        or not all(line.startswith("?? ") for line in git.get("porcelain_lines", ()))
    ):
        raise ValueError("tracked source was dirty or porcelain binding differs")
    return {
        "commit": commit,
        "branch": git["branch"],
        "tracked_sources_clean_at_launch": True,
        "porcelain_sha256": git["porcelain_sha256"],
    }


def verify_report(
    report_path: Path,
    expected_report_sha256: str,
    *,
    workers: int,
    enforce_source_boundary: bool = True,
) -> dict:
    report_path = report_path.resolve()
    report_hash = sha256(report_path)
    if report_hash != expected_report_sha256:
        raise ValueError("report hash differs")
    report = json.loads(report_path.read_text(encoding="utf-8"))
    graphs, indices, current = validate_upstream_and_report(report)
    provenance = (
        validate_source_boundary(report)
        if enforce_source_boundary
        else {"status": "SKIPPED_FOR_TEST"}
    )
    current_by_index = {
        int(record["index"]): record for record in current["records"]
    }
    selected = set(indices)
    prior = current_independent.base_independent.read_failure_keys(
        current_independent.base_independent.PRIOR_CERTIFICATES,
        selected,
        "dual_failure_witnesses",
    )
    tetrad = current_independent.base_independent.read_failure_keys(
        current_independent.base_independent.TETRAD_CERTIFICATES,
        selected,
        "tetrad_failure_witnesses",
    )
    tetrad_rows = current_independent.base_independent.read_tetrad_rows(selected)
    payloads = [
        (
            graph,
            (
                graph,
                [
                    [list(seed), zmask]
                    for seed, zmask in sorted(
                        prior.get(int(graph["index"]), set())
                    )
                ],
                [
                    [list(seed), zmask]
                    for seed, zmask in sorted(
                        tetrad.get(int(graph["index"]), set())
                    )
                ],
                int(
                    tetrad_rows[int(graph["index"])][
                        "tetrad_passing_covers"
                    ]
                ),
            ),
            current_by_index[int(graph["index"])],
        )
        for graph in graphs
    ]
    started = time.monotonic()
    if workers == 1:
        rebuilt = list(map(evaluate_graph_from_independent_current, payloads))
    else:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            rebuilt = list(
                executor.map(
                    evaluate_graph_from_independent_current,
                    payloads,
                    chunksize=1,
                )
            )
    records = [item["star_record"] for item in rebuilt]
    if records != report.get("records"):
        raise ValueError("independent family/sign replay differs from report")

    current_totals: Counter[str] = Counter()
    for item in rebuilt:
        current_totals.update(item["current_layer_counts"])

    totals: Counter[str] = Counter()
    for record in records:
        totals.update(record["counts"])
    rejected = [
        record["index"] for record in records if record["decision"] == "REJECTED"
    ]
    rejected_set = set(rejected)
    survivors = [index for index in indices if index not in rejected_set]
    summary = report["summary"]
    if (
        totals["current_passing_families"] != 88
        or totals["star_failed_families"] != 4
        or rejected != EXPECTED_REJECTIONS
        or summary.get("totals") != dict(sorted(totals.items()))
        or summary.get("graphs_rejected") != len(rejected)
        or summary.get("rejected_indices") != rejected
        or summary.get("rejected_indices_sha256") != stable_hash(rejected)
        or summary.get("graphs_surviving") != len(survivors)
        or summary.get("ordered_survivor_indices") != survivors
        or summary.get("ordered_survivor_indices_sha256") != stable_hash(survivors)
    ):
        raise ValueError("independent summary differs")

    positive = tuple(map(int, lower_bound_18_graph()))
    reference.validate_graph(positive)
    positive_k7 = sum(1 for _ in reference.clique_masks(positive, 7))
    optional = line_system_decision(
        ((ONE, ZERO, ZERO),), ONE, q(-1)
    )
    expected_controls = {
        "known_realizable_18": {
            "vertices": 18,
            "K7_seeds": 0,
            "status": "PASS_NOT_APPLICABLE_NO_K7",
        },
        "optional_zero_line_hyperbola": {
            "identity": "(x-1)(y-1)=-1 with required line x=0",
            "witness": {"x": 0, "y": 2},
            "classification": "rank1_fixed_coordinate",
            "status": "PASS_FEASIBLE",
        },
    }
    if positive_k7 or not optional["feasible"] or report.get(
        "positive_controls"
    ) != expected_controls:
        raise ValueError("positive semantics control failed")

    return {
        "schema": 1,
        "kind": "d6_k7_one_two_star_increment_verification",
        "status": "PASS",
        "claim": (
            "Every current passing family and every correlated sign "
            "assignment was independently reconstructed without importing "
            "the producer or production star kernel."
        ),
        "report": {"path": str(report_path), "sha256": report_hash},
        "verifier_source_sha256": sha256(Path(__file__).resolve()),
        "upstream_artifact_sha256": dict(sorted(EXPECTED_UPSTREAM.items())),
        "provenance": provenance,
        "checked": {
            "graphs": len(indices),
            "K7_seeds": current_totals["seeds"],
            "eligible_covers": current_totals["eligible_covers"],
            "inherited_current_covers": current_totals["current_covers"],
            "labeled_support_families": current_totals["labeled_z_families"],
            "current_passing_families": totals["current_passing_families"],
            "star_failed_families": totals["star_failed_families"],
            "certificates": sum(len(record["certificates"]) for record in records),
            "rejected_graphs": len(rejected),
            "surviving_graphs": len(survivors),
            "positive_18_control": True,
            "optional_zero_control": True,
        },
        "rejected_indices": rejected,
        "rejected_indices_sha256": stable_hash(rejected),
        "ordered_survivor_indices": survivors,
        "ordered_survivor_indices_sha256": stable_hash(survivors),
        "semantics": report["semantics"],
        "execution": {
            "workers": workers,
            "started_utc": datetime.now(UTC).isoformat(),
            "elapsed_seconds": time.monotonic() - started,
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python_version": sys.version,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--report-sha256", required=True)
    parser.add_argument("--workers", type=int, default=9)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("workers must be positive")
    result = verify_report(
        args.report,
        args.report_sha256,
        workers=args.workers,
    )
    atomic_json(args.output.resolve(), result)
    print(
        json.dumps(
            {
                "output": str(args.output),
                "sha256": sha256(args.output.resolve()),
                "status": result["status"],
                "checked": result["checked"],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
