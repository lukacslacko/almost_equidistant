#!/usr/bin/env python3
"""Independent verifier for the exact K7 one-free-edge conjunction.

This checker imports neither the production conjunction builder nor either
production one-free-edge locator.  It independently replays every graph,
K7 seed, inherited-passing cover, and labeled support family in the frozen
24-graph residue.  The singleton and correlated certificates are rebuilt by
a separate exact ``Q(sqrt(7))`` implementation with an independently
transcribed coordinate-component/sign-correlation enumerator.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import sys
import time
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from pathlib import Path
from typing import Sequence

import d6_k7_rank_reference as reference
import verify_d6_k7_double_pin_conjunction as base_independent
import verify_d6_k7_full_pin_increment as full_pin_independent
import verify_d6_k7_sparse_value_full as sparse_independent
import verify_d6_k7_support_full as support_independent


ROOT = Path(__file__).resolve().parent
REPORT = ROOT / "d6_k7_one_free_conjunction_report.json"
BASE_REPORT = ROOT / "d6_k7_double_pin_conjunction_report.json"
BASE_VERIFICATION = ROOT / "d6_k7_double_pin_conjunction_verification.json"
FULL_PIN_REPORT = ROOT / "d6_k7_full_pin_increment_report.json"
FULL_PIN_VERIFICATION = ROOT / "d6_k7_full_pin_increment_verification.json"
PIPELINES = ("old", "singleton", "correlated", "combined")
Quadratic = tuple[Fraction, Fraction]
PinComponent = tuple[int, tuple[int, ...], tuple[int, int]]
ComponentKey = tuple[int, tuple[int, ...]]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_hash(value: object) -> str:
    return hashlib.sha256(
        json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
        ).encode("ascii")
    ).hexdigest()


def atomic_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def bits(mask: int):
    """Yield set-bit positions in increasing order."""

    while mask:
        bit = mask & -mask
        yield bit.bit_length() - 1
        mask ^= bit


def validate_branch(
    graph: Sequence[int], masks: Sequence[int]
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    adjacency = tuple(map(int, graph))
    supports = tuple(map(int, masks))
    if len(adjacency) != len(supports):
        raise ValueError("graph/mask size mismatch")
    full = (1 << len(adjacency)) - 1
    for vertex, row in enumerate(adjacency):
        if row < 0 or row & ~full or row & (1 << vertex):
            raise ValueError("invalid branch graph")
        for neighbour in bits(row):
            if not adjacency[neighbour] & (1 << vertex):
                raise ValueError("asymmetric branch graph")
    if any(mask <= 0 or mask >= 128 for mask in supports):
        raise ValueError("branch masks must be nonempty seven-bit masks")
    return adjacency, supports


def coordinate_graph(
    required: Sequence[int], masks: Sequence[int], coordinate: int
) -> tuple[int, ...]:
    """Required edges whose propagated-mask intersection is one bit."""

    bit = 1 << coordinate
    local = [0] * len(required)
    for high in range(len(required)):
        for low in range(high):
            if (
                required[high] & (1 << low)
                and masks[high] & masks[low] == bit
            ):
                local[high] |= 1 << low
                local[low] |= 1 << high
    return tuple(local)


def nonbipartite_components(
    graph: Sequence[int],
) -> tuple[tuple[tuple[int, ...], tuple[int, int]], ...]:
    """Deterministically expose every odd component and a color conflict."""

    unseen = set(range(len(graph)))
    output = []
    while unseen:
        root = min(unseen)
        colors = {root: 0}
        queue = [root]
        unseen.remove(root)
        conflict = None
        for vertex in queue:
            for neighbour in bits(graph[vertex]):
                if neighbour not in colors:
                    colors[neighbour] = colors[vertex] ^ 1
                    unseen.discard(neighbour)
                    queue.append(neighbour)
                elif (
                    colors[neighbour] == colors[vertex]
                    and conflict is None
                ):
                    conflict = (min(vertex, neighbour), max(vertex, neighbour))
        if conflict is not None:
            output.append((tuple(sorted(colors)), conflict))
    return tuple(output)


def pinning_data(
    required: Sequence[int], masks: Sequence[int]
) -> tuple[tuple[int, ...], dict[tuple[int, int], PinComponent]]:
    pinned = [0] * len(required)
    data: dict[tuple[int, int], PinComponent] = {}
    for coordinate in range(7):
        local = coordinate_graph(required, masks, coordinate)
        for component, conflict in nonbipartite_components(local):
            for vertex in component:
                pinned[vertex] |= 1 << coordinate
                data[(vertex, coordinate)] = (
                    coordinate,
                    component,
                    conflict,
                )
    return tuple(pinned), data


def sign_sums(mask_size: int) -> tuple[int, ...]:
    if not 1 <= mask_size <= 7:
        raise ValueError("mask size lies outside [1,7]")
    return tuple(range(1 - mask_size, mask_size, 2))


def one_free_value(mask_size: int, sign_sum: int) -> Quadratic:
    """Return the forced free component in exact ``Q(sqrt(7))`` form."""

    if sign_sum not in sign_sums(mask_size):
        raise ValueError("inadmissible pinned-sign sum")
    denominator = 2 * (7 - sign_sum * sign_sum)
    return (
        Fraction(
            sign_sum * (sign_sum * sign_sum - mask_size - 5),
            denominator,
        ),
        Fraction(
            9 - mask_size - sign_sum * sign_sum,
            denominator,
        ),
    )


def quadratic_product(first: Quadratic, second: Quadratic) -> Quadratic:
    first_rational, first_radical = first
    second_rational, second_radical = second
    return (
        first_rational * second_rational
        + 7 * first_radical * second_radical,
        first_rational * second_radical
        + first_radical * second_rational,
    )


def reciprocal_possible(first_size: int, second_size: int) -> bool:
    target = (Fraction(1), Fraction(0))
    for first_sum in sign_sums(first_size):
        first_value = one_free_value(first_size, first_sum)
        for second_sum in sign_sums(second_size):
            if quadratic_product(
                first_value, one_free_value(second_size, second_sum)
            ) == target:
                return True
    return False


def pins_json(entries: Sequence[PinComponent]) -> list[dict]:
    return [
        {
            "coordinate": coordinate,
            "component": list(component),
            "conflict_edge": list(conflict),
        }
        for coordinate, component, conflict in entries
    ]


def independent_singleton_locator(
    graph: Sequence[int], masks: Sequence[int]
) -> dict | None:
    required, supports = validate_branch(graph, masks)
    pinned, data = pinning_data(required, supports)
    for second in range(len(required)):
        for first in bits(required[second] & ((1 << second) - 1)):
            intersection = supports[first] & supports[second]
            if intersection.bit_count() != 1:
                continue
            shared = intersection.bit_length() - 1
            first_needed = supports[first] ^ intersection
            second_needed = supports[second] ^ intersection
            if pinned[first] & supports[first] != first_needed:
                continue
            if pinned[second] & supports[second] != second_needed:
                continue
            if reciprocal_possible(
                supports[first].bit_count(), supports[second].bit_count()
            ):
                continue
            first_pins = tuple(
                data[(first, coordinate)] for coordinate in bits(first_needed)
            )
            second_pins = tuple(
                data[(second, coordinate)] for coordinate in bits(second_needed)
            )
            certificate = {
                "first": first,
                "second": second,
                "shared_coordinate": shared,
                "first_mask": supports[first],
                "second_mask": supports[second],
                "first_pin_components": pins_json(first_pins),
                "second_pin_components": pins_json(second_pins),
            }
            verify_singleton_certificate(required, supports, certificate)
            return certificate
    return None


def pin_entries_from_json(entries: Sequence[dict]) -> tuple[PinComponent, ...]:
    return tuple(
        (
            int(entry["coordinate"]),
            tuple(map(int, entry["component"])),
            tuple(map(int, entry["conflict_edge"])),
        )
        for entry in entries
    )


def verify_pin_entries(
    required: Sequence[int],
    supports: Sequence[int],
    vertex: int,
    free_bit: int,
    entries: tuple[PinComponent, ...],
) -> None:
    expected = tuple(bits(supports[vertex] ^ free_bit))
    if tuple(entry[0] for entry in entries) != expected:
        raise ValueError("pin entries do not cover exactly the nonfree mask")
    for coordinate, component, conflict in entries:
        components = dict(
            nonbipartite_components(
                coordinate_graph(required, supports, coordinate)
            )
        )
        if vertex not in component or components.get(component) != conflict:
            raise ValueError("invalid odd-component pin certificate")


def verify_singleton_certificate(
    graph: Sequence[int], masks: Sequence[int], certificate: dict
) -> None:
    required, supports = validate_branch(graph, masks)
    first = int(certificate["first"])
    second = int(certificate["second"])
    if not 0 <= first < second < len(required):
        raise ValueError("invalid singleton endpoints")
    if not required[first] & (1 << second):
        raise ValueError("singleton pair is not a required edge")
    if (
        supports[first] != int(certificate["first_mask"])
        or supports[second] != int(certificate["second_mask"])
    ):
        raise ValueError("singleton certificate masks differ")
    intersection = supports[first] & supports[second]
    if (
        intersection.bit_count() != 1
        or intersection.bit_length() - 1
        != int(certificate["shared_coordinate"])
    ):
        raise ValueError("singleton intersection differs")
    if reciprocal_possible(
        supports[first].bit_count(), supports[second].bit_count()
    ):
        raise ValueError("singleton mask-size pair admits a reciprocal")
    verify_pin_entries(
        required,
        supports,
        first,
        intersection,
        pin_entries_from_json(certificate["first_pin_components"]),
    )
    verify_pin_entries(
        required,
        supports,
        second,
        intersection,
        pin_entries_from_json(certificate["second_pin_components"]),
    )


def component_key(entry: PinComponent) -> ComponentKey:
    return entry[0], entry[1]


def independent_compatibility(
    first_mask: int,
    second_mask: int,
    free_coordinate: int,
    first_pins: Sequence[PinComponent],
    second_pins: Sequence[PinComponent],
) -> tuple[bool, int, int]:
    """Enumerate correlated odd-component signs with exact edge arithmetic."""

    free_bit = 1 << free_coordinate
    if (
        first_mask <= 0
        or second_mask <= 0
        or first_mask >= 128
        or second_mask >= 128
        or not first_mask & free_bit
        or not second_mask & free_bit
    ):
        raise ValueError("invalid correlated masks/free coordinate")
    if tuple(entry[0] for entry in first_pins) != tuple(
        bits(first_mask ^ free_bit)
    ):
        raise ValueError("first correlated pins are incomplete")
    if tuple(entry[0] for entry in second_pins) != tuple(
        bits(second_mask ^ free_bit)
    ):
        raise ValueError("second correlated pins are incomplete")

    entries = tuple(first_pins) + tuple(second_pins)
    keys = tuple(sorted({component_key(entry) for entry in entries}))
    positions = {key: position for position, key in enumerate(keys)}
    first_variables = {
        entry[0]: positions[component_key(entry)] for entry in first_pins
    }
    second_variables = {
        entry[0]: positions[component_key(entry)] for entry in second_pins
    }
    shared_pinned = (first_mask & second_mask) ^ free_bit
    if not shared_pinned:
        raise ValueError("correlated layer has no pinned overlap")
    assignments = 1 << len(keys)
    for assignment in range(assignments):
        signs = tuple(
            1 if assignment & (1 << position) else -1
            for position in range(len(keys))
        )
        first_sum = sum(signs[position] for position in first_variables.values())
        second_sum = sum(
            signs[position] for position in second_variables.values()
        )
        shared_sum = sum(
            signs[first_variables[coordinate]]
            * signs[second_variables[coordinate]]
            for coordinate in bits(shared_pinned)
        )
        rational, radical = quadratic_product(
            one_free_value(first_mask.bit_count(), first_sum),
            one_free_value(second_mask.bit_count(), second_sum),
        )
        if rational + shared_sum == 1 and radical == 0:
            return True, len(keys), assignment + 1
    return False, len(keys), assignments


def independent_correlated_locator(
    graph: Sequence[int], masks: Sequence[int]
) -> dict | None:
    required, supports = validate_branch(graph, masks)
    pinned, data = pinning_data(required, supports)
    for second in range(len(required)):
        for first in bits(required[second] & ((1 << second) - 1)):
            first_free = supports[first] & ~pinned[first]
            second_free = supports[second] & ~pinned[second]
            if first_free != second_free or first_free.bit_count() != 1:
                continue
            if not (supports[first] & supports[second] & first_free):
                continue
            shared_pinned = (
                supports[first] & supports[second]
            ) ^ first_free
            if not shared_pinned:
                continue
            free_coordinate = first_free.bit_length() - 1
            first_pins = tuple(
                data[(first, coordinate)]
                for coordinate in bits(supports[first] ^ first_free)
            )
            second_pins = tuple(
                data[(second, coordinate)]
                for coordinate in bits(supports[second] ^ second_free)
            )
            feasible, variables, assignments = independent_compatibility(
                supports[first],
                supports[second],
                free_coordinate,
                first_pins,
                second_pins,
            )
            if feasible:
                continue
            certificate = {
                "first": first,
                "second": second,
                "free_coordinate": free_coordinate,
                "first_mask": supports[first],
                "second_mask": supports[second],
                "first_pin_components": pins_json(first_pins),
                "second_pin_components": pins_json(second_pins),
                "sign_variables": variables,
                "assignments_checked": assignments,
            }
            verify_correlated_certificate(required, supports, certificate)
            return certificate
    return None


def verify_correlated_certificate(
    graph: Sequence[int], masks: Sequence[int], certificate: dict
) -> None:
    required, supports = validate_branch(graph, masks)
    first = int(certificate["first"])
    second = int(certificate["second"])
    if not 0 <= first < second < len(required):
        raise ValueError("invalid correlated endpoints")
    if not required[first] & (1 << second):
        raise ValueError("correlated pair is not a required edge")
    if (
        supports[first] != int(certificate["first_mask"])
        or supports[second] != int(certificate["second_mask"])
    ):
        raise ValueError("correlated certificate masks differ")
    coordinate = int(certificate["free_coordinate"])
    free_bit = 1 << coordinate
    shared = supports[first] & supports[second]
    if not shared & free_bit or not (shared ^ free_bit):
        raise ValueError("bad correlated free coordinate/overlap")
    first_pins = pin_entries_from_json(certificate["first_pin_components"])
    second_pins = pin_entries_from_json(certificate["second_pin_components"])
    verify_pin_entries(required, supports, first, free_bit, first_pins)
    verify_pin_entries(required, supports, second, free_bit, second_pins)
    feasible, variables, assignments = independent_compatibility(
        supports[first],
        supports[second],
        coordinate,
        first_pins,
        second_pins,
    )
    if feasible:
        raise ValueError("correlated certificate has a compatible assignment")
    if (
        variables != int(certificate["sign_variables"])
        or assignments != int(certificate["assignments_checked"])
    ):
        raise ValueError("correlated enumeration counts differ")


def classify_cover(
    graph_n: tuple[int, ...],
    z_allowed: tuple[int, ...],
    n_allowed: tuple[int, ...],
) -> dict:
    counts: Counter = Counter()
    first_witness = {pipeline: None for pipeline in PIPELINES}
    singleton_certificates = []
    correlated_certificates = []

    for fixed in support_independent.independent_labeled_support_families(
        z_allowed
    ):
        counts["labeled_z_families"] += 1
        masks, _ = sparse_independent.propagated_masks(n_allowed, fixed)
        failure = sparse_independent.simple_propagation_failure(graph_n, masks)
        if failure is not None:
            counts[f"propagation_failure:{failure}"] += 1
            continue
        sparse = sparse_independent.independent_small_support_check(
            graph_n, masks
        )
        if not sparse.feasible:
            counts["sparse_value_failure"] += 1
            continue
        counts["pre_pinning_families"] += 1

        double_pin = base_independent.independent_double_pin(graph_n, masks)
        full_pin = None
        if double_pin is not None:
            counts["old_double_pin_failure"] += 1
        else:
            full_pin = full_pin_independent.independent_full_pin(
                graph_n, masks
            )
            if full_pin is not None:
                counts["old_full_pin_failure"] += 1
        if double_pin is not None or full_pin is not None:
            continue
        counts["pre_new_families"] += 1

        singleton = independent_singleton_locator(graph_n, masks)
        correlated = independent_correlated_locator(graph_n, masks)
        if singleton is not None:
            counts["singleton_failure"] += 1
            singleton_certificates.append(
                {
                    "z_supports": list(fixed),
                    "propagated_masks": list(masks),
                    "certificate": singleton,
                }
            )
        if correlated is not None:
            counts["correlated_failure"] += 1
            correlated_certificates.append(
                {
                    "z_supports": list(fixed),
                    "propagated_masks": list(masks),
                    "certificate": correlated,
                }
            )

        witness = {
            "z_supports": list(fixed),
            "propagated_masks": list(masks),
        }
        survives = {
            "old": True,
            "singleton": singleton is None,
            "correlated": correlated is None,
            "combined": singleton is None and correlated is None,
        }
        for pipeline, passing in survives.items():
            if passing:
                counts[f"{pipeline}_passing_families"] += 1
                if first_witness[pipeline] is None:
                    first_witness[pipeline] = witness
            else:
                counts[f"{pipeline}_infeasible_families"] += 1

    status_by_pipeline = {
        pipeline: (
            "PASSING" if first_witness[pipeline] is not None else "INFEASIBLE"
        )
        for pipeline in PIPELINES
    }
    return {
        "status_by_pipeline": status_by_pipeline,
        "family_counts": dict(sorted(counts.items())),
        "singleton_certificates": singleton_certificates,
        "correlated_certificates": correlated_certificates,
        "first_passing_witness_by_pipeline": first_witness,
    }


def verify_graph(payload: tuple) -> dict:
    graph, serialized_prior, serialized_tetrad, expected_current = payload
    prior = {(tuple(seed), int(zmask)) for seed, zmask in serialized_prior}
    tetrad = {(tuple(seed), int(zmask)) for seed, zmask in serialized_tetrad}
    index = int(graph["index"])
    adj = tuple(map(int, graph["adjacency"]))
    reference.validate_graph(adj)
    support_solver = reference.SupportSolver()
    zero_forcing = reference.ZeroForcingSolver()
    clique_solver = reference.CliqueStructureSolver()
    totals: Counter = Counter()
    seed_records = []
    observed_current_keys = set()

    for seed_mask in reference.clique_masks(adj, 7):
        totals["seeds"] += 1
        outside, defects, ladj, eligible = (
            support_independent.independent_seed_instance(adj, seed_mask)
        )
        outside = tuple(map(int, outside))
        defects = tuple(map(int, defects))
        seed = tuple(support_independent.bit_positions(seed_mask))
        covers = base_independent.independently_all_eligible_covers(
            ladj, eligible
        )
        if covers != tuple(reference.eligible_covers(ladj, eligible)):
            raise ValueError("independent eligible-cover enumeration differs")
        total_term_rank = reference.matching_size(defects)
        cover_records = []
        for zmask in covers:
            totals["eligible_covers"] += 1
            baseline = reference.analyze_cover(
                adj,
                outside,
                defects,
                zmask,
                support_solver,
                zero_forcing,
                total_term_rank,
                clique_solver,
            )
            inherited_status = base_independent.inherited_cover_status(
                adj,
                outside,
                defects,
                seed,
                zmask,
                baseline,
                prior,
                tetrad,
            )
            totals[f"inherited_cover_{inherited_status}"] += 1
            if inherited_status != "passing":
                continue
            observed_current_keys.add((seed, zmask))
            totals["current_covers"] += 1
            zvertices = tuple(support_independent.bit_positions(zmask))
            nlocal = tuple(
                position
                for position in range(len(outside))
                if not zmask & (1 << position)
            )
            graph_n = support_independent.independent_induced_graph(
                adj, tuple(outside[position] for position in nlocal)
            )
            classified = classify_cover(
                graph_n,
                tuple(defects[position] for position in zvertices),
                tuple(defects[position] for position in nlocal),
            )
            for pipeline, status in classified["status_by_pipeline"].items():
                totals[f"{pipeline}_{status.lower()}_covers"] += 1
            totals.update(classified["family_counts"])
            cover_records.append(
                {
                    "seed": list(seed),
                    "zmask": zmask,
                    **classified,
                }
            )

        if not cover_records:
            raise ValueError("frozen survivor has no inherited-passing cover")
        status_by_pipeline = {
            pipeline: (
                "PASSING"
                if any(
                    cover["status_by_pipeline"][pipeline] == "PASSING"
                    for cover in cover_records
                )
                else "INFEASIBLE"
            )
            for pipeline in PIPELINES
        }
        seed_records.append(
            {
                "seed": list(seed),
                "seed_mask": seed_mask,
                "status_by_pipeline": status_by_pipeline,
                "current_covers": cover_records,
            }
        )

    if len(observed_current_keys) != int(expected_current):
        raise ValueError("independent inherited current-cover count mismatch")
    decision_by_pipeline = {}
    first_rejecting_seed_mask_by_pipeline = {}
    for pipeline in PIPELINES:
        rejecting = next(
            (
                seed_record
                for seed_record in seed_records
                if seed_record["status_by_pipeline"][pipeline] == "INFEASIBLE"
            ),
            None,
        )
        decision_by_pipeline[pipeline] = (
            "REJECTED" if rejecting is not None else "SURVIVOR"
        )
        first_rejecting_seed_mask_by_pipeline[pipeline] = (
            int(rejecting["seed_mask"]) if rejecting is not None else 0
        )
    return {
        "index": index,
        "decision_by_pipeline": decision_by_pipeline,
        "first_rejecting_seed_mask_by_pipeline": (
            first_rejecting_seed_mask_by_pipeline
        ),
        "counts": dict(sorted(totals.items())),
        "seeds": seed_records,
    }


def load_and_validate_artifacts(
    report: dict, artifact_hashes: dict[str, str]
) -> tuple[dict, dict, dict, dict, list[int]]:
    paths = {
        BASE_REPORT.name: BASE_REPORT,
        BASE_VERIFICATION.name: BASE_VERIFICATION,
        FULL_PIN_REPORT.name: FULL_PIN_REPORT,
        FULL_PIN_VERIFICATION.name: FULL_PIN_VERIFICATION,
    }
    for name, path in paths.items():
        if sha256(path) != artifact_hashes[name]:
            raise ValueError(f"upstream artifact hash mismatch: {name}")
    for name, expected in base_independent.EXPECTED_HASHES.items():
        if sha256(ROOT / name) != expected:
            raise ValueError(f"frozen root hash mismatch: {name}")

    base_report = json.loads(BASE_REPORT.read_text(encoding="utf-8"))
    base_verification = json.loads(
        BASE_VERIFICATION.read_text(encoding="utf-8")
    )
    full_report = json.loads(FULL_PIN_REPORT.read_text(encoding="utf-8"))
    full_verification = json.loads(
        FULL_PIN_VERIFICATION.read_text(encoding="utf-8")
    )
    expected_base = {
        BASE_REPORT.name: artifact_hashes[BASE_REPORT.name],
        BASE_VERIFICATION.name: artifact_hashes[BASE_VERIFICATION.name],
    }
    if (
        base_report.get("schema") != 1
        or base_report.get("kind")
        != "d6_k7_two_defect_double_pin_seed_conjunction"
        or base_report.get("status") != "COMPLETE"
        or base_report.get("upstream_sha256")
        != dict(sorted(base_independent.EXPECTED_HASHES.items()))
        or base_verification.get("schema") != 1
        or base_verification.get("kind")
        != "d6_k7_two_defect_double_pin_seed_conjunction_verification"
        or base_verification.get("status") != "PASS"
        or base_verification.get("report", {}).get("sha256")
        != artifact_hashes[BASE_REPORT.name]
        or full_report.get("schema") != 1
        or full_report.get("kind") != "d6_k7_full_pin_odd_cycle_increment"
        or full_report.get("status") != "COMPLETE"
        or full_report.get("base_sha256") != dict(sorted(expected_base.items()))
        or full_verification.get("schema") != 1
        or full_verification.get("kind")
        != "d6_k7_full_pin_odd_cycle_increment_verification"
        or full_verification.get("status") != "PASS"
        or full_verification.get("report", {}).get("sha256")
        != artifact_hashes[FULL_PIN_REPORT.name]
        or full_verification.get("base_sha256")
        != dict(sorted(expected_base.items()))
    ):
        raise ValueError("upstream semantic/hash-binding gate failed")

    selected = list(map(int, full_report.get("exact_survivors", ())))
    if (
        len(selected) != 24
        or stable_hash(selected) != full_report.get("exact_survivors_sha256")
        or full_verification.get("exact_survivors_sha256")
        != stable_hash(selected)
    ):
        raise ValueError("upstream full-pin residue is not a bound 24-list")
    expected_artifacts = dict(sorted(artifact_hashes.items()))
    if (
        report.get("schema") != 1
        or report.get("kind") != "d6_k7_one_free_edge_seed_conjunction"
        or report.get("status") != "COMPLETE"
        or report.get("upstream_artifact_sha256") != expected_artifacts
        or report.get("input_indices") != selected
        or report.get("input_indices_sha256") != stable_hash(selected)
    ):
        raise ValueError("one-free report schema, roots, or input binding differ")
    expected_semantics = {
        "candidate_nonedges_optional": True,
        "only_required_edges_enter_rejection": True,
        "propagated_masks_are_support_supersets": True,
        "floating_point_enters_rejection": False,
        "graph_rejected_if_any_k7_seed_is_infeasible": True,
    }
    if report.get("semantics") != expected_semantics:
        raise ValueError("one-free report semantics differ")
    for name, expected in report.get("source_sha256", {}).items():
        if Path(name).name != name or sha256(ROOT / name) != expected:
            raise ValueError(f"reported production source changed: {name}")
    if not report.get("source_sha256"):
        raise ValueError("one-free report has no production source bindings")
    return (
        base_report,
        base_verification,
        full_report,
        full_verification,
        selected,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, default=REPORT)
    parser.add_argument("--report-sha256", required=True)
    parser.add_argument("--base-report-sha256", required=True)
    parser.add_argument("--base-verification-sha256", required=True)
    parser.add_argument("--full-pin-report-sha256", required=True)
    parser.add_argument("--full-pin-verification-sha256", required=True)
    parser.add_argument(
        "--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1)
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "d6_k7_one_free_conjunction_verification.json",
    )
    args = parser.parse_args()
    if args.workers <= 0:
        parser.error("workers must be positive")
    report_hash = sha256(args.report)
    if report_hash != args.report_sha256:
        raise ValueError("report differs from explicit SHA-256 pin")
    report = json.loads(args.report.read_text(encoding="utf-8"))
    artifact_hashes = {
        BASE_REPORT.name: args.base_report_sha256,
        BASE_VERIFICATION.name: args.base_verification_sha256,
        FULL_PIN_REPORT.name: args.full_pin_report_sha256,
        FULL_PIN_VERIFICATION.name: args.full_pin_verification_sha256,
    }
    _, _, full_report, _, selected_indices = load_and_validate_artifacts(
        report, artifact_hashes
    )

    residue_verification = json.loads(
        base_independent.RESIDUE_VERIFICATION.read_text(encoding="utf-8")
    )
    if residue_verification.get("status") != "PASS":
        raise ValueError("v5 residue verification is not PASS")
    residue = json.loads(
        base_independent.RESIDUE.read_text(encoding="utf-8")
    )
    all_graphs = residue["classes"]["K7"]["graphs"]
    all_indices = [int(graph["index"]) for graph in all_graphs]
    if (
        len(all_graphs) != 155
        or all_indices != residue["classes"]["K7"]["indices"]
    ):
        raise ValueError("v5 K7 corpus population/order mismatch")
    graph_by_index = {int(graph["index"]): graph for graph in all_graphs}
    if any(index not in graph_by_index for index in selected_indices):
        raise ValueError("upstream survivor is absent from the K7 corpus")
    graphs = [graph_by_index[index] for index in selected_indices]
    selected = set(selected_indices)
    prior = base_independent.read_failure_keys(
        base_independent.PRIOR_CERTIFICATES,
        selected,
        "dual_failure_witnesses",
    )
    tetrad = base_independent.read_failure_keys(
        base_independent.TETRAD_CERTIFICATES,
        selected,
        "tetrad_failure_witnesses",
    )
    tetrad_rows = base_independent.read_tetrad_rows(selected)
    payloads = [
        (
            graph,
            [
                [list(seed), zmask]
                for seed, zmask in sorted(prior.get(index, set()))
            ],
            [
                [list(seed), zmask]
                for seed, zmask in sorted(tetrad.get(index, set()))
            ],
            int(tetrad_rows[index]["tetrad_passing_covers"]),
        )
        for graph, index in zip(graphs, selected_indices, strict=True)
    ]

    started = time.monotonic()
    if args.workers == 1:
        records = list(map(verify_graph, payloads))
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            records = list(executor.map(verify_graph, payloads, chunksize=1))
    if [record["index"] for record in records] != selected_indices:
        raise ValueError("independent worker output order changed")
    if records != report.get("records"):
        raise ValueError("independent full quantifier records differ from report")

    rejected_by_pipeline = {
        pipeline: [
            record["index"]
            for record in records
            if record["decision_by_pipeline"][pipeline] == "REJECTED"
        ]
        for pipeline in PIPELINES
    }
    survivors_by_pipeline = {
        pipeline: [
            record["index"]
            for record in records
            if record["decision_by_pipeline"][pipeline] == "SURVIVOR"
        ]
        for pipeline in PIPELINES
    }
    if (
        rejected_by_pipeline["old"]
        or survivors_by_pipeline["old"] != selected_indices
    ):
        raise ValueError("independent old-pinning replay differs from residue")
    singleton = set(rejected_by_pipeline["singleton"])
    correlated = set(rejected_by_pipeline["correlated"])
    combined = set(rejected_by_pipeline["combined"])
    if not singleton <= combined or not correlated <= combined:
        raise ValueError("combined pipeline lost an individual rejection")
    totals: Counter = Counter()
    for record in records:
        totals.update(record["counts"])
    expected_summary = {
        "input_graphs": len(selected_indices),
        "rejected_by_pipeline": rejected_by_pipeline,
        "rejected_by_pipeline_sha256": {
            pipeline: stable_hash(indices)
            for pipeline, indices in rejected_by_pipeline.items()
        },
        "survivors_by_pipeline": survivors_by_pipeline,
        "survivors_by_pipeline_sha256": {
            pipeline: stable_hash(indices)
            for pipeline, indices in survivors_by_pipeline.items()
        },
        "singleton_correlated_graph_overlap": sorted(singleton & correlated),
        "combined_only_synergy_rejections": sorted(
            combined - singleton - correlated
        ),
        "totals": dict(sorted(totals.items())),
    }
    if report.get("summary") != expected_summary:
        raise ValueError("independent aggregate summary differs from report")

    output = {
        "schema": 1,
        "kind": "d6_k7_one_free_edge_seed_conjunction_verification",
        "status": "PASS",
        "report": {"path": str(args.report), "sha256": report_hash},
        "upstream_artifact_sha256": dict(sorted(artifact_hashes.items())),
        "checked": {
            "graphs": len(records),
            "K7_seeds": totals["seeds"],
            "eligible_covers": totals["eligible_covers"],
            "current_covers": totals["current_covers"],
            "labeled_z_families": totals["labeled_z_families"],
            "pre_pinning_families": totals["pre_pinning_families"],
            "pre_new_families": totals["pre_new_families"],
            "singleton_certificates": totals["singleton_failure"],
            "correlated_certificates": totals["correlated_failure"],
            "combined_rejected_graphs": len(
                rejected_by_pipeline["combined"]
            ),
            "combined_survivors": len(survivors_by_pipeline["combined"]),
        },
        "rejected_by_pipeline_sha256": {
            pipeline: stable_hash(indices)
            for pipeline, indices in rejected_by_pipeline.items()
        },
        "survivors_by_pipeline_sha256": {
            pipeline: stable_hash(indices)
            for pipeline, indices in survivors_by_pipeline.items()
        },
        "exact_survivors_sha256": stable_hash(
            survivors_by_pipeline["combined"]
        ),
        "input_indices_sha256": stable_hash(selected_indices),
        "root_sha256": dict(sorted(base_independent.EXPECTED_HASHES.items())),
        "source_sha256": {
            Path(__file__).name: sha256(Path(__file__)),
            "verify_d6_k7_double_pin_conjunction.py": sha256(
                ROOT / "verify_d6_k7_double_pin_conjunction.py"
            ),
            "verify_d6_k7_full_pin_increment.py": sha256(
                ROOT / "verify_d6_k7_full_pin_increment.py"
            ),
            "verify_d6_k7_sparse_value_full.py": sha256(
                ROOT / "verify_d6_k7_sparse_value_full.py"
            ),
            "verify_d6_k7_support_full.py": sha256(
                ROOT / "verify_d6_k7_support_full.py"
            ),
        },
        "execution": {
            "workers": args.workers,
            "logical_cpus": os.cpu_count(),
            "elapsed_seconds": time.monotonic() - started,
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python_version": sys.version,
        },
        "claim": (
            "A source-independent exact replay reconstructed every frozen "
            "24-residue graph, seed, cover, support family, singleton "
            "certificate, correlated sign-component certificate, and all "
            "four pipeline partitions."
        ),
    }
    atomic_json(args.output, output)
    print(
        json.dumps(
            {
                "output": str(args.output),
                "sha256": sha256(args.output),
                "status": "PASS",
                "checked": output["checked"],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
