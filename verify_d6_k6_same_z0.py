#!/usr/bin/env python3
"""Independent exact verifier for the dimension-six K6 same-Z0 layer.

This checker does not import the production evaluator, its block-support
prototype, or the earlier joint-support implementation.  It independently
reconstructs graphs, K6 seeds, defect masks, Lorentz components, every Z0,
the block-Hall inequalities, light-ray colorings, and the actual-support CSP.
It also changes the two central algorithms: inertia comes from exact SymPy
characteristic polynomials with Sturm root counts, and matching is checked by
direct Hall-subset enumeration rather than augmenting paths.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import sys
import time
from collections import Counter
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Iterator, Sequence

import sympy as sp

from test_d6_k6_lorentz import lower_bound_18_graph


ROOT = Path(__file__).resolve().parent
COORDINATES = 6
VARIABLE = sp.symbols("lambda")
REPORT_SCHEMA = "d6-k6-same-z0-v1"
VERIFICATION_SCHEMA = "d6-k6-same-z0-verification-v1"
EXPECTED_INPUT = 990
EXPECTED_INDICES_SHA256 = (
    "16cb562473ebb8c496018891362642108113a16f007929d413770ccdb665acb9"
)
EXPECTED = {
    "d6_current_residue_manifest.json": (
        "d06cdb23b03b1f8523fbe3c3ede5bc297ace21a844282830c5cf51a6cee3512d"
    ),
    "d6_k6_bipartite_rank_input.json": (
        "dfee979e6ca2759e80cf803c3054c8954d98543875d467644b91ee434c22e845"
    ),
    "d6_k6_same_z0.py": (
        "6d6dd948a245d5d8245f6ee10ce58c725aba29b9d3a9b325c8be11b59d519b5f"
    ),
    "d6_k6_same_z0_report.json": (
        "cea11ab95f24aa5fe79f06cc4d1a719b775ac5647551e0fb737b85835c443ca1"
    ),
    "test_d6_k6_lorentz.py": (
        "f4fc8544a4311f742267808e6c785f82139accb9da99b496d206283b52831a47"
    ),
}


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


def atomic_json(path: Path, value: object) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def bits(mask: int) -> list[int]:
    answer = []
    while mask:
        bit = mask & -mask
        mask ^= bit
        answer.append(bit.bit_length() - 1)
    return answer


def validate_graph(adj: Sequence[int]) -> None:
    size = len(adj)
    full = (1 << size) - 1
    for vertex, row in enumerate(adj):
        if row & ~full or row & (1 << vertex):
            raise ValueError("independent verifier received invalid adjacency")
        for other in range(vertex):
            if bool(row & (1 << other)) != bool(adj[other] & (1 << vertex)):
                raise ValueError("independent verifier received asymmetric graph")
    for first in range(size):
        later = full & ~adj[first] & ~((1 << (first + 1)) - 1)
        remaining = later
        while remaining:
            bit = remaining & -remaining
            remaining ^= bit
            second = bit.bit_length() - 1
            if remaining & ~adj[second]:
                raise ValueError("independent triple violates K6 rule hypothesis")


def clique_masks(adj: Sequence[int], size: int) -> Iterator[int]:
    def visit(candidates: int, need: int, selected: int) -> Iterator[int]:
        if need == 0:
            yield selected
            return
        while candidates.bit_count() >= need:
            bit = candidates & -candidates
            candidates ^= bit
            vertex = bit.bit_length() - 1
            yield from visit(candidates & adj[vertex], need - 1, selected | bit)

    yield from visit((1 << len(adj)) - 1, size, 0)


@dataclass(frozen=True)
class Instance:
    seed: tuple[int, ...]
    outside: tuple[int, ...]
    defects: tuple[int, ...]
    l_adj: tuple[int, ...]
    eligible_z0: int


def build_instance(adj: Sequence[int], seed_mask: int) -> Instance:
    seed = tuple(bits(seed_mask))
    if len(seed) != COORDINATES:
        raise ValueError("independent verifier seed is not a K6")
    outside = tuple(
        vertex for vertex in range(len(adj)) if not seed_mask & (1 << vertex)
    )
    defects = []
    for vertex in outside:
        allowed = 0
        for coordinate, seed_vertex in enumerate(seed):
            if not adj[vertex] & (1 << seed_vertex):
                allowed |= 1 << coordinate
        defects.append(allowed)
    l_adj = [0] * len(outside)
    for left, right in combinations(range(len(outside)), 2):
        if (
            adj[outside[left]] & (1 << outside[right])
            and not defects[left] & defects[right]
        ):
            l_adj[left] |= 1 << right
            l_adj[right] |= 1 << left
    eligible = sum(
        1 << local
        for local, allowed in enumerate(defects)
        if allowed.bit_count() >= 3
    )
    return Instance(seed, outside, tuple(defects), tuple(l_adj), eligible)


def z0_subsets(eligible: int) -> Iterator[int]:
    candidates = bits(eligible)
    for size in range(min(COORDINATES, len(candidates)) + 1):
        for selected in combinations(candidates, size):
            yield sum(1 << local for local in selected)


def hall_matchable(selected: int, masks: Sequence[int]) -> bool:
    """Check structural full row rank by Hall's theorem, subset by subset."""

    chosen = bits(selected)
    if len(chosen) > COORDINATES:
        return False
    for encoded in range(1, 1 << len(chosen)):
        union = 0
        count = 0
        for position, local in enumerate(chosen):
            if encoded & (1 << position):
                union |= masks[local]
                count += 1
        if union.bit_count() < count:
            return False
    return True


@dataclass(frozen=True)
class Component:
    vertices: int
    bipartite: bool
    side_a: int
    side_b: int


def components(instance: Instance, deleted: int) -> tuple[Component, ...]:
    remaining = ((1 << len(instance.outside)) - 1) & ~deleted
    answer = []
    while remaining:
        root = remaining & -remaining
        reached = root
        frontier = root
        side_b = 0
        bipartite = True
        while frontier:
            bit = frontier & -frontier
            frontier ^= bit
            local = bit.bit_length() - 1
            neighbors = instance.l_adj[local] & remaining
            if side_b & bit:
                if neighbors & side_b:
                    bipartite = False
            elif neighbors & (reached & ~side_b):
                bipartite = False
            new = neighbors & ~reached
            if not side_b & bit:
                side_b |= new
            reached |= new
            frontier |= new
        side_b &= reached
        answer.append(Component(reached, bipartite, reached & ~side_b, side_b))
        remaining &= ~reached
    return tuple(answer)


def side_pattern(
    adj: Sequence[int], instance: Instance, selected: int
) -> tuple[int, ...]:
    absolute = tuple(instance.outside[local] for local in bits(selected))
    return tuple(
        (1 << row)
        | sum(
            1 << column
            for column, other in enumerate(absolute)
            if adj[vertex] & (1 << other)
        )
        for row, vertex in enumerate(absolute)
    )


class SympyInertiaCache:
    """Independent exact inertia via characteristic polynomials and Sturm."""

    def __init__(self) -> None:
        self.values: dict[tuple[int, ...], tuple[int, int, int]] = {}

    def solve(self, rows: tuple[int, ...]) -> tuple[int, int, int]:
        if rows in self.values:
            return self.values[rows]
        size = len(rows)
        matrix = sp.Matrix(
            size,
            size,
            lambda row, column: (rows[row] >> column) & 1,
        )
        polynomial = matrix.charpoly(VARIABLE).as_poly()
        positive = 0
        negative = 0
        zero = 0
        _, factors = polynomial.sqf_list()
        for factor, multiplicity in factors:
            current = factor
            if current.eval(0) == 0:
                zero += multiplicity
                current = current.exquo(sp.Poly(VARIABLE, VARIABLE))
            if current.degree() > 0:
                negative += multiplicity * int(current.count_roots(-sp.oo, 0))
                positive += multiplicity * int(current.count_roots(0, sp.oo))
        if positive + negative + zero != size:
            raise AssertionError("independent exact inertia count does not sum")
        answer = positive, negative, zero
        self.values[rows] = answer
        return answer


def allowed_union(instance: Instance, selected: int) -> int:
    answer = 0
    for local in bits(selected):
        answer |= instance.defects[local]
    return answer


def block_subsets_pass(blocks: Sequence[tuple[int, int]]) -> bool:
    for selected in range(1, 1 << len(blocks)):
        rank = 0
        allowed = 0
        for position, (rank_lower, coordinate_mask) in enumerate(blocks):
            if selected & (1 << position):
                rank += rank_lower
                allowed |= coordinate_mask
        if rank > allowed.bit_count():
            return False
    return True


def component_block_passes(
    adj: Sequence[int],
    instance: Instance,
    z0: int,
    component: Component,
    inertia: SympyInertiaCache,
) -> bool:
    pattern_a = side_pattern(adj, instance, component.side_a)
    pattern_b = side_pattern(adj, instance, component.side_b)
    positive_a, negative_a, _ = inertia.solve(pattern_a)
    positive_b, negative_b, _ = inertia.solve(pattern_b)
    size_a = component.side_a.bit_count()
    size_b = component.side_b.bit_count()
    orientations = (
        (size_a - negative_a, size_b - positive_b),
        (size_a - positive_a, size_b - negative_b),
    )
    for rank_a, rank_b in orientations:
        blocks = [
            (rank_a, allowed_union(instance, component.side_a)),
            (rank_b, allowed_union(instance, component.side_b)),
        ]
        blocks.extend((1, instance.defects[local]) for local in bits(z0))
        if block_subsets_pass(blocks):
            return True

    light = z0 | component.vertices
    return (
        light.bit_count() <= COORDINATES
        and hall_matchable(light, instance.defects)
    )


def actual_support_domains(allowed: int, zero_factor: bool) -> tuple[int, ...]:
    minimum = 3 if zero_factor else 1
    choices = []
    current = allowed
    while current:
        if current.bit_count() >= minimum:
            choices.append(current)
        current = (current - 1) & allowed
    # Deliberately use a different tie order from production.
    choices.sort(key=lambda mask: (mask.bit_count(), -mask))
    return tuple(choices)


def assigned_hall_passes(
    selected: int, assignments: dict[int, int], outside_size: int
) -> bool:
    assigned = selected & sum(1 << local for local in assignments)
    masks = [0] * outside_size
    for local, support in assignments.items():
        masks[local] = support
    return hall_matchable(assigned, masks)


def joint_actual_support_exists(
    instance: Instance, z0: int, bins: tuple[int, int]
) -> bool:
    if bins[0] & bins[1] != z0:
        raise ValueError("independent support bins intersect outside Z0")
    involved = bins[0] | bins[1]
    domains = {
        local: actual_support_domains(
            instance.defects[local], bool(z0 & (1 << local))
        )
        for local in bits(involved)
    }
    peers = {}
    for local in bits(involved):
        peer_mask = 0
        for bin_mask in bins:
            if bin_mask & (1 << local):
                peer_mask |= bin_mask & ~(1 << local)
        peers[local] = peer_mask

    assignments: dict[int, int] = {}

    def compatible_choices(local: int) -> list[int]:
        return [
            support
            for support in domains[local]
            if all(
                (support & other_support).bit_count() != 1
                for other, other_support in assignments.items()
                if peers[local] & (1 << other)
            )
        ]

    def visit(unassigned: int) -> bool:
        if not unassigned:
            return True
        best_local = -1
        best_choices = None
        for local in bits(unassigned):
            choices = compatible_choices(local)
            key = (len(choices), local)
            if best_choices is None or key < (len(best_choices), best_local):
                best_local = local
                best_choices = choices
        assert best_choices is not None
        if not best_choices:
            return False
        bit = 1 << best_local
        for support in best_choices:
            assignments[best_local] = support
            if all(
                assigned_hall_passes(bin_mask, assignments, len(instance.outside))
                for bin_mask in bins
            ) and visit(unassigned ^ bit):
                return True
            assignments.pop(best_local, None)
        return False

    return visit(involved)


def nonbipartite_system_passes(
    instance: Instance, z0: int, odd: Sequence[int]
) -> bool:
    if not odd:
        bins = (z0, z0)
        return (
            hall_matchable(z0, instance.defects)
            and joint_actual_support_exists(instance, z0, bins)
        )
    for encoded in range(1 << (len(odd) - 1)):
        bins = [z0 | odd[0], z0]
        for position, component in enumerate(odd[1:]):
            color = (encoded >> position) & 1
            bins[color] |= component
        pair = (bins[0], bins[1])
        if not all(hall_matchable(bin_mask, instance.defects) for bin_mask in pair):
            continue
        if joint_actual_support_exists(instance, z0, pair):
            return True
    return False


def independent_cross_classification(
    adj: Sequence[int], instance: Instance
) -> tuple[tuple[tuple[int, ...], bool, bool, str], ...]:
    """Evaluate both existential subsystems independently for every Z0."""

    inertia = SympyInertiaCache()
    rows = []
    for z0 in z0_subsets(instance.eligible_z0):
        absolute_z0 = tuple(instance.outside[local] for local in bits(z0))
        if not hall_matchable(z0, instance.defects):
            block_passed = False
            nonbipartite_passed = False
        else:
            local_components = components(instance, z0)
            block_passed = all(
                component_block_passes(adj, instance, z0, component, inertia)
                for component in local_components
                if component.bipartite
            )
            odd = [
                component.vertices
                for component in local_components
                if not component.bipartite
            ]
            nonbipartite_passed = nonbipartite_system_passes(
                instance, z0, odd
            )
        if block_passed and nonbipartite_passed:
            category = "both"
        elif block_passed:
            category = "block_only"
        elif nonbipartite_passed:
            category = "nonbipartite_only"
        else:
            category = "neither"
        rows.append(
            (absolute_z0, block_passed, nonbipartite_passed, category)
        )
    return tuple(rows)


CORE_FIELDS = (
    "rejected",
    "seeds_checked",
    "impossible_seeds",
    "z0_considered",
    "z0_matchable",
    "z0_block_passed",
    "z0_block_failed",
    "z0_nonbipartite_failed",
)


@dataclass
class SeedResult:
    feasible: bool
    z0_considered: int
    z0_matchable: int
    z0_block_passed: int
    z0_block_failed: int
    z0_nonbipartite_failed: int
    reasons: tuple[tuple[tuple[int, ...], str], ...]


def solve_seed(
    adj: Sequence[int], instance: Instance, inertia: SympyInertiaCache
) -> SeedResult:
    considered = 0
    matchable = 0
    block_passed = 0
    block_failed = 0
    nonbip_failed = 0
    reasons = []
    for z0 in z0_subsets(instance.eligible_z0):
        considered += 1
        absolute_z0 = tuple(instance.outside[local] for local in bits(z0))
        if not hall_matchable(z0, instance.defects):
            reasons.append((absolute_z0, "Z0_allowed_mask_matching"))
            continue
        matchable += 1
        local_components = components(instance, z0)
        all_blocks_pass = True
        for component in local_components:
            if component.bipartite and not component_block_passes(
                adj, instance, z0, component, inertia
            ):
                all_blocks_pass = False
        if not all_blocks_pass:
            block_failed += 1
            reasons.append((absolute_z0, "bipartite_block_support"))
            continue
        block_passed += 1
        odd = [
            component.vertices
            for component in local_components
            if not component.bipartite
        ]
        if nonbipartite_system_passes(instance, z0, odd):
            return SeedResult(
                True,
                considered,
                matchable,
                block_passed,
                block_failed,
                nonbip_failed,
                tuple(reasons),
            )
        nonbip_failed += 1
        reasons.append((absolute_z0, "nonbipartite_joint_actual_support"))
    return SeedResult(
        False,
        considered,
        matchable,
        block_passed,
        block_failed,
        nonbip_failed,
        tuple(reasons),
    )


def evaluate_graph(adj: Sequence[int]) -> dict:
    validate_graph(adj)
    if next(clique_masks(adj, 7), 0):
        raise ValueError("independent same-Z0 target unexpectedly contains K7")
    inertia = SympyInertiaCache()
    totals = {
        "z0_considered": 0,
        "z0_matchable": 0,
        "z0_block_passed": 0,
        "z0_block_failed": 0,
        "z0_nonbipartite_failed": 0,
    }
    seeds_checked = 0
    first_seed = None
    first_reasons = None
    first_cross_classification = None
    for seed_mask in clique_masks(adj, COORDINATES):
        seeds_checked += 1
        instance = build_instance(adj, seed_mask)
        decision = solve_seed(adj, instance, inertia)
        for name in totals:
            totals[name] += getattr(decision, name)
        if not decision.feasible:
            first_seed = list(instance.seed)
            first_reasons = decision.reasons
            first_cross_classification = independent_cross_classification(
                adj, instance
            )
            break
    rejected = first_seed is not None
    return {
        "rejected": rejected,
        "seeds_checked": seeds_checked,
        "impossible_seeds": int(rejected),
        **totals,
        "first_impossible_seed": first_seed,
        "first_failure_reasons": first_reasons,
        "first_cross_classification": first_cross_classification,
    }


def verify_report_internal_aggregates(report: dict) -> None:
    aggregate_fields = (
        "seeds_checked",
        "impossible_seeds",
        "z0_considered",
        "z0_matchable",
        "z0_block_passed",
        "z0_block_failed",
        "z0_nonbipartite_failed",
        "bipartite_components_checked",
        "bipartite_components_failed",
        "generic_orientation_cases_checked",
        "generic_orientation_support_new_failures",
        "block_subsets_checked",
        "pure_colorings_considered",
        "pure_colorings_matchable",
        "joint_support_searches",
        "joint_support_dfs_nodes",
        "inertia_cache_entries",
        "inertia_cache_hits",
    )
    rejected = [
        item["index"]
        for item in report["graph_results"]
        if item["decision"]["rejected"]
    ]
    if rejected != report["rejected_indices"]:
        raise AssertionError("production report rejection list is not self-consistent")
    if report["graphs_rejected"] != len(rejected):
        raise AssertionError("production report rejected count is not self-consistent")
    for name in aggregate_fields:
        observed = sum(
            item["decision"][name] for item in report["graph_results"]
        )
        if observed != report[name]:
            raise AssertionError(f"production report aggregate mismatch for {name}")


def verify(output: Path) -> dict:
    observed = {name: sha256(ROOT / name) for name in EXPECTED}
    if observed != EXPECTED:
        raise ValueError(f"independent verifier source/input hash mismatch: {observed}")
    report = json.loads(
        (ROOT / "d6_k6_same_z0_report.json").read_text(encoding="utf-8")
    )
    if report.get("schema") != REPORT_SCHEMA or report.get("status") != "COMPLETE":
        raise ValueError("production same-Z0 report is not complete expected schema")
    verify_report_internal_aggregates(report)

    manifest = json.loads(
        (ROOT / "d6_current_residue_manifest.json").read_text(encoding="utf-8")
    )
    indices = manifest["classes"]["K6_only"]["final_indices"]
    if len(indices) != EXPECTED_INPUT or stable_hash(indices) != EXPECTED_INDICES_SHA256:
        raise ValueError("independent same-Z0 selection changed")
    payload = json.loads(
        (ROOT / "d6_k6_bipartite_rank_input.json").read_text(encoding="utf-8")
    )
    by_index = {record["index"]: record for record in payload["graphs"]}
    production_by_index = {
        item["index"]: item["decision"] for item in report["graph_results"]
    }
    if list(production_by_index) != indices:
        raise AssertionError("production graph order differs from exact manifest")

    started = time.perf_counter()
    independently_rejected = []
    rejected_certificates = []
    for position, index in enumerate(indices, 1):
        independent = evaluate_graph(by_index[index]["adjacency"])
        production = production_by_index[index]
        for name in CORE_FIELDS:
            if independent[name] != production[name]:
                raise AssertionError(
                    f"graph {index} field {name}: independent={independent[name]} "
                    f"production={production[name]}"
                )
        if independent["first_impossible_seed"] != production["first_impossible_seed"]:
            raise AssertionError(f"graph {index} first impossible seed differs")
        if not independent["rejected"]:
            continue
        independently_rejected.append(index)
        certificate = production["certificate"]
        if certificate is None:
            raise AssertionError(f"rejected graph {index} lacks production certificate")
        production_reasons = tuple(
            (tuple(item["Z0"]), item["reason"])
            for item in certificate["Z0_failures"]
        )
        if production_reasons != independent["first_failure_reasons"]:
            raise AssertionError(f"graph {index} per-Z0 failure partition differs")
        production_cross = tuple(
            (
                tuple(item["Z0"]),
                item["block_support_passed"],
                item["nonbipartite_joint_support_passed"],
                item["category"],
            )
            for item in certificate["same_Z0_cross_classification"]["choices"]
        )
        if production_cross != independent["first_cross_classification"]:
            raise AssertionError(f"graph {index} full Z0 cross-classification differs")
        histogram = Counter(reason for _, reason in production_reasons)
        cross_histogram = Counter(row[3] for row in production_cross)
        rejected_certificates.append({
            "index": index,
            "first_impossible_seed": independent["first_impossible_seed"],
            "eligible_Z0_choices": len(production_reasons),
            "failure_reason_histogram": dict(sorted(histogram.items())),
            "cross_classification_histogram": dict(sorted(cross_histogram.items())),
        })

    if independently_rejected != report["rejected_indices"]:
        raise AssertionError("independent full rejection set differs")

    positive = evaluate_graph(lower_bound_18_graph())
    if positive["rejected"] or positive["seeds_checked"] != 32:
        raise AssertionError("independent known realizable control failed")

    result = {
        "schema": VERIFICATION_SCHEMA,
        "status": "PASS",
        "method": (
            "independent graph/CSP reconstruction; SymPy exact characteristic-"
            "polynomial Sturm inertia; Hall-subset matching; separately "
            "implemented actual-support DFS"
        ),
        "graphs_recomputed": len(indices),
        "graphs_rejected": len(independently_rejected),
        "graphs_surviving": len(indices) - len(independently_rejected),
        "rejected_indices": independently_rejected,
        "all_graph_core_fields_matched": list(CORE_FIELDS),
        "all_first_impossible_seeds_matched": True,
        "all_rejected_seed_Z0_failure_partitions_matched": True,
        "all_rejected_seed_Z0_cross_classifications_matched": True,
        "rejected_certificates": rejected_certificates,
        "positive_control": {
            "name": "known realizable 18-point construction",
            "passed": not positive["rejected"],
            "K6_seeds": positive["seeds_checked"],
        },
        "inputs": {
            **observed,
            Path(__file__).name: sha256(Path(__file__).resolve()),
        },
        "runtime": {
            "command": " ".join(sys.argv),
            "workers": 1,
            "wall_seconds": time.perf_counter() - started,
            "python": sys.version,
            "platform": platform.platform(),
        },
    }
    atomic_json(output, result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path, default=ROOT / "d6_k6_same_z0_verification.json"
    )
    args = parser.parse_args()
    result = verify(args.output.resolve())
    print(json.dumps({
        "status": result["status"],
        "graphs_recomputed": result["graphs_recomputed"],
        "graphs_rejected": result["graphs_rejected"],
        "rejected_indices": result["rejected_indices"],
        "output": str(args.output.resolve()),
    }, sort_keys=True))


if __name__ == "__main__":
    main()
