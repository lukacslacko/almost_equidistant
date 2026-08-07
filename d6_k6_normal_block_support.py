#!/usr/bin/env python3
"""Prototype support-Hall refinement of the K6 normal-inertia rule.

For one bipartite Lorentz component, the A span, B span, and every individual
zero-factor vector are mutually orthogonal blocks in R^6.  Each block is
contained in the coordinate span of the union of its allowed defect masks.
Thus every subset of the at-most-eight blocks must have total inertia rank
lower bound no larger than its coordinate-support union.  Both generic
Lorentz sign orientations are tested; the lightlike case is retained as a
separate conservative alternative.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import sys
import time
from dataclasses import asdict, dataclass
from itertools import combinations
from pathlib import Path
from typing import Iterator, Sequence

from d6_k6_bipartite_rank_reference import lorentz_components
from d6_k6_lorentz_reference import (
    COORDINATES,
    K6LorentzInstance,
    build_instance,
    clique_masks,
    support_matching,
    validate_graph,
    vertices,
)
from d6_k6_normal_inertia import InertiaCache, side_unit_pattern


ROOT = Path(__file__).resolve().parent
EXPECTED = {
    "d6_current_residue_manifest.json": (
        "d06cdb23b03b1f8523fbe3c3ede5bc297ace21a844282830c5cf51a6cee3512d"
    ),
    "d6_k6_bipartite_rank_input.json": (
        "dfee979e6ca2759e80cf803c3054c8954d98543875d467644b91ee434c22e845"
    ),
    "d6_k6_lorentz_reference.py": (
        "ef6a0877a5987119c92ec1bd9271a6a47774fc9e4992c10e1f28b6d5373592c2"
    ),
    "d6_k6_support_reference.py": (
        "d6481137ef49d88154882285660dd755eb7cb652ab87280f05c3792de9742577"
    ),
    "d6_k6_bipartite_rank_reference.py": (
        "e760fef74c423a89f1e1f9d08909fa1223e0389900cf3d45b0140ee296dea922"
    ),
    "d6_k6_normal_inertia.py": (
        "2d00ab40eb97aa79ffec3b8134b83c4fd5f703cbd691f094a25bf3ff9f3ee023"
    ),
}
EXPECTED_INPUT = 990
EXPECTED_INDICES_SHA256 = (
    "16cb562473ebb8c496018891362642108113a16f007929d413770ccdb665acb9"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_hash(value: object) -> str:
    raw = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(raw.encode("ascii")).hexdigest()


def z0_subsets(eligible: int) -> Iterator[int]:
    candidates = vertices(eligible)
    for size in range(min(COORDINATES, len(candidates)) + 1):
        for selected in combinations(candidates, size):
            yield sum(1 << vertex for vertex in selected)


def allowed_union(instance: K6LorentzInstance, selected: int) -> int:
    union = 0
    for vertex in vertices(selected):
        union |= instance.defects[vertex]
    return union


@dataclass(frozen=True)
class OrthogonalBlock:
    label: str
    rank_lower: int
    allowed_coordinates: int


@dataclass(frozen=True)
class BlockSubsetResult:
    passed: bool
    subsets_checked: int
    first_failure: dict | None


def check_block_subsets(blocks: Sequence[OrthogonalBlock]) -> BlockSubsetResult:
    """Apply every subset dimension/support-union inequality exactly."""

    if len(blocks) > 8:
        raise ValueError("K6 generic component has at most eight blocks")
    checked = 0
    for selected in range(1, 1 << len(blocks)):
        checked += 1
        rank = 0
        coordinates = 0
        labels = []
        for position, block in enumerate(blocks):
            if selected & (1 << position):
                rank += block.rank_lower
                coordinates |= block.allowed_coordinates
                labels.append(block.label)
        capacity = coordinates.bit_count()
        if rank > capacity:
            return BlockSubsetResult(
                False,
                checked,
                {
                    "blocks": labels,
                    "rank_lower": rank,
                    "coordinate_capacity": capacity,
                    "allowed_coordinates": vertices(coordinates),
                },
            )
    return BlockSubsetResult(True, checked, None)


def generic_blocks(
    instance: K6LorentzInstance,
    z0: int,
    side_a: int,
    side_b: int,
    rank_a: int,
    rank_b: int,
) -> tuple[OrthogonalBlock, ...]:
    blocks = [
        OrthogonalBlock("A", rank_a, allowed_union(instance, side_a)),
        OrthogonalBlock("B", rank_b, allowed_union(instance, side_b)),
    ]
    for local in vertices(z0):
        blocks.append(
            OrthogonalBlock(
                f"Z0:{instance.outside[local]}", 1, instance.defects[local]
            )
        )
    return tuple(blocks)


@dataclass(frozen=True)
class ComponentDecision:
    passed: bool
    old_ambient_rule_passed: bool
    old_passed_new_failed: bool
    generic_cases: tuple[dict, dict]
    lightlike_case: dict
    subsets_checked: int


def check_component(
    adj: Sequence[int],
    instance: K6LorentzInstance,
    z0: int,
    component,
    cache: InertiaCache,
) -> ComponentDecision:
    if not component.bipartite:
        raise ValueError("support-Hall component checker expects bipartite input")
    _, pattern_a = side_unit_pattern(adj, instance, component.side_a)
    _, pattern_b = side_unit_pattern(adj, instance, component.side_b)
    inertia_a = cache.solve(pattern_a)
    inertia_b = cache.solve(pattern_b)
    size_a = component.side_a.bit_count()
    size_b = component.side_b.bit_count()
    orientations = (
        (
            "A_positive",
            size_a - inertia_a[1],
            size_b - inertia_b[0],
        ),
        (
            "A_negative",
            size_a - inertia_a[0],
            size_b - inertia_b[1],
        ),
    )
    generic = []
    subsets_checked = 0
    for name, rank_a, rank_b in orientations:
        blocks = generic_blocks(
            instance, z0, component.side_a, component.side_b, rank_a, rank_b
        )
        result = check_block_subsets(blocks)
        subsets_checked += result.subsets_checked
        ambient_dimension = z0.bit_count() + rank_a + rank_b
        generic.append({
            "case": name,
            "rank_lower_A": rank_a,
            "rank_lower_B": rank_b,
            "ambient_dimension_lower": ambient_dimension,
            "ambient_dimension_passed": ambient_dimension <= COORDINATES,
            "blocks": [asdict(block) for block in blocks],
            "passed": result.passed,
            "subsets_checked": result.subsets_checked,
            "first_failure": result.first_failure,
        })

    # When the Lorentz direction is lightlike, all component and Z0 defect
    # vectors are individually orthonormal.  A realization therefore needs
    # at most six such vectors and a full matching into their allowed masks.
    light_selected = z0 | component.component
    light_dimension = light_selected.bit_count()
    light_matching = (
        light_dimension <= COORDINATES
        and support_matching(light_selected, instance.defects) is not None
    )
    light = {
        "case": "lightlike",
        "orthonormal_vectors": light_dimension,
        "dimension_passed": light_dimension <= COORDINATES,
        "allowed_mask_matching_passed": light_matching,
        "passed": light_matching,
    }
    old_passed = bool(
        generic[0]["ambient_dimension_passed"]
        or generic[1]["ambient_dimension_passed"]
    )
    new_passed = bool(
        generic[0]["passed"] or generic[1]["passed"] or light["passed"]
    )
    return ComponentDecision(
        new_passed,
        old_passed,
        old_passed and not new_passed,
        (generic[0], generic[1]),
        light,
        subsets_checked,
    )


@dataclass
class SeedDecision:
    feasible: bool
    z0_considered: int
    z0_matchable: int
    z0_failed: int
    bipartite_components_checked: int
    generic_cases_checked: int
    generic_cases_ambient_passed: int
    generic_cases_failed: int
    generic_cases_support_new_failures: int
    block_subsets_checked: int
    lightlike_cases_passed: int
    components_old_passed_new_failed: int
    z0_old_passed_new_failed: int
    chosen_z0: list[int] | None
    witness: dict | None


def solve_seed(
    adj: Sequence[int], instance: K6LorentzInstance, cache: InertiaCache
) -> SeedDecision:
    counters = {
        "z0_considered": 0,
        "z0_matchable": 0,
        "z0_failed": 0,
        "components": 0,
        "generic": 0,
        "generic_ambient_passed": 0,
        "generic_failed": 0,
        "generic_support_new_failures": 0,
        "subsets": 0,
        "light_passed": 0,
        "components_old_passed_new_failed": 0,
        "z0_old_passed_new_failed": 0,
    }
    failed_examples = []
    for z0 in z0_subsets(instance.eligible_z0_mask):
        counters["z0_considered"] += 1
        if support_matching(z0, instance.defects) is None:
            continue
        counters["z0_matchable"] += 1
        failed = []
        old_z0_passed = True
        for component in lorentz_components(instance, z0):
            if not component.bipartite:
                continue
            decision = check_component(adj, instance, z0, component, cache)
            counters["components"] += 1
            counters["generic"] += 2
            counters["generic_ambient_passed"] += sum(
                case["ambient_dimension_passed"]
                for case in decision.generic_cases
            )
            counters["generic_failed"] += sum(
                not case["passed"] for case in decision.generic_cases
            )
            counters["generic_support_new_failures"] += sum(
                case["ambient_dimension_passed"] and not case["passed"]
                for case in decision.generic_cases
            )
            counters["subsets"] += decision.subsets_checked
            counters["light_passed"] += decision.lightlike_case["passed"]
            counters["components_old_passed_new_failed"] += (
                decision.old_passed_new_failed
            )
            old_z0_passed = old_z0_passed and decision.old_ambient_rule_passed
            if not decision.passed:
                failed.append({
                    "component": [
                        instance.outside[local]
                        for local in vertices(component.component)
                    ],
                    "generic_cases": decision.generic_cases,
                    "lightlike_case": decision.lightlike_case,
                })
        if not failed:
            return SeedDecision(
                True,
                counters["z0_considered"],
                counters["z0_matchable"],
                counters["z0_failed"],
                counters["components"],
                counters["generic"],
                counters["generic_ambient_passed"],
                counters["generic_failed"],
                counters["generic_support_new_failures"],
                counters["subsets"],
                counters["light_passed"],
                counters["components_old_passed_new_failed"],
                counters["z0_old_passed_new_failed"],
                [instance.outside[local] for local in vertices(z0)],
                None,
            )
        counters["z0_failed"] += 1
        counters["z0_old_passed_new_failed"] += old_z0_passed
        if len(failed_examples) < 8:
            failed_examples.append({
                "Z0": [instance.outside[local] for local in vertices(z0)],
                "failed_components": failed,
            })
    return SeedDecision(
        False,
        counters["z0_considered"],
        counters["z0_matchable"],
        counters["z0_failed"],
        counters["components"],
        counters["generic"],
        counters["generic_ambient_passed"],
        counters["generic_failed"],
        counters["generic_support_new_failures"],
        counters["subsets"],
        counters["light_passed"],
        counters["components_old_passed_new_failed"],
        counters["z0_old_passed_new_failed"],
        None,
        {
            "seed": list(instance.seed),
            "outside": list(instance.outside),
            "allowed_defects": {
                str(instance.outside[local]): vertices(mask)
                for local, mask in enumerate(instance.defects)
            },
            "failed_examples": failed_examples,
        },
    )


def evaluate_graph(adj: Sequence[int]) -> dict:
    # The disjoint-support consequence for candidate nonedges uses
    # alpha(G) <= 2, so keep that hypothesis at every public graph entry.
    validate_graph(adj, require_alpha_two=True)
    totals = {
        "z0_considered": 0,
        "z0_matchable": 0,
        "z0_failed": 0,
        "bipartite_components_checked": 0,
        "generic_cases_checked": 0,
        "generic_cases_ambient_passed": 0,
        "generic_cases_failed": 0,
        "generic_cases_support_new_failures": 0,
        "block_subsets_checked": 0,
        "lightlike_cases_passed": 0,
        "components_old_passed_new_failed": 0,
        "z0_old_passed_new_failed": 0,
    }
    cache = InertiaCache()
    seeds = 0
    impossible = 0
    first_witness = None
    for seed_mask in clique_masks(adj, COORDINATES):
        seeds += 1
        instance = build_instance(adj, vertices(seed_mask))
        decision = solve_seed(adj, instance, cache)
        for key in totals:
            totals[key] += getattr(decision, key)
        if not decision.feasible:
            impossible += 1
            if first_witness is None:
                first_witness = decision.witness
    return {
        "rejected": impossible > 0,
        "seeds_checked": seeds,
        "impossible_seeds": impossible,
        **totals,
        "inertia_cache_entries": len(cache.values),
        "inertia_cache_hits": cache.hits,
        "first_witness": first_witness,
    }


def atomic_json(path: Path, value: object) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def run(output: Path) -> dict:
    observed = {name: sha256(ROOT / name) for name in EXPECTED}
    if observed != EXPECTED:
        raise ValueError("support-Hall prototype source hash mismatch")
    manifest = json.loads((ROOT / "d6_current_residue_manifest.json").read_text())
    selected = manifest["classes"]["K6_only"]["final_indices"]
    if len(selected) != EXPECTED_INPUT or stable_hash(selected) != EXPECTED_INDICES_SHA256:
        raise ValueError("current K6 residue selection mismatch")
    payload = json.loads((ROOT / "d6_k6_bipartite_rank_input.json").read_text())
    by_index = {record["index"]: record["adjacency"] for record in payload["graphs"]}
    if any(index not in by_index for index in selected):
        raise ValueError("a selected K6 index is absent from adjacency source")
    started = time.time()
    graph_results = []
    aggregate = {
        "seeds_checked": 0,
        "impossible_seeds": 0,
        "z0_considered": 0,
        "z0_matchable": 0,
        "z0_failed": 0,
        "bipartite_components_checked": 0,
        "generic_cases_checked": 0,
        "generic_cases_ambient_passed": 0,
        "generic_cases_failed": 0,
        "generic_cases_support_new_failures": 0,
        "block_subsets_checked": 0,
        "lightlike_cases_passed": 0,
        "components_old_passed_new_failed": 0,
        "z0_old_passed_new_failed": 0,
    }
    rejected = []
    for index in selected:
        decision = evaluate_graph(by_index[index])
        graph_results.append({"index": index, "decision": decision})
        if decision["rejected"]:
            rejected.append(index)
        for key in aggregate:
            aggregate[key] += decision[key]
    report = {
        "schema": "d6-k6-normal-block-support-v1",
        "status": "COMPLETE",
        "description": (
            "Exact block-subset support-Hall refinement for both generic "
            "Lorentz orientations, with a separate conservative lightlike case."
        ),
        "input_graphs": len(selected),
        "input_indices_sha256": stable_hash(selected),
        "graphs_rejected": len(rejected),
        "graphs_surviving": len(selected) - len(rejected),
        "rejected_indices": rejected,
        **aggregate,
        "graph_results": graph_results,
        "sources": {
            **observed,
            Path(__file__).name: sha256(Path(__file__).resolve()),
        },
        "runtime": {
            "workers": 1,
            "wall_seconds": time.time() - started,
            "python": sys.version,
            "platform": platform.platform(),
            "command": " ".join(sys.argv),
        },
        "semantics": {
            "candidate_nonedges": "unconstrained and may be unit",
            "allowed_masks": "upper bounds on actual supports; coordinates may be zero",
            "survivor": "filter non-rejection only, not a realizability claim",
        },
    }
    atomic_json(output, report)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path,
        default=ROOT / "d6_k6_normal_block_support_report.json",
    )
    args = parser.parse_args()
    report = run(args.output.resolve())
    print(json.dumps({
        "status": report["status"],
        "input_graphs": report["input_graphs"],
        "graphs_rejected": report["graphs_rejected"],
        "graphs_surviving": report["graphs_surviving"],
        "output": str(args.output.resolve()),
    }, sort_keys=True))


if __name__ == "__main__":
    main()
