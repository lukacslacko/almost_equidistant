#!/usr/bin/env python3
"""Independent verifier for the K6 fused side-rank layer.

This checker imports no production evaluator or K6 production kernel.  It
uses the frozen independent graph/block/nonbipartite reconstruction and the
independent simultaneous-closure zero-forcing solver, then implements the
sidewise max fusion and every block-subset inequality afresh.  It checks all
990 graphs and every archived fused-component certificate.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import verify_d6_k6_all_same_z0 as prior
from test_d6_k6_lorentz import lower_bound_18_graph


independent = prior.base
ROOT = Path(__file__).resolve().parent
COORDINATES = 6
REPORT_SCHEMA = "d6-k6-fused-side-rank-v1"
CERTIFICATE_SCHEMA = "d6-k6-fused-side-rank-certificates-v1"
VERIFICATION_SCHEMA = "d6-k6-fused-side-rank-verification-v1"
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
    "d6_k6_fused_side_rank.py": (
        "5fc6416feca1898f38b1bb3c2f06a1289fcdd6f179654a4ab746e6f3cff869c2"
    ),
    "d6_k6_fused_side_rank_report.json": (
        "60911844102dfeb494170f6b7aa5d5142b8f545fcc968e43078eb07c7ed102ec"
    ),
    "d6_k6_fused_side_rank_certificates.json": (
        "1d497df83b329948c2e09702ca8a94ea94d2d959f7a866a532084d9d920ea664"
    ),
    "verify_d6_k6_all_same_z0.py": (
        "9ad84985cffbc6f53562730f1aa536a6d7046709c21af0dfd1bfbab907a32162"
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


def atomic_json(path: Path, value: object) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


@dataclass(frozen=True)
class BlockResult:
    passed: bool
    first_failure: dict | None


def block_result(blocks: Sequence[tuple[str, int, int]]) -> BlockResult:
    for selected in range(1, 1 << len(blocks)):
        rank = 0
        allowed = 0
        labels = []
        for position, (label, rank_lower, coordinate_mask) in enumerate(blocks):
            if selected & (1 << position):
                labels.append(label)
                rank += rank_lower
                allowed |= coordinate_mask
        if rank > allowed.bit_count():
            return BlockResult(False, {
                "blocks": labels,
                "rank_lower": rank,
                "coordinate_capacity": allowed.bit_count(),
                "allowed_coordinates": independent.bits(allowed),
            })
    return BlockResult(True, None)


def allowed_union(instance: independent.Instance, selected: int) -> int:
    answer = 0
    for local in independent.bits(selected):
        answer |= instance.defects[local]
    return answer


def make_blocks(
    instance: independent.Instance,
    z0: int,
    component: independent.Component,
    rank_a: int,
    rank_b: int,
) -> tuple[tuple[str, int, int], ...]:
    blocks = [
        ("A", rank_a, allowed_union(instance, component.side_a)),
        ("B", rank_b, allowed_union(instance, component.side_b)),
    ]
    blocks.extend(
        (
            f"Z0:{instance.outside[local]}",
            1,
            instance.defects[local],
        )
        for local in independent.bits(z0)
    )
    return tuple(blocks)


def fused_component(
    adj: Sequence[int],
    instance: independent.Instance,
    z0: int,
    component: independent.Component,
    inertia: independent.SympyInertiaCache,
    zero_forcing: prior.IndependentZeroForcing,
) -> tuple[bool, dict]:
    graph_a, absolute_a = prior.induced_adjacency(
        adj, instance, component.side_a
    )
    graph_b, absolute_b = prior.induced_adjacency(
        adj, instance, component.side_b
    )
    zf_a = zero_forcing.number(graph_a)
    zf_b = zero_forcing.number(graph_b)
    lower_zf_a = len(absolute_a) - zf_a
    lower_zf_b = len(absolute_b) - zf_b
    zf_total = z0.bit_count() + lower_zf_a + lower_zf_b
    zf_passed = zf_total <= COORDINATES

    inertia_a = inertia.solve(independent.side_pattern(adj, instance, component.side_a))
    inertia_b = inertia.solve(independent.side_pattern(adj, instance, component.side_b))
    positive_a, negative_a, zero_a = inertia_a
    positive_b, negative_b, zero_b = inertia_b
    size_a = component.side_a.bit_count()
    size_b = component.side_b.bit_count()
    orientations = (
        ("A_positive", size_a - negative_a, size_b - positive_b),
        ("A_negative", size_a - positive_a, size_b - negative_b),
    )
    details = []
    for name, rank_a, rank_b in orientations:
        old = block_result(make_blocks(instance, z0, component, rank_a, rank_b))
        fused_a = max(lower_zf_a, rank_a)
        fused_b = max(lower_zf_b, rank_b)
        fused = block_result(
            make_blocks(instance, z0, component, fused_a, fused_b)
        )
        fusion_only = zf_passed and old.passed and not fused.passed
        details.append({
            "case": name,
            "inertia_rank_lower_A": rank_a,
            "inertia_rank_lower_B": rank_b,
            "zero_forcing_rank_lower_A": lower_zf_a,
            "zero_forcing_rank_lower_B": lower_zf_b,
            "fused_rank_lower_A": fused_a,
            "fused_rank_lower_B": fused_b,
            "old_inertia_support_passed": old.passed,
            "old_inertia_first_failure": old.first_failure,
            "fused_support_passed": fused.passed,
            "fused_first_failure": fused.first_failure,
            "old_separate_tests_passed_but_fused_failed": fusion_only,
        })

    light_selected = z0 | component.vertices
    light_dimension = light_selected.bit_count()
    light_matching = (
        light_dimension <= COORDINATES
        and independent.hall_matchable(light_selected, instance.defects)
    )
    light = {
        "orthonormal_vectors": light_dimension,
        "dimension_passed": light_dimension <= COORDINATES,
        "allowed_mask_matching_passed": light_matching,
        "passed": light_matching,
    }
    old_block_passed = bool(
        any(case["old_inertia_support_passed"] for case in details)
        or light_matching
    )
    old_both = zf_passed and old_block_passed
    fused_passed = bool(
        any(case["fused_support_passed"] for case in details) or light_matching
    )
    detail = {
        "component": [
            instance.outside[local]
            for local in independent.bits(component.vertices)
        ],
        "side_A": list(absolute_a),
        "side_B": list(absolute_b),
        "zero_forcing": {
            "number_A": zf_a,
            "number_B": zf_b,
            "rank_lower_A": lower_zf_a,
            "rank_lower_B": lower_zf_b,
            "zero_factor_dimension": z0.bit_count(),
            "required_dimension": zf_total,
            "passed": zf_passed,
        },
        "inertia_A": [positive_a, negative_a, zero_a],
        "inertia_B": [positive_b, negative_b, zero_b],
        "orientations": details,
        "lightlike": light,
        "old_separate_systems_passed": old_both,
        "fused_passed": fused_passed,
    }
    return fused_passed, detail


@dataclass
class SeedDecision:
    feasible: bool
    z0_considered: int
    z0_matchable: int
    z0_fused_passed: int
    z0_fused_failed: int
    z0_nonbipartite_failed: int


def solve_seed(
    adj: Sequence[int],
    instance: independent.Instance,
    inertia: independent.SympyInertiaCache,
    zero_forcing: prior.IndependentZeroForcing,
) -> SeedDecision:
    counts = {
        "z0_considered": 0,
        "z0_matchable": 0,
        "z0_fused_passed": 0,
        "z0_fused_failed": 0,
        "z0_nonbipartite_failed": 0,
    }
    for z0 in independent.z0_subsets(instance.eligible_z0):
        counts["z0_considered"] += 1
        if not independent.hall_matchable(z0, instance.defects):
            continue
        counts["z0_matchable"] += 1
        local_components = independent.components(instance, z0)
        fused_passed = all(
            fused_component(
                adj, instance, z0, component, inertia, zero_forcing
            )[0]
            for component in local_components
            if component.bipartite
        )
        if not fused_passed:
            counts["z0_fused_failed"] += 1
            continue
        counts["z0_fused_passed"] += 1
        odd = [
            component.vertices
            for component in local_components
            if not component.bipartite
        ]
        if independent.nonbipartite_system_passes(instance, z0, odd):
            return SeedDecision(True, **counts)
        counts["z0_nonbipartite_failed"] += 1
    return SeedDecision(False, **counts)


CORE_FIELDS = (
    "rejected",
    "seeds_checked",
    "impossible_seeds",
    "z0_considered",
    "z0_matchable",
    "z0_fused_passed",
    "z0_fused_failed",
    "z0_nonbipartite_failed",
)


def evaluate_graph(adj: Sequence[int]) -> dict:
    independent.validate_graph(adj)
    if next(independent.clique_masks(adj, 7), 0):
        raise ValueError("independent fused target unexpectedly contains K7")
    inertia = independent.SympyInertiaCache()
    zero_forcing = prior.IndependentZeroForcing()
    totals = {name: 0 for name in CORE_FIELDS[3:]}
    seeds_checked = 0
    first_seed = None
    first_instance = None
    for seed_mask in independent.clique_masks(adj, COORDINATES):
        seeds_checked += 1
        instance = independent.build_instance(adj, seed_mask)
        decision = solve_seed(adj, instance, inertia, zero_forcing)
        for name in totals:
            totals[name] += getattr(decision, name)
        if not decision.feasible:
            first_seed = list(instance.seed)
            first_instance = instance
            break
    rejected = first_seed is not None
    return {
        "rejected": rejected,
        "seeds_checked": seeds_checked,
        "impossible_seeds": int(rejected),
        **totals,
        "first_impossible_seed": first_seed,
        "first_instance": first_instance,
    }


def certificate_rows(
    adj: Sequence[int], instance: independent.Instance
) -> tuple[dict, ...]:
    inertia = independent.SympyInertiaCache()
    zero_forcing = prior.IndependentZeroForcing()
    rows = []
    for z0 in independent.z0_subsets(instance.eligible_z0):
        z0_absolute = tuple(instance.outside[local] for local in independent.bits(z0))
        matchable = independent.hall_matchable(z0, instance.defects)
        if matchable:
            local_components = independent.components(instance, z0)
            decisions = [
                fused_component(
                    adj, instance, z0, component, inertia, zero_forcing
                )
                for component in local_components
                if component.bipartite
            ]
            fused_passed = all(item[0] for item in decisions)
            failed = tuple(item[1] for item in decisions if not item[0])
            odd = [
                component.vertices
                for component in local_components
                if not component.bipartite
            ]
            nonbipartite_passed = independent.nonbipartite_system_passes(
                instance, z0, odd
            )
        else:
            fused_passed = False
            failed = ()
            nonbipartite_passed = False
        if fused_passed and nonbipartite_passed:
            category = "fused+nonbipartite"
        elif fused_passed:
            category = "fused_only"
        elif nonbipartite_passed:
            category = "nonbipartite_only"
        else:
            category = "neither"
        rows.append({
            "Z0": z0_absolute,
            "matchable": matchable,
            "fused": fused_passed,
            "nonbipartite": nonbipartite_passed,
            "category": category,
            "failed_fused_components": failed,
        })
    return tuple(rows)


def verify_internal(report: dict, archive: dict) -> None:
    if report.get("schema") != REPORT_SCHEMA or report.get("status") != "COMPLETE":
        raise ValueError("fused production report has wrong schema/status")
    if archive.get("schema") != CERTIFICATE_SCHEMA or archive.get("status") != "COMPLETE":
        raise ValueError("fused production archive has wrong schema/status")
    if report["certificate_archive"]["sha256"] != sha256(
        ROOT / report["certificate_archive"]["path"]
    ):
        raise AssertionError("fused certificate archive hash differs")
    rejected = [
        item["index"]
        for item in report["graph_results"]
        if item["decision"]["rejected"]
    ]
    if rejected != report["rejected_indices"] or rejected != archive["rejected_indices"]:
        raise AssertionError("fused production rejection lists disagree")
    numeric = tuple(
        name
        for name in report["graph_results"][0]["decision"]
        if name not in ("rejected", "first_impossible_seed")
    )
    for name in numeric:
        observed = sum(
            item["decision"][name] for item in report["graph_results"]
        )
        if observed != report[name]:
            raise AssertionError(f"fused report aggregate differs for {name}")


def verify(output: Path) -> dict:
    observed = {name: sha256(ROOT / name) for name in EXPECTED}
    if observed != EXPECTED:
        raise ValueError(f"fused verifier dependency hash mismatch: {observed}")
    report = json.loads(
        (ROOT / "d6_k6_fused_side_rank_report.json").read_text(encoding="utf-8")
    )
    archive = json.loads(
        (ROOT / "d6_k6_fused_side_rank_certificates.json").read_text(
            encoding="utf-8"
        )
    )
    verify_internal(report, archive)
    manifest = json.loads(
        (ROOT / "d6_current_residue_manifest.json").read_text(encoding="utf-8")
    )
    indices = manifest["classes"]["K6_only"]["final_indices"]
    if len(indices) != EXPECTED_INPUT or independent.stable_hash(indices) != EXPECTED_INDICES_SHA256:
        raise ValueError("independent fused selection differs")
    payload = json.loads(
        (ROOT / "d6_k6_bipartite_rank_input.json").read_text(encoding="utf-8")
    )
    by_index = {record["index"]: record for record in payload["graphs"]}
    production = {
        item["index"]: item["decision"] for item in report["graph_results"]
    }
    certificates = {
        item["index"]: item["certificate"] for item in archive["certificates"]
    }

    started = time.perf_counter()
    rejected = []
    summaries = []
    for index in indices:
        adjacency = by_index[index]["adjacency"]
        result = evaluate_graph(adjacency)
        expected = production[index]
        for name in CORE_FIELDS:
            if result[name] != expected[name]:
                raise AssertionError(
                    f"graph {index} {name}: independent={result[name]} "
                    f"production={expected[name]}"
                )
        if result["first_impossible_seed"] != expected["first_impossible_seed"]:
            raise AssertionError(f"graph {index} first impossible seed differs")
        if not result["rejected"]:
            continue
        rejected.append(index)
        certificate = certificates[index]
        independent_rows = certificate_rows(adjacency, result["first_instance"])
        production_rows = certificate["choices"]
        if len(independent_rows) != len(production_rows):
            raise AssertionError(f"graph {index} certificate length differs")
        histogram: dict[str, int] = {}
        for left, right in zip(independent_rows, production_rows):
            observed_row = {
                "Z0": tuple(right["Z0"]),
                "matchable": right["Z0_allowed_mask_matchable"],
                "fused": right["fused_side_rank_passed"],
                "nonbipartite": right["nonbipartite_joint_support_passed"],
                "category": right["category"],
                "failed_fused_components": tuple(right["failed_fused_components"]),
            }
            if left != observed_row:
                raise AssertionError(f"graph {index} per-Z0 fused certificate differs")
            histogram[left["category"]] = histogram.get(left["category"], 0) + 1
        if histogram != certificate["histogram"]:
            raise AssertionError(f"graph {index} fused histogram differs")
        if histogram.get("fused+nonbipartite", 0):
            raise AssertionError(f"graph {index} archive contains a common pass")
        summaries.append({
            "index": index,
            "first_impossible_seed": result["first_impossible_seed"],
            "Z0_choices": len(independent_rows),
            "histogram": histogram,
        })

    if rejected != report["rejected_indices"]:
        raise AssertionError("independent fused rejection set differs")
    positive = evaluate_graph(lower_bound_18_graph())
    if positive["rejected"] or positive["seeds_checked"] != 32:
        raise AssertionError("known realizable 18-point fused control failed")

    result = {
        "schema": VERIFICATION_SCHEMA,
        "status": "PASS",
        "method": (
            "independent SymPy/Sturm inertia, independent simultaneous-closure "
            "zero forcing, fresh sidewise-max block enumeration, and frozen "
            "independent nonbipartite support reconstruction"
        ),
        "graphs_recomputed": len(indices),
        "graphs_rejected": len(rejected),
        "graphs_surviving": len(indices) - len(rejected),
        "rejected_indices": rejected,
        "new_rejected_indices": report["new_rejected_indices"],
        "all_graph_core_fields_matched": list(CORE_FIELDS),
        "all_first_impossible_seeds_matched": True,
        "all_archived_fused_component_certificates_matched": True,
        "all_archived_nonbipartite_booleans_matched": True,
        "certificate_summaries": summaries,
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
        "--output",
        type=Path,
        default=ROOT / "d6_k6_fused_side_rank_verification.json",
    )
    args = parser.parse_args()
    result = verify(args.output.resolve())
    print(json.dumps({
        "status": result["status"],
        "graphs_recomputed": result["graphs_recomputed"],
        "graphs_rejected": result["graphs_rejected"],
        "rejected_indices": result["rejected_indices"],
        "new_rejected_indices": result["new_rejected_indices"],
        "output": str(args.output.resolve()),
    }, sort_keys=True))


if __name__ == "__main__":
    main()
