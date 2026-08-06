#!/usr/bin/env python3
"""Exact sidewise fusion of K6 zero-forcing and inertia rank bounds.

For one generic bipartite Lorentz component C=A union B, ordinary zero
forcing and orientation-dependent inertia both lower-bound the dimensions of
the same actual Euclidean spans U_A and U_B.  Therefore each side has the
lower bound max(zero-forcing, inertia), separately, before the orthogonal
block-support inequalities are applied.  This is stronger than checking the
two total-dimension rules independently.

Every eligible Z0 is exhausted, the lightlike component alternative is kept
separate, and surviving fused choices are coupled to the non-bipartite joint
support system on the same Z0.  All decisions are exact.  Candidate nonedges
remain unconstrained and allowed defect coordinates may be zero.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

from d6_k6_all_same_z0 import (
    absolute_vertices,
    check_nonbipartite_system,
)
from d6_k6_bipartite_rank_reference import (
    ZeroForcingSolver,
    induced_required_graph,
    lorentz_components,
)
from d6_k6_lorentz_reference import (
    COORDINATES,
    K6LorentzInstance,
    build_instance,
    clique_masks,
    find_clique_mask,
    support_matching,
    validate_graph,
    vertices,
)
from d6_k6_normal_block_support import (
    OrthogonalBlock,
    allowed_union,
    check_block_subsets,
    z0_subsets,
)
from d6_k6_normal_inertia import InertiaCache, side_unit_pattern
from test_d6_k6_lorentz import lower_bound_18_graph


ROOT = Path(__file__).resolve().parent
REPORT_SCHEMA = "d6-k6-fused-side-rank-v1"
CERTIFICATE_SCHEMA = "d6-k6-fused-side-rank-certificates-v1"
EXPECTED_INPUT = 990
EXPECTED_INDICES_SHA256 = (
    "16cb562473ebb8c496018891362642108113a16f007929d413770ccdb665acb9"
)
EXPECTED_PRIOR_REJECTIONS = [652_900, 2_301_548, 2_842_523, 3_289_061]
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
    "d6_k6_normal_block_support.py": (
        "14faafc6b2f9c3c34c954bb9cbaa21e9a7f646cfb76fcdfb30dfe7af742cf9cf"
    ),
    "d6_k6_all_same_z0.py": (
        "d15fe376993393a5106e5d01f8c0b2316f43945c5fd9ee01b1ebdabfb4fe2b94"
    ),
    "d6_k6_all_same_z0_report.json": (
        "8fa616455d930a1ac3c2255ecb123d8095f7731780f6b6d0bec1ba71c173edeb"
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


def make_blocks(
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
    blocks.extend(
        OrthogonalBlock(
            f"Z0:{instance.outside[local]}", 1, instance.defects[local]
        )
        for local in vertices(z0)
    )
    return tuple(blocks)


@dataclass(frozen=True)
class FusedComponentDecision:
    passed: bool
    old_separate_systems_passed: bool
    old_passed_fused_failed: bool
    detail: dict
    old_subsets_checked: int
    fused_subsets_checked: int
    orientation_fusion_new_failures: int


def check_fused_component(
    adj: Sequence[int],
    instance: K6LorentzInstance,
    z0: int,
    component,
    inertia: InertiaCache,
    zero_forcing: ZeroForcingSolver,
) -> FusedComponentDecision:
    if not component.bipartite:
        raise ValueError("fused side-rank checker requires a bipartite component")

    graph_a, absolute_a = induced_required_graph(
        adj, instance, component.side_a
    )
    graph_b, absolute_b = induced_required_graph(
        adj, instance, component.side_b
    )
    zf_a = zero_forcing.solve(graph_a).number
    zf_b = zero_forcing.solve(graph_b).number
    lower_zf_a = len(absolute_a) - zf_a
    lower_zf_b = len(absolute_b) - zf_b
    zf_total = z0.bit_count() + lower_zf_a + lower_zf_b
    zf_passed = zf_total <= COORDINATES

    _, pattern_a = side_unit_pattern(adj, instance, component.side_a)
    _, pattern_b = side_unit_pattern(adj, instance, component.side_b)
    positive_a, negative_a, zero_a = inertia.solve(pattern_a)
    positive_b, negative_b, zero_b = inertia.solve(pattern_b)
    size_a = component.side_a.bit_count()
    size_b = component.side_b.bit_count()
    orientations = (
        ("A_positive", size_a - negative_a, size_b - positive_b),
        ("A_negative", size_a - positive_a, size_b - negative_b),
    )

    orientation_details = []
    old_checked = 0
    fused_checked = 0
    new_failures = 0
    for name, inertia_a, inertia_b in orientations:
        old_blocks = make_blocks(
            instance,
            z0,
            component.side_a,
            component.side_b,
            inertia_a,
            inertia_b,
        )
        fused_a = max(lower_zf_a, inertia_a)
        fused_b = max(lower_zf_b, inertia_b)
        fused_blocks = make_blocks(
            instance,
            z0,
            component.side_a,
            component.side_b,
            fused_a,
            fused_b,
        )
        old = check_block_subsets(old_blocks)
        fused = check_block_subsets(fused_blocks)
        old_checked += old.subsets_checked
        fused_checked += fused.subsets_checked
        fusion_only = zf_passed and old.passed and not fused.passed
        new_failures += fusion_only
        orientation_details.append({
            "case": name,
            "inertia_rank_lower_A": inertia_a,
            "inertia_rank_lower_B": inertia_b,
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

    light_selected = z0 | component.component
    light_dimension = light_selected.bit_count()
    light_matching = (
        light_dimension <= COORDINATES
        and support_matching(light_selected, instance.defects) is not None
    )
    light = {
        "orthonormal_vectors": light_dimension,
        "dimension_passed": light_dimension <= COORDINATES,
        "allowed_mask_matching_passed": light_matching,
        "passed": light_matching,
    }
    old_block_passed = bool(
        any(case["old_inertia_support_passed"] for case in orientation_details)
        or light_matching
    )
    old_both_passed = zf_passed and old_block_passed
    fused_passed = bool(
        any(case["fused_support_passed"] for case in orientation_details)
        or light_matching
    )
    detail = {
        "component": absolute_vertices(instance, component.component),
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
        "orientations": orientation_details,
        "lightlike": light,
        "old_separate_systems_passed": old_both_passed,
        "fused_passed": fused_passed,
    }
    return FusedComponentDecision(
        fused_passed,
        old_both_passed,
        old_both_passed and not fused_passed,
        detail,
        old_checked,
        fused_checked,
        new_failures,
    )


COUNTER_FIELDS = (
    "z0_considered",
    "z0_matchable",
    "z0_fused_passed",
    "z0_fused_failed",
    "z0_nonbipartite_failed",
    "bipartite_components_checked",
    "bipartite_components_failed",
    "components_old_separate_passed_fused_failed",
    "generic_orientation_cases_checked",
    "orientation_fusion_new_failures",
    "old_block_subsets_checked",
    "fused_block_subsets_checked",
    "pure_colorings_considered",
    "pure_colorings_matchable",
    "joint_support_searches",
    "joint_support_dfs_nodes",
)


def empty_counts() -> dict[str, int]:
    return {name: 0 for name in COUNTER_FIELDS}


def add_counts(target: dict[str, int], source: dict[str, int]) -> None:
    for name in COUNTER_FIELDS:
        target[name] += source[name]


def fused_system(
    adj: Sequence[int],
    instance: K6LorentzInstance,
    z0: int,
    components,
    inertia: InertiaCache,
    zero_forcing: ZeroForcingSolver,
) -> tuple[bool, list[dict], dict[str, int]]:
    failures = []
    counts = {
        "components": 0,
        "failed": 0,
        "old_passed_fused_failed": 0,
        "generic": 0,
        "orientation_new": 0,
        "old_subsets": 0,
        "fused_subsets": 0,
    }
    for component in components:
        if not component.bipartite:
            continue
        decision = check_fused_component(
            adj, instance, z0, component, inertia, zero_forcing
        )
        counts["components"] += 1
        counts["generic"] += 2
        counts["orientation_new"] += decision.orientation_fusion_new_failures
        counts["old_subsets"] += decision.old_subsets_checked
        counts["fused_subsets"] += decision.fused_subsets_checked
        counts["old_passed_fused_failed"] += decision.old_passed_fused_failed
        if not decision.passed:
            counts["failed"] += 1
            failures.append(decision.detail)
    return not failures, failures, counts


def nonbipartite_json(decision) -> dict:
    return asdict(decision)


def exhaustive_seed_certificate(
    adj: Sequence[int], instance: K6LorentzInstance
) -> dict:
    inertia = InertiaCache()
    zero_forcing = ZeroForcingSolver()
    choices = []
    histogram: dict[str, int] = {}
    representatives: dict[str, list[int]] = {}
    for z0 in z0_subsets(instance.eligible_z0_mask):
        z0_absolute = absolute_vertices(instance, z0)
        matchable = support_matching(z0, instance.defects) is not None
        if matchable:
            local_components = lorentz_components(instance, z0)
            fused_passed, fused_failures, _ = fused_system(
                adj,
                instance,
                z0,
                local_components,
                inertia,
                zero_forcing,
            )
            nonbipartite = check_nonbipartite_system(
                instance, z0, local_components
            )
        else:
            fused_passed = False
            fused_failures = []
            nonbipartite = None
        nonbipartite_passed = bool(
            nonbipartite is not None and nonbipartite.passed
        )
        if fused_passed and nonbipartite_passed:
            category = "fused+nonbipartite"
        elif fused_passed:
            category = "fused_only"
        elif nonbipartite_passed:
            category = "nonbipartite_only"
        else:
            category = "neither"
        histogram[category] = histogram.get(category, 0) + 1
        representatives.setdefault(category, z0_absolute)
        choices.append({
            "Z0": z0_absolute,
            "Z0_allowed_mask_matchable": matchable,
            "fused_side_rank_passed": fused_passed,
            "nonbipartite_joint_support_passed": nonbipartite_passed,
            "category": category,
            "failed_fused_components": fused_failures,
            "nonbipartite": (
                nonbipartite_json(nonbipartite)
                if nonbipartite is not None
                else None
            ),
        })
    return {
        "seed": list(instance.seed),
        "outside": list(instance.outside),
        "eligible_Z0": absolute_vertices(instance, instance.eligible_z0_mask),
        "allowed_defects": {
            str(instance.outside[local]): [
                instance.seed[coordinate] for coordinate in vertices(mask)
            ]
            for local, mask in enumerate(instance.defects)
        },
        "histogram": histogram,
        "representative_Z0": representatives,
        "common_pass_count": histogram.get("fused+nonbipartite", 0),
        "choices": choices,
    }


@dataclass
class FusedSeedDecision:
    feasible: bool
    counts: dict[str, int]
    chosen_z0: list[int] | None
    certificate: dict | None


def solve_seed(
    adj: Sequence[int],
    instance: K6LorentzInstance,
    inertia: InertiaCache,
    zero_forcing: ZeroForcingSolver,
) -> FusedSeedDecision:
    counts = empty_counts()
    for z0 in z0_subsets(instance.eligible_z0_mask):
        counts["z0_considered"] += 1
        if support_matching(z0, instance.defects) is None:
            continue
        counts["z0_matchable"] += 1
        local_components = lorentz_components(instance, z0)
        fused_passed, fused_failures, local = fused_system(
            adj, instance, z0, local_components, inertia, zero_forcing
        )
        counts["bipartite_components_checked"] += local["components"]
        counts["bipartite_components_failed"] += local["failed"]
        counts["components_old_separate_passed_fused_failed"] += local[
            "old_passed_fused_failed"
        ]
        counts["generic_orientation_cases_checked"] += local["generic"]
        counts["orientation_fusion_new_failures"] += local["orientation_new"]
        counts["old_block_subsets_checked"] += local["old_subsets"]
        counts["fused_block_subsets_checked"] += local["fused_subsets"]
        if not fused_passed:
            counts["z0_fused_failed"] += 1
            continue
        counts["z0_fused_passed"] += 1
        nonbipartite = check_nonbipartite_system(instance, z0, local_components)
        counts["pure_colorings_considered"] += nonbipartite.raw_colorings
        counts["pure_colorings_matchable"] += nonbipartite.matchable_colorings
        counts["joint_support_searches"] += nonbipartite.support_searches
        counts["joint_support_dfs_nodes"] += nonbipartite.support_dfs_nodes
        if nonbipartite.passed:
            return FusedSeedDecision(
                True, counts, absolute_vertices(instance, z0), None
            )
        counts["z0_nonbipartite_failed"] += 1

    certificate = exhaustive_seed_certificate(adj, instance)
    if certificate["common_pass_count"]:
        raise AssertionError("fused main search and certificate disagree")
    return FusedSeedDecision(False, counts, None, certificate)


def evaluate_graph(adj: Sequence[int], stop_on_rejection: bool = True) -> dict:
    validate_graph(adj, require_alpha_two=True)
    if find_clique_mask(adj, 7):
        raise ValueError("fused side-rank production target must be K6-only")
    inertia = InertiaCache()
    zero_forcing = ZeroForcingSolver()
    totals = empty_counts()
    seeds_checked = 0
    impossible_seeds = 0
    first_seed = None
    first_certificate = None
    for seed_mask in clique_masks(adj, COORDINATES):
        seeds_checked += 1
        instance = build_instance(adj, vertices(seed_mask))
        decision = solve_seed(adj, instance, inertia, zero_forcing)
        add_counts(totals, decision.counts)
        if not decision.feasible:
            impossible_seeds += 1
            if first_seed is None:
                first_seed = list(instance.seed)
                first_certificate = decision.certificate
            if stop_on_rejection:
                break
    return {
        "rejected": impossible_seeds > 0,
        "seeds_checked": seeds_checked,
        "impossible_seeds": impossible_seeds,
        **totals,
        "inertia_cache_entries": len(inertia.values),
        "inertia_cache_hits": inertia.hits,
        "zero_forcing_cache_entries": len(zero_forcing.cache),
        "zero_forcing_cache_hits": zero_forcing.cache_hits,
        "zero_forcing_initial_sets_checked": zero_forcing.initial_sets_checked,
        "first_impossible_seed": first_seed,
        "certificate": first_certificate,
    }


def evaluate_record(record: dict) -> dict:
    return {
        "index": record["index"],
        "decision": evaluate_graph(record["adjacency"], stop_on_rejection=True),
    }


def verify_inputs() -> tuple[dict[str, str], list[dict], list[int]]:
    observed = {name: sha256(ROOT / name) for name in EXPECTED}
    if observed != EXPECTED:
        raise ValueError(f"fused side-rank dependency hash mismatch: {observed}")
    prior = json.loads(
        (ROOT / "d6_k6_all_same_z0_report.json").read_text(encoding="utf-8")
    )
    if prior.get("rejected_indices") != EXPECTED_PRIOR_REJECTIONS:
        raise ValueError("prior four-index all-same-Z0 result changed")
    manifest = json.loads(
        (ROOT / "d6_current_residue_manifest.json").read_text(encoding="utf-8")
    )
    indices = manifest["classes"]["K6_only"]["final_indices"]
    if len(indices) != EXPECTED_INPUT or stable_hash(indices) != EXPECTED_INDICES_SHA256:
        raise ValueError("fused side-rank K6 residue selection changed")
    payload = json.loads(
        (ROOT / "d6_k6_bipartite_rank_input.json").read_text(encoding="utf-8")
    )
    by_index = {record["index"]: record for record in payload["graphs"]}
    if any(index not in by_index for index in indices):
        raise ValueError("fused selected index absent from adjacency input")
    return observed, [by_index[index] for index in indices], indices


def run(report_path: Path, certificate_path: Path, workers: int) -> dict:
    observed, records, indices = verify_inputs()
    started = time.perf_counter()
    if workers == 1:
        raw_results = [evaluate_record(record) for record in records]
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            raw_results = list(pool.map(evaluate_record, records, chunksize=1))
    wall = time.perf_counter() - started
    rejected = [item for item in raw_results if item["decision"]["rejected"]]

    certificate_archive = {
        "schema": CERTIFICATE_SCHEMA,
        "status": "COMPLETE",
        "input_indices_sha256": stable_hash(indices),
        "production_source_sha256": sha256(Path(__file__).resolve()),
        "rejected_indices": [item["index"] for item in rejected],
        "certificates": [
            {
                "index": item["index"],
                "certificate": item["decision"]["certificate"],
            }
            for item in rejected
        ],
    }
    atomic_json(certificate_path, certificate_archive)

    graph_results = []
    for item in raw_results:
        decision = item["decision"].copy()
        decision.pop("certificate")
        graph_results.append({"index": item["index"], "decision": decision})
    aggregate_fields = (
        "seeds_checked",
        "impossible_seeds",
        *COUNTER_FIELDS,
        "inertia_cache_entries",
        "inertia_cache_hits",
        "zero_forcing_cache_entries",
        "zero_forcing_cache_hits",
        "zero_forcing_initial_sets_checked",
    )
    totals = {
        name: sum(item["decision"][name] for item in raw_results)
        for name in aggregate_fields
    }

    positive = evaluate_graph(lower_bound_18_graph(), stop_on_rejection=True)
    if positive["rejected"] or positive["seeds_checked"] != 32:
        raise AssertionError("known realizable 18-point positive control failed")
    rejected_indices = [item["index"] for item in rejected]
    if not set(EXPECTED_PRIOR_REJECTIONS).issubset(rejected_indices):
        raise AssertionError("fused result lost a prior exact rejection")
    report = {
        "schema": REPORT_SCHEMA,
        "status": "COMPLETE",
        "description": (
            "Exact sidewise max fusion of zero-forcing and orientation inertia "
            "inside block support, coupled to nonbipartite support on one Z0."
        ),
        "input_graphs": len(graph_results),
        "input_indices_sha256": stable_hash(indices),
        "graphs_rejected": len(rejected),
        "graphs_surviving": len(graph_results) - len(rejected),
        "rejected_indices": rejected_indices,
        "prior_all_same_Z0_rejected_indices": EXPECTED_PRIOR_REJECTIONS,
        "new_rejected_indices": sorted(
            set(rejected_indices) - set(EXPECTED_PRIOR_REJECTIONS)
        ),
        **totals,
        "graph_results": graph_results,
        "certificate_archive": {
            "path": certificate_path.name,
            "sha256": sha256(certificate_path),
            "bytes": certificate_path.stat().st_size,
            "certificates": len(rejected),
        },
        "positive_control": {
            "name": "known realizable 18-point construction",
            "passed": not positive["rejected"],
            "K6_seeds": positive["seeds_checked"],
            "impossible_seeds": positive["impossible_seeds"],
        },
        "sources": {
            **observed,
            Path(__file__).name: sha256(Path(__file__).resolve()),
        },
        "runtime": {
            "command": " ".join(sys.argv),
            "workers": workers,
            "wall_seconds": wall,
            "python": sys.version,
            "platform": platform.platform(),
        },
        "quantifiers": {
            "same_spans": (
                "for each generic orientation, ZF and inertia lower-bound the "
                "same U_A and U_B, so each side uses their maximum"
            ),
            "component": (
                "every bipartite component passes a fused generic orientation "
                "or the exact lightlike alternative"
            ),
            "Z0": (
                "some one eligible Z0 passes both fused bipartite and "
                "nonbipartite joint-support systems"
            ),
            "seed": "one impossible required K6 seed rejects a graph",
        },
        "semantics": {
            "candidate_nonedges": "unconstrained and may be unit",
            "allowed_defects": "upper bounds; allowed coordinates may be zero",
            "survivor": "filter non-rejection only, not realizability",
            "arithmetic": "exact integer/rational decisions; no tolerance",
        },
    }
    atomic_json(report_path, report)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path, default=ROOT / "d6_k6_fused_side_rank_report.json"
    )
    parser.add_argument(
        "--certificates",
        type=Path,
        default=ROOT / "d6_k6_fused_side_rank_certificates.json",
    )
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("workers must be positive")
    report = run(
        args.output.resolve(), args.certificates.resolve(), args.workers
    )
    print(json.dumps({
        "status": report["status"],
        "input_graphs": report["input_graphs"],
        "graphs_rejected": report["graphs_rejected"],
        "rejected_indices": report["rejected_indices"],
        "new_rejected_indices": report["new_rejected_indices"],
        "output": str(args.output.resolve()),
        "certificates": str(args.certificates.resolve()),
    }, sort_keys=True))


if __name__ == "__main__":
    main()
