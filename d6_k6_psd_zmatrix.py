#!/usr/bin/env python3
"""Exact PSD Z-matrix rank refinement for the K6-only residue.

For a generic Lorentz-component side, the actual Euclidean Gram matrix is

    I + kappa D (I + Adj(F)) D,

where ``D`` is nonsingular.  On every connected component of ``F``, a
negative ``kappa`` gives a congruent irreducible positive-semidefinite
Z-matrix.  Its nullity is at most one by Perron--Frobenius.  A positive
``kappa`` has the same property when that component is bipartite, after a
signature switch.  The resulting component rank lower bound is fused with
ordinary zero forcing and exact inertia on that same component, then fed to
the frozen orthogonal block-support and same-Z0 nonbipartite systems.

The lightlike Lorentz case remains a separate exact alternative.  All proof
decisions are integer or exact rational computations; floating point is used
only for elapsed-time reporting.
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

import d6_k6_fused_side_rank as frozen
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
from d6_k6_normal_block_support import check_block_subsets, z0_subsets
from d6_k6_normal_inertia import InertiaCache
from test_d6_k6_lorentz import lower_bound_18_graph


ROOT = Path(__file__).resolve().parent
REPORT_SCHEMA = "d6-k6-psd-zmatrix-v1"
CERTIFICATE_SCHEMA = "d6-k6-psd-zmatrix-certificates-v1"
EXPECTED_INPUT = 977
EXPECTED_INPUT_SHA256 = (
    "27e435506ecaecb936acaef63e5873c3538771be17d1a52f12df82df8dedd947"
)
EXPECTED_PARENT_REJECTIONS = [
    122_871,
    129_414,
    652_900,
    665_962,
    719_091,
    1_401_532,
    1_948_946,
    2_301_548,
    2_842_523,
    3_274_548,
    3_289_061,
    3_959_774,
    3_962_868,
]
EXPECTED = {
    "d6_k6_fused_side_rank.py": (
        "5fc6416feca1898f38b1bb3c2f06a1289fcdd6f179654a4ab746e6f3cff869c2"
    ),
    "d6_k6_fused_side_rank_report.json": (
        "60911844102dfeb494170f6b7aa5d5142b8f545fcc968e43078eb07c7ed102ec"
    ),
    "d6_k6_fused_side_rank_certificates.json": (
        "1d497df83b329948c2e09702ca8a94ea94d2d959f7a866a532084d9d920ea664"
    ),
    "d6_k6_bipartite_rank_input.json": (
        "dfee979e6ca2759e80cf803c3054c8954d98543875d467644b91ee434c22e845"
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


def graph_components(rows: Sequence[int]) -> tuple[tuple[int, ...], ...]:
    """Connected components of a loopless bit-row graph."""

    unseen = (1 << len(rows)) - 1
    answer = []
    while unseen:
        start = unseen & -unseen
        unseen ^= start
        component = start
        frontier = start
        while frontier:
            bit = frontier & -frontier
            frontier ^= bit
            vertex = bit.bit_length() - 1
            new = rows[vertex] & unseen
            unseen &= ~new
            component |= new
            frontier |= new
        answer.append(tuple(vertices(component)))
    return tuple(answer)


def induced_rows(
    rows: Sequence[int], selected: Sequence[int]
) -> tuple[int, ...]:
    position = {vertex: local for local, vertex in enumerate(selected)}
    return tuple(
        sum(
            1 << position[other]
            for other in selected
            if rows[vertex] & (1 << other)
        )
        for vertex in selected
    )


def is_bipartite(rows: Sequence[int]) -> bool:
    colors: dict[int, int] = {}
    for start in range(len(rows)):
        if start in colors:
            continue
        colors[start] = 0
        stack = [start]
        while stack:
            vertex = stack.pop()
            for other in vertices(rows[vertex]):
                wanted = 1 - colors[vertex]
                if other in colors:
                    if colors[other] != wanted:
                        return False
                else:
                    colors[other] = wanted
                    stack.append(other)
    return True


@dataclass(frozen=True)
class SideRank:
    rank_lower: int
    prior_fused_rank_lower: int
    zero_forcing_rank_lower: int
    inertia_rank_lower: int
    components: tuple[dict, ...]


def side_rank_lower(
    graph_rows: Sequence[int],
    absolute: Sequence[int],
    lorentz_sign: str,
    inertia: InertiaCache,
    zero_forcing: ZeroForcingSolver,
) -> SideRank:
    """Fuse all exact bounds inside each connected side-graph block."""

    if lorentz_sign not in ("positive", "negative"):
        raise ValueError("Lorentz sign must be positive or negative")
    total = 0
    total_zf = 0
    total_inertia = 0
    details = []
    for selected in graph_components(graph_rows):
        local = induced_rows(graph_rows, selected)
        size = len(local)
        pattern = tuple(row | (1 << index) for index, row in enumerate(local))
        positive, negative, zero = inertia.solve(pattern)
        zf_number = zero_forcing.solve(local).number
        zf_rank = size - zf_number
        inertia_rank = (
            size - negative if lorentz_sign == "positive" else size - positive
        )
        bipartite = is_bipartite(local)
        zmatrix_applicable = lorentz_sign == "negative" or bipartite
        zmatrix_rank = size - 1 if zmatrix_applicable else 0
        fused = max(zf_rank, inertia_rank, zmatrix_rank)
        total += fused
        total_zf += zf_rank
        total_inertia += inertia_rank
        details.append({
            "vertices": [absolute[local_vertex] for local_vertex in selected],
            "adjacency_rows": list(local),
            "size": size,
            "bipartite": bipartite,
            "inertia_I_plus_Adj": [positive, negative, zero],
            "ordinary_zero_forcing_number": zf_number,
            "ordinary_zero_forcing_rank_lower": zf_rank,
            "inertia_rank_lower": inertia_rank,
            "psd_zmatrix_applicable": zmatrix_applicable,
            "psd_zmatrix_rank_lower": zmatrix_rank,
            "componentwise_fused_rank_lower": fused,
            "psd_zmatrix_strict_improvement": (
                zmatrix_rank > max(zf_rank, inertia_rank)
            ),
        })
    return SideRank(
        rank_lower=total,
        prior_fused_rank_lower=max(total_zf, total_inertia),
        zero_forcing_rank_lower=total_zf,
        inertia_rank_lower=total_inertia,
        components=tuple(details),
    )


@dataclass(frozen=True)
class PsdComponentDecision:
    passed: bool
    prior_fused_passed: bool
    detail: dict
    old_subsets_checked: int
    new_subsets_checked: int
    side_blocks_checked: int
    zmatrix_blocks_applicable: int
    zmatrix_strict_improvements: int
    strict_orientation_failures: int


def check_psd_component(
    adj: Sequence[int],
    instance: K6LorentzInstance,
    z0: int,
    component,
    inertia: InertiaCache,
    zero_forcing: ZeroForcingSolver,
) -> PsdComponentDecision:
    if not component.bipartite:
        raise ValueError("PSD Z-matrix generic checker requires bipartite L component")

    graph_a, absolute_a = induced_required_graph(adj, instance, component.side_a)
    graph_b, absolute_b = induced_required_graph(adj, instance, component.side_b)
    orientation_details = []
    old_checked = 0
    new_checked = 0
    side_blocks = 0
    applicable = 0
    improvements = 0
    strict_failures = 0
    for name, sign_a, sign_b in (
        ("A_positive", "positive", "negative"),
        ("A_negative", "negative", "positive"),
    ):
        side_a = side_rank_lower(
            graph_a, absolute_a, sign_a, inertia, zero_forcing
        )
        side_b = side_rank_lower(
            graph_b, absolute_b, sign_b, inertia, zero_forcing
        )
        old = check_block_subsets(frozen.make_blocks(
            instance,
            z0,
            component.side_a,
            component.side_b,
            side_a.prior_fused_rank_lower,
            side_b.prior_fused_rank_lower,
        ))
        new = check_block_subsets(frozen.make_blocks(
            instance,
            z0,
            component.side_a,
            component.side_b,
            side_a.rank_lower,
            side_b.rank_lower,
        ))
        all_blocks = side_a.components + side_b.components
        side_blocks += len(all_blocks)
        applicable += sum(item["psd_zmatrix_applicable"] for item in all_blocks)
        improvements += sum(
            item["psd_zmatrix_strict_improvement"] for item in all_blocks
        )
        old_checked += old.subsets_checked
        new_checked += new.subsets_checked
        strict = old.passed and not new.passed
        strict_failures += strict
        orientation_details.append({
            "case": name,
            "side_A_sign": sign_a,
            "side_B_sign": sign_b,
            "side_A": asdict(side_a),
            "side_B": asdict(side_b),
            "prior_fused_support_passed": old.passed,
            "prior_fused_first_failure": old.first_failure,
            "psd_zmatrix_support_passed": new.passed,
            "psd_zmatrix_first_failure": new.first_failure,
            "prior_fused_passed_psd_zmatrix_failed": strict,
        })

    light_selected = z0 | component.component
    light_dimension = light_selected.bit_count()
    light_matching = (
        light_dimension <= COORDINATES
        and support_matching(light_selected, instance.defects) is not None
    )
    lightlike = {
        "orthonormal_vectors": light_dimension,
        "dimension_passed": light_dimension <= COORDINATES,
        "allowed_mask_matching_passed": light_matching,
        "passed": light_matching,
    }
    prior_passed = bool(
        any(item["prior_fused_support_passed"] for item in orientation_details)
        or light_matching
    )
    passed = bool(
        any(item["psd_zmatrix_support_passed"] for item in orientation_details)
        or light_matching
    )
    if prior_passed is False and passed is True:
        raise AssertionError("PSD Z-matrix refinement cannot revive a component")
    detail = {
        "component": frozen.absolute_vertices(instance, component.component),
        "side_A": list(absolute_a),
        "side_B": list(absolute_b),
        "orientations": orientation_details,
        "lightlike": lightlike,
        "prior_fused_passed": prior_passed,
        "psd_zmatrix_passed": passed,
    }
    return PsdComponentDecision(
        passed=passed,
        prior_fused_passed=prior_passed,
        detail=detail,
        old_subsets_checked=old_checked,
        new_subsets_checked=new_checked,
        side_blocks_checked=side_blocks,
        zmatrix_blocks_applicable=applicable,
        zmatrix_strict_improvements=improvements,
        strict_orientation_failures=strict_failures,
    )


COUNTER_FIELDS = (
    "z0_considered",
    "z0_matchable",
    "z0_psd_zmatrix_passed",
    "z0_psd_zmatrix_failed",
    "z0_nonbipartite_failed",
    "bipartite_components_checked",
    "bipartite_components_failed",
    "prior_fused_components_newly_failed",
    "generic_orientation_cases_checked",
    "prior_fused_orientations_newly_failed",
    "side_graph_connected_blocks_checked",
    "psd_zmatrix_blocks_applicable",
    "psd_zmatrix_strict_component_rank_improvements",
    "prior_fused_block_subsets_checked",
    "psd_zmatrix_block_subsets_checked",
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


def psd_system(
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
        "prior_new": 0,
        "orientations": 0,
        "orientation_new": 0,
        "side_blocks": 0,
        "applicable": 0,
        "improvements": 0,
        "old_subsets": 0,
        "new_subsets": 0,
    }
    for component in components:
        if not component.bipartite:
            continue
        decision = check_psd_component(
            adj, instance, z0, component, inertia, zero_forcing
        )
        counts["components"] += 1
        counts["orientations"] += 2
        counts["orientation_new"] += decision.strict_orientation_failures
        counts["side_blocks"] += decision.side_blocks_checked
        counts["applicable"] += decision.zmatrix_blocks_applicable
        counts["improvements"] += decision.zmatrix_strict_improvements
        counts["old_subsets"] += decision.old_subsets_checked
        counts["new_subsets"] += decision.new_subsets_checked
        counts["prior_new"] += decision.prior_fused_passed and not decision.passed
        if not decision.passed:
            counts["failed"] += 1
            failures.append(decision.detail)
    return not failures, failures, counts


def exhaustive_seed_certificate(
    adj: Sequence[int], instance: K6LorentzInstance
) -> dict:
    inertia = InertiaCache()
    zero_forcing = ZeroForcingSolver()
    choices = []
    histogram: dict[str, int] = {}
    representatives: dict[str, list[int]] = {}
    for z0 in z0_subsets(instance.eligible_z0_mask):
        z0_absolute = frozen.absolute_vertices(instance, z0)
        matchable = support_matching(z0, instance.defects) is not None
        if matchable:
            local_components = lorentz_components(instance, z0)
            psd_passed, failures, _ = psd_system(
                adj, instance, z0, local_components, inertia, zero_forcing
            )
            nonbipartite = frozen.check_nonbipartite_system(
                instance, z0, local_components
            )
        else:
            psd_passed = False
            failures = []
            nonbipartite = None
        nonbipartite_passed = bool(
            nonbipartite is not None and nonbipartite.passed
        )
        if psd_passed and nonbipartite_passed:
            category = "psd_zmatrix+nonbipartite"
        elif psd_passed:
            category = "psd_zmatrix_only"
        elif nonbipartite_passed:
            category = "nonbipartite_only"
        else:
            category = "neither"
        histogram[category] = histogram.get(category, 0) + 1
        representatives.setdefault(category, z0_absolute)
        choices.append({
            "Z0": z0_absolute,
            "Z0_allowed_mask_matchable": matchable,
            "psd_zmatrix_side_rank_passed": psd_passed,
            "nonbipartite_joint_support_passed": nonbipartite_passed,
            "category": category,
            "first_failed_psd_zmatrix_component": (
                failures[0] if failures else None
            ),
            "nonbipartite_odd_components": (
                nonbipartite.odd_components
                if nonbipartite is not None
                else None
            ),
        })
    return {
        "seed": list(instance.seed),
        "outside": list(instance.outside),
        "eligible_Z0": frozen.absolute_vertices(
            instance, instance.eligible_z0_mask
        ),
        "allowed_defects": {
            str(instance.outside[local]): [
                instance.seed[coordinate] for coordinate in vertices(mask)
            ]
            for local, mask in enumerate(instance.defects)
        },
        "histogram": histogram,
        "representative_Z0": representatives,
        "common_pass_count": histogram.get("psd_zmatrix+nonbipartite", 0),
        "choices": choices,
    }


@dataclass
class SeedDecision:
    feasible: bool
    counts: dict[str, int]
    chosen_z0: list[int] | None
    certificate: dict | None


def solve_seed(
    adj: Sequence[int],
    instance: K6LorentzInstance,
    inertia: InertiaCache,
    zero_forcing: ZeroForcingSolver,
) -> SeedDecision:
    counts = empty_counts()
    for z0 in z0_subsets(instance.eligible_z0_mask):
        counts["z0_considered"] += 1
        if support_matching(z0, instance.defects) is None:
            continue
        counts["z0_matchable"] += 1
        local_components = lorentz_components(instance, z0)
        psd_passed, _, local = psd_system(
            adj, instance, z0, local_components, inertia, zero_forcing
        )
        counts["bipartite_components_checked"] += local["components"]
        counts["bipartite_components_failed"] += local["failed"]
        counts["prior_fused_components_newly_failed"] += local["prior_new"]
        counts["generic_orientation_cases_checked"] += local["orientations"]
        counts["prior_fused_orientations_newly_failed"] += local[
            "orientation_new"
        ]
        counts["side_graph_connected_blocks_checked"] += local["side_blocks"]
        counts["psd_zmatrix_blocks_applicable"] += local["applicable"]
        counts["psd_zmatrix_strict_component_rank_improvements"] += local[
            "improvements"
        ]
        counts["prior_fused_block_subsets_checked"] += local["old_subsets"]
        counts["psd_zmatrix_block_subsets_checked"] += local["new_subsets"]
        if not psd_passed:
            counts["z0_psd_zmatrix_failed"] += 1
            continue
        counts["z0_psd_zmatrix_passed"] += 1
        nonbipartite = frozen.check_nonbipartite_system(
            instance, z0, local_components
        )
        counts["pure_colorings_considered"] += nonbipartite.raw_colorings
        counts["pure_colorings_matchable"] += nonbipartite.matchable_colorings
        counts["joint_support_searches"] += nonbipartite.support_searches
        counts["joint_support_dfs_nodes"] += nonbipartite.support_dfs_nodes
        if nonbipartite.passed:
            return SeedDecision(
                True, counts, frozen.absolute_vertices(instance, z0), None
            )
        counts["z0_nonbipartite_failed"] += 1

    certificate = exhaustive_seed_certificate(adj, instance)
    if certificate["common_pass_count"]:
        raise AssertionError("main search and exhaustive certificate disagree")
    return SeedDecision(False, counts, None, certificate)


def evaluate_graph(adj: Sequence[int], stop_on_rejection: bool = True) -> dict:
    validate_graph(adj, require_alpha_two=True)
    if find_clique_mask(adj, 7):
        raise ValueError("PSD Z-matrix production target must be K6-only")
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
        raise ValueError(f"PSD Z-matrix dependency hash mismatch: {observed}")
    _, records, indices = frozen.verify_inputs()
    parent = json.loads(
        (ROOT / "d6_k6_fused_side_rank_report.json").read_text(encoding="utf-8")
    )
    if parent.get("rejected_indices") != EXPECTED_PARENT_REJECTIONS:
        raise ValueError("parent fused rejection set changed")
    survivor_indices = [
        item["index"]
        for item in parent["graph_results"]
        if not item["decision"]["rejected"]
    ]
    if (
        len(survivor_indices) != EXPECTED_INPUT
        or stable_hash(survivor_indices) != EXPECTED_INPUT_SHA256
    ):
        raise ValueError("parent fused survivor boundary changed")
    by_index = {record["index"]: record for record in records}
    if any(index not in by_index for index in survivor_indices):
        raise ValueError("fused survivor absent from pinned adjacency input")
    return observed, [by_index[index] for index in survivor_indices], survivor_indices


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
    report = {
        "schema": REPORT_SCHEMA,
        "status": "COMPLETE",
        "description": (
            "Exact componentwise fusion of PSD Z-matrix/Perron-Frobenius, "
            "ordinary zero-forcing, and inertia side-rank bounds inside "
            "block support, coupled to nonbipartite support on one Z0."
        ),
        "input_graphs": len(graph_results),
        "input_indices_sha256": stable_hash(indices),
        "parent_fused_rejected_indices": EXPECTED_PARENT_REJECTIONS,
        "graphs_rejected": len(rejected),
        "graphs_surviving": len(graph_results) - len(rejected),
        "rejected_indices": rejected_indices,
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
            "side_graph_component": (
                "ordinary zero forcing, exact inertia, and every applicable "
                "PSD Z-matrix n-1 bound are fused on the same connected F block"
            ),
            "lorentz_component": (
                "each bipartite L-Z0 component passes a generic orientation "
                "with componentwise fused side ranks or the lightlike alternative"
            ),
            "Z0": (
                "some eligible Z0 passes both PSD Z-matrix side-rank and "
                "nonbipartite joint-support systems"
            ),
            "seed": "one impossible required K6 seed rejects a graph",
        },
        "semantics": {
            "candidate_nonedges": "unconstrained and may be unit",
            "allowed_defects": "upper bounds; allowed coordinates may be zero",
            "lightlike": "checked separately, never subjected to generic signs",
            "survivor": "filter non-rejection only, not realizability",
            "arithmetic": "exact integer/rational decisions; no tolerance",
        },
    }
    atomic_json(report_path, report)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path, default=ROOT / "d6_k6_psd_zmatrix_report.json"
    )
    parser.add_argument(
        "--certificates",
        type=Path,
        default=ROOT / "d6_k6_psd_zmatrix_certificates.json",
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
        "graphs_surviving": report["graphs_surviving"],
        "rejected_indices": report["rejected_indices"],
        "output": str(args.output.resolve()),
        "certificates": str(args.certificates.resolve()),
    }, sort_keys=True))


if __name__ == "__main__":
    main()
