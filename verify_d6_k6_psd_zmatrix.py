#!/usr/bin/env python3
"""Independent verifier for the K6 PSD Z-matrix exact layer.

This checker imports neither the new production evaluator nor its kernel.  It
uses the frozen independent K6 reconstruction, exact SymPy characteristic
polynomials with Sturm root counts, an independent simultaneous ordinary
zero-forcing solver, fresh side-graph components and bipartitions, and a
fresh implementation of the Perron--Frobenius rank bound and block subsets.
It recomputes every one of the 977 graph decisions and every archived Z0 row.
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

import verify_d6_k6_fused_side_rank as prior
import verify_d6_k6_all_same_z0 as zf_independent
from test_d6_k6_lorentz import lower_bound_18_graph


independent = prior.independent
ROOT = Path(__file__).resolve().parent
COORDINATES = 6
REPORT_SCHEMA = "d6-k6-psd-zmatrix-v1"
CERTIFICATE_SCHEMA = "d6-k6-psd-zmatrix-certificates-v1"
VERIFICATION_SCHEMA = "d6-k6-psd-zmatrix-verification-v1"
EXPECTED_INPUT = 977
EXPECTED_INPUT_SHA256 = (
    "27e435506ecaecb936acaef63e5873c3538771be17d1a52f12df82df8dedd947"
)
EXPECTED = {
    "d6_k6_fused_side_rank.py": (
        "5fc6416feca1898f38b1bb3c2f06a1289fcdd6f179654a4ab746e6f3cff869c2"
    ),
    "d6_k6_fused_side_rank_report.json": (
        "60911844102dfeb494170f6b7aa5d5142b8f545fcc968e43078eb07c7ed102ec"
    ),
    "verify_d6_k6_fused_side_rank.py": (
        "6621f9a18e6ad2ca19c93776040c837ff37e6021be6e05d2f42f7dcf616f09ff"
    ),
    "d6_k6_bipartite_rank_input.json": (
        "dfee979e6ca2759e80cf803c3054c8954d98543875d467644b91ee434c22e845"
    ),
    "test_d6_k6_lorentz.py": (
        "f4fc8544a4311f742267808e6c785f82139accb9da99b496d206283b52831a47"
    ),
    "d6_k6_psd_zmatrix.py": (
        "24c2b3abfce5074560183d17f77778f96298bef6e38ece5bc937ccb2a4cba736"
    ),
    "d6_k6_psd_zmatrix_report.json": (
        "abb920fbd9eb848502cca37e32c2aa86bf220ab1773e6006ea1fbee4bc2f4e30"
    ),
    "d6_k6_psd_zmatrix_certificates.json": (
        "62fc7170d97056941bcaf553c45cfa1cabdf20117ffc69fbd8429b3764ec999c"
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


def graph_components(rows: Sequence[int]) -> tuple[tuple[int, ...], ...]:
    unseen = set(range(len(rows)))
    answer = []
    while unseen:
        start = min(unseen)
        unseen.remove(start)
        reached = {start}
        frontier = [start]
        while frontier:
            vertex = frontier.pop()
            neighbors = [
                other
                for other in tuple(unseen)
                if rows[vertex] & (1 << other)
            ]
            for other in neighbors:
                unseen.remove(other)
                reached.add(other)
                frontier.append(other)
        answer.append(tuple(sorted(reached)))
    return tuple(answer)


def induced_rows(
    rows: Sequence[int], selected: Sequence[int]
) -> tuple[int, ...]:
    return tuple(
        sum(
            1 << column
            for column, other in enumerate(selected)
            if rows[vertex] & (1 << other)
        )
        for vertex in selected
    )


def is_bipartite(rows: Sequence[int]) -> bool:
    color: list[int | None] = [None] * len(rows)
    for root in range(len(rows)):
        if color[root] is not None:
            continue
        color[root] = 0
        frontier = [root]
        while frontier:
            vertex = frontier.pop(0)
            for other in independent.bits(rows[vertex]):
                wanted = 1 - int(color[vertex])
                if color[other] is None:
                    color[other] = wanted
                    frontier.append(other)
                elif color[other] != wanted:
                    return False
    return True


def pf_nullity_upper(rows: Sequence[int], lorentz_sign: str) -> int | None:
    """Fresh reconstruction of the irreducible PSD Z-matrix conclusion."""

    if len(graph_components(rows)) != 1:
        raise ValueError("PF bound must be called on one connected block")
    if lorentz_sign == "negative":
        return 1
    if lorentz_sign == "positive" and is_bipartite(rows):
        return 1
    if lorentz_sign != "positive":
        raise ValueError("unknown Lorentz sign")
    return None


def side_rank_lower(
    graph_rows: Sequence[int],
    absolute: Sequence[int],
    lorentz_sign: str,
    inertia: independent.SympyInertiaCache,
    zero_forcing: zf_independent.IndependentZeroForcing,
) -> dict:
    total = 0
    total_zf = 0
    total_inertia = 0
    details = []
    for selected in graph_components(graph_rows):
        local = induced_rows(graph_rows, selected)
        size = len(local)
        pattern = tuple(row | (1 << index) for index, row in enumerate(local))
        positive, negative, zero = inertia.solve(pattern)
        zf_number = zero_forcing.number(local)
        zf_rank = size - zf_number
        inertia_rank = (
            size - negative if lorentz_sign == "positive" else size - positive
        )
        pf_nullity = pf_nullity_upper(local, lorentz_sign)
        zmatrix_applicable = pf_nullity is not None
        zmatrix_rank = size - pf_nullity if pf_nullity is not None else 0
        fused = max(zf_rank, inertia_rank, zmatrix_rank)
        total += fused
        total_zf += zf_rank
        total_inertia += inertia_rank
        details.append({
            "vertices": [absolute[local_vertex] for local_vertex in selected],
            "adjacency_rows": list(local),
            "size": size,
            "bipartite": is_bipartite(local),
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
    return {
        "rank_lower": total,
        "prior_fused_rank_lower": max(total_zf, total_inertia),
        "zero_forcing_rank_lower": total_zf,
        "inertia_rank_lower": total_inertia,
        "components": details,
    }


@dataclass(frozen=True)
class BlockResult:
    passed: bool
    first_failure: dict | None


def block_result(blocks: Sequence[tuple[str, int, int]]) -> BlockResult:
    for selected in range(1, 1 << len(blocks)):
        rank = 0
        allowed = 0
        labels = []
        for position, (label, lower, mask) in enumerate(blocks):
            if selected & (1 << position):
                labels.append(label)
                rank += lower
                allowed |= mask
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
        (f"Z0:{instance.outside[local]}", 1, instance.defects[local])
        for local in independent.bits(z0)
    )
    return tuple(blocks)


def check_component(
    adj: Sequence[int],
    instance: independent.Instance,
    z0: int,
    component: independent.Component,
    inertia: independent.SympyInertiaCache,
    zero_forcing: zf_independent.IndependentZeroForcing,
) -> tuple[bool, dict]:
    graph_a, absolute_a = zf_independent.induced_adjacency(
        adj, instance, component.side_a
    )
    graph_b, absolute_b = zf_independent.induced_adjacency(
        adj, instance, component.side_b
    )
    orientations = []
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
        old = block_result(make_blocks(
            instance,
            z0,
            component,
            side_a["prior_fused_rank_lower"],
            side_b["prior_fused_rank_lower"],
        ))
        new = block_result(make_blocks(
            instance,
            z0,
            component,
            side_a["rank_lower"],
            side_b["rank_lower"],
        ))
        orientations.append({
            "case": name,
            "side_A_sign": sign_a,
            "side_B_sign": sign_b,
            "side_A": side_a,
            "side_B": side_b,
            "prior_fused_support_passed": old.passed,
            "prior_fused_first_failure": old.first_failure,
            "psd_zmatrix_support_passed": new.passed,
            "psd_zmatrix_first_failure": new.first_failure,
            "prior_fused_passed_psd_zmatrix_failed": (
                old.passed and not new.passed
            ),
        })

    light = z0 | component.vertices
    light_dimension = light.bit_count()
    light_matching = (
        light_dimension <= COORDINATES
        and independent.hall_matchable(light, instance.defects)
    )
    lightlike = {
        "orthonormal_vectors": light_dimension,
        "dimension_passed": light_dimension <= COORDINATES,
        "allowed_mask_matching_passed": light_matching,
        "passed": light_matching,
    }
    prior_passed = bool(
        any(item["prior_fused_support_passed"] for item in orientations)
        or light_matching
    )
    passed = bool(
        any(item["psd_zmatrix_support_passed"] for item in orientations)
        or light_matching
    )
    detail = {
        "component": [
            instance.outside[local]
            for local in independent.bits(component.vertices)
        ],
        "side_A": list(absolute_a),
        "side_B": list(absolute_b),
        "orientations": orientations,
        "lightlike": lightlike,
        "prior_fused_passed": prior_passed,
        "psd_zmatrix_passed": passed,
    }
    return passed, detail


def psd_system(
    adj: Sequence[int],
    instance: independent.Instance,
    z0: int,
    components: Sequence[independent.Component],
    inertia: independent.SympyInertiaCache,
    zero_forcing: zf_independent.IndependentZeroForcing,
) -> tuple[bool, list[dict]]:
    decisions = [
        check_component(adj, instance, z0, component, inertia, zero_forcing)
        for component in components
        if component.bipartite
    ]
    return all(item[0] for item in decisions), [
        item[1] for item in decisions if not item[0]
    ]


@dataclass
class SeedDecision:
    feasible: bool
    z0_considered: int
    z0_matchable: int
    z0_psd_zmatrix_passed: int
    z0_psd_zmatrix_failed: int
    z0_nonbipartite_failed: int


def solve_seed(
    adj: Sequence[int],
    instance: independent.Instance,
    inertia: independent.SympyInertiaCache,
    zero_forcing: zf_independent.IndependentZeroForcing,
) -> SeedDecision:
    counts = {
        "z0_considered": 0,
        "z0_matchable": 0,
        "z0_psd_zmatrix_passed": 0,
        "z0_psd_zmatrix_failed": 0,
        "z0_nonbipartite_failed": 0,
    }
    for z0 in independent.z0_subsets(instance.eligible_z0):
        counts["z0_considered"] += 1
        if not independent.hall_matchable(z0, instance.defects):
            continue
        counts["z0_matchable"] += 1
        local_components = independent.components(instance, z0)
        passed, _ = psd_system(
            adj, instance, z0, local_components, inertia, zero_forcing
        )
        if not passed:
            counts["z0_psd_zmatrix_failed"] += 1
            continue
        counts["z0_psd_zmatrix_passed"] += 1
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
    "z0_psd_zmatrix_passed",
    "z0_psd_zmatrix_failed",
    "z0_nonbipartite_failed",
)


def evaluate_graph(adj: Sequence[int]) -> dict:
    independent.validate_graph(adj)
    if next(independent.clique_masks(adj, 7), 0):
        raise ValueError("independent PSD Z-matrix target unexpectedly contains K7")
    inertia = independent.SympyInertiaCache()
    zero_forcing = zf_independent.IndependentZeroForcing()
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


def exhaustive_seed_certificate(
    adj: Sequence[int], instance: independent.Instance
) -> dict:
    inertia = independent.SympyInertiaCache()
    zero_forcing = zf_independent.IndependentZeroForcing()
    choices = []
    histogram: dict[str, int] = {}
    representatives: dict[str, list[int]] = {}
    for z0 in independent.z0_subsets(instance.eligible_z0):
        z0_absolute = [
            instance.outside[local] for local in independent.bits(z0)
        ]
        matchable = independent.hall_matchable(z0, instance.defects)
        if matchable:
            local_components = independent.components(instance, z0)
            psd_passed, failures = psd_system(
                adj, instance, z0, local_components, inertia, zero_forcing
            )
            odd = [
                component.vertices
                for component in local_components
                if not component.bipartite
            ]
            odd_absolute = [
                [
                    instance.outside[local]
                    for local in independent.bits(component)
                ]
                for component in odd
            ]
            nonbipartite_passed = independent.nonbipartite_system_passes(
                instance, z0, odd
            )
        else:
            psd_passed = False
            failures = []
            odd_absolute = None
            nonbipartite_passed = False
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
            "nonbipartite_odd_components": odd_absolute,
        })
    return {
        "seed": list(instance.seed),
        "outside": list(instance.outside),
        "eligible_Z0": [
            instance.outside[local]
            for local in independent.bits(instance.eligible_z0)
        ],
        "allowed_defects": {
            str(instance.outside[local]): [
                instance.seed[coordinate]
                for coordinate in independent.bits(mask)
            ]
            for local, mask in enumerate(instance.defects)
        },
        "histogram": histogram,
        "representative_Z0": representatives,
        "common_pass_count": histogram.get("psd_zmatrix+nonbipartite", 0),
        "choices": choices,
    }


def kernel_controls() -> dict:
    cases = {
        "K1": (0,),
        "K2": (2, 1),
        "P3": (2, 5, 2),
        "C5": (18, 5, 10, 20, 9),
    }
    expected = {
        "K1": {"positive": 1, "negative": 0},
        "K2": {"positive": 2, "negative": 1},
        "P3": {"positive": 2, "negative": 2},
        "C5": {"positive": 3, "negative": 4},
    }
    observed = {}
    for name, rows in cases.items():
        observed[name] = {}
        for sign in ("positive", "negative"):
            side = side_rank_lower(
                rows,
                tuple(range(len(rows))),
                sign,
                independent.SympyInertiaCache(),
                zf_independent.IndependentZeroForcing(),
            )
            observed[name][sign] = side["rank_lower"]
    if observed != expected:
        raise AssertionError(f"independent small-graph controls differ: {observed}")
    return observed


def verify_internal(report: dict, archive: dict) -> None:
    if report.get("schema") != REPORT_SCHEMA or report.get("status") != "COMPLETE":
        raise ValueError("PSD Z-matrix production report has wrong schema/status")
    if archive.get("schema") != CERTIFICATE_SCHEMA or archive.get("status") != "COMPLETE":
        raise ValueError("PSD Z-matrix archive has wrong schema/status")
    if report["certificate_archive"]["sha256"] != sha256(
        ROOT / report["certificate_archive"]["path"]
    ):
        raise AssertionError("PSD Z-matrix certificate archive hash differs")
    rejected = [
        item["index"]
        for item in report["graph_results"]
        if item["decision"]["rejected"]
    ]
    if rejected != report["rejected_indices"] or rejected != archive["rejected_indices"]:
        raise AssertionError("PSD Z-matrix production rejection lists disagree")
    numeric = tuple(
        name
        for name, value in report["graph_results"][0]["decision"].items()
        if isinstance(value, int) and not isinstance(value, bool)
    )
    for name in numeric:
        observed = sum(item["decision"][name] for item in report["graph_results"])
        if observed != report[name]:
            raise AssertionError(f"PSD Z-matrix aggregate differs for {name}")


def verify(output: Path) -> dict:
    observed_hashes = {name: sha256(ROOT / name) for name in EXPECTED}
    if observed_hashes != EXPECTED:
        raise ValueError(f"PSD Z-matrix verifier dependency mismatch: {observed_hashes}")
    report = json.loads(
        (ROOT / "d6_k6_psd_zmatrix_report.json").read_text(encoding="utf-8")
    )
    archive = json.loads(
        (ROOT / "d6_k6_psd_zmatrix_certificates.json").read_text(
            encoding="utf-8"
        )
    )
    verify_internal(report, archive)
    parent = json.loads(
        (ROOT / "d6_k6_fused_side_rank_report.json").read_text(encoding="utf-8")
    )
    indices = [
        item["index"]
        for item in parent["graph_results"]
        if not item["decision"]["rejected"]
    ]
    if len(indices) != EXPECTED_INPUT or independent.stable_hash(indices) != EXPECTED_INPUT_SHA256:
        raise ValueError("independent PSD Z-matrix input boundary differs")
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

    controls = kernel_controls()
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
        certificate = exhaustive_seed_certificate(
            adjacency, result["first_instance"]
        )
        if certificate != certificates[index]:
            raise AssertionError(f"graph {index} exhaustive certificate differs")
        if certificate["common_pass_count"]:
            raise AssertionError(f"graph {index} certificate contains common pass")
        summaries.append({
            "index": index,
            "first_impossible_seed": result["first_impossible_seed"],
            "Z0_choices": len(certificate["choices"]),
            "histogram": certificate["histogram"],
        })

    if rejected != report["rejected_indices"]:
        raise AssertionError("independent PSD Z-matrix rejection set differs")
    positive = evaluate_graph(lower_bound_18_graph())
    if positive["rejected"] or positive["seeds_checked"] != 32:
        raise AssertionError("known realizable 18-point PSD Z-matrix control failed")

    result = {
        "schema": VERIFICATION_SCHEMA,
        "status": "PASS",
        "method": (
            "fresh graph components and PF nullity reconstruction, independent "
            "SymPy/Sturm inertia and simultaneous ordinary zero forcing, fresh "
            "componentwise fusion/block subsets, and frozen independent "
            "same-Z0 nonbipartite support reconstruction"
        ),
        "graphs_recomputed": len(indices),
        "graphs_rejected": len(rejected),
        "graphs_surviving": len(indices) - len(rejected),
        "rejected_indices": rejected,
        "all_graph_core_fields_matched": list(CORE_FIELDS),
        "all_first_impossible_seeds_matched": True,
        "all_exhaustive_Z0_certificates_matched": True,
        "certificate_summaries": summaries,
        "small_graph_controls": controls,
        "positive_control": {
            "name": "known realizable 18-point construction",
            "passed": not positive["rejected"],
            "K6_seeds": positive["seeds_checked"],
        },
        "inputs": {
            **observed_hashes,
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
        default=ROOT / "d6_k6_psd_zmatrix_verification.json",
    )
    args = parser.parse_args()
    result = verify(args.output.resolve())
    print(json.dumps({
        "status": result["status"],
        "graphs_recomputed": result["graphs_recomputed"],
        "graphs_rejected": result["graphs_rejected"],
        "graphs_surviving": result["graphs_surviving"],
        "rejected_indices": result["rejected_indices"],
        "output": str(args.output.resolve()),
    }, sort_keys=True))


if __name__ == "__main__":
    main()
