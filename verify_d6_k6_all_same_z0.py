#!/usr/bin/env python3
"""Independent verifier for the all-K6-systems same-Z0 layer.

The previously frozen independent K6 same-Z0 checker supplies its separate
block-Hall and nonbipartite actual-support reconstructions.  This verifier
adds a new, independently written ordinary zero-forcing solver, rebuilds the
three-way same-Z0 quantifier, and checks every graph and every archived Z0
certificate.  It imports no production evaluator or K6 production engine.
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
from itertools import combinations
from pathlib import Path
from typing import Sequence

import verify_d6_k6_same_z0 as base
from test_d6_k6_lorentz import lower_bound_18_graph


ROOT = Path(__file__).resolve().parent
COORDINATES = 6
REPORT_SCHEMA = "d6-k6-all-same-z0-v1"
CERTIFICATE_SCHEMA = "d6-k6-all-same-z0-certificates-v1"
VERIFICATION_SCHEMA = "d6-k6-all-same-z0-verification-v1"
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
    "d6_k6_all_same_z0.py": (
        "d15fe376993393a5106e5d01f8c0b2316f43945c5fd9ee01b1ebdabfb4fe2b94"
    ),
    "d6_k6_all_same_z0_report.json": (
        "8fa616455d930a1ac3c2255ecb123d8095f7731780f6b6d0bec1ba71c173edeb"
    ),
    "d6_k6_all_same_z0_certificates.json": (
        "bce608e936cc6826962ecae8ea26baf3b9162072a861f568fc15570350e3f401"
    ),
    "verify_d6_k6_same_z0.py": (
        "4db05ff06910a11395e4bbdd2fd13161d6c19f58f81cc99851b8c37558d7000f"
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


class IndependentZeroForcing:
    """Exact ordinary zero forcing by subset enumeration and parallel closure."""

    def __init__(self) -> None:
        self.cache: dict[tuple[int, ...], int] = {}

    @staticmethod
    def closure(adj: Sequence[int], initial: int) -> int:
        black = initial
        while True:
            forced = 0
            for vertex, neighbors in enumerate(adj):
                if not black & (1 << vertex):
                    continue
                white = neighbors & ~black
                if white and not white & (white - 1):
                    forced |= white
            new = forced & ~black
            if not new:
                return black
            black |= new

    def number(self, adj: Sequence[int]) -> int:
        key = tuple(adj)
        if key in self.cache:
            return self.cache[key]
        size = len(key)
        full = (1 << size) - 1
        for count in range(size + 1):
            for chosen in combinations(range(size), count):
                initial = sum(1 << vertex for vertex in chosen)
                if self.closure(key, initial) == full:
                    self.cache[key] = count
                    return count
        raise AssertionError("full vertex set is always zero forcing")


def induced_adjacency(
    adj: Sequence[int], instance: base.Instance, selected: int
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    absolute = tuple(instance.outside[local] for local in base.bits(selected))
    rows = tuple(
        sum(
            1 << column
            for column, other in enumerate(absolute)
            if adj[vertex] & (1 << other)
        )
        for vertex in absolute
    )
    return rows, absolute


def zero_forcing_system(
    adj: Sequence[int],
    instance: base.Instance,
    z0: int,
    local_components: Sequence[base.Component],
    solver: IndependentZeroForcing,
) -> tuple[bool, tuple[dict, ...]]:
    checks = []
    for component in local_components:
        if not component.bipartite:
            continue
        graph_a, absolute_a = induced_adjacency(
            adj, instance, component.side_a
        )
        graph_b, absolute_b = induced_adjacency(
            adj, instance, component.side_b
        )
        zero_a = solver.number(graph_a)
        zero_b = solver.number(graph_b)
        lower_a = len(absolute_a) - zero_a
        lower_b = len(absolute_b) - zero_b
        required = z0.bit_count() + lower_a + lower_b
        checks.append({
            "component": [
                instance.outside[local]
                for local in base.bits(component.vertices)
            ],
            "side_a": list(absolute_a),
            "side_b": list(absolute_b),
            "zero_forcing_a": zero_a,
            "zero_forcing_b": zero_b,
            "rank_lower_a": lower_a,
            "rank_lower_b": lower_b,
            "zero_factor_dimension": z0.bit_count(),
            "required_dimension": required,
            "passed": required <= COORDINATES,
        })
    return all(check["passed"] for check in checks), tuple(checks)


def passed_category(block: bool, zero_forcing: bool, nonbipartite: bool) -> str:
    names = []
    if block:
        names.append("block")
    if zero_forcing:
        names.append("zero_forcing")
    if nonbipartite:
        names.append("nonbipartite")
    return "+".join(names) if names else "none"


@dataclass
class SeedDecision:
    feasible: bool
    z0_considered: int
    z0_matchable: int
    z0_zero_forcing_passed: int
    z0_zero_forcing_failed: int
    z0_block_passed_after_zero_forcing: int
    z0_block_failed_after_zero_forcing: int
    z0_nonbipartite_failed: int


def solve_seed(
    adj: Sequence[int],
    instance: base.Instance,
    inertia: base.SympyInertiaCache,
    zero_forcing: IndependentZeroForcing,
) -> SeedDecision:
    counts = {
        "z0_considered": 0,
        "z0_matchable": 0,
        "z0_zero_forcing_passed": 0,
        "z0_zero_forcing_failed": 0,
        "z0_block_passed_after_zero_forcing": 0,
        "z0_block_failed_after_zero_forcing": 0,
        "z0_nonbipartite_failed": 0,
    }
    for z0 in base.z0_subsets(instance.eligible_z0):
        counts["z0_considered"] += 1
        if not base.hall_matchable(z0, instance.defects):
            continue
        counts["z0_matchable"] += 1
        local_components = base.components(instance, z0)
        rank_passed, _ = zero_forcing_system(
            adj, instance, z0, local_components, zero_forcing
        )
        if not rank_passed:
            counts["z0_zero_forcing_failed"] += 1
            continue
        counts["z0_zero_forcing_passed"] += 1
        block_passed = all(
            base.component_block_passes(
                adj, instance, z0, component, inertia
            )
            for component in local_components
            if component.bipartite
        )
        if not block_passed:
            counts["z0_block_failed_after_zero_forcing"] += 1
            continue
        counts["z0_block_passed_after_zero_forcing"] += 1
        odd = [
            component.vertices
            for component in local_components
            if not component.bipartite
        ]
        if base.nonbipartite_system_passes(instance, z0, odd):
            return SeedDecision(True, **counts)
        counts["z0_nonbipartite_failed"] += 1
    return SeedDecision(False, **counts)


CORE_FIELDS = (
    "rejected",
    "seeds_checked",
    "impossible_seeds",
    "z0_considered",
    "z0_matchable",
    "z0_zero_forcing_passed",
    "z0_zero_forcing_failed",
    "z0_block_passed_after_zero_forcing",
    "z0_block_failed_after_zero_forcing",
    "z0_nonbipartite_failed",
)


def evaluate_graph(adj: Sequence[int]) -> dict:
    base.validate_graph(adj)
    if next(base.clique_masks(adj, 7), 0):
        raise ValueError("independent all-same-Z0 target unexpectedly has K7")
    inertia = base.SympyInertiaCache()
    zero_forcing = IndependentZeroForcing()
    totals = {name: 0 for name in CORE_FIELDS[3:]}
    seeds_checked = 0
    first_seed = None
    first_instance = None
    for seed_mask in base.clique_masks(adj, COORDINATES):
        seeds_checked += 1
        instance = base.build_instance(adj, seed_mask)
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


def independent_certificate_rows(
    adj: Sequence[int], instance: base.Instance
) -> tuple[dict, ...]:
    inertia = base.SympyInertiaCache()
    zero_forcing = IndependentZeroForcing()
    rows = []
    for z0 in base.z0_subsets(instance.eligible_z0):
        absolute_z0 = tuple(instance.outside[local] for local in base.bits(z0))
        matchable = base.hall_matchable(z0, instance.defects)
        if matchable:
            local_components = base.components(instance, z0)
            rank_passed, rank_checks = zero_forcing_system(
                adj, instance, z0, local_components, zero_forcing
            )
            block_passed = all(
                base.component_block_passes(
                    adj, instance, z0, component, inertia
                )
                for component in local_components
                if component.bipartite
            )
            odd = [
                component.vertices
                for component in local_components
                if not component.bipartite
            ]
            nonbipartite_passed = base.nonbipartite_system_passes(
                instance, z0, odd
            )
        else:
            rank_passed = False
            rank_checks = ()
            block_passed = False
            nonbipartite_passed = False
        rows.append({
            "Z0": absolute_z0,
            "matchable": matchable,
            "block": block_passed,
            "zero_forcing": rank_passed,
            "nonbipartite": nonbipartite_passed,
            "category": passed_category(
                block_passed, rank_passed, nonbipartite_passed
            ),
            "failed_zero_forcing_components": tuple(
                check for check in rank_checks if not check["passed"]
            ),
        })
    return tuple(rows)


def verify_internal_report(report: dict, archive: dict) -> None:
    if report.get("schema") != REPORT_SCHEMA or report.get("status") != "COMPLETE":
        raise ValueError("production all-same-Z0 report has wrong schema/status")
    if archive.get("schema") != CERTIFICATE_SCHEMA or archive.get("status") != "COMPLETE":
        raise ValueError("production certificate archive has wrong schema/status")
    if report["certificate_archive"]["sha256"] != sha256(
        ROOT / report["certificate_archive"]["path"]
    ):
        raise AssertionError("production certificate archive hash differs")
    rejected = [
        item["index"]
        for item in report["graph_results"]
        if item["decision"]["rejected"]
    ]
    if rejected != report["rejected_indices"] or rejected != archive["rejected_indices"]:
        raise AssertionError("production rejection lists disagree")
    aggregate_fields = tuple(
        name
        for name in report["graph_results"][0]["decision"]
        if name not in ("rejected", "first_impossible_seed")
    )
    for name in aggregate_fields:
        observed = sum(
            item["decision"][name] for item in report["graph_results"]
        )
        if observed != report[name]:
            raise AssertionError(f"production aggregate differs for {name}")


def verify(output: Path) -> dict:
    observed = {name: sha256(ROOT / name) for name in EXPECTED}
    if observed != EXPECTED:
        raise ValueError(f"all-same-Z0 verifier dependency hash mismatch: {observed}")
    report = json.loads(
        (ROOT / "d6_k6_all_same_z0_report.json").read_text(encoding="utf-8")
    )
    archive = json.loads(
        (ROOT / "d6_k6_all_same_z0_certificates.json").read_text(encoding="utf-8")
    )
    verify_internal_report(report, archive)

    manifest = json.loads(
        (ROOT / "d6_current_residue_manifest.json").read_text(encoding="utf-8")
    )
    indices = manifest["classes"]["K6_only"]["final_indices"]
    if len(indices) != EXPECTED_INPUT or base.stable_hash(indices) != EXPECTED_INDICES_SHA256:
        raise ValueError("independent all-same-Z0 selection differs")
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
    if list(production) != indices:
        raise AssertionError("production graph order differs from manifest")

    started = time.perf_counter()
    rejected = []
    certificate_summaries = []
    for index in indices:
        adjacency = by_index[index]["adjacency"]
        independent = evaluate_graph(adjacency)
        expected = production[index]
        for name in CORE_FIELDS:
            if independent[name] != expected[name]:
                raise AssertionError(
                    f"graph {index} {name}: independent={independent[name]} "
                    f"production={expected[name]}"
                )
        if independent["first_impossible_seed"] != expected["first_impossible_seed"]:
            raise AssertionError(f"graph {index} first impossible seed differs")
        if not independent["rejected"]:
            continue
        rejected.append(index)
        certificate = certificates[index]
        rows = independent_certificate_rows(
            adjacency, independent["first_instance"]
        )
        production_rows = certificate["choices"]
        if len(rows) != len(production_rows):
            raise AssertionError(f"graph {index} certificate row count differs")
        histogram: dict[str, int] = {}
        for independent_row, production_row in zip(rows, production_rows):
            fields = {
                "Z0": tuple(production_row["Z0"]),
                "matchable": production_row["Z0_allowed_mask_matchable"],
                "block": production_row["block_support_passed"],
                "zero_forcing": production_row["zero_forcing_rank_passed"],
                "nonbipartite": production_row[
                    "nonbipartite_joint_support_passed"
                ],
                "category": production_row["category"],
                "failed_zero_forcing_components": tuple(
                    production_row["failed_zero_forcing_components"]
                ),
            }
            if fields != independent_row:
                raise AssertionError(f"graph {index} per-Z0 certificate differs")
            category = independent_row["category"]
            histogram[category] = histogram.get(category, 0) + 1
        if histogram != certificate["histogram"]:
            raise AssertionError(f"graph {index} certificate histogram differs")
        if histogram.get("block+zero_forcing+nonbipartite", 0):
            raise AssertionError(f"graph {index} certificate has common Z0")
        certificate_summaries.append({
            "index": index,
            "first_impossible_seed": independent["first_impossible_seed"],
            "Z0_choices": len(rows),
            "histogram": histogram,
        })

    if rejected != report["rejected_indices"]:
        raise AssertionError("independent four-index rejection set differs")

    positive = evaluate_graph(lower_bound_18_graph())
    if positive["rejected"] or positive["seeds_checked"] != 32:
        raise AssertionError("known realizable 18-point control failed")

    result = {
        "schema": VERIFICATION_SCHEMA,
        "status": "PASS",
        "method": (
            "frozen independent block/nonbipartite reconstruction plus a new "
            "independent ordinary-zero-forcing subset/parallel-closure solver"
        ),
        "graphs_recomputed": len(indices),
        "graphs_rejected": len(rejected),
        "graphs_surviving": len(indices) - len(rejected),
        "rejected_indices": rejected,
        "new_rejected_indices": report["new_rejected_indices"],
        "all_graph_core_fields_matched": list(CORE_FIELDS),
        "all_first_impossible_seeds_matched": True,
        "all_archived_Z0_system_booleans_matched": True,
        "all_archived_zero_forcing_component_certificates_matched": True,
        "certificate_summaries": certificate_summaries,
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
        default=ROOT / "d6_k6_all_same_z0_verification.json",
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
