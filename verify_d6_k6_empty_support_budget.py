#!/usr/bin/env python3
"""Independent full replay of the K6 empty-defect support-budget result.

The checker never imports ``d6_k6_empty_support_budget``.  It independently
rebuilds the singleton alternatives, runs its own 64-state Hall DP, replays
all 822 parent-residue graphs in a fresh ordered process pool, and exactly
compares every archived Z0/component certificate row.  The caller must supply
the expected SHA-256 of the production report as an explicit review boundary.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
import os
import platform
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Sequence

import d6_k6_arbitrary_subset_hall as parent
import d6_k6_psd_z_hereditary as groups_parent
import d6_k6_psd_zmatrix as ranks
from d6_k6_bipartite_rank_reference import (
    ZeroForcingSolver,
    induced_required_graph,
    lorentz_components,
)
from d6_k6_lorentz_reference import (
    COORDINATES,
    build_instance,
    clique_masks,
    support_matching,
    vertices,
)
from d6_k6_normal_inertia import InertiaCache
from test_d6_k6_lorentz import lower_bound_18_graph


ROOT = Path(__file__).resolve().parent
EXPECTED_INPUT = 822
EXPECTED_INPUT_SHA256 = (
    "cf94aac41eba28904ee6254f9c50996e73d38dcc15316fb6a411b29ff5cd05d8"
)
EXPECTED_PRODUCTION_SOURCE_SHA256 = (
    "284b8c3bda4d2a43581665ade3d52d9295b29d460418454e3a760897cabc5fac"
)
REPORT_SCHEMA = "d6-k6-empty-support-budget-v1"
CERTIFICATE_SCHEMA = "d6-k6-empty-support-budget-certificates-v1"
VERIFICATION_SCHEMA = "d6-k6-empty-support-budget-verification-v1"
GLOBAL_EMPTY_BUDGET = 2
DEFAULT_WORKERS = 11

EXPECTED_DEPENDENCIES = {
    "d6_k6_arbitrary_subset_hall.py": (
        "cadc0048fcc7a77299bc2950e2eca57bfed87d54ed6ffd59886d155b2c9a88ff"
    ),
    "d6_k6_arbitrary_subset_hall_report.json": (
        "bcc4baecf28a95443a5b13e87a796cbba0a97cd42fa76573119bfd81ef3735fc"
    ),
    "d6_k6_arbitrary_subset_hall_certificates.json": (
        "b7fb04f769ddaf5d36630372f5134f2e4a32e572f620b26059667e9cf65f4a23"
    ),
    "d6_k6_arbitrary_subset_hall_verification.json": (
        "e6042edfff3683043903f9be544322379e85673207e35ee87349637ff4c73058"
    ),
    "d6_k6_psd_z_hereditary.py": (
        "7b869895cb9b4f7f3bbf1a1d1b11ec71507b6355c8ec0b7d7f8eaec2a51b180d"
    ),
    "d6_k6_psd_zmatrix.py": (
        "24c2b3abfce5074560183d17f77778f96298bef6e38ece5bc937ccb2a4cba736"
    ),
    "d6_k6_bipartite_rank_reference.py": (
        "e760fef74c423a89f1e1f9d08909fa1223e0389900cf3d45b0140ee296dea922"
    ),
    "d6_k6_lorentz_reference.py": (
        "ef6a0877a5987119c92ec1bd9271a6a47774fc9e4992c10e1f28b6d5373592c2"
    ),
    "d6_k6_normal_inertia.py": (
        "2d00ab40eb97aa79ffec3b8134b83c4fd5f703cbd691f094a25bf3ff9f3ee023"
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


def worker_initializer() -> None:
    for name in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[name] = "1"


def dependency_hashes() -> dict[str, str]:
    observed = {
        name: sha256(ROOT / name) for name in EXPECTED_DEPENDENCIES
    }
    if observed != EXPECTED_DEPENDENCIES:
        raise ValueError(
            "independent empty-support dependency mismatch: "
            f"observed={observed}"
        )
    return observed


def independent_hall(groups: Sequence[groups_parent.Group]) -> bool:
    """Exhaust one selectable rank/support alternative from every group."""

    states = {0: 0}
    for group in groups:
        new: dict[int, int] = {}
        for coordinates, rank in states.items():
            for choice in group.choices:
                union = coordinates | choice.coordinates
                total = rank + choice.rank
                if total > union.bit_count():
                    return False
                if total > new.get(union, -1):
                    new[union] = total
        states = new
    return True


def singleton_data(instance, side) -> dict[int, tuple[int, int]]:
    local = {vertex: i for i, vertex in enumerate(instance.outside)}
    answer = {}
    for group_index, detail in enumerate(side.components):
        absolute = tuple(detail["vertices"])
        rank = int(detail["componentwise_fused_rank_lower"])
        if len(absolute) == 1 and rank == 0:
            vertex = absolute[0]
            answer[group_index] = (vertex, instance.defects[local[vertex]])
        elif len(absolute) == 1 and rank != 1:
            raise AssertionError("unexpected singleton rank lower")
    return answer


def impose_singleton_status(
    groups: Sequence[groups_parent.Group],
    zero: dict[int, tuple[int, int]],
    empty: int | None,
) -> tuple[groups_parent.Group, ...]:
    answer = []
    for index, group in enumerate(groups):
        if index not in zero:
            answer.append(group)
            continue
        vertex, mask = zero[index]
        omit = groups_parent.Choice(0, 0, (), "independent_omit")
        if index == empty:
            choices = (omit,)
        else:
            choices = (
                omit,
                groups_parent.Choice(
                    1, mask, (vertex,), "independent_forced_nonempty"
                ),
            )
        answer.append(groups_parent.Group(group.label, choices, True))
    return tuple(answer)


def orientation_decision(
    instance, z0, side_a, side_b, rank_solver, negative: str
) -> dict:
    groups_a, _ = parent.make_side_groups(
        instance, "A", side_a, rank_solver
    )
    groups_b, _ = parent.make_side_groups(
        instance, "B", side_b, rank_solver
    )
    zero = singleton_data(instance, side_a if negative == "A" else side_b)
    zgroups = tuple(
        groups_parent.singleton_group(
            f"Z0:{instance.outside[i]}",
            instance.outside[i],
            instance.defects[i],
        )
        for i in vertices(z0)
    )
    for empty in (None, *zero.keys()):
        if negative == "A":
            test_a = impose_singleton_status(groups_a, zero, empty)
            test_b = groups_b
        else:
            test_a = groups_a
            test_b = impose_singleton_status(groups_b, zero, empty)
        if independent_hall(test_a + test_b + zgroups):
            return {
                "possible": True,
                "minimum_empty": int(empty is not None),
                "negative_singletons": [value[0] for value in zero.values()],
                "chosen_empty": None if empty is None else zero[empty][0],
            }
    return {
        "possible": False,
        "minimum_empty": None,
        "negative_singletons": [value[0] for value in zero.values()],
        "chosen_empty": None,
    }


def component_decision(
    adjacency, instance, z0, component, inertia, forcing, subset
) -> dict:
    graph_a, absolute_a = induced_required_graph(
        adjacency, instance, component.side_a
    )
    graph_b, absolute_b = induced_required_graph(
        adjacency, instance, component.side_b
    )
    alternatives = []
    for kind, sign_a, sign_b, negative in (
        ("A_positive", "positive", "negative", "B"),
        ("A_negative", "negative", "positive", "A"),
    ):
        side_a = ranks.side_rank_lower(
            graph_a, absolute_a, sign_a, inertia, forcing
        )
        side_b = ranks.side_rank_lower(
            graph_b, absolute_b, sign_b, inertia, forcing
        )
        decision = orientation_decision(
            instance, z0, side_a, side_b, subset, negative
        )
        alternatives.append({"kind": kind, **decision})
    light = z0 | component.component
    light_passed = (
        light.bit_count() <= COORDINATES
        and support_matching(light, instance.defects) is not None
    )
    alternatives.append({
        "kind": "lightlike",
        "possible": light_passed,
        "minimum_empty": 0 if light_passed else None,
        "orthonormal_vectors": light.bit_count(),
    })
    costs = [
        row["minimum_empty"] for row in alternatives if row["possible"]
    ]
    return {
        "possible": bool(costs),
        "minimum_empty": min(costs) if costs else None,
        "alternatives": alternatives,
    }


def absolute_vertices(instance, mask: int) -> list[int]:
    return [instance.outside[i] for i in vertices(mask)]


def seed_decision(adjacency, instance) -> dict:
    inertia = InertiaCache()
    forcing = ZeroForcingSolver()
    subset = parent.PositiveSubsetRank(inertia, forcing)
    rows = []
    counts = {"z0": 0, "components": 0, "orientations": 0}
    for z0 in ranks.z0_subsets(instance.eligible_z0_mask):
        counts["z0"] += 1
        if support_matching(z0, instance.defects) is None:
            rows.append({
                "Z0": absolute_vertices(instance, z0),
                "matchable": False,
                "terminal": "Z0_ALLOWED_SUPPORT_UNMATCHABLE",
            })
            continue
        total_empty = 0
        component_costs: list[int | None] = []
        component_rows = []
        possible = True
        terminal = None
        for component in lorentz_components(instance, z0):
            if not component.bipartite:
                continue
            counts["components"] += 1
            counts["orientations"] += 2
            decision = component_decision(
                adjacency, instance, z0, component,
                inertia, forcing, subset,
            )
            component_rows.append({
                "vertices": absolute_vertices(instance, component.component),
                **decision,
            })
            if not decision["possible"]:
                component_costs.append(None)
                possible = False
                terminal = "COMPONENT_IMPOSSIBLE"
                break
            cost = decision["minimum_empty"]
            assert cost is not None
            component_costs.append(cost)
            total_empty += cost
            if (
                any(item is None for item in component_costs)
                or sum(int(item) for item in component_costs if item is not None)
                > GLOBAL_EMPTY_BUDGET
            ):
                possible = False
                terminal = "GLOBAL_EMPTY_BUDGET_EXCEEDED"
                break
        if possible:
            return {
                "feasible": True,
                "counts": counts,
                "chosen_Z0": absolute_vertices(instance, z0),
                "minimum_empty_used": total_empty,
                "certificate_rows": None,
            }
        rows.append({
            "Z0": absolute_vertices(instance, z0),
            "matchable": True,
            "terminal": terminal,
            "minimum_empty_before_failure": total_empty,
            "global_empty_budget": GLOBAL_EMPTY_BUDGET,
            "components": component_rows,
        })
    return {
        "feasible": False,
        "counts": counts,
        "chosen_Z0": None,
        "minimum_empty_used": None,
        "certificate_rows": rows,
    }


def evaluate(adjacency) -> dict:
    totals = {"z0": 0, "components": 0, "orientations": 0}
    seeds = 0
    for seed_mask in clique_masks(adjacency, COORDINATES):
        seeds += 1
        instance = build_instance(adjacency, vertices(seed_mask))
        decision = seed_decision(adjacency, instance)
        for name in totals:
            totals[name] += decision["counts"][name]
        if not decision["feasible"]:
            return {
                "rejected": True,
                "seeds_checked": seeds,
                **totals,
                "first_impossible_seed": list(instance.seed),
                "certificate": {
                    "seed": list(instance.seed),
                    "outside": list(instance.outside),
                    "allowed_defects": {
                        str(instance.outside[local]): [
                            instance.seed[coordinate]
                            for coordinate in vertices(mask)
                        ]
                        for local, mask in enumerate(instance.defects)
                    },
                    "choices": decision["certificate_rows"],
                },
            }
    return {
        "rejected": False,
        "seeds_checked": seeds,
        **totals,
        "first_impossible_seed": None,
        "certificate": None,
    }


def evaluate_task(task: tuple[int, Sequence[int]]) -> dict:
    index, adjacency = task
    return {"index": index, "decision": evaluate(adjacency)}


def positive_seed_task(task: tuple[Sequence[int], int]) -> dict:
    adjacency, seed_mask = task
    instance = build_instance(adjacency, vertices(seed_mask))
    decision = seed_decision(adjacency, instance)
    return {
        "seed": list(instance.seed),
        "passed": bool(decision["feasible"]),
        "z0": decision["counts"]["z0"],
        "components": decision["counts"]["components"],
        "orientations": decision["counts"]["orientations"],
    }


def parent_records() -> tuple[dict[str, str], dict[int, dict], list[int], list[int]]:
    dependencies = dependency_hashes()
    _, records, indices = parent.verify_inputs()
    report = json.loads(
        (ROOT / "d6_k6_arbitrary_subset_hall_report.json").read_text(
            encoding="utf-8"
        )
    )
    archive = json.loads(
        (ROOT / "d6_k6_arbitrary_subset_hall_certificates.json").read_text(
            encoding="utf-8"
        )
    )
    verification = json.loads(
        (ROOT / "d6_k6_arbitrary_subset_hall_verification.json").read_text(
            encoding="utf-8"
        )
    )
    if report.get("status") != "COMPLETE" or verification.get("status") != "PASS":
        raise ValueError("pinned parent report/verification is not complete")
    rejected = list(report["rejected_indices"])
    if rejected != archive.get("rejected_indices") or rejected != verification.get("rejected_indices"):
        raise ValueError("pinned parent result manifests disagree")
    rejected_set = set(rejected)
    survivors = [index for index in indices if index not in rejected_set]
    if (
        len(survivors) != EXPECTED_INPUT
        or stable_hash(survivors) != EXPECTED_INPUT_SHA256
    ):
        raise ValueError("independent 822-parent reconstruction failed")
    by_index = {record["index"]: record for record in records}
    return dependencies, by_index, survivors, rejected


def assert_import_independence() -> None:
    tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    imports = {
        node.module
        for node in ast.walk(tree)
        if isinstance(node, ast.ImportFrom) and node.module
    }
    imports.update(
        alias.name
        for node in ast.walk(tree)
        if isinstance(node, ast.Import)
        for alias in node.names
    )
    forbidden = {
        "d6_k6_empty_support_budget",
        "probe_d6_k6_empty_support_budget",
    }
    if imports & forbidden:
        raise AssertionError(
            f"independent verifier imports production kernel: {imports & forbidden}"
        )


class VerificationInfrastructureAbort(RuntimeError):
    pass


def verify(
    report_path: Path,
    certificates_path: Path,
    expected_report_sha256: str,
    output: Path,
    workers: int,
) -> dict:
    if (
        len(expected_report_sha256) != 64
        or any(character not in "0123456789abcdef" for character in expected_report_sha256)
    ):
        raise ValueError("--expected-report-sha256 must be 64 lowercase hex digits")
    actual_report_sha256 = sha256(report_path)
    if actual_report_sha256 != expected_report_sha256:
        raise ValueError(
            "explicit report SHA-256 boundary mismatch: "
            f"actual={actual_report_sha256}"
        )
    assert_import_independence()
    dependencies, by_index, indices, inherited = parent_records()
    report = json.loads(report_path.read_text(encoding="utf-8"))
    archive = json.loads(certificates_path.read_text(encoding="utf-8"))
    if report.get("schema") != REPORT_SCHEMA or report.get("status") != "COMPLETE":
        raise ValueError("unexpected production report schema/status")
    if (
        archive.get("schema") != CERTIFICATE_SCHEMA
        or archive.get("status") != "COMPLETE"
    ):
        raise ValueError("unexpected certificate archive schema/status")
    production_source = sha256(ROOT / "d6_k6_empty_support_budget.py")
    if (
        production_source != EXPECTED_PRODUCTION_SOURCE_SHA256
        or report.get("sources", {}).get("d6_k6_empty_support_budget.py")
        != production_source
        or archive.get("production_source_sha256") != production_source
    ):
        raise ValueError("production source boundary mismatch")
    if report.get("sources", {}) | dependencies != report.get("sources", {}):
        raise ValueError("production report omits or changes a pinned dependency")
    if archive.get("dependency_hashes") != dependencies:
        raise ValueError("certificate archive dependency boundary mismatch")
    if (
        report.get("input_graphs") != EXPECTED_INPUT
        or report.get("input_indices_sha256") != EXPECTED_INPUT_SHA256
        or report.get("ordered_input_indices") != indices
        or report.get("inherited_parent_rejections") != inherited
    ):
        raise ValueError("production ordered 822 input boundary mismatch")
    if (
        report.get("certificate_archive", {}).get("sha256")
        != sha256(certificates_path)
        or Path(report["certificate_archive"]["path"]).name
        != certificates_path.name
    ):
        raise ValueError("production report certificate binding mismatch")
    graph_rows = report.get("graph_results")
    if (
        not isinstance(graph_rows, list)
        or [row.get("index") for row in graph_rows] != indices
    ):
        raise ValueError("production graph results are not the ordered residue input")
    archive_by_index = {
        row["index"]: row["certificate"]
        for row in archive.get("certificates", [])
    }
    if (
        archive.get("rejected_indices") != report.get("rejected_indices")
        or list(archive_by_index) != report.get("rejected_indices")
    ):
        raise ValueError("production report/archive rejection ordering mismatch")

    started = time.perf_counter()
    tasks = [(index, by_index[index]["adjacency"]) for index in indices]
    try:
        if workers == 1:
            recomputed_rows = [evaluate_task(task) for task in tasks]
        else:
            with ProcessPoolExecutor(
                max_workers=workers, initializer=worker_initializer
            ) as pool:
                recomputed_rows = list(
                    pool.map(evaluate_task, tasks, chunksize=1)
                )
    except BaseException as error:
        raise VerificationInfrastructureAbort(
            "fresh 822-graph verifier pool aborted; no rejection conclusion"
        ) from error

    if [row["index"] for row in recomputed_rows] != indices:
        raise AssertionError("independent process pool changed ordered input")
    production_by_index = {
        row["index"]: row["decision"] for row in graph_rows
    }
    rejected = []
    z0_rows_validated = 0
    component_witnesses_validated = 0
    for row in recomputed_rows:
        index = row["index"]
        decision = row["decision"]
        certificate = decision.pop("certificate")
        if decision != production_by_index[index]:
            raise AssertionError(
                f"independent full decision mismatch at graph {index}"
            )
        if not decision["rejected"]:
            if index in archive_by_index:
                raise AssertionError(f"survivor {index} has an archived rejection")
            continue
        rejected.append(index)
        archived = archive_by_index.get(index)
        if archived != certificate:
            raise AssertionError(
                f"independent archived Z0/component witness mismatch at {index}"
            )
        choices = archived["choices"]
        z0_rows_validated += len(choices)
        component_witnesses_validated += sum(
            len(choice.get("components", [])) for choice in choices
        )
        for choice in choices:
            if choice["matchable"]:
                if choice["terminal"] not in {
                    "COMPONENT_IMPOSSIBLE",
                    "GLOBAL_EMPTY_BUDGET_EXCEEDED",
                }:
                    raise AssertionError("archived matchable Z0 has no exact terminal")
                for component in choice["components"]:
                    kinds = [row["kind"] for row in component["alternatives"]]
                    if kinds != ["A_positive", "A_negative", "lightlike"]:
                        raise AssertionError("archived component omits an alternative")
            elif choice["terminal"] != "Z0_ALLOWED_SUPPORT_UNMATCHABLE":
                raise AssertionError("archived unmatchable Z0 reason differs")

    if rejected != report["rejected_indices"]:
        raise AssertionError("independent rejection ordering differs")
    residue = [index for index in indices if index not in set(rejected)]
    if (
        report.get("ordered_residue_indices") != residue
        or report.get("ordered_residue_indices_sha256") != stable_hash(residue)
        or report.get("graphs_rejected") != len(rejected)
        or report.get("graphs_surviving") != len(residue)
    ):
        raise AssertionError("production ordered post-filter residue differs")

    positive_graph = lower_bound_18_graph()
    positive_seeds = list(clique_masks(positive_graph, COORDINATES))
    positive_rows = [
        positive_seed_task((positive_graph, seed)) for seed in positive_seeds
    ]
    if len(positive_rows) != 32 or not all(row["passed"] for row in positive_rows):
        raise AssertionError("known realizable control failed one of 32 K6 seeds")
    if report.get("positive_control", {}).get("seed_results") != positive_rows:
        raise AssertionError("production positive-control seed rows differ")

    result = {
        "schema": VERIFICATION_SCHEMA,
        "status": "PASS",
        "method": (
            "fresh ordered 822-graph process pool; independently rebuilt "
            "singleton alternatives and 64-state Hall DP; exact comparison "
            "of every archived Z0/component row; all 32 positive K6 seeds"
        ),
        "report_sha256": actual_report_sha256,
        "certificate_archive_sha256": sha256(certificates_path),
        "production_source_sha256": production_source,
        "verifier_source_sha256": sha256(Path(__file__).resolve()),
        "graphs_recomputed": len(indices),
        "graphs_rejected": len(rejected),
        "graphs_surviving": len(residue),
        "rejected_indices": rejected,
        "ordered_residue_indices_sha256": stable_hash(residue),
        "archived_z0_rows_validated": z0_rows_validated,
        "archived_component_witnesses_validated": component_witnesses_validated,
        "positive_control": {
            "passed": True,
            "K6_seeds": len(positive_rows),
            "seed_results": positive_rows,
        },
        "dependency_hashes": dependencies,
        "import_independence": {
            "production_kernel_imported": False,
            "checked_by_AST": True,
        },
        "runtime": {
            "command": list(sys.argv),
            "workers": workers,
            "wall_seconds": time.perf_counter() - started,
            "python": sys.version,
            "platform": platform.platform(),
        },
    }
    atomic_json(output, result)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--report", type=Path,
        default=ROOT / "d6_k6_empty_support_budget_report.json",
    )
    parser.add_argument(
        "--certificates", type=Path,
        default=ROOT / "d6_k6_empty_support_budget_certificates.json",
    )
    parser.add_argument("--expected-report-sha256", required=True)
    parser.add_argument(
        "--output", type=Path,
        default=ROOT / "d6_k6_empty_support_budget_verification.json",
    )
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    args = parser.parse_args()
    if args.workers < 1 or args.workers > DEFAULT_WORKERS:
        parser.error(f"workers must lie in 1..{DEFAULT_WORKERS}")
    try:
        result = verify(
            args.report.resolve(), args.certificates.resolve(),
            args.expected_report_sha256, args.output.resolve(), args.workers,
        )
    except VerificationInfrastructureAbort as error:
        payload = {
            "schema": VERIFICATION_SCHEMA,
            "status": "INFRA_ABORT",
            "meaning": "no mathematical rejection is inferred",
            "message": str(error),
            "report_sha256": sha256(args.report.resolve()),
            "verifier_source_sha256": sha256(Path(__file__).resolve()),
        }
        atomic_json(args.output.resolve(), payload)
        print(str(error), file=sys.stderr)
        return 2
    print(json.dumps({
        "status": result["status"],
        "graphs_recomputed": result["graphs_recomputed"],
        "graphs_rejected": result["graphs_rejected"],
        "archived_z0_rows_validated": result["archived_z0_rows_validated"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
