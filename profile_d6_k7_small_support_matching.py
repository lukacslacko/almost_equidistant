#!/usr/bin/env python3
"""Replay the K7 sample to measure the disjoint-support rule's marginal.

The proof predicts zero marginal coverage in the current N-layer because
``alpha(graph_n)<=2`` and the frozen support-intersection CSP already rejects
every required edge whose assigned supports are disjoint.  This profiler
replays the 61 frozen sparse-value sample survivors and compares the frozen
and refined CSP result on every propagation family reached by the search.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path

import d6_k7_rank_reference as reference
import d6_k7_small_support_matching as matching
import d6_k7_small_support_value as frozen_value
import d6_k7_support_propagation as propagation


ROOT = Path(__file__).resolve().parent
EXPECTED_SAMPLE_INPUT_SHA256 = (
    "bf699f36a1c90b6d6751498e52fab73e6b58db322c5b71b04accf32623c2c0dd"
)
EXPECTED_FROZEN_REPORT_SHA256 = (
    "9beeb04e709fcc0e4d755bb8a481b7d320a2c9470221f4b66dad5b7f7f0ddf07"
)
EXPECTED_FROZEN_SHA256 = {
    "d6_k7_rank_reference.py": (
        "e94e2fd92ad03a921b26440827b8365cb4740d9a58399379860f7106fec37db0"
    ),
    "d6_k7_support_propagation.py": (
        "578103f70d2906fbb3d449b7b7b5c1b6c86e5ae17afb0c341602a543bea69f0c"
    ),
    "d6_k7_small_support_value.py": (
        "bed5987c89fbef2ef36762051ae7fa3414fd381a6246032b06857a557f13743f"
    ),
    "d6_k7_small_support_value.md": (
        "e8ec6a012f9ea824e9cddb76510840537dba0e603e207cd3220c1b06866cb42a"
    ),
    "profile_d6_k7_small_support_value.py": (
        "c8231ff47263b348e3daba2b32cd3eb693ff70f3db5b10aa808b1b5296589573"
    ),
    "test_d6_k7_small_support_value.py": (
        "6a2ee7a97844b5922a5467e1c3a8cf4c8a49240014b4154f45307335dfea5c44"
    ),
}
EXPECTED_FROZEN_SAMPLE_SELECTED = 80
EXPECTED_FROZEN_SAMPLE_REJECTED = 19
EXPECTED_SELECTED = 61


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def git_provenance() -> dict:
    def run(*arguments: str) -> str:
        return subprocess.check_output(
            ("git", *arguments), cwd=ROOT, text=True, encoding="utf-8"
        ).strip()

    status = run("status", "--porcelain=v1", "--untracked-files=all")
    return {
        "commit": run("rev-parse", "HEAD"),
        "branch": run("branch", "--show-current"),
        "dirty": bool(status),
        "status_porcelain": status.splitlines(),
        "status_porcelain_sha256": hashlib.sha256(
            status.encode("utf-8")
        ).hexdigest(),
    }


def verify_and_select(
    input_path: Path, frozen_report_path: Path
) -> tuple[list[dict], dict]:
    input_hash = sha256(input_path)
    report_hash = sha256(frozen_report_path)
    if input_hash != EXPECTED_SAMPLE_INPUT_SHA256:
        raise ValueError(f"sample input hash mismatch: {input_hash}")
    if report_hash != EXPECTED_FROZEN_REPORT_SHA256:
        raise ValueError(f"frozen sample report hash mismatch: {report_hash}")
    observed = {name: sha256(ROOT / name) for name in EXPECTED_FROZEN_SHA256}
    if observed != EXPECTED_FROZEN_SHA256:
        raise ValueError(
            f"frozen dependency hash mismatch: {observed} != "
            f"{EXPECTED_FROZEN_SHA256}"
        )

    payload = json.loads(input_path.read_text(encoding="utf-8"))
    graphs = payload.get("graphs")
    if not isinstance(graphs, list):
        raise ValueError("sample input has no graph list")
    graph_by_index = {int(graph["index"]): graph for graph in graphs}
    if len(graph_by_index) != len(graphs):
        raise ValueError("sample input has duplicate graph indices")

    frozen_report = json.loads(frozen_report_path.read_text(encoding="utf-8"))
    if (
        frozen_report.get("input_sha256") != input_hash
        or int(frozen_report.get("selected_graphs", -1))
        != EXPECTED_FROZEN_SAMPLE_SELECTED
        or int(frozen_report.get("rejected_graphs", -1))
        != EXPECTED_FROZEN_SAMPLE_REJECTED
        or int(frozen_report.get("surviving_graphs", -1)) != EXPECTED_SELECTED
        or frozen_report.get("dependency_sha256")
        != {name: observed[name] for name in (
            "d6_k7_rank_reference.py", "d6_k7_support_propagation.py"
        )}
        or frozen_report.get("value_source_sha256")
        != observed["d6_k7_small_support_value.py"]
        or frozen_report.get("profiler_source_sha256")
        != observed["profile_d6_k7_small_support_value.py"]
        or frozen_report.get("proof_sha256")
        != observed["d6_k7_small_support_value.md"]
        or frozen_report.get("test_source_sha256")
        != observed["test_d6_k7_small_support_value.py"]
    ):
        raise ValueError("frozen sample report metadata is inconsistent")

    selected_indices = [
        int(entry["index"])
        for entry in frozen_report.get("per_graph", ())
        if entry.get("decision") == "SURVIVOR"
    ]
    if len(selected_indices) != EXPECTED_SELECTED:
        raise ValueError("unexpected frozen sample survivor count")
    if missing := set(selected_indices) - set(graph_by_index):
        raise ValueError(f"sample survivors absent from input: {sorted(missing)}")
    return [graph_by_index[index] for index in selected_indices], {
        "input_sha256": input_hash,
        "frozen_report_sha256": report_hash,
        "frozen_dependency_sha256": observed,
        "selected_indices": selected_indices,
        "selected_indices_sha256": hashlib.sha256(
            json.dumps(selected_indices, separators=(",", ":")).encode("ascii")
        ).hexdigest(),
    }


def analyze_graph(graph: dict) -> dict:
    adjacency = tuple(int(row) for row in graph["adjacency"])
    index = int(graph["index"])
    reference.validate_graph(adjacency)
    support_solver = reference.SupportSolver()
    zero_forcing = reference.ZeroForcingSolver()
    clique_solver = reference.CliqueStructureSolver()
    counts: Counter[str] = Counter()
    first_failing_seed = 0

    for seed_mask in reference.clique_masks(adjacency, 7):
        counts["seeds"] += 1
        _, outside, defects, ladj, eligible = reference.seed_instance(
            adjacency, seed_mask
        )
        total_term_rank = reference.matching_size(defects)
        frozen_seed_passes = False
        refined_seed_passes = False
        for zmask in reference.eligible_covers(ladj, eligible, cap=3):
            counts["covers"] += 1
            baseline = reference.analyze_cover(
                adjacency,
                outside,
                defects,
                zmask,
                support_solver,
                zero_forcing,
                total_term_rank,
                clique_solver,
            )
            if baseline.enhanced_joint_failed:
                continue
            counts["baseline_passing_covers"] += 1
            zvertices = tuple(reference.bits(zmask))
            nvertices = tuple(
                vertex for vertex in range(len(outside))
                if not (zmask & (1 << vertex))
            )
            graph_n = reference.induced_graph(
                adjacency, [outside[vertex] for vertex in nvertices]
            )
            # This is the premise making the new rule redundant.
            reference.validate_graph(graph_n)
            counts["alpha_two_graph_n_checks"] += 1
            z_allowed = tuple(defects[vertex] for vertex in zvertices)
            n_allowed = tuple(defects[vertex] for vertex in nvertices)
            cover_passes = False
            for z_supports in propagation.labeled_support_families(z_allowed):
                counts["labeled_z_families"] += 1
                propagated = propagation.analyze_support_assignment(
                    graph_n,
                    z_supports,
                    n_allowed,
                    zero_forcing,
                    clique_solver,
                )
                if propagated.failure is not None:
                    counts[f"propagation_failure:{propagated.failure}"] += 1
                    continue
                counts["propagation_passing_families"] += 1
                support_only = frozen_value.check_small_support_masks(
                    graph_n,
                    propagated.propagated_masks,
                    apply_value_constraints=False,
                )
                if not support_only.feasible:
                    counts["small_support_intersection_failures"] += 1
                    continue
                counts["small_support_intersection_passes"] += 1

                frozen = frozen_value.check_small_support_masks(
                    graph_n, propagated.propagated_masks
                )
                refined = matching.check_small_support_masks(
                    graph_n, propagated.propagated_masks
                )
                counts["frozen_value_passes" if frozen.feasible
                       else "frozen_value_failures"] += 1
                counts["refined_value_passes" if refined.feasible
                       else "refined_value_failures"] += 1
                packing_hits = dict(refined.failures).get(
                    matching.PACKING_FAILURE, 0
                )
                counts["packing_failure_branches"] += packing_hits
                if refined != frozen:
                    counts["frozen_refined_result_disagreements"] += 1
                if frozen.feasible:
                    frozen_seed_passes = True
                if refined.feasible:
                    refined_seed_passes = True
                    cover_passes = True
                    break
            if cover_passes:
                counts["refined_passing_covers"] += 1
                break
        if not frozen_seed_passes:
            raise AssertionError(
                f"frozen survivor graph {index} unexpectedly has no old "
                f"witness at seed {seed_mask}"
            )
        if not refined_seed_passes:
            first_failing_seed = seed_mask
            break

    return {
        "index": index,
        "decision": "REJECTED" if first_failing_seed else "SURVIVOR",
        "first_failing_seed": first_failing_seed,
        "counts": dict(sorted(counts.items())),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input", type=Path, default=ROOT / "d6_k7_rank_sample.json"
    )
    parser.add_argument(
        "--frozen-report",
        type=Path,
        default=ROOT / "d6_k7_small_support_value_report.json",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "d6_k7_small_support_matching_report.json",
    )
    parser.add_argument("--limit", type=int)
    args = parser.parse_args()
    if args.limit is not None and args.limit < 1:
        parser.error("--limit must be positive")

    selected, selection = verify_and_select(args.input, args.frozen_report)
    if args.limit is not None:
        selected = selected[:args.limit]
    started_unix = time.time()
    started = time.perf_counter()
    totals: Counter[str] = Counter()
    per_graph = []
    for position, graph in enumerate(selected, 1):
        result = analyze_graph(graph)
        per_graph.append(result)
        totals.update(result["counts"])
        print(
            f"{position}/{len(selected)} graph {result['index']}: "
            f"{result['decision']}",
            flush=True,
        )

    rejected = [
        entry["index"] for entry in per_graph
        if entry["decision"] == "REJECTED"
    ]
    report = {
        "schema": 1,
        "description": (
            "Exact replay of the three-disjoint-small-support refinement on "
            "the frozen 61-graph sparse-value sample residue."
        ),
        "mathematical_conclusion": (
            "The set-packing rule is exact but branchwise redundant in this "
            "N-layer because alpha(graph_n)<=2 and required-edge support "
            "intersection is already enforced."
        ),
        "full_campaign_recommendation": (
            "Do not run: the proved redundancy forces zero marginal coverage."
        ),
        "input": str(args.input),
        "frozen_report": str(args.frozen_report),
        "selection": selection,
        "source_sha256": sha256(ROOT / "d6_k7_small_support_matching.py"),
        "proof_sha256": sha256(ROOT / "d6_k7_small_support_matching.md"),
        "test_sha256": sha256(
            ROOT / "test_d6_k7_small_support_matching.py"
        ),
        "profiler_sha256": sha256(Path(__file__)),
        "git": git_provenance(),
        "command": [sys.executable, *sys.argv],
        "started_unix": started_unix,
        "elapsed_seconds": time.perf_counter() - started,
        "platform": platform.platform(),
        "machine": platform.machine(),
        "python_version": sys.version,
        "selected_graphs": len(selected),
        "rejected_graphs": len(rejected),
        "surviving_graphs": len(selected) - len(rejected),
        "rejected_indices": rejected,
        "packing_failure_branches": totals["packing_failure_branches"],
        "frozen_refined_result_disagreements": totals[
            "frozen_refined_result_disagreements"
        ],
        "totals": dict(sorted(totals.items())),
        "per_graph": per_graph,
    }
    propagation.atomic_json(args.output, report)
    print(
        f"wrote {args.output}; {len(rejected)} marginal graph rejections, "
        f"{totals['packing_failure_branches']} packing branch hits; "
        f"SHA-256 {sha256(args.output)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
