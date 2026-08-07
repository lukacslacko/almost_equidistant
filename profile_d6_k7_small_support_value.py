#!/usr/bin/env python3
"""Profile exact one/two-defect value constraints on a small K7 sample.

This is intentionally a single-process research profiler, not a production
full-corpus runner.  It replays the frozen K7 cover analysis and labeled
support propagation, then existentially checks the sparse-value constraints
from :mod:`d6_k7_small_support_value` on every propagation-passing family.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path

import d6_k7_rank_reference as reference
import d6_k7_small_support_value as value
import d6_k7_support_propagation as propagation


ROOT = Path(__file__).resolve().parent
REFERENCE_PATH = ROOT / "d6_k7_rank_reference.py"
PROPAGATION_PATH = ROOT / "d6_k7_support_propagation.py"
EXPECTED_REFERENCE_SHA256 = (
    "e94e2fd92ad03a921b26440827b8365cb4740d9a58399379860f7106fec37db0"
)
EXPECTED_PROPAGATION_SHA256 = (
    "578103f70d2906fbb3d449b7b7b5c1b6c86e5ae17afb0c341602a543bea69f0c"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def git_provenance() -> dict:
    """Capture the exact local Git state at profiler start."""

    def run(*arguments: str) -> str:
        return subprocess.check_output(
            ("git", *arguments),
            cwd=ROOT,
            text=True,
            encoding="utf-8",
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


def analyze_graph(graph: dict) -> dict:
    adj = tuple(int(row) for row in graph["adjacency"])
    index = int(graph["index"])
    reference.validate_graph(adj)
    support_solver = reference.SupportSolver()
    zero_forcing = reference.ZeroForcingSolver()
    clique_solver = reference.CliqueStructureSolver()
    counts: Counter[str] = Counter()
    rejected_seed = 0

    for seed_mask in reference.clique_masks(adj, 7):
        counts["seeds"] += 1
        _, outside, defects, ladj, eligible = reference.seed_instance(
            adj, seed_mask
        )
        total_term_rank = reference.matching_size(defects)
        propagation_seed_passes = False
        value_seed_passes = False
        for zmask in reference.eligible_covers(ladj, eligible, cap=3):
            counts["covers"] += 1
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
            if baseline.enhanced_joint_failed:
                continue
            counts["baseline_passing_covers"] += 1
            zvertices = tuple(reference.bits(zmask))
            nvertices = tuple(
                vertex for vertex in range(len(outside))
                if not (zmask & (1 << vertex))
            )
            graph_n = reference.induced_graph(
                adj, [outside[vertex] for vertex in nvertices]
            )
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
                propagation_seed_passes = True
                counts["propagation_passing_families"] += 1
                support_only = value.check_small_support_masks(
                    graph_n,
                    propagated.propagated_masks,
                    apply_value_constraints=False,
                )
                if not support_only.feasible:
                    counts["small_support_intersection_failures"] += 1
                    continue
                counts["small_support_intersection_passes"] += 1
                sparse_value = value.check_small_support_masks(
                    graph_n, propagated.propagated_masks
                )
                if not sparse_value.feasible:
                    counts["sparse_value_failures"] += 1
                    for reason, number in sparse_value.failures:
                        counts[f"sparse_value_branch:{reason}"] += number
                    continue
                counts["sparse_value_passes"] += 1
                cover_passes = True
                value_seed_passes = True
                break
            if cover_passes:
                counts["value_passing_covers"] += 1
                break
        if not propagation_seed_passes:
            raise AssertionError(
                f"input graph {index} was not a propagation survivor at "
                f"seed {seed_mask}"
            )
        if not value_seed_passes:
            rejected_seed = seed_mask
            break

    return {
        "index": index,
        "decision": "REJECTED" if rejected_seed else "SURVIVOR",
        "first_failing_seed": rejected_seed,
        "counts": dict(sorted(counts.items())),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input", type=Path, default=ROOT / "d6_k7_rank_sample.json"
    )
    parser.add_argument(
        "--support-decisions",
        type=Path,
        default=ROOT / "d6_k7_support_propagation_sample_decisions.tsv",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "d6_k7_small_support_value_report.json",
    )
    parser.add_argument(
        "--strict-h-report",
        type=Path,
        default=ROOT / "d6_k7_strict_h_report.json",
        help="optional independent sample report used only for set overlap",
    )
    parser.add_argument("--limit", type=int)
    parser.add_argument(
        "--indices",
        help="optional comma-separated graph-index override",
    )
    args = parser.parse_args()

    dependency_hashes = {
        "d6_k7_rank_reference.py": sha256(REFERENCE_PATH),
        "d6_k7_support_propagation.py": sha256(PROPAGATION_PATH),
    }
    expected_dependency_hashes = {
        "d6_k7_rank_reference.py": EXPECTED_REFERENCE_SHA256,
        "d6_k7_support_propagation.py": EXPECTED_PROPAGATION_SHA256,
    }
    if dependency_hashes != expected_dependency_hashes:
        raise SystemExit(
            "imported exact dependency hash mismatch: "
            f"observed {dependency_hashes}, expected "
            f"{expected_dependency_hashes}"
        )
    git_state = git_provenance()

    payload = json.loads(args.input.read_text(encoding="utf-8"))
    graph_by_index = {
        int(graph["index"]): graph for graph in payload["graphs"]
    }
    if args.indices:
        selected = [int(item) for item in args.indices.split(",") if item]
    else:
        with args.support_decisions.open(newline="", encoding="ascii") as stream:
            selected = [
                int(row["index"])
                for row in csv.DictReader(stream, delimiter="\t")
                if row["refined_decision"] == "SURVIVOR"
            ]
    if args.limit is not None:
        selected = selected[: args.limit]
    missing = set(selected) - set(graph_by_index)
    if missing:
        raise SystemExit(f"selected indices absent from input: {sorted(missing)}")

    started_unix = time.time()
    started = time.monotonic()
    per_graph = []
    totals: Counter[str] = Counter()
    for position, index in enumerate(selected, 1):
        result = analyze_graph(graph_by_index[index])
        per_graph.append(result)
        totals.update(result["counts"])
        print(
            f"{position}/{len(selected)} graph {index}: {result['decision']}",
            flush=True,
        )
    rejected = [entry["index"] for entry in per_graph
                if entry["decision"] == "REJECTED"]
    strict_rejected: set[int] = set()
    if args.strict_h_report is not None:
        strict_payload = json.loads(
            args.strict_h_report.read_text(encoding="utf-8")
        )
        strict_rejected = {
            int(index)
            for index in strict_payload["graph_decisions"]["combined_rejected"]
        }
    elapsed = time.monotonic() - started
    report = {
        "schema": 1,
        "description": (
            "Exact sparse one/two-defect value profile on selected labeled-"
            "support-propagation survivors; single-process research sample."
        ),
        "input": str(args.input),
        "input_sha256": sha256(args.input),
        "support_decisions": str(args.support_decisions),
        "support_decisions_sha256": sha256(args.support_decisions),
        "strict_h_report": str(args.strict_h_report),
        "strict_h_report_sha256": sha256(args.strict_h_report),
        "dependency_sha256": dependency_hashes,
        "expected_dependency_sha256": expected_dependency_hashes,
        "value_source_sha256": sha256(ROOT / "d6_k7_small_support_value.py"),
        "profiler_source_sha256": sha256(Path(__file__)),
        "proof_sha256": sha256(ROOT / "d6_k7_small_support_value.md"),
        "test_source_sha256": sha256(
            ROOT / "test_d6_k7_small_support_value.py"
        ),
        "git": git_state,
        "command": [sys.executable, *sys.argv],
        "started_unix": started_unix,
        "platform": platform.platform(),
        "machine": platform.machine(),
        "python_version": sys.version,
        "selected_graphs": len(selected),
        "rejected_graphs": len(rejected),
        "surviving_graphs": len(selected) - len(rejected),
        "rejected_indices": rejected,
        "strict_h_rejected_overlap": sorted(set(rejected) & strict_rejected),
        "incremental_beyond_strict_h": sorted(
            set(rejected) - strict_rejected
        ),
        "elapsed_seconds": elapsed,
        "totals": dict(sorted(totals.items())),
        "per_graph": per_graph,
    }
    propagation.atomic_json(args.output, report)
    print(
        f"wrote {args.output}; {len(rejected)} exact rejections; "
        f"SHA-256 {sha256(args.output)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
