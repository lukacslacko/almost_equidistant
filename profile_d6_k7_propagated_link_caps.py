#!/usr/bin/env python3
"""Measure propagated clique-link caps beyond the K7 sparse-value sample."""

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

import d6_k7_propagated_link_caps as link_caps
import d6_k7_rank_reference as reference
import d6_k7_small_support_value as sparse_value
import d6_k7_support_propagation as propagation


ROOT = Path(__file__).resolve().parent
EXPECTED_SHA256 = {
    "d6_k7_rank_reference.py": (
        "e94e2fd92ad03a921b26440827b8365cb4740d9a58399379860f7106fec37db0"
    ),
    "d6_k7_support_propagation.py": (
        "578103f70d2906fbb3d449b7b7b5c1b6c86e5ae17afb0c341602a543bea69f0c"
    ),
    "d6_k7_small_support_value.py": (
        "bed5987c89fbef2ef36762051ae7fa3414fd381a6246032b06857a557f13743f"
    ),
    "d6_k7_small_support_value_report.json": (
        "9beeb04e709fcc0e4d755bb8a481b7d320a2c9470221f4b66dad5b7f7f0ddf07"
    ),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def git_provenance() -> dict:
    def git(*arguments: str) -> str:
        return subprocess.check_output(
            ("git", *arguments), cwd=ROOT, text=True, encoding="utf-8"
        ).strip()

    status = git("status", "--porcelain=v1", "--untracked-files=all")
    return {
        "commit": git("rev-parse", "HEAD"),
        "branch": git("branch", "--show-current"),
        "dirty": bool(status),
        "porcelain_sha256": hashlib.sha256(status.encode("utf-8")).hexdigest(),
        "porcelain_lines": status.splitlines(),
    }


def analyze_graph(graph: dict) -> dict:
    adj = tuple(int(row) for row in graph["adjacency"])
    index = int(graph["index"])
    reference.validate_graph(adj)
    support_solver = reference.SupportSolver()
    zero_forcing = reference.ZeroForcingSolver()
    clique_solver = reference.CliqueStructureSolver()
    counts: Counter[str] = Counter()
    rejecting_seed = 0

    for seed_mask in reference.clique_masks(adj, 7):
        counts["seeds"] += 1
        _, outside, defects, ladj, eligible = reference.seed_instance(
            adj, seed_mask
        )
        total_term_rank = reference.matching_size(defects)
        prior_seed_passes = False
        combined_seed_passes = False
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
                counts["propagation_passing_families"] += 1
                value_result = sparse_value.check_small_support_masks(
                    graph_n, propagated.propagated_masks
                )
                link_failure = link_caps.first_link_cap_failure(
                    z_supports, propagated.propagated_masks
                )
                if link_failure is not None:
                    counts["link_cap_failures"] += 1
                    counts[f"link_cap_failure:k{link_failure.subset_size}"] += 1
                else:
                    counts["link_cap_passes"] += 1
                if not value_result.feasible:
                    counts["prior_value_failures"] += 1
                    continue
                prior_seed_passes = True
                counts["prior_value_passes"] += 1
                if link_failure is not None:
                    counts["marginal_link_cap_failures"] += 1
                    counts[
                        f"marginal_link_cap_failure:k{link_failure.subset_size}"
                    ] += 1
                    continue
                combined_seed_passes = True
                cover_passes = True
                counts["combined_passing_families"] += 1
                break
            if cover_passes:
                counts["combined_passing_covers"] += 1
                break
        if not prior_seed_passes:
            raise AssertionError(
                f"selected prior survivor {index} fails seed {seed_mask}"
            )
        if not combined_seed_passes:
            rejecting_seed = seed_mask
            break

    return {
        "index": index,
        "decision": "REJECTED" if rejecting_seed else "SURVIVOR",
        "first_failing_seed": rejecting_seed,
        "counts": dict(sorted(counts.items())),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input", type=Path, default=ROOT / "d6_k7_rank_sample.json"
    )
    parser.add_argument(
        "--value-report", type=Path,
        default=ROOT / "d6_k7_small_support_value_report.json",
    )
    parser.add_argument(
        "--output", type=Path,
        default=ROOT / "d6_k7_propagated_link_caps_report.json",
    )
    args = parser.parse_args()

    observed = {name: sha256(ROOT / name) for name in EXPECTED_SHA256}
    if observed != EXPECTED_SHA256:
        raise SystemExit(f"dependency hash mismatch: {observed} != {EXPECTED_SHA256}")
    input_hash = sha256(args.input)
    if input_hash != "bf699f36a1c90b6d6751498e52fab73e6b58db322c5b71b04accf32623c2c0dd":
        raise SystemExit(f"unexpected sample input hash {input_hash}")
    value_report = json.loads(args.value_report.read_text(encoding="utf-8"))
    if sha256(args.value_report) != EXPECTED_SHA256[
        "d6_k7_small_support_value_report.json"
    ]:
        raise SystemExit("value report path is not the pinned sample report")
    survivor_indices = [
        int(entry["index"])
        for entry in value_report["per_graph"]
        if entry["decision"] == "SURVIVOR"
    ]
    if len(survivor_indices) != 61:
        raise SystemExit(f"expected 61 sparse-value survivors, got {len(survivor_indices)}")
    payload = json.loads(args.input.read_text(encoding="utf-8"))
    graph_by_index = {
        int(graph["index"]): graph for graph in payload["graphs"]
    }

    started = time.perf_counter()
    started_unix = time.time()
    per_graph = []
    totals: Counter[str] = Counter()
    for position, index in enumerate(survivor_indices, 1):
        result = analyze_graph(graph_by_index[index])
        per_graph.append(result)
        totals.update(result["counts"])
        print(
            f"{position}/{len(survivor_indices)} graph {index}: "
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
            "Exact propagated clique-link-cap profile on the 61 sparse-value "
            "sample survivors."
        ),
        "input": str(args.input),
        "input_sha256": input_hash,
        "value_report": str(args.value_report),
        "value_report_sha256": sha256(args.value_report),
        "dependencies_sha256": observed,
        "link_source_sha256": sha256(ROOT / "d6_k7_propagated_link_caps.py"),
        "proof_sha256": sha256(ROOT / "d6_k7_propagated_link_caps.md"),
        "test_sha256": sha256(ROOT / "test_d6_k7_propagated_link_caps.py"),
        "profiler_sha256": sha256(Path(__file__)),
        "command": [sys.executable, *sys.argv],
        "started_unix": started_unix,
        "elapsed_seconds": time.perf_counter() - started,
        "platform": platform.platform(),
        "machine": platform.machine(),
        "python_version": sys.version,
        "git": git_provenance(),
        "selected_graphs": len(survivor_indices),
        "rejected_graphs": len(rejected),
        "surviving_graphs": len(survivor_indices) - len(rejected),
        "rejected_indices": rejected,
        "totals": dict(sorted(totals.items())),
        "per_graph": per_graph,
    }
    propagation.atomic_json(args.output, report)
    print(
        f"wrote {args.output}; {len(rejected)} marginal exact rejections; "
        f"SHA-256 {sha256(args.output)}",
        flush=True,
    )


if __name__ == "__main__":
    main()

