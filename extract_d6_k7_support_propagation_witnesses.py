#!/usr/bin/env python3
"""Reconstruct compact witnesses for new K7 support-propagation rejections."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import d6_k7_rank_reference as reference
import d6_k7_support_propagation as propagation


ROOT = Path(__file__).resolve().parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--sample", type=Path, default=ROOT / "d6_k7_rank_sample.json"
    )
    parser.add_argument(
        "--decisions",
        type=Path,
        default=ROOT / "d6_k7_support_propagation_sample_decisions.tsv",
    )
    parser.add_argument(
        "--report",
        type=Path,
        default=ROOT / "d6_k7_support_propagation_report.json",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "d6_k7_support_propagation_witnesses.json",
    )
    args = parser.parse_args()

    report = json.loads(args.report.read_text(encoding="utf-8"))
    expected = {
        "input_sha256": propagation.sha256(args.sample),
        "decisions_sha256": propagation.sha256(args.decisions),
        "source_sha256": propagation.sha256(
            ROOT / "d6_k7_support_propagation.py"
        ),
        "reference_sha256": propagation.sha256(
            ROOT / "d6_k7_rank_reference.py"
        ),
    }
    for field, observed in expected.items():
        if report[field] != observed:
            raise SystemExit(
                f"report {field} is {report[field]}, current artifact is "
                f"{observed}"
            )

    sample = json.loads(args.sample.read_text(encoding="utf-8"))
    graph_by_index = {
        int(graph["index"]): graph for graph in sample["graphs"]
    }
    with args.decisions.open(newline="", encoding="ascii") as stream:
        decision_by_index = {
            int(row["index"]): row
            for row in csv.DictReader(stream, delimiter="\t")
        }

    witnesses = []
    for index in report["incremental_rejected_indices"]:
        graph = graph_by_index[index]
        decision = decision_by_index[index]
        if (
            decision["baseline_decision"] != "SURVIVOR"
            or decision["refined_decision"] != "REJECTED"
        ):
            raise AssertionError(f"bad incremental decision row for {index}")
        adj = tuple(int(row) for row in graph["adjacency"])
        seed_mask = int(decision["first_refined_failing_seed"])
        seed, outside, defects, ladj, eligible = reference.seed_instance(
            adj, seed_mask
        )
        support_solver = reference.SupportSolver()
        zero_forcing = reference.ZeroForcingSolver()
        clique_solver = reference.CliqueStructureSolver()
        total_term_rank = reference.matching_size(defects)
        eligible_covers = reference.eligible_covers(ladj, eligible, cap=3)
        cover_witnesses = []
        for zmask in eligible_covers:
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
            zvertices = tuple(reference.bits(zmask))
            nvertices = tuple(
                vertex for vertex in range(len(outside))
                if not (zmask & (1 << vertex))
            )
            graph_n = reference.induced_graph(
                adj, [outside[vertex] for vertex in nvertices]
            )
            refined = propagation.refine_baseline_passing_cover(
                graph_n,
                tuple(defects[vertex] for vertex in zvertices),
                tuple(defects[vertex] for vertex in nvertices),
                zero_forcing,
                clique_solver,
            )
            if refined.passes:
                raise AssertionError(
                    f"recorded failing seed {seed_mask} passes cover {zmask}"
                )
            cover_witnesses.append(
                {
                    "zmask_outside_local": zmask,
                    "z_global_vertices": [outside[v] for v in zvertices],
                    "z_allowed_masks": [defects[v] for v in zvertices],
                    "labeled_support_families_exhausted": (
                        refined.families_checked
                    ),
                    "failure_counts": dict(refined.failures),
                }
            )
        if not cover_witnesses:
            raise AssertionError(
                f"new rejection {index} had no baseline-passing cover"
            )
        witnesses.append(
            {
                "index": index,
                "failing_seed_mask_global": seed_mask,
                "failing_seed_vertices": seed,
                "eligible_cap3_covers": len(eligible_covers),
                "baseline_passing_covers": len(cover_witnesses),
                "refined_passing_covers": 0,
                "cover_exhaustions": cover_witnesses,
            }
        )

    output = {
        "schema": 1,
        "description": (
            "Compact exact failing-seed and formerly-passing-cover audit for "
            "the incremental K7 labeled-support sample rejections."
        ),
        "inputs": {
            "sample": str(args.sample),
            "sample_sha256": expected["input_sha256"],
            "decisions": str(args.decisions),
            "decisions_sha256": expected["decisions_sha256"],
            "report": str(args.report),
            "report_sha256": propagation.sha256(args.report),
            "propagation_source_sha256": expected["source_sha256"],
            "frozen_reference_sha256": expected["reference_sha256"],
        },
        "incremental_rejections": len(witnesses),
        "witnesses": witnesses,
    }
    propagation.atomic_json(args.output, output)
    print(
        f"wrote {len(witnesses)} exact rejection audits to {args.output}; "
        f"SHA-256 {propagation.sha256(args.output)}"
    )


if __name__ == "__main__":
    main()
