#!/usr/bin/env python3
"""Extract the complete 1,098-graph input to the K6 bipartite-rank filter.

The parent residue contains all 1,106 survivors of the pure K6 Lorentz
filter.  The independently recorded actual-support refinement rejects eight
of them.  This script joins the two tracked artifacts by corpus index, checks
their complete hashes and decision coverage, and preserves the original
adjacency rows of precisely the 1,098 joint survivors.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


EXPECTED_PARENT_GRAPHS = 1_106
EXPECTED_GRAPHS = 1_098
EXPECTED_PARENT_SHA256 = (
    "3d207db88e66f022258edc79a4684757119b59ea9fda17a00c6e44ab6917825b"
)
EXPECTED_SUPPORT_REPORT_SHA256 = (
    "8f344566b9fa84d042a05712bf48892f013510740141471ac16e2f186fb95c14"
)
EXPECTED_REJECTED = {
    132_876,
    550_637,
    1_323_686,
    2_640_615,
    2_642_489,
    3_107_353,
    3_819_932,
    3_952_058,
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def build_input(parent_path: Path, report_path: Path) -> dict:
    parent_hash = sha256(parent_path)
    if parent_hash != EXPECTED_PARENT_SHA256:
        raise ValueError(f"unexpected parent residue SHA-256 {parent_hash}")
    report_hash = sha256(report_path)
    if report_hash != EXPECTED_SUPPORT_REPORT_SHA256:
        raise ValueError(f"unexpected support report SHA-256 {report_hash}")

    with parent_path.open(encoding="utf-8") as stream:
        parent = json.load(stream)
    with report_path.open(encoding="utf-8") as stream:
        report = json.load(stream)

    graphs = parent["graphs"]
    if len(graphs) != EXPECTED_PARENT_GRAPHS:
        raise ValueError(f"parent residue has {len(graphs)} graphs")
    parent_indices = [record["index"] for record in graphs]
    if parent_indices != sorted(parent_indices) or len(set(parent_indices)) != len(
        parent_indices
    ):
        raise ValueError("parent indices are not unique and increasing")

    decisions = report["graph_results"]
    if len(decisions) != EXPECTED_PARENT_GRAPHS:
        raise ValueError(f"support report has {len(decisions)} graph decisions")
    decision_indices = [record["index"] for record in decisions]
    if decision_indices != parent_indices:
        raise ValueError("support report does not cover the parent in exact order")
    rejected = {
        record["index"]
        for record in decisions
        if record["decision"]["rejected"]
    }
    if rejected != EXPECTED_REJECTED:
        raise ValueError(f"unexpected support rejection set {sorted(rejected)}")
    if report["graphs_rejected"] != len(rejected):
        raise ValueError("support aggregate rejection count disagrees with decisions")

    survivors = [record for record in graphs if record["index"] not in rejected]
    if len(survivors) != EXPECTED_GRAPHS:
        raise ValueError(f"extracted {len(survivors)} survivors")
    return {
        "schema": 1,
        "description": (
            "Complete 1,098-graph K6-only residue after the pure Lorentz "
            "filter and the joint actual-support intersection refinement."
        ),
        "sources": {
            "parent_K6_residue": {
                "file": parent_path.name,
                "sha256": parent_hash,
                "graphs": EXPECTED_PARENT_GRAPHS,
            },
            "actual_support_report": {
                "file": report_path.name,
                "sha256": report_hash,
                "graphs_rejected": len(rejected),
                "graphs_surviving": EXPECTED_GRAPHS,
            },
            **parent["sources"],
        },
        "selection": {
            "rule": "not rejected by K6_joint_actual_support_intersection_CSP",
            "complete_not_sampled": True,
            "indices_in_increasing_corpus_order": True,
            "removed_indices": sorted(rejected),
        },
        "preconditions": {
            **parent["preconditions"],
            "pure_K6_Lorentz_CSP_passed": True,
            "joint_actual_support_intersection_CSP_passed": True,
        },
        "graphs": survivors,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--parent", type=Path, default=Path("d6_k6_support_residue.json")
    )
    parser.add_argument(
        "--support-report", type=Path, default=Path("d6_k6_support_report.json")
    )
    parser.add_argument(
        "--output", type=Path, default=Path("d6_k6_bipartite_rank_input.json")
    )
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = build_input(args.parent, args.support_report)
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.check:
        if args.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit(f"FAIL: {args.output} is not reproducible")
        print(f"PASS: {args.output} is reproducible")
    else:
        args.output.write_text(rendered, encoding="utf-8")
        print(f"wrote {args.output}: {len(result['graphs'])} complete survivors")


if __name__ == "__main__":
    main()
