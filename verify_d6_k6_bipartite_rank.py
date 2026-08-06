#!/usr/bin/env python3
"""Verify provenance, aggregate arithmetic, and exact K6 rank decisions."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from d6_k6_bipartite_rank_reference import evaluate_rank_graph


ROOT = Path(__file__).resolve().parent
EXPECTED_REPORT_SHA256 = (
    "ede04a71a7f9b6901bbc1e4c198a117797bfd90d4d1366912ca9efd1ca559bdc"
)
EXPECTED_INPUT_SHA256 = (
    "dfee979e6ca2759e80cf803c3054c8954d98543875d467644b91ee434c22e845"
)
EXPECTED_REJECTED = [461_363]
COUNT_KEYS = (
    "seeds_checked",
    "impossible_seeds",
    "z0_subsets_considered",
    "z0_support_matchable",
    "z0_bipartite_rank_passed",
    "bipartite_components_checked",
    "bipartite_rank_failures",
    "pure_colorings_considered",
    "pure_colorings_matchable",
    "support_searches",
    "support_dfs_nodes",
    "zero_forcing_initial_sets_checked",
    "zero_forcing_cache_entries",
    "zero_forcing_cache_hits",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_file_hashes(report_path: Path, input_path: Path, report: dict) -> None:
    if sha256(report_path) != EXPECTED_REPORT_SHA256:
        raise AssertionError("tracked report SHA-256 changed")
    if sha256(input_path) != EXPECTED_INPUT_SHA256:
        raise AssertionError("tracked input SHA-256 changed")
    if report["inputs"]["input"]["sha256"] != EXPECTED_INPUT_SHA256:
        raise AssertionError("report does not bind the tracked input")

    bound = [report["inputs"]["source"], *report["inputs"]["dependencies"]]
    for record in bound:
        path = ROOT / record["file"]
        if sha256(path) != record["sha256"]:
            raise AssertionError(f"source/dependency hash changed: {path.name}")


def verify_aggregates(sample: dict, report: dict) -> dict[int, list[int]]:
    graphs = sample["graphs"]
    results = report["graph_results"]
    input_indices = [record["index"] for record in graphs]
    result_indices = [record["index"] for record in results]
    if input_indices != result_indices:
        raise AssertionError("report decisions do not cover the input in order")
    if len(graphs) != 1_098 or report["input_graphs"] != 1_098:
        raise AssertionError("unexpected complete-input size")

    rejected = [
        record["index"]
        for record in results
        if record["decision"]["rejected"]
    ]
    if rejected != EXPECTED_REJECTED or rejected != report["rejected_indices"]:
        raise AssertionError("unexpected rejection decision set")
    if report["graphs_rejected"] != 1 or report["graphs_surviving"] != 1_097:
        raise AssertionError("unexpected rejection aggregates")
    for key in COUNT_KEYS:
        total = sum(record["decision"][key] for record in results)
        if total != report[key]:
            raise AssertionError(f"aggregate {key} is {report[key]}, sum is {total}")

    return {record["index"]: record["adjacency"] for record in graphs}


def verify_decisions(
    adjacency: dict[int, list[int]], report: dict, full: bool
) -> None:
    recorded = {
        record["index"]: record["decision"] for record in report["graph_results"]
    }
    indices = list(adjacency) if full else EXPECTED_REJECTED
    for index in indices:
        actual = evaluate_rank_graph(adjacency[index], scan_all=True).jsonable()
        witness = actual.pop("first_witness")
        actual["first_impossible_seed"] = (
            None if witness is None else witness["seed"]
        )
        if actual != recorded[index]:
            raise AssertionError(f"recomputed decision differs at index {index}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--report", type=Path, default=ROOT / "d6_k6_bipartite_rank_report.json"
    )
    parser.add_argument(
        "--input", type=Path, default=ROOT / "d6_k6_bipartite_rank_input.json"
    )
    parser.add_argument(
        "--full",
        action="store_true",
        help="recompute all 1,098 graph decisions rather than the rejection witness",
    )
    args = parser.parse_args()
    with args.report.open(encoding="utf-8") as stream:
        report = json.load(stream)
    with args.input.open(encoding="utf-8") as stream:
        sample = json.load(stream)
    verify_file_hashes(args.report, args.input, report)
    adjacency = verify_aggregates(sample, report)
    verify_decisions(adjacency, report, args.full)
    scope = "all 1,098 graph decisions" if args.full else "the rejection witness"
    print(f"PASS: K6 bipartite-rank provenance, aggregates, and {scope}")


if __name__ == "__main__":
    main()
