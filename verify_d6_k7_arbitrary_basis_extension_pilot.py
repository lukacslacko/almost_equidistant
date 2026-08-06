#!/usr/bin/env python3
"""Independent exact verifier for the fixed-clique basis-extension pilot."""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

import d6_k7_arbitrary_basis_psd_dual as dual
import run_d6_k7_arbitrary_basis_extension_pilot as pilot


ROOT = Path(__file__).resolve().parent


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_report(report_path: Path, expected_sha256: str) -> dict:
    observed_sha256 = sha256(report_path)
    if observed_sha256 != expected_sha256:
        raise ValueError(
            f"report hash {observed_sha256} != explicit pin {expected_sha256}"
        )
    report = json.loads(report_path.read_text(encoding="utf-8"))
    if (
        report.get("schema") != 1
        or report.get("kind")
        != "d6_k7_arbitrary_basis_fixed_maximum_clique_extension_pilot"
        or report.get("status") != "COMPLETE"
    ):
        raise ValueError("unexpected report schema/kind/status")

    provenance = report["provenance"]
    rank_input = ROOT / provenance["rank_input"]["path"]
    selection = ROOT / provenance["selection"]["path"]
    union = ROOT / provenance["union"]["path"]
    graph_by_index, graph_indices, witnesses, replayed_provenance = (
        pilot.load_campaign_inputs(rank_input, selection, union)
    )
    for name in ("rank_input", "selection", "prior_certificates", "union"):
        if replayed_provenance[name]["sha256"] != provenance[name]["sha256"]:
            raise ValueError(f"provenance mismatch for {name}")

    all_targets = pilot.extract_no_near_covers(
        graph_by_index, graph_indices, witnesses
    )
    if len(all_targets) != report["population"]["no_near_covers"]:
        raise ValueError("no-near cover population mismatch")
    selected_count = int(report["configuration"]["selected_cover_prefix"])
    expected_targets = all_targets[:selected_count]
    observed_targets = report["covers"]
    if len(observed_targets) != len(expected_targets):
        raise ValueError("selected cover count mismatch")

    aggregate = Counter()
    direct_certificates = 0
    dual_certificates = 0
    unresolved_bases = 0
    for ordinal, (expected, observed) in enumerate(
        zip(expected_targets, observed_targets)
    ):
        for field in (
            "graph_index", "seed", "zmask", "nvertices", "graph_n",
            "rank_upper", "maximum_clique_size", "fixed_maximum_clique",
        ):
            if observed[field] != expected[field]:
                raise ValueError(
                    f"cover {ordinal} replay mismatch in field {field}"
                )
        if int(observed["ordinal"]) != ordinal:
            raise ValueError("cover ordinal mismatch")
        graph_n = tuple(int(row) for row in observed["graph_n"])
        clique = tuple(int(vertex) for vertex in observed["fixed_maximum_clique"])
        clique_mask = sum(1 << vertex for vertex in clique)
        if any(
            not (graph_n[first] & (1 << second))
            for first in clique for second in clique if first < second
        ):
            raise ValueError("fixed maximum clique is not a required clique")
        # Maximality is independently exhaustive at n=12.
        larger = next(
            pilot.prior.clique_masks(graph_n, len(clique) + 1), None
        )
        if larger is not None:
            raise ValueError("fixed clique is not maximum")
        if clique_mask != next(
            pilot.prior.clique_masks(graph_n, len(clique))
        ):
            raise ValueError("fixed clique is not the deterministic first maximum clique")

        expected_cores = pilot.extension_cores(expected)
        if int(observed["extension_core_count"]) != len(expected_cores):
            raise ValueError("extension core count mismatch")
        if observed["extension_cores_sha256"] != pilot.stable_hash(
            [list(core) for core in expected_cores]
        ):
            raise ValueError("extension core universe hash mismatch")
        records = observed["basis_records"]
        observed_cores = [tuple(int(vertex) for vertex in item["core"]) for item in records]
        if observed_cores != list(expected_cores):
            raise ValueError("basis records do not exhaust the extension universe")

        cover_counts = Counter()
        for core, record in zip(expected_cores, records):
            if int(record["cover_ordinal"]) != ordinal:
                raise ValueError("basis cover ordinal mismatch")
            if int(record["rank"]) != len(core):
                raise ValueError("basis rank mismatch")
            system = dual.basis_affine_system(graph_n, core)
            outcome = record["outcome"]
            certificate = record["certificate"]
            if outcome == "EASY_REJECTED":
                if certificate is None:
                    raise ValueError("missing easy certificate")
                pilot.verify_easy_basis_certificate(system, certificate)
                # The search uses the first direct contradiction in a fixed order.
                if pilot.easy_basis_certificate(system) != certificate:
                    raise ValueError("easy certificate is not canonical")
                direct_certificates += 1
            elif outcome == "DUAL_REJECTED":
                if certificate is None:
                    raise ValueError("missing dual certificate")
                if pilot.easy_basis_certificate(system) is not None:
                    raise ValueError("dual certificate masks an easy contradiction")
                dual.verify_certificate(system, certificate)
                dual_certificates += 1
            elif outcome == "UNRESOLVED":
                if certificate is not None:
                    raise ValueError("unresolved basis carries a certificate")
                if pilot.easy_basis_certificate(system) is not None:
                    raise ValueError("unresolved basis has a direct contradiction")
                unresolved_bases += 1
            else:
                raise ValueError(f"unknown basis outcome {outcome!r}")
            cover_counts[outcome] += 1
            aggregate[outcome] += 1

        if observed["basis_outcome_counts"] != dict(cover_counts):
            raise ValueError("cover outcome counts mismatch")
        rejected = cover_counts["UNRESOLVED"] == 0
        if bool(observed["cover_rejected"]) != rejected:
            raise ValueError("cover rejection flag mismatch")
        expected_by_rank = {}
        for rank in range(len(clique), int(observed["rank_upper"]) + 1):
            expected_by_rank[str(rank)] = dict(Counter(
                item["outcome"] for item in records if int(item["rank"]) == rank
            ))
        if observed["basis_outcome_counts_by_rank"] != expected_by_rank:
            raise ValueError("cover rank-stratified counts mismatch")

    summary = report["summary"]
    if int(summary["basis_candidates"]) != sum(aggregate.values()):
        raise ValueError("aggregate basis count mismatch")
    if summary["basis_outcome_counts"] != dict(aggregate):
        raise ValueError("aggregate outcome counts mismatch")
    rejected_covers = sum(bool(item["cover_rejected"]) for item in observed_targets)
    if int(summary["covers_rejected"]) != rejected_covers:
        raise ValueError("aggregate rejected-cover count mismatch")
    if int(summary["covers_unresolved"]) != len(observed_targets) - rejected_covers:
        raise ValueError("aggregate unresolved-cover count mismatch")

    return {
        "schema": 1,
        "kind": "d6_k7_arbitrary_basis_extension_pilot_verification",
        "status": "PASS",
        "report": {"path": str(report_path), "sha256": observed_sha256},
        "checked": {
            "full_no_near_cover_population": len(all_targets),
            "selected_covers": len(observed_targets),
            "basis_candidates": sum(aggregate.values()),
            "direct_certificates": direct_certificates,
            "rational_psd_atom_certificates": dual_certificates,
            "unresolved_bases": unresolved_bases,
            "covers_rejected": rejected_covers,
        },
        "claim": (
            "Every rejection certificate was rechecked exactly. Unresolved bases "
            "and covers carry no non-realizability or realizability conclusion."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--report", type=Path,
        default=Path("d6_k7_arbitrary_basis_extension_pilot_report.json"),
    )
    parser.add_argument("--report-sha256", required=True)
    parser.add_argument(
        "--output", type=Path,
        default=Path("d6_k7_arbitrary_basis_extension_pilot_verification.json"),
    )
    args = parser.parse_args()
    result = verify_report(args.report, args.report_sha256)
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "output": str(args.output),
        "sha256": sha256(args.output),
        "status": result["status"],
        "checked": result["checked"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
