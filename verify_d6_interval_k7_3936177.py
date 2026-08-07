#!/usr/bin/env python3
"""Independently verify the focused interval replay of K7 graph 3936177.

The checker imports neither launch wrapper.  It reconstructs the target
directly from the pinned v6 manifest, audits the committed-source boundary,
checks the checkpoint/report/TSV accounting, and replays all 24 winning
slices through the frozen outward-rounded kernel.  The report passes only if
the single target is ``KILLED``; no exploratory witness is trusted.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import time
from collections import Counter
from pathlib import Path

import cdriver6
import verify_d6_interval_benchmarks as common


ROOT = Path(__file__).resolve().parent
TARGET_INDEX = 3_936_177
TARGET_INDICES_SHA256 = (
    "26e5f40303a16a85deb1f94cadaf3af5826c55d87cd42a01cde5e897a23e0106"
)
EXPECTED_MANIFEST_SHA256 = (
    "c513bcc40c037af23340dce3cc91d593a101dd1d073f93feab82abd556efadc9"
)
EXPECTED_VERIFICATION_SHA256 = (
    "3d11cbd20cb59b42a8c8fd5252fb3ef3dea649576982ea8017622fad3279a7bc"
)
EXPECTED_COMMON_SHA256 = (
    "84ec7ae7e682de660329fd1fb02bd88d774bde542032ece95a7eb7d3c54ece4a"
)
EXPECTED_SOURCES = {
    "focused_wrapper": (
        "run_d6_interval_k7_3936177.py",
        "bb129161b88efe8609ec24c42666578cd74a67fa11bb12774266b178015c1309",
    ),
    "v6_wrapper": (
        "run_d6_interval_v6.py",
        "509ec8055542bb60863756d1ce64a24cccfe57ccc9f56f13925edf3290ac38be",
    ),
    "interval_runner": (
        "run_d6_interval_residue.py",
        "58fbbbc548ec3cc6b9d5b231a6790249b02a75d7bfd82f1afcd7fc77fce59404",
    ),
    "cdriver6.py": (
        "cdriver6.py",
        "8f65948c8bc95242d0fde50a8cbcb798cb0e225b78a8ce626382d00135b95173",
    ),
    "ckernel6.c": (
        "ckernel6.c",
        "383ef7a17c328869c338f64a16f2b90875ff968fe2fc20e9a5c3094c830e0874",
    ),
    "ival.py": (
        "ival.py",
        "524e41e0d0637a5352c59ec998009b968f2c9f5e7f31b7f7927e8bc36c2e6fa0",
    ),
}
EXPECTED_INPUTS = {
    "v6_manifest": (
        "d6_current_residue_manifest_v6.json",
        EXPECTED_MANIFEST_SHA256,
    ),
    "v6_verification": (
        "d6_current_residue_manifest_v6_verification.json",
        EXPECTED_VERIFICATION_SHA256,
    ),
}
EXPECTED_COUNTS = {"K7": 19, "K6_only": 625}
EXPECTED_CLASS_HASHES = {
    "K7": "af25b842ca9eaa6e9721f51d0400be2436d930b72cc33c72bc366f0ac5e41290",
    "K6_only": "04d9475217e294efed0c6eac34a7fad88616c90ad7d37a2525776a713e6dd0e1",
}
EXPECTED_COMBINED_ORDER_SHA256 = (
    "907b6fe460846ab613e7d35a94dbfc6eb1eb94322a36f10b4500a8ce4fcc6ccf"
)
EXPECTED_TRUST = {
    "basic_arithmetic": (
        "IEEE-754 binary64 basic operations and sqrt are correctly rounded; "
        "each interval endpoint is expanded with nextafter."
    ),
    "transcendentals": (
        "macOS libm cos endpoint values are assumed within 8 ulps; the "
        "kernel pads both directions by 8 nextafter steps and includes all "
        "interior extrema."
    ),
    "candidate_nonedges": "unconstrained and allowed to be unit distance",
}
EXPECTED_SEARCH = {
    "orders": 4,
    "cap": 2_000_000,
    "slices": 24,
    "zero_circle_only": False,
    "include_bulk_order": False,
}
EXPECTED_SELECTOR = {
    "kind": "exact_v6_K7_corpus_index",
    "index": TARGET_INDEX,
    "population": "K7",
}
RESULT_STATUSES = common.RESULT_STATUSES


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def file_sha256(path: Path) -> str:
    return common.sha256(path)


def object_sha256(value: object) -> str:
    return common.stable_hash(value)


def load_json(path: Path) -> dict:
    return common.load_json(path)


def reconstruct_target(configuration: dict) -> list[dict]:
    """Rebuild all v6 records, then select the exact target by corpus index."""

    manifest_path = ROOT / EXPECTED_INPUTS["v6_manifest"][0]
    verification_path = ROOT / EXPECTED_INPUTS["v6_verification"][0]
    require(file_sha256(manifest_path) == EXPECTED_MANIFEST_SHA256, "manifest hash")
    require(
        file_sha256(verification_path) == EXPECTED_VERIFICATION_SHA256,
        "manifest verification hash",
    )
    manifest = load_json(manifest_path)
    verification = load_json(verification_path)
    require(
        manifest.get("schema") == "d6-current-certified-residue-v6"
        and manifest.get("status") == "COMPLETE_MIXED_CERTIFICATE_UNION",
        "manifest status/schema",
    )
    require(manifest.get("class_order") == ["K7", "K6_only"], "class order")
    require(
        verification.get("schema")
        == "d6-current-certified-residue-v6-verification-v1"
        and verification.get("status") == "PASS",
        "verification status/schema",
    )
    require(
        verification.get("manifest", {}).get("sha256")
        == EXPECTED_MANIFEST_SHA256,
        "verification binding",
    )
    require(
        verification.get("counts", {}).get("combined") == 644
        and verification.get("counts", {}).get("exact_rejections_from_v5") == 264
        and verification.get("counts", {}).get("interval_incremental_rejections")
        == 3,
        "verification counts",
    )

    selection_layer = configuration["selection"]["selection_layer"]
    require(
        selection_layer.get("kind") == "mixed_certificate_d6_current_residue_v6",
        "selection kind",
    )
    require(
        selection_layer.get("manifest", {}).get("sha256")
        == EXPECTED_MANIFEST_SHA256,
        "selection manifest hash",
    )
    require(
        selection_layer.get("verification", {}).get("sha256")
        == EXPECTED_VERIFICATION_SHA256
        and selection_layer.get("verification", {}).get("status") == "PASS",
        "selection verification root",
    )
    require(selection_layer.get("counts") == EXPECTED_COUNTS, "selection counts")
    require(selection_layer.get("combined") == 644, "selection combined count")
    require(
        selection_layer.get("ordered_indices_sha256")
        == EXPECTED_COMBINED_ORDER_SHA256,
        "selection combined hash",
    )
    accounting = selection_layer.get("certificate_accounting", {})
    require(
        accounting.get("exact_rejections_from_v5") == 264
        and accounting.get("interval_incremental_rejections") == 3
        and accounting.get("interval_non_KILLED_used") is False
        and accounting.get("interval_trust_assumptions") == EXPECTED_TRUST,
        "selection trust tiers",
    )
    manifest_accounting = manifest.get("certificate_accounting", {})
    require(
        manifest_accounting.get("exact_algebra_and_graph_logic", {}).get(
            "floating_point_enters_rejection"
        )
        is False
        and manifest_accounting.get("interval", {}).get(
            "non_KILLED_statuses_used_for_rejection"
        )
        is False
        and manifest_accounting.get("interval", {}).get("trust_assumptions")
        == EXPECTED_TRUST,
        "manifest trust tiers",
    )

    all_graphs: list[dict] = []
    ordered_indices: list[int] = []
    seen: set[int] = set()
    for class_name, population in (("K7", "K7"), ("K6_only", "K6")):
        block = manifest["classes"][class_name]
        records = block["graphs"]
        indices = block["indices"]
        require(
            len(records) == len(indices) == EXPECTED_COUNTS[class_name],
            "class count",
        )
        require(indices == [record["index"] for record in records], "class order")
        require(object_sha256(indices) == block["indices_sha256"], "class index hash")
        require(
            block["indices_sha256"] == EXPECTED_CLASS_HASHES[class_name],
            "frozen class hash",
        )
        require(object_sha256(records) == block["graphs_sha256"], "class graph hash")
        for record in records:
            index = int(record["index"])
            require(index not in seen, "duplicate index")
            seen.add(index)
            adjacency = tuple(map(int, record["adjacency"]))
            common.validate_adjacency(adjacency)
            all_graphs.append(
                {"index": index, "population": population, "adjacency": adjacency}
            )
            ordered_indices.append(index)
    require(len(all_graphs) == 644, "base graph count")
    require(
        object_sha256(ordered_indices)
        == manifest["combined"]["ordered_indices_sha256"]
        == EXPECTED_COMBINED_ORDER_SHA256,
        "base order hash",
    )

    require(configuration.get("selector") == EXPECTED_SELECTOR, "exact selector")
    graphs = [graph for graph in all_graphs if graph["index"] == TARGET_INDEX]
    require(len(graphs) == 1, "target multiplicity")
    require(graphs[0]["population"] == "K7", "target K7 class")
    graphs[0]["ordinal"] = 0
    indices = [graph["index"] for graph in graphs]
    selection = configuration["selection"]
    require(selection.get("graphs") == 1, "selected graph count")
    require(selection.get("target_index") == TARGET_INDEX, "selected target index")
    require(selection.get("indices_sha256") == object_sha256(indices), "target hash")
    require(selection.get("indices_sha256") == TARGET_INDICES_SHA256, "frozen target hash")
    require(selection.get("population_counts") == {"K7": 1}, "target population")
    return graphs


def validate_source_boundary(configuration: dict) -> dict:
    require(
        file_sha256(ROOT / "verify_d6_interval_benchmarks.py")
        == EXPECTED_COMMON_SHA256,
        "common verifier helper changed",
    )
    sources = configuration["sources"]
    for key, (filename, expected) in EXPECTED_SOURCES.items():
        require(sources.get(key) == expected, f"recorded source hash: {key}")
        require(file_sha256(ROOT / filename) == expected, f"working source hash: {key}")
    require(
        sources.get(cdriver6.LIBNAME) == file_sha256(ROOT / cdriver6.LIBNAME),
        "kernel binary hash",
    )
    require(configuration.get("trust_assumptions") == EXPECTED_TRUST, "trust assumptions")

    inputs = configuration.get("inputs")
    require(isinstance(inputs, dict), "missing input boundary")
    for key, (filename, expected) in EXPECTED_INPUTS.items():
        record = inputs.get(key, {})
        require(record.get("sha256") == expected, f"recorded input hash: {key}")
        recorded_path = record.get("path")
        require(
            isinstance(recorded_path, str)
            and Path(recorded_path).resolve() == (ROOT / filename).resolve(),
            f"recorded input path: {key}",
        )
        require(file_sha256(ROOT / filename) == expected, f"working input hash: {key}")

    git = configuration["git"]
    require(git.get("available") is True, "git provenance unavailable")
    require(git.get("branch") == "codex/dimension6", "launch branch")
    commit = git.get("commit")
    require(isinstance(commit, str) and len(commit) == 40, "launch commit")
    checker_hash = file_sha256(Path(__file__).resolve())
    committed_files = [
        *EXPECTED_SOURCES.values(),
        *EXPECTED_INPUTS.values(),
        ("verify_d6_interval_benchmarks.py", EXPECTED_COMMON_SHA256),
        (Path(__file__).name, checker_hash),
    ]
    for filename, expected in committed_files:
        blob = subprocess.run(
            ["git", "show", f"{commit}:{filename}"],
            cwd=ROOT,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        ).stdout
        require(
            hashlib.sha256(blob).hexdigest() == expected,
            f"commit blob: {filename}",
        )
    lines = git.get("porcelain_lines")
    require(isinstance(lines, list), "porcelain lines")
    porcelain = "\n".join(lines)
    require(
        hashlib.sha256(porcelain.encode("utf-8")).hexdigest()
        == git.get("porcelain_sha256"),
        "porcelain hash",
    )
    require(git.get("dirty") == bool(lines), "dirty flag")
    require(
        all(isinstance(line, str) and line.startswith("?? ") for line in lines),
        "tracked file was dirty at launch",
    )
    boundary = configuration.get("source_boundary")
    require(
        boundary
        == {
            "status": "PASS",
            "commit": commit,
            "branch": "codex/dimension6",
            "tracked_files_clean": True,
            "committed_source_and_input_blobs_match": True,
            "porcelain_sha256": git["porcelain_sha256"],
        },
        "runner source boundary",
    )
    return {
        "commit": commit,
        "branch": git["branch"],
        "tracked_sources_clean_at_launch": True,
        "committed_source_and_input_blobs_match": True,
        "porcelain_sha256": git["porcelain_sha256"],
        "verifier_commit_blob_sha256": checker_hash,
        "common_helper_commit_blob_sha256": EXPECTED_COMMON_SHA256,
    }


def verify_report(report_path: Path, expected_hash: str | None, replay: bool) -> dict:
    report_path = report_path.resolve()
    report_hash = file_sha256(report_path)
    if expected_hash is not None:
        require(report_hash == expected_hash, "report hash")
    report = load_json(report_path)
    require(report.get("schema") == 2, "report schema")
    configuration = report["configuration"]
    require(configuration.get("schema") == 1, "configuration schema")
    require(object_sha256(configuration) == report["config_sha256"], "config hash")
    require(configuration.get("search") == EXPECTED_SEARCH, "search parameters")
    require(configuration.get("selector") == EXPECTED_SELECTOR, "selector parameters")
    launch = configuration["launch"]
    require(
        launch.get("outer_launch_command_user_supplied") is True,
        "outer launch command",
    )
    outer_launch = launch.get("outer_launch_command")
    require(
        isinstance(outer_launch, str) and "caffeinate" in outer_launch,
        "sleep prevention not recorded",
    )

    graphs = reconstruct_target(configuration)
    require(report.get("total_graphs") == report.get("completed_graphs") == 1, "coverage")
    require(
        report["session"].get("config_sha256") == report["config_sha256"],
        "session config hash",
    )
    require(report["session"].get("workers") == 1, "worker count")
    require(len(report.get("results", [])) == 1, "result count")
    provenance = validate_source_boundary(configuration)
    control_audit = common.validate_controls(
        configuration["kernel_controls"], replay=replay
    )

    result = report["results"][0]
    graph = graphs[0]
    common.validate_result_structure(result, graph, configuration)
    require(result["status"] == "KILLED", "target was not certified KILLED")
    statuses = Counter(row["status"] for row in report["results"])
    status_counts = {name: statuses[name] for name in RESULT_STATUSES}
    require(report["status_counts"] == status_counts, "status counts")
    require(status_counts["KILLED"] == 1, "focused kill count")
    require(report["certified_killed"] == 1, "certified kill count")
    require(report["unresolved_total"] == 0, "focused residue count")
    require(report["killed_by_population"] == {"K7": 1, "K6": 0}, "kill partition")
    kernel_totals = {
        "KILLED": result["kernel_killed"],
        "SURVIVORS": result["kernel_survivors"],
        "ABORT": result["kernel_aborts"],
    }
    require(report["kernel_status_totals"] == kernel_totals, "kernel totals")

    replayed_slices = 0
    replayed_nodes = 0
    if replay:
        for record in result["winning_slice_records"]:
            observed = cdriver6.decide6(
                graph["adjacency"],
                common.N,
                seed=result["winning_seed"],
                order=result["winning_placement_order"],
                th0=(record["lo"], record["hi"]),
                max_nodes=EXPECTED_SEARCH["cap"],
            )
            expected = (
                record["status"],
                record["nodes"],
                record["unresolved_cells"],
            )
            require(observed == expected, "winning slice replay")
            replayed_slices += 1
            replayed_nodes += observed[1]
    checkpoint = common.validate_checkpoint_directory(report, report_path, graphs)
    witness = {
        key: result[key]
        for key in (
            "index",
            "winning_order",
            "winning_circles",
            "winning_seed",
            "winning_placement_order",
            "winning_slice_records",
        )
    }
    return {
        "schema": "d6-interval-k7-3936177-verification-v1",
        "status": "PASS",
        "claim": (
            "K7 residue graph 3936177 is certified non-realizable in R6 by "
            "the source-bound outward-rounded interval replay."
        ),
        "report": {"path": str(report_path), "sha256": report_hash},
        "verifier_source_sha256": file_sha256(Path(__file__).resolve()),
        "common_helper_sha256": EXPECTED_COMMON_SHA256,
        "provenance": provenance,
        "selection": {
            "graphs": 1,
            "index": TARGET_INDEX,
            "indices_sha256": TARGET_INDICES_SHA256,
            "population_counts": {"K7": 1},
            "v6_manifest_sha256": EXPECTED_MANIFEST_SHA256,
            "v6_verification_sha256": EXPECTED_VERIFICATION_SHA256,
        },
        "search": {**EXPECTED_SEARCH, "workers": 1},
        "status_counts": status_counts,
        "certified_killed": 1,
        "winning_witness": witness,
        "winning_witness_sha256": object_sha256(witness),
        "witness_replay": {
            "status": "PASS" if replay else "SKIPPED",
            "slices": replayed_slices,
            "nodes": replayed_nodes,
            "exact_status_node_cell_match": replay,
        },
        "kernel_controls": control_audit,
        "checkpoint": checkpoint,
        "trust_assumptions": EXPECTED_TRUST,
    }


def atomic_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--expected-report-sha256")
    parser.add_argument("--no-replay", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = verify_report(
        args.report,
        args.expected_report_sha256,
        replay=not args.no_replay,
    )
    atomic_json(args.output.resolve(), result)
    print(
        json.dumps(
            {
                "status": "PASS",
                "index": TARGET_INDEX,
                "certified_killed": result["certified_killed"],
                "output": str(args.output),
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
