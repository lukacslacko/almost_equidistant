#!/usr/bin/env python3
"""Independent, kernel-free audit of the completed d=6 n=18 campaign.

The checker imports neither production runner.  It rebuilds the 181-class
selection from the frozen v7 residue and deletion manifest, checks every
atomic result against that selection, reconstructs the decisions TSV, and
binds the immutable campaign/report/decision hashes.  It deliberately does
not interpret ABORT or UNRESOLVED as a rejection and does not replay the
interval kernel.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import os
import subprocess
import time
from collections import Counter, defaultdict
from datetime import UTC, datetime
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parent
RUN_DIR = ROOT / ".runs/d6_interval_18_cover_v7_cap100000_w3"
RAW_CAMPAIGN = RUN_DIR / "campaign.json"
RAW_REPORT = ROOT / ".runs/d6_interval_18_cover_v7_cap100000_w3_report.json"
RAW_DECISIONS = ROOT / ".runs/d6_interval_18_cover_v7_cap100000_w3_decisions.tsv"

RAW_CAMPAIGN_SHA256 = (
    "70ba451094eb4eb5d039fc059ab8e3c351ebeed157e2eae39cc35fe4dcb3e348"
)
RAW_REPORT_SHA256 = (
    "0a8dc0d875f3f3f547c2feec3992152d2b841c4a4ba4d6b1e2dad97397ab0a53"
)
RAW_DECISIONS_SHA256 = (
    "f7416adec9ad442d4749eef3f5afd3e1a4ec25972267d432302dc8cf68a5aab8"
)
CONFIG_SHA256 = "82c8d332af4ac4e3984949b98cde969f5b239eadb640e9633f704f8ec88fbcba"
CHECKPOINT_INDEX_SHA256 = (
    "cf5807d8408a76dda996ca5bbfaccbb76260d45c2311cce9f0ac3c9989060b5c"
)
CAMPAIGN_SOURCE_COMMIT = "22be775b1dafde30122ab579f66a7a539581ba73"

CAMPAIGN_SOURCES = {
    "cdriver6.py": "8f65948c8bc95242d0fde50a8cbcb798cb0e225b78a8ce626382d00135b95173",
    "ckernel6.c": "383ef7a17c328869c338f64a16f2b90875ff968fe2fc20e9a5c3094c830e0874",
    "ival.py": "524e41e0d0637a5352c59ec998009b968f2c9f5e7f31b7f7927e8bc36c2e6fa0",
    "run_d6_interval_18_cover_v7.py": (
        "5eeca35b328d930c94fa33d96c714ef756f4a777fcce8c8bffccbe472d744019"
    ),
    "run_d6_interval_residue.py": (
        "58fbbbc548ec3cc6b9d5b231a6790249b02a75d7bfd82f1afcd7fc77fce59404"
    ),
}
KERNEL_BINARY_SHA256 = (
    "ec93d78302d6fcf0f3a21ebcdb65333004eb7a406ab4ada8445231536fed07bd"
)
INPUTS = {
    "d6_current_residue_manifest_v7.json": (
        "1ab0948780d73cdd6aca2925107fc240dcae5df51d7cf5393accd8de3ae79f95"
    ),
    "d6_current_residue_manifest_v7_verification.json": (
        "ff02a9ca5bfe769ce77bbf67dedbc8e928561fb4f291a78d2ea577e06e784255"
    ),
    "d6_residue_18_deletions_v2.json": (
        "7a13d7a204f866c2d7c421db8a46c5b042bbad8bad6a2cade4a94eb4d340521d"
    ),
    "d6_residue_18_deletions_v2_verification.json": (
        "99927f8b41b48fb2c0151f3f6d04eca37701861cfef9fce4d3bba3f5dab91090"
    ),
    "d6_standard18_geometry_report.json": (
        "97904946b382eb1e6cc325083e1c1ab20e5e8c426fcbaaed1b79357646c846ef"
    ),
    "d6_standard18_geometry_verification.json": (
        "12e3e1df3919dd465a54b1c470f524a3907db2a19110a34876137eeedeb926d6"
    ),
    "d6_18_exact_reconstruction.json": (
        "69fa1c0be5f46e285773f13cc05270219d606673703fd16b3e36a9551cb9a990"
    ),
    "d6_18_exact_reconstruction_verification.json": (
        "e86f8b538f4f78c8a0ace71a13d2e029be054e2135d19ed710faf4fc21aca887"
    ),
}

V7_COUNTS = {"K7": 12, "K6_only": 251}
V7_INDEX_HASHES = {
    "K7": "e545ea5a8fc9691ed98e885919e149f84dc7e8b4796d0fe9e6cb3a519df6fd8c",
    "K6_only": (
        "a23da7bd9ae21eb7c078a08a7bce92f22d397d5c71f1e4f51d9cf8e59bbfa907"
    ),
}
V7_COMBINED_INDEX_HASH = (
    "a50baf54bfff642e7c1711fe7041b5184df12f90de2cf52557cd086f10875127"
)
GREEDY_INDEX_HASH = "ef4ac7e73ea924d3d9c6382f590dc73cdf2e65b357e86ddc69719372f717d6fe"
PRODUCTION_INDEX_HASH = (
    "8c624c969b517bb5663dc32bdf577ab3addbe8a108458b2803443c17f294d8e6"
)
KNOWN_POSITIVE_CLASS_IDS = {
    "u18-761add09c775ad5ba966a0eef97848e076faacc647744c4079b7bee48a65d92f",
    "u18-84a0a88e385b1aa4673cc149e0c2458f37b7cc6aef3def3cf2589f2075be4841",
}
EXPECTED_BACKUPS = [126, 5_673]
EXPECTED_STATUS_COUNTS = {
    "KILLED": 1,
    "ABORT": 179,
    "UNRESOLVED": 1,
    "INFRA_ERROR": 0,
}
EXPECTED_KERNEL_TOTALS = {"KILLED": 161, "SURVIVORS": 0, "ABORT": 368}
KILLED_CLASS = 2_100
NO_ORDER_CLASS = 7_259
RESULT_STATUSES = frozenset(EXPECTED_STATUS_COUNTS)
KERNEL_STATUSES = frozenset({"KILLED", "SURVIVORS", "ABORT"})
DECISION_FIELDS = (
    "ordinal",
    "index",
    "population",
    "status",
    "winning_order",
    "winning_circles",
    "orders_available",
    "orders_tried",
    "kernel_calls",
    "kernel_killed",
    "kernel_survivors",
    "kernel_aborts",
    "total_nodes",
    "unresolved_cells",
    "elapsed_seconds",
    "error_type",
    "error_message",
)
PACKAGE_FILES = (
    "verify_d6_interval_18_cover_v7_completion.py",
    "test_d6_interval_18_cover_v7_completion.py",
    "d6_interval_18_cover_v7_completion.md",
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def bytes_sha256(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def stable_hash(value: object) -> str:
    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")
    return bytes_sha256(encoded)


def load_json(path: Path) -> dict:
    value = json.loads(path.read_text(encoding="utf-8"))
    require(isinstance(value, dict), f"{path}: expected a JSON object")
    return value


def git_blob_sha256(commit: str, name: str) -> str:
    value = subprocess.run(
        ["git", "show", f"{commit}:{name}"],
        cwd=ROOT,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    ).stdout
    return bytes_sha256(value)


def validate_adjacency(value: object, order: int, label: str) -> tuple[int, ...]:
    require(
        isinstance(value, list)
        and len(value) == order
        and all(type(row) is int and row >= 0 for row in value),
        f"{label}: adjacency shape",
    )
    rows = tuple(value)
    full = (1 << order) - 1
    for vertex, row in enumerate(rows):
        require(not row & ~full and not row & (1 << vertex), f"{label}: range/loop")
        for other in range(vertex):
            require(
                bool(row & (1 << other)) == bool(rows[other] & (1 << vertex)),
                f"{label}: asymmetric adjacency",
            )
    return rows


def parent_sort_key(parent: tuple[str, int]) -> tuple[int, int]:
    return ({"K7": 0, "K6_only": 1}[parent[0]], parent[1])


def tagged(parent: tuple[str, int]) -> str:
    return f"{parent[0]}:{parent[1]}"


def known_positive(record: dict) -> bool:
    return bool(record.get("standard18_compatible")) or record.get(
        "class_id"
    ) in KNOWN_POSITIVE_CLASS_IDS


def load_v7_parents() -> list[tuple[str, int]]:
    manifest = load_json(ROOT / "d6_current_residue_manifest_v7.json")
    verification = load_json(
        ROOT / "d6_current_residue_manifest_v7_verification.json"
    )
    require(
        manifest.get("schema") == "d6-current-certified-residue-v7"
        and manifest.get("status") == "COMPLETE_MIXED_CERTIFICATE_UNION"
        and verification.get("status") == "PASS"
        and verification.get("manifest", {}).get("sha256")
        == INPUTS["d6_current_residue_manifest_v7.json"]
        and all(verification.get("checks", {}).values()),
        "v7 manifest boundary",
    )
    parents: list[tuple[str, int]] = []
    ordered_indices: list[int] = []
    for class_name in ("K7", "K6_only"):
        block = manifest.get("classes", {}).get(class_name, {})
        indices = block.get("indices")
        graphs = block.get("graphs")
        require(
            isinstance(indices, list)
            and isinstance(graphs, list)
            and len(indices) == len(graphs) == V7_COUNTS[class_name]
            and indices == [record.get("index") for record in graphs]
            and stable_hash(indices)
            == block.get("indices_sha256")
            == V7_INDEX_HASHES[class_name],
            f"v7 {class_name} roots",
        )
        for index, graph in zip(indices, graphs, strict=True):
            validate_adjacency(graph.get("adjacency"), 19, f"v7 {index}")
            parents.append((class_name, int(index)))
            ordered_indices.append(int(index))
    require(
        stable_hash(ordered_indices)
        == manifest.get("combined", {}).get("ordered_indices_sha256")
        == V7_COMBINED_INDEX_HASH,
        "v7 combined index root",
    )
    return parents


def rebuild_selection(configuration: dict) -> tuple[list[dict], dict]:
    """Reimplement the dense/greedy/backup selection from frozen inputs."""

    parents = load_v7_parents()
    parent_set = set(parents)
    require(len(parents) == len(parent_set) == 263, "v7 parent set")
    deletion = load_json(ROOT / "d6_residue_18_deletions_v2.json")
    deletion_check = load_json(ROOT / "d6_residue_18_deletions_v2_verification.json")
    records = deletion.get("unique_deletions")
    require(
        deletion.get("schema") == "d6-residue-18-deletion-manifest-v2"
        and isinstance(records, list)
        and len(records) == 11_975
        and deletion_check.get("status") == "PASS"
        and deletion_check.get("manifest", {}).get("sha256")
        == INPUTS["d6_residue_18_deletions_v2.json"]
        and all(deletion_check.get("checks", {}).values()),
        "deletion manifest boundary",
    )
    class_ids = set()
    for index, record in enumerate(records):
        class_id = record.get("class_id")
        require(
            isinstance(class_id, str) and class_id not in class_ids,
            f"deletion {index}: class id",
        )
        class_ids.add(class_id)
        rows = validate_adjacency(record.get("adjacency"), 18, f"deletion {index}")
        require(
            sum(row.bit_count() for row in rows) // 2 == record.get("edges"),
            f"deletion {index}: edge count",
        )

    occurrences: dict[tuple[str, int], list[tuple[int, int, int]]] = defaultdict(list)
    for class_index, record in enumerate(records):
        for occurrence in record.get("occurrences", ()):
            parent = (occurrence.get("parent_class"), occurrence.get("parent_index"))
            if parent in parent_set:
                occurrences[parent].append(
                    (
                        class_index,
                        int(record["edges"]),
                        int(occurrence["deleted_vertex"]),
                    )
                )
    require(set(occurrences) == parent_set, "deletion coverage of v7 parents")

    coverage: dict[int, set[tuple[str, int]]] = defaultdict(set)
    dense_by_parent: dict[tuple[str, int], list[int]] = {}
    for parent in sorted(parent_set, key=parent_sort_key):
        options = occurrences[parent]
        require(
            len(options) == 19 and {item[2] for item in options} == set(range(19)),
            f"rooted deletion partition {tagged(parent)}",
        )
        maximum = max(item[1] for item in options)
        dense = sorted({item[0] for item in options if item[1] == maximum})
        dense_by_parent[parent] = dense
        for class_index in dense:
            coverage[class_index].add(parent)

    uncovered = set(parent_set)
    greedy: list[int] = []
    while uncovered:
        candidates = [
            index for index, covered in coverage.items() if covered & uncovered
        ]
        require(bool(candidates), "greedy cover stalled")
        chosen = max(
            candidates,
            key=lambda index: (
                len(coverage[index] & uncovered),
                records[index]["class_id"],
            ),
        )
        greedy.append(chosen)
        uncovered.difference_update(coverage[chosen])
    require(
        len(greedy) == len(set(greedy)) == 179
        and stable_hash(greedy) == GREEDY_INDEX_HASH,
        "greedy selection root",
    )

    selected_coverers: dict[tuple[str, int], list[int]] = defaultdict(list)
    for index in greedy:
        for parent in coverage[index]:
            selected_coverers[parent].append(index)
    blocked = sorted(
        [
            parent
            for parent, coverers in selected_coverers.items()
            if coverers and all(known_positive(records[index]) for index in coverers)
        ],
        key=parent_sort_key,
    )
    require(
        [tagged(parent) for parent in blocked]
        == ["K6_only:364827", "K6_only:3335955"],
        "blocked positive-only parents",
    )
    backups = []
    for parent in blocked:
        choices = []
        for index in dense_by_parent[parent]:
            if known_positive(records[index]):
                continue
            choices.append((len(coverage[index]), records[index]["class_id"], index))
        require(bool(choices), f"backup for {tagged(parent)}")
        backups.append(max(choices)[2])
    require(backups == EXPECTED_BACKUPS, "backup selection")

    production_indices = sorted(greedy + backups)
    require(
        len(production_indices) == len(set(production_indices)) == 181
        and stable_hash(production_indices) == PRODUCTION_INDEX_HASH,
        "production selection root",
    )
    production = []
    for ordinal, index in enumerate(production_indices):
        active_classes = {parent[0] for parent in coverage[index]}
        require(len(active_classes) == 1, f"class {index}: mixed population")
        population = "K7" if active_classes == {"K7"} else "K6"
        production.append(
            {"ordinal": ordinal, "index": index, "population": population}
        )
    require(
        Counter(item["population"] for item in production)
        == Counter({"K6": 169, "K7": 12}),
        "production population counts",
    )

    selection = configuration.get("selection", {})
    run = selection.get("run_selection", {})
    augmentation = selection.get("actionable_augmentation", {})
    dense_cover = selection.get("full_dense_cover", {})
    profile = selection.get("placement_order_profile", {})
    require(
        run.get("graphs") == 181
        and run.get("class_indices") == production_indices
        and run.get("class_indices_sha256") == PRODUCTION_INDEX_HASH
        and run.get("population_counts") == {"K6": 169, "K7": 12},
        "raw run selection",
    )
    require(
        augmentation.get("production_class_indices") == production_indices
        and augmentation.get("production_class_indices_sha256")
        == PRODUCTION_INDEX_HASH
        and [item.get("class_index") for item in augmentation.get("backup_records", ())]
        == backups,
        "raw actionable augmentation",
    )
    require(
        dense_cover.get("selected_class_indices_greedy_order") == greedy
        and dense_cover.get("selected_class_indices_sha256") == GREEDY_INDEX_HASH,
        "raw greedy cover",
    )
    require(
        profile.get("graphs_profiled") == 181
        and profile.get("graphs_with_orders") == 180
        and profile.get("no_order_class_indices") == [NO_ORDER_CLASS],
        "placement-order profile",
    )
    return production, {
        "v7_parents": len(parents),
        "deletion_classes": len(records),
        "eligible_dense_classes": len(coverage),
        "greedy_classes": len(greedy),
        "greedy_indices_sha256": stable_hash(greedy),
        "blocked_positive_only_parents": [tagged(parent) for parent in blocked],
        "backup_classes": backups,
        "production_classes": len(production),
        "production_population_counts": dict(
            sorted(Counter(item["population"] for item in production).items())
        ),
        "production_indices_sha256": stable_hash(production_indices),
    }


def verify_sources_and_inputs(configuration: dict) -> dict:
    for name, expected in INPUTS.items():
        require(sha256(ROOT / name) == expected, f"input hash: {name}")
    git = configuration.get("git", {})
    sources = configuration.get("sources", {})
    require(
        git.get("branch") == "codex/dimension6"
        and git.get("source_commit") == CAMPAIGN_SOURCE_COMMIT
        and git.get("sources") == dict(sorted(CAMPAIGN_SOURCES.items())),
        "campaign git boundary",
    )
    for name, expected in CAMPAIGN_SOURCES.items():
        require(sources.get(name) == expected, f"recorded source hash: {name}")
        require(sha256(ROOT / name) == expected, f"working source hash: {name}")
        require(
            git_blob_sha256(CAMPAIGN_SOURCE_COMMIT, name) == expected,
            f"committed source hash: {name}",
        )
    require(
        sources.get("ckernel6.dylib") == KERNEL_BINARY_SHA256
        and sha256(ROOT / "ckernel6.dylib") == KERNEL_BINARY_SHA256,
        "kernel binary hash",
    )
    recorded_inputs = configuration.get("inputs", {})
    for key, name in (
        ("v7_manifest", "d6_current_residue_manifest_v7.json"),
        ("v7_verification", "d6_current_residue_manifest_v7_verification.json"),
        ("deletion_manifest", "d6_residue_18_deletions_v2.json"),
        ("deletion_verification", "d6_residue_18_deletions_v2_verification.json"),
    ):
        value = recorded_inputs.get(key, {})
        require(
            Path(str(value.get("path"))).name == name
            and value.get("sha256") == INPUTS[name],
            f"recorded input boundary: {key}",
        )
    positive = configuration.get("kernel_controls", {}).get(
        "positive_source_hashes"
    )
    positive_names = {
        name: INPUTS[name]
        for name in sorted(INPUTS)
        if not name.startswith("d6_current_residue_manifest_v7")
    }
    require(positive == positive_names, "recorded positive-control input hashes")
    require(
        configuration.get("search")
        == {
            "cap": 100_000,
            "include_bulk_order": True,
            "orders": 4,
            "slices": 24,
            "zero_circle_only": False,
        },
        "search configuration",
    )
    return {
        "source_commit": CAMPAIGN_SOURCE_COMMIT,
        "campaign_sources_checked": len(CAMPAIGN_SOURCES),
        "kernel_binary_sha256": KERNEL_BINARY_SHA256,
        "inputs_checked": len(INPUTS),
    }


def validate_slice_records(value: object) -> list[dict]:
    require(isinstance(value, list) and len(value) == 24, "winning slice records")
    records = []
    for part, record in enumerate(value):
        require(isinstance(record, dict), f"slice {part}: record")
        lo = 2.0 * math.pi * part / 24
        hi = 2.0 * math.pi * (part + 1) / 24
        require(
            record.get("part") == part
            and record.get("lo") == lo
            and record.get("hi") == hi
            and record.get("status") == "KILLED"
            and type(record.get("nodes")) is int
            and record["nodes"] >= 0
            and type(record.get("unresolved_cells")) is int
            and record["unresolved_cells"] >= 0,
            f"slice {part}: fields",
        )
        records.append(record)
    return records


def validate_result(result: dict, expected: dict) -> None:
    require(
        result.get("ordinal") == expected["ordinal"]
        and result.get("index") == expected["index"]
        and result.get("population") == expected["population"],
        f"result identity at ordinal {expected['ordinal']}",
    )
    status = result.get("status")
    require(status in RESULT_STATUSES, f"class {expected['index']}: status")
    integer_fields = (
        "winning_order",
        "winning_circles",
        "orders_available",
        "orders_tried",
        "kernel_calls",
        "kernel_killed",
        "kernel_survivors",
        "kernel_aborts",
        "total_nodes",
        "unresolved_cells",
    )
    require(
        all(type(result.get(field)) is int for field in integer_fields),
        f"class {expected['index']}: integer counters",
    )
    require(
        result["orders_available"] >= result["orders_tried"] >= 0
        and result["kernel_calls"]
        == result["kernel_killed"]
        + result["kernel_survivors"]
        + result["kernel_aborts"]
        and all(
            result[field] >= 0
            for field in (
                "kernel_calls",
                "kernel_killed",
                "kernel_survivors",
                "kernel_aborts",
                "total_nodes",
                "unresolved_cells",
            )
        )
        and isinstance(result.get("elapsed_seconds"), (int, float))
        and math.isfinite(result["elapsed_seconds"])
        and result["elapsed_seconds"] >= 0,
        f"class {expected['index']}: counters",
    )
    if status == "KILLED":
        slices = validate_slice_records(result.get("winning_slice_records"))
        require(
            result["winning_order"] >= 0
            and result["winning_circles"] >= 0
            and result["orders_tried"] == result["winning_order"] + 1
            and result["kernel_killed"] >= len(slices)
            and isinstance(result.get("winning_seed"), list)
            and len(result["winning_seed"]) == 7
            and isinstance(result.get("winning_placement_order"), list)
            and len(result["winning_placement_order"]) == 11,
            f"class {expected['index']}: killed witness",
        )
    else:
        require(
            result.get("winning_order") == -1
            and result.get("winning_circles") == -1
            and result.get("winning_seed") is None
            and result.get("winning_placement_order") is None
            and result.get("winning_slice_records") is None,
            f"class {expected['index']}: non-killed witness fields",
        )
    if status == "INFRA_ERROR":
        require(result.get("error_type"), f"class {expected['index']}: infra detail")
    else:
        require(
            result.get("error_type") is None
            and result.get("error_message") is None
            and result.get("error_traceback") is None,
            f"class {expected['index']}: unexpected error detail",
        )


def decision_text(results: Iterable[dict]) -> str:
    lines = ["\t".join(DECISION_FIELDS) + "\n"]
    for result in results:
        values = []
        for field in DECISION_FIELDS:
            value = result.get(field)
            if field == "elapsed_seconds":
                value = f"{float(value):.9f}"
            if value is None:
                value = ""
            values.append(str(value).replace("\t", " ").replace("\n", " "))
        lines.append("\t".join(values) + "\n")
    return "".join(lines)


def verify_results(
    report: dict, production: list[dict], checkpoint_dir: Path, decisions_path: Path
) -> dict:
    report_results = report.get("results")
    require(
        isinstance(report_results, list) and len(report_results) == len(production) == 181,
        "report result count",
    )
    result_dir = checkpoint_dir / "results"
    actual_files = sorted(result_dir.glob("result_*.json"))
    require(len(actual_files) == 181, "checkpoint result-file count")
    checkpoint_entries = []
    checkpoint_results = []
    for expected, report_result in zip(production, report_results, strict=True):
        name = f"result_{expected['ordinal']:06d}_{expected['index']:07d}.json"
        path = result_dir / name
        require(path in actual_files, f"missing checkpoint {name}")
        checkpoint = load_json(path)
        require(
            checkpoint.get("schema") == 1
            and checkpoint.get("config_sha256") == CONFIG_SHA256
            and checkpoint.get("search_slices") == 24
            and isinstance(checkpoint.get("completed_utc"), str)
            and checkpoint.get("result") == report_result,
            f"checkpoint envelope {name}",
        )
        validate_result(report_result, expected)
        checkpoint_results.append(report_result)
        checkpoint_entries.append([name, sha256(path)])
    require(
        stable_hash(checkpoint_entries)
        == report.get("checkpoint_result_index_sha256")
        == CHECKPOINT_INDEX_SHA256,
        "checkpoint index root",
    )

    statuses = Counter(result["status"] for result in checkpoint_results)
    status_counts = {status: statuses[status] for status in EXPECTED_STATUS_COUNTS}
    require(status_counts == EXPECTED_STATUS_COUNTS, "result status accounting")
    killed = [result for result in checkpoint_results if result["status"] == "KILLED"]
    unresolved = [
        result for result in checkpoint_results if result["status"] == "UNRESOLVED"
    ]
    require(
        len(killed) == 1
        and killed[0]["index"] == KILLED_CLASS
        and killed[0]["population"] == "K7",
        "sole certified KILLED result",
    )
    require(
        len(unresolved) == 1
        and unresolved[0]["index"] == NO_ORDER_CLASS
        and unresolved[0]["orders_available"] == 0
        and unresolved[0]["kernel_calls"] == 0,
        "sole no-order unresolved result",
    )
    kernel_totals = {
        "KILLED": sum(result["kernel_killed"] for result in checkpoint_results),
        "SURVIVORS": sum(result["kernel_survivors"] for result in checkpoint_results),
        "ABORT": sum(result["kernel_aborts"] for result in checkpoint_results),
    }
    require(kernel_totals == EXPECTED_KERNEL_TOTALS, "kernel status accounting")
    require(
        report.get("total_graphs") == 181
        and report.get("completed_graphs") == 181
        and report.get("status_counts") == EXPECTED_STATUS_COUNTS
        and report.get("certified_killed") == 1
        and report.get("killed_by_population") == {"K6": 0, "K7": 1}
        and report.get("kernel_status_totals") == EXPECTED_KERNEL_TOTALS,
        "report summary accounting",
    )
    runtime = report.get("runtime", {})
    require(
        runtime.get("sum_graph_seconds")
        == sum(result["elapsed_seconds"] for result in checkpoint_results)
        and runtime.get("maximum_graph_seconds")
        == max(result["elapsed_seconds"] for result in checkpoint_results),
        "report runtime accounting",
    )

    decisions_bytes = decisions_path.read_bytes()
    require(
        decisions_bytes == decision_text(checkpoint_results).encode("utf-8"),
        "decisions TSV reconstruction",
    )
    with io.StringIO(decisions_bytes.decode("utf-8")) as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    require(
        len(rows) == 181
        and [int(row["index"]) for row in rows]
        == [item["index"] for item in production],
        "decisions TSV identities",
    )
    require(
        report.get("decisions", {}).get("rows_excluding_header") == 181
        and report.get("decisions", {}).get("sha256") == sha256(decisions_path),
        "report decisions boundary",
    )
    return {
        "checkpoint_files": len(checkpoint_entries),
        "checkpoint_result_index_sha256": stable_hash(checkpoint_entries),
        "status_counts": status_counts,
        "kernel_status_totals": kernel_totals,
        "sole_certified_killed_class": KILLED_CLASS,
        "sole_unresolved_no_order_class": NO_ORDER_CLASS,
        "nonrejecting_results": statuses["ABORT"] + statuses["UNRESOLVED"],
        "infra_errors": statuses["INFRA_ERROR"],
        "decisions_rows": len(rows),
        "decisions_exactly_reconstructed": True,
    }


def verify_package_source(source_commit: str | None) -> dict:
    hashes = {name: sha256(ROOT / name) for name in PACKAGE_FILES}
    if source_commit is None:
        return {
            "status": "SOURCE_ONLY_UNCOMMITTED_BOUNDARY",
            "files_sha256": hashes,
        }
    require(len(source_commit) == 40, "package source commit")
    for name, expected in hashes.items():
        require(
            git_blob_sha256(source_commit, name) == expected,
            f"package committed source: {name}",
        )
    return {
        "status": "PASS",
        "commit": source_commit,
        "files_sha256": hashes,
    }


def verify_completion(
    campaign_path: Path = RAW_CAMPAIGN,
    report_path: Path = RAW_REPORT,
    decisions_path: Path = RAW_DECISIONS,
    checkpoint_dir: Path = RUN_DIR,
    *,
    enforce_raw_hashes: bool = True,
    source_commit: str | None = None,
) -> dict:
    if enforce_raw_hashes:
        require(sha256(campaign_path) == RAW_CAMPAIGN_SHA256, "raw campaign hash")
        require(sha256(report_path) == RAW_REPORT_SHA256, "raw report hash")
        require(sha256(decisions_path) == RAW_DECISIONS_SHA256, "raw decisions hash")
    campaign = load_json(campaign_path)
    report = load_json(report_path)
    configuration = campaign.get("configuration")
    require(isinstance(configuration, dict), "campaign configuration")
    require(
        campaign.get("schema") == 1
        and campaign.get("config_sha256") == CONFIG_SHA256
        and stable_hash(configuration) == CONFIG_SHA256,
        "campaign configuration identity",
    )
    require(
        report.get("schema") == 2
        and report.get("config_sha256") == CONFIG_SHA256
        and report.get("configuration") == configuration
        and report.get("claim")
        == (
            "Only KILLED is a certified non-realizability result. ABORT, "
            "UNRESOLVED, and INFRA_ERROR make no mathematical claim."
        ),
        "report configuration/semantics",
    )
    sources = verify_sources_and_inputs(configuration)
    production, selection = rebuild_selection(configuration)
    results = verify_results(report, production, checkpoint_dir, decisions_path)
    if enforce_raw_hashes:
        require(
            report.get("decisions", {}).get("sha256") == RAW_DECISIONS_SHA256,
            "pinned decisions hash",
        )
    run_state_path = checkpoint_dir / "run_state.json"
    run_state = load_json(run_state_path)
    require(
        run_state.get("status") == "COMPLETE"
        and run_state.get("completed_graphs") == 181
        and run_state.get("status_counts") == EXPECTED_STATUS_COUNTS
        and run_state.get("report_sha256") == sha256(report_path),
        "final run state",
    )
    return {
        "schema": "d6-interval-18-cover-v7-completion-verification-v1",
        "status": "PASS",
        "claim": (
            "The raw 181-class campaign is complete. Its sole certified "
            "KILLED class is 2100. All 179 ABORT and one UNRESOLVED results "
            "make no non-realizability claim."
        ),
        "raw_artifacts": {
            "campaign": {
                "path": str(campaign_path),
                "sha256": sha256(campaign_path),
            },
            "report": {"path": str(report_path), "sha256": sha256(report_path)},
            "decisions": {
                "path": str(decisions_path),
                "sha256": sha256(decisions_path),
            },
            "config_sha256": CONFIG_SHA256,
        },
        "sources_and_inputs": sources,
        "selection_reconstruction": selection,
        "result_audit": results,
        "semantics": {
            "certifying_status": "KILLED",
            "noncertifying_statuses": ["ABORT", "UNRESOLVED", "INFRA_ERROR"],
            "candidate_nonedges": "unconstrained and allowed to be unit distance",
            "kernel_replayed_by_this_checker": False,
        },
        "package": verify_package_source(source_commit),
    }


def atomic_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}.{time.time_ns()}")
    with temporary.open("w", encoding="utf-8", newline="") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", type=Path, default=RAW_CAMPAIGN)
    parser.add_argument("--report", type=Path, default=RAW_REPORT)
    parser.add_argument("--decisions", type=Path, default=RAW_DECISIONS)
    parser.add_argument("--checkpoint-dir", type=Path, default=RUN_DIR)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "d6_interval_18_cover_v7_completion_verification.json",
    )
    parser.add_argument(
        "--source-commit",
        help="require all three package sources to equal blobs in this commit",
    )
    args = parser.parse_args()
    result = verify_completion(
        args.campaign,
        args.report,
        args.decisions,
        args.checkpoint_dir,
        source_commit=args.source_commit,
    )
    result["execution"] = {
        "completed_utc": datetime.now(UTC).isoformat(),
        "command": " ".join(os.sys.argv),
        "python": os.sys.version.split()[0],
    }
    atomic_json(args.output, result)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
