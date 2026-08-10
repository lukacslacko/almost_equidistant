#!/usr/bin/env python3
"""Package the first certified kill from the live d=6 n=18 campaign.

The live campaign is intentionally allowed to continue.  This producer reads
only its immutable ``campaign.json`` and the atomic result for deletion class
2100.  It writes a deterministic gzip evidence object containing the exact
raw bytes and a compact theorem-level report.  It never reads ``progress`` or
``run_state`` and never modifies the checkpoint directory.

The independent checker does not import this module or the launch wrapper.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import platform
import subprocess
import time
from datetime import UTC, datetime
from pathlib import Path
from typing import Sequence

import run_d6_interval_18_cover_v7 as campaign_runner


ROOT = Path(__file__).resolve().parent
TARGET_CLASS_INDEX = 2_100
TARGET_ORDINAL = 32
TARGET_PARENT_INDEX = 3_945_564
TARGET_CLASS_ID = (
    "u18-b61a16db4126f818c181e6b5ef3340229711dec6162bceb4e08b7414ee70f0a1"
)
TARGET_GRAPH6 = "QTm}BxB{nTzVnwX}i}y^^^^]}~o"
TARGET_RAW_CAMPAIGN_SHA256 = (
    "70ba451094eb4eb5d039fc059ab8e3c351ebeed157e2eae39cc35fe4dcb3e348"
)
TARGET_RAW_RESULT_SHA256 = (
    "0e246dc1e2a84d5fafb3746366e5f1b06fb3b2368c0ad1b61ac12446d835815a"
)
TARGET_CONFIG_SHA256 = (
    "82c8d332af4ac4e3984949b98cde969f5b239eadb640e9633f704f8ec88fbcba"
)
CAMPAIGN_SOURCE_COMMIT = "22be775b1dafde30122ab579f66a7a539581ba73"

PACKAGE_FILES = (
    "build_d6_interval_18_class2100_increment.py",
    "verify_d6_interval_18_class2100_increment.py",
    "test_d6_interval_18_class2100_increment.py",
    "d6_interval_18_class2100_increment.md",
)
CAMPAIGN_SOURCES = {
    "run_d6_interval_18_cover_v7.py": (
        "5eeca35b328d930c94fa33d96c714ef756f4a777fcce8c8bffccbe472d744019"
    ),
    "run_d6_interval_residue.py": (
        "58fbbbc548ec3cc6b9d5b231a6790249b02a75d7bfd82f1afcd7fc77fce59404"
    ),
    "cdriver6.py": "8f65948c8bc95242d0fde50a8cbcb798cb0e225b78a8ce626382d00135b95173",
    "ckernel6.c": "383ef7a17c328869c338f64a16f2b90875ff968fe2fc20e9a5c3094c830e0874",
    "ival.py": "524e41e0d0637a5352c59ec998009b968f2c9f5e7f31b7f7927e8bc36c2e6fa0",
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
    "d6_current_residue_manifest_v8.json": (
        "9ea10a7794f033e66152c477a130f11c6c2862e87b4bdd705c14021521c18285"
    ),
    "d6_current_residue_manifest_v8_verification.json": (
        "117ac6833a24eb69cec7a514445bb4c2913c3f3ac7c413aa81fad19fe305602e"
    ),
}
EXPECTED_SEARCH = {
    "orders": 4,
    "cap": 100_000,
    "slices": 24,
    "zero_circle_only": False,
    "include_bulk_order": True,
}
TRUST_ASSUMPTIONS = {
    "basic_arithmetic": (
        "IEEE-754 binary64 basic operations and sqrt are correctly rounded; "
        "each interval endpoint is expanded with nextafter."
    ),
    "candidate_nonedges": "unconstrained and allowed to be unit distance",
    "transcendentals": (
        "macOS libm cos endpoint values are assumed within 8 ulps; the "
        "kernel pads both directions by 8 nextafter steps and includes all "
        "interior extrema."
    ),
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def bytes_sha256(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_hash(value: object) -> str:
    raw = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")
    return bytes_sha256(raw)


def load_json_bytes(raw: bytes, label: str) -> dict:
    value = json.loads(raw.decode("utf-8"))
    require(isinstance(value, dict), f"{label} is not a JSON object")
    return value


def validate_rows(value: object, order: int, label: str) -> tuple[int, ...]:
    require(
        isinstance(value, list)
        and len(value) == order
        and all(type(row) is int and row >= 0 for row in value),
        f"{label}: malformed adjacency",
    )
    rows = tuple(value)
    full = (1 << order) - 1
    for vertex, row in enumerate(rows):
        require(not row & ~full and not row & (1 << vertex), f"{label}: loop/range")
        for other in range(vertex):
            require(
                bool(row & (1 << other))
                == bool(rows[other] & (1 << vertex)),
                f"{label}: asymmetry",
            )
    return rows


def decode_graph6(text: str) -> tuple[int, ...]:
    require(text and 63 <= ord(text[0]) <= 125, "unsupported graph6 header")
    order = ord(text[0]) - 63
    require(order == 18, "target graph6 order")
    bits: list[int] = []
    for character in text[1:]:
        value = ord(character) - 63
        require(0 <= value < 64, "invalid graph6 character")
        bits.extend((value >> shift) & 1 for shift in range(5, -1, -1))
    require(len(bits) >= order * (order - 1) // 2, "short graph6 payload")
    rows = [0] * order
    position = 0
    for right in range(1, order):
        for left in range(right):
            if bits[position]:
                rows[left] |= 1 << right
                rows[right] |= 1 << left
            position += 1
    require(not any(bits[position:]), "nonzero graph6 padding")
    return tuple(rows)


def delete_vertex(rows: Sequence[int], deleted: int) -> tuple[int, ...]:
    old_vertices = [vertex for vertex in range(len(rows)) if vertex != deleted]
    return tuple(
        sum(
            1 << new_other
            for new_other, old_other in enumerate(old_vertices)
            if rows[old_vertex] & (1 << old_other)
        )
        for old_vertex in old_vertices
    )


def complement(rows: Sequence[int]) -> tuple[int, ...]:
    full = (1 << len(rows)) - 1
    return tuple(full ^ (1 << vertex) ^ row for vertex, row in enumerate(rows))


def spanning_subgraph(pattern: Sequence[int], target: Sequence[int]) -> bool:
    """Exact backtracking test: some relabeling maps every pattern edge to target."""

    require(len(pattern) == len(target), "subgraph orders differ")
    order = len(pattern)
    full = (1 << order) - 1
    pattern_degree = [row.bit_count() for row in pattern]
    target_degree = [row.bit_count() for row in target]
    domains = [
        sum(
            1 << vertex
            for vertex in range(order)
            if target_degree[vertex] >= pattern_degree[source]
        )
        for source in range(order)
    ]

    def search(current: list[int], unassigned: int) -> bool:
        if not unassigned:
            return True
        sources = [
            source for source in range(order) if unassigned & (1 << source)
        ]
        source = min(
            sources,
            key=lambda item: (current[item].bit_count(), -pattern_degree[item]),
        )
        choices = current[source]
        remaining = unassigned ^ (1 << source)
        while choices:
            choice = choices & -choices
            choices ^= choice
            vertex = choice.bit_length() - 1
            updated = current.copy()
            possible = True
            for other in sources:
                if other == source:
                    continue
                domain = updated[other] & ~choice
                if pattern[source] & (1 << other):
                    domain &= target[vertex]
                if not domain:
                    possible = False
                    break
                updated[other] = domain
            if possible and search(updated, remaining):
                return True
        return False

    return search(domains, full)


def load_v8_graphs() -> list[dict]:
    path = ROOT / "d6_current_residue_manifest_v8.json"
    check_path = ROOT / "d6_current_residue_manifest_v8_verification.json"
    require(sha256(path) == INPUTS[path.name], "v8 manifest hash")
    require(sha256(check_path) == INPUTS[check_path.name], "v8 verification hash")
    manifest = load_json_bytes(path.read_bytes(), path.name)
    verification = load_json_bytes(check_path.read_bytes(), check_path.name)
    require(
        manifest.get("schema") == "d6-current-certified-residue-v8"
        and manifest.get("status") == "COMPLETE_MIXED_CERTIFICATE_UNION"
        and verification.get("status") == "PASS"
        and verification.get("manifest", {}).get("sha256") == INPUTS[path.name],
        "v8 boundary",
    )
    graphs: list[dict] = []
    for class_name in ("K7", "K6_only"):
        block = manifest.get("classes", {}).get(class_name, {})
        records = block.get("graphs")
        require(isinstance(records, list), f"v8 {class_name} records")
        require(
            block.get("indices") == [record.get("index") for record in records]
            and stable_hash(block["indices"]) == block.get("indices_sha256")
            and stable_hash(records) == block.get("graphs_sha256"),
            f"v8 {class_name} roots",
        )
        for record in records:
            rows = validate_rows(record.get("adjacency"), 19, "v8 parent")
            graphs.append(
                {
                    "class": class_name,
                    "index": int(record["index"]),
                    "adjacency": rows,
                }
            )
    require(len(graphs) == 261, "v8 graph count")
    return graphs


def scan_v8_containment(target_rows: Sequence[int]) -> dict:
    """Scan every current parent deletion for the 106-edge obstruction.

    In complement language, a target unit-edge subgraph occurs exactly when
    the candidate complement is a spanning subgraph of the target complement.
    This direct domain-propagation implementation is deliberately different
    from the independent checker's deleted-edge enumeration.
    """

    target_edges = sum(row.bit_count() for row in target_rows) // 2
    target_complement = complement(target_rows)
    eligible = 0
    maximum_edges = 0
    occurrences = []
    for parent in load_v8_graphs():
        parent_rows = parent["adjacency"]
        parent_edges = sum(row.bit_count() for row in parent_rows) // 2
        for deleted in range(19):
            deletion_edges = parent_edges - parent_rows[deleted].bit_count()
            maximum_edges = max(maximum_edges, deletion_edges)
            if deletion_edges < target_edges:
                continue
            eligible += 1
            deletion = delete_vertex(parent_rows, deleted)
            if spanning_subgraph(complement(deletion), target_complement):
                occurrences.append(
                    {
                        "parent_class": parent["class"],
                        "parent_index": parent["index"],
                        "deleted_vertex": deleted,
                        "deletion_edges": deletion_edges,
                    }
                )
    require(
        occurrences
        == [
            {
                "parent_class": "K7",
                "parent_index": TARGET_PARENT_INDEX,
                "deleted_vertex": 4,
                "deletion_edges": 106,
            },
            {
                "parent_class": "K7",
                "parent_index": TARGET_PARENT_INDEX,
                "deleted_vertex": 15,
                "deletion_edges": 106,
            },
        ],
        "unexpected v8 containment result",
    )
    return {
        "parent_boundary": "d6_current_residue_manifest_v8.json",
        "parents_scanned": 261,
        "deletions_scanned": 261 * 19,
        "minimum_required_edges": target_edges,
        "eligible_deletions": eligible,
        "maximum_deletion_edges": maximum_edges,
        "occurrences": occurrences,
        "occurrences_sha256": stable_hash(occurrences),
        "rejected_parent_indices": [TARGET_PARENT_INDEX],
        "method": (
            "exact spanning-subgraph backtracking on complements with degree "
            "domains and edge-domain propagation"
        ),
    }


def git_blob_sha256(commit: str, name: str) -> str:
    blob = subprocess.run(
        ["git", "show", f"{commit}:{name}"],
        cwd=ROOT,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    ).stdout
    return bytes_sha256(blob)


def committed_package_boundary() -> dict:
    commit = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=ROOT,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    ).stdout.strip()
    branch = subprocess.run(
        ["git", "branch", "--show-current"],
        cwd=ROOT,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    ).stdout.strip()
    require(branch == "codex/dimension6", "official package branch")
    sources = {name: sha256(ROOT / name) for name in PACKAGE_FILES}
    for name, expected in sources.items():
        require(git_blob_sha256(commit, name) == expected, f"uncommitted package: {name}")
    return {
        "status": "PASS",
        "commit": commit,
        "branch": branch,
        "package_sources": dict(sorted(sources.items())),
        "package_sources_equal_committed_blobs": True,
    }


def validate_campaign_sources(configuration: dict) -> dict:
    sources = configuration.get("sources")
    require(isinstance(sources, dict), "campaign sources")
    git = configuration.get("git", {})
    require(
        git.get("branch") == "codex/dimension6"
        and git.get("source_commit") == CAMPAIGN_SOURCE_COMMIT
        and git.get("sources") == dict(sorted(CAMPAIGN_SOURCES.items()))
        and git.get("all_sources_equal_head_and_source_commit_blobs") is True,
        "campaign Git boundary",
    )
    for name, expected in CAMPAIGN_SOURCES.items():
        require(sources.get(name) == expected, f"recorded campaign source: {name}")
        require(sha256(ROOT / name) == expected, f"working campaign source: {name}")
        require(
            git_blob_sha256(CAMPAIGN_SOURCE_COMMIT, name) == expected,
            f"committed campaign source: {name}",
        )
    require(
        sources.get("ckernel6.dylib") == KERNEL_BINARY_SHA256
        and sha256(ROOT / "ckernel6.dylib") == KERNEL_BINARY_SHA256,
        "campaign kernel binary",
    )
    for name, expected in INPUTS.items():
        require(sha256(ROOT / name) == expected, f"input artifact hash: {name}")
    for key in ("v7_manifest", "v7_verification", "deletion_manifest", "deletion_verification"):
        record = configuration.get("inputs", {}).get(key, {})
        require(
            isinstance(record.get("path"), str)
            and Path(record["path"]).name in INPUTS
            and record.get("sha256") == INPUTS[Path(record["path"]).name],
            f"campaign input record: {key}",
        )
    return {
        "source_commit": CAMPAIGN_SOURCE_COMMIT,
        "source_sha256": dict(sorted(CAMPAIGN_SOURCES.items())),
        "kernel_binary_sha256": KERNEL_BINARY_SHA256,
        "input_sha256": dict(sorted(INPUTS.items())),
    }


def reconstruct_selection(configuration: dict) -> tuple[dict, dict, list[dict]]:
    parents, parent_boundary = campaign_runner.load_v7_parents(
        ROOT / "d6_current_residue_manifest_v7.json",
        ROOT / "d6_current_residue_manifest_v7_verification.json",
    )
    deletions, deletion_boundary = campaign_runner.load_deletion_boundary(
        ROOT / "d6_residue_18_deletions_v2.json",
        ROOT / "d6_residue_18_deletions_v2_verification.json",
    )
    base, cover = campaign_runner.build_dense_deletion_cover(parents, deletions)
    production, augmentation = campaign_runner.build_actionable_campaign(
        base, cover, deletions
    )
    profile = campaign_runner.profile_placement_orders(production, 4)
    graphs = campaign_runner.select_run_graphs(
        production,
        shards=1,
        shard=0,
        sample=None,
        sample_seed=618_263_179,
        indices=None,
    )
    run_selection = {
        "graphs": len(graphs),
        "class_indices": [graph["index"] for graph in graphs],
        "class_indices_sha256": stable_hash([graph["index"] for graph in graphs]),
        "population_counts": {"K7": 12, "K6": 169},
        "sampling": {
            "sample": None,
            "sample_seed": 618_263_179,
            "shards": 1,
            "shard": 0,
            "requested_indices": [],
        },
    }
    expected = {
        "parent_boundary": parent_boundary,
        "deletion_boundary": deletion_boundary,
        "full_dense_cover": cover,
        "actionable_augmentation": augmentation,
        "placement_order_profile": profile,
        "run_selection": run_selection,
    }
    require(configuration.get("selection") == expected, "campaign selection reconstruction")
    require(graphs[TARGET_ORDINAL]["index"] == TARGET_CLASS_INDEX, "target ordinal")
    return deletions, cover, graphs


def validate_target_checkpoint(
    checkpoint: dict, configuration: dict, graph: dict
) -> dict:
    with campaign_runner.engine_order_18():
        result = campaign_runner.engine.validate_checkpoint(
            checkpoint,
            TARGET_CONFIG_SHA256,
            graph,
            EXPECTED_SEARCH["slices"],
        )
    require(
        result.get("status") == "KILLED"
        and result.get("ordinal") == TARGET_ORDINAL
        and result.get("index") == TARGET_CLASS_INDEX
        and result.get("population") == "K7"
        and result.get("orders_available") == 2
        and result.get("orders_tried") == 2
        and result.get("winning_order") == 1
        and result.get("winning_circles") == 1
        and result.get("kernel_calls") == 31
        and result.get("kernel_killed") == 30
        and result.get("kernel_survivors") == 0
        and result.get("kernel_aborts") == 1
        and result.get("unresolved_cells") == 0,
        "target checkpoint counters",
    )
    records = result.get("winning_slice_records")
    require(
        isinstance(records, list)
        and len(records) == 24
        and all(record.get("status") == "KILLED" for record in records),
        "target winning slices",
    )
    require(configuration.get("search") == EXPECTED_SEARCH, "campaign search")
    return result


def make_evidence(campaign_raw: bytes, result_raw: bytes) -> bytes:
    evidence = {
        "schema": "d6-interval-18-class2100-raw-evidence-v1",
        "campaign": {
            "name": "campaign.json",
            "sha256": bytes_sha256(campaign_raw),
            "raw_utf8": campaign_raw.decode("utf-8"),
        },
        "result": {
            "name": "result_000032_0002100.json",
            "sha256": bytes_sha256(result_raw),
            "raw_utf8": result_raw.decode("utf-8"),
        },
    }
    payload = json.dumps(
        evidence, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")
    return gzip.compress(payload, compresslevel=9, mtime=0)


def atomic_bytes(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    with temporary.open("wb") as stream:
        stream.write(payload)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def atomic_json(path: Path, value: object) -> None:
    payload = (json.dumps(value, indent=2, sort_keys=True) + "\n").encode("utf-8")
    atomic_bytes(path, payload)


def build_increment(
    campaign_path: Path,
    result_path: Path,
    evidence_path: Path,
    output_path: Path,
    *,
    source_boundary: dict | None = None,
) -> dict:
    campaign_raw = campaign_path.resolve().read_bytes()
    result_raw = result_path.resolve().read_bytes()
    require(bytes_sha256(campaign_raw) == TARGET_RAW_CAMPAIGN_SHA256, "raw campaign hash")
    require(bytes_sha256(result_raw) == TARGET_RAW_RESULT_SHA256, "raw result hash")
    campaign = load_json_bytes(campaign_raw, "raw campaign")
    checkpoint = load_json_bytes(result_raw, "raw result")
    configuration = campaign.get("configuration")
    require(
        campaign.get("schema") == 1
        and isinstance(configuration, dict)
        and campaign.get("config_sha256") == TARGET_CONFIG_SHA256
        and stable_hash(configuration) == TARGET_CONFIG_SHA256,
        "raw campaign/configuration boundary",
    )
    require(checkpoint.get("config_sha256") == TARGET_CONFIG_SHA256, "checkpoint campaign binding")
    require(configuration.get("trust_assumptions") == TRUST_ASSUMPTIONS, "trust assumptions")
    campaign_sources = validate_campaign_sources(configuration)
    deletions, cover, graphs = reconstruct_selection(configuration)
    target_record = deletions["unique_deletions"][TARGET_CLASS_INDEX]
    target_rows = validate_rows(target_record.get("adjacency"), 18, "target class")
    require(
        target_record.get("class_id") == TARGET_CLASS_ID
        and target_record.get("canonical_graph6") == TARGET_GRAPH6
        and hashlib.sha256(TARGET_GRAPH6.encode("ascii")).hexdigest()
        == TARGET_CLASS_ID.removeprefix("u18-")
        and decode_graph6(TARGET_GRAPH6) == target_rows
        and target_record.get("edges") == 106
        and target_record.get("clique_number") == 7
        and target_record.get("standard18_compatible") is False,
        "target deletion-class boundary",
    )
    target_graph = graphs[TARGET_ORDINAL]
    require(tuple(target_graph["adjacency"]) == target_rows, "target run adjacency")
    result = validate_target_checkpoint(checkpoint, configuration, target_graph)
    selected_record = next(
        record
        for record in cover["selected_records"]
        if record["class_index"] == TARGET_CLASS_INDEX
    )
    require(
        selected_record.get("active_parents") == [f"K7:{TARGET_PARENT_INDEX}"]
        and selected_record.get("known_realizable_positive_control") is False,
        "target cover record",
    )
    raw_occurrences = target_record.get("occurrences")
    require(
        raw_occurrences
        == [
            {
                "attachment_degree": 11,
                "attachment_mask": 117_613,
                "deleted_vertex": 4,
                "parent_class": "K7",
                "parent_index": TARGET_PARENT_INDEX,
            },
            {
                "attachment_degree": 11,
                "attachment_mask": 95_925,
                "deleted_vertex": 15,
                "parent_class": "K7",
                "parent_index": TARGET_PARENT_INDEX,
            },
        ],
        "target occurrence records",
    )
    controls = configuration.get("kernel_controls", {})
    require(
        controls.get("status") == "PASS"
        and controls.get("graph_order") == 18
        and controls.get("flexible_circle_and_K8", {}).get("status") == "PASS"
        and all(
            record.get("status") != "KILLED"
            for record in controls.get("known_realizable_exact18", ())
        ),
        "recorded controls",
    )
    containment = scan_v8_containment(target_rows)
    evidence_bytes = make_evidence(campaign_raw, result_raw)
    atomic_bytes(evidence_path.resolve(), evidence_bytes)
    package_boundary = source_boundary or committed_package_boundary()
    report = {
        "schema": "d6-interval-18-class2100-increment-v1",
        "status": "COMPLETE_INTERVAL_REJECTION",
        "claim": (
            "The 18-vertex required-edge graph in deletion class 2100 has no "
            "realization by distinct points in R6 under the recorded interval "
            "trust assumptions.  Hence current parent 3945564 is impossible."
        ),
        "semantics": {
            "required_edges": "distance exactly one",
            "candidate_nonedges": "unconstrained and may also have distance one",
            "distinct_points_required": True,
            "only_KILLED_is_a_rejection": True,
            "subgraph_rule": (
                "a realization of a parent would restrict to a realization of "
                "every required-edge subgraph on any vertex deletion"
            ),
            "arithmetic": "outward-rounded binary64 interval branch-and-bound",
        },
        "trust_assumptions": TRUST_ASSUMPTIONS,
        "raw_evidence": {
            "path": evidence_path.name,
            "sha256": bytes_sha256(evidence_bytes),
            "campaign_sha256": TARGET_RAW_CAMPAIGN_SHA256,
            "result_sha256": TARGET_RAW_RESULT_SHA256,
            "config_sha256": TARGET_CONFIG_SHA256,
            "evidence_contains_exact_raw_utf8": True,
        },
        "campaign_source_boundary": campaign_sources,
        "search": {**EXPECTED_SEARCH, "workers": 3},
        "deletion_class": {
            "class_index": TARGET_CLASS_INDEX,
            "ordinal_in_production_run": TARGET_ORDINAL,
            "class_id": TARGET_CLASS_ID,
            "canonical_graph6": TARGET_GRAPH6,
            "adjacency": list(target_rows),
            "adjacency_sha256": stable_hash(list(target_rows)),
            "edges": 106,
            "clique_number": 7,
            "parent_occurrences": raw_occurrences,
            "parent_occurrences_sha256": stable_hash(raw_occurrences),
        },
        "interval_certificate": {
            "status": "KILLED",
            "winning_order": result["winning_order"],
            "winning_circles": result["winning_circles"],
            "winning_seed": result["winning_seed"],
            "winning_placement_order": result["winning_placement_order"],
            "winning_slice_records": result["winning_slice_records"],
            "winning_slice_records_sha256": stable_hash(
                result["winning_slice_records"]
            ),
            "raw_result": result,
            "raw_result_sha256": stable_hash(result),
        },
        "controls": {
            "recorded_status": controls["status"],
            "recorded_controls_sha256": stable_hash(controls),
            "flexible_circle_slices": 24,
            "K8_slices": 24,
            "exact_positive_graphs": 4,
            "independent_replay_required": True,
        },
        "current_v8_containment": containment,
        "conclusion": {
            "rejected_parent_indices": [TARGET_PARENT_INDEX],
            "rejected_parent_indices_sha256": stable_hash([TARGET_PARENT_INDEX]),
            "new_interval_rejections": 1,
        },
        "inputs_sha256": dict(sorted(INPUTS.items())),
        "execution": {
            "created_utc": datetime.now(UTC).isoformat(),
            "python": platform.python_version(),
            "platform": platform.platform(),
            "package": package_boundary,
        },
    }
    atomic_json(output_path.resolve(), report)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", type=Path, required=True)
    parser.add_argument("--result", type=Path, required=True)
    parser.add_argument(
        "--evidence",
        type=Path,
        default=ROOT / "d6_interval_18_class2100_increment_evidence.json.gz",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "d6_interval_18_class2100_increment_report.json",
    )
    args = parser.parse_args()
    report = build_increment(
        args.campaign, args.result, args.evidence, args.output
    )
    print(
        json.dumps(
            {
                "status": report["status"],
                "class_index": TARGET_CLASS_INDEX,
                "rejected_parent_indices": [TARGET_PARENT_INDEX],
                "evidence": str(args.evidence),
                "output": str(args.output),
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
