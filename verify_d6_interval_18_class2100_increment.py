#!/usr/bin/env python3
"""Independent checker for the partial d=6 n=18 class-2100 kill.

This file imports neither the producer, the n=18 launch wrapper, nor the
generic campaign manager.  It independently rebuilds the dense cover from
the frozen v7/deletion corpora, implements the placement-order generator,
reconstructs both parent deletions by exact graph isomorphism, replays the
complete two-order interval search and all controls through the pinned low-
level kernel, and performs an independent all-v8 containment scan.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import itertools
import json
import math
import os
import platform
import subprocess
import time
from collections import Counter, defaultdict
from datetime import UTC, datetime
from pathlib import Path
from typing import Callable, Sequence

import cdriver6


ROOT = Path(__file__).resolve().parent
N = 18
D = 6
TARGET_CLASS_INDEX = 2_100
TARGET_ORDINAL = 32
TARGET_PARENT_INDEX = 3_945_564
TARGET_CLASS_ID = (
    "u18-b61a16db4126f818c181e6b5ef3340229711dec6162bceb4e08b7414ee70f0a1"
)
TARGET_GRAPH6 = "QTm}BxB{nTzVnwX}i}y^^^^]}~o"
RAW_CAMPAIGN_SHA256 = (
    "70ba451094eb4eb5d039fc059ab8e3c351ebeed157e2eae39cc35fe4dcb3e348"
)
RAW_RESULT_SHA256 = (
    "0e246dc1e2a84d5fafb3746366e5f1b06fb3b2368c0ad1b61ac12446d835815a"
)
CONFIG_SHA256 = "82c8d332af4ac4e3984949b98cde969f5b239eadb640e9633f704f8ec88fbcba"
CAMPAIGN_SOURCE_COMMIT = "22be775b1dafde30122ab579f66a7a539581ba73"
PACKAGE_FILES = {
    "build_d6_interval_18_class2100_increment.py",
    "verify_d6_interval_18_class2100_increment.py",
    "test_d6_interval_18_class2100_increment.py",
    "d6_interval_18_class2100_increment.md",
}
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
EXPECTED_V7_COUNTS = {"K7": 12, "K6_only": 251}
EXPECTED_V7_CLASS_HASHES = {
    "K7": "e545ea5a8fc9691ed98e885919e149f84dc7e8b4796d0fe9e6cb3a519df6fd8c",
    "K6_only": "a23da7bd9ae21eb7c078a08a7bce92f22d397d5c71f1e4f51d9cf8e59bbfa907",
}
EXPECTED_V7_COMBINED_HASH = (
    "a50baf54bfff642e7c1711fe7041b5184df12f90de2cf52557cd086f10875127"
)
EXPECTED_V8_COUNTS = {"K7": 12, "K6_only": 249}
EXPECTED_V8_CLASS_HASHES = {
    "K7": "e545ea5a8fc9691ed98e885919e149f84dc7e8b4796d0fe9e6cb3a519df6fd8c",
    "K6_only": "f274aab3f2d408fce34c7cd7eb4518621fa9d634a8943e3ee1a38d76df9f3d31",
}
EXPECTED_V8_COMBINED_HASH = (
    "2391a93a3629517363106603bdad00be9b6f960966d9089ef37d93be2988c213"
)
EXPECTED_PRODUCTION_HASH = (
    "8c624c969b517bb5663dc32bdf577ab3addbe8a108458b2803443c17f294d8e6"
)
EXPECTED_GREEDY_HASH = (
    "ef4ac7e73ea924d3d9c6382f590dc73cdf2e65b357e86ddc69719372f717d6fe"
)
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
        raise AssertionError(message)


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


def load_json(path: Path) -> dict:
    value = json.loads(path.read_text(encoding="utf-8"))
    require(isinstance(value, dict), f"{path.name}: not a JSON object")
    return value


def parse_json_raw(raw: bytes, label: str) -> dict:
    value = json.loads(raw.decode("utf-8"))
    require(isinstance(value, dict), f"{label}: not a JSON object")
    return value


def validate_rows(value: object, order: int, label: str) -> tuple[int, ...]:
    require(
        isinstance(value, list)
        and len(value) == order
        and all(type(row) is int and row >= 0 for row in value),
        f"{label}: adjacency structure",
    )
    rows = tuple(value)
    full = (1 << order) - 1
    for vertex, row in enumerate(rows):
        require(not row & ~full and not row & (1 << vertex), f"{label}: range/loop")
        for other in range(vertex):
            require(
                bool(row & (1 << other))
                == bool(rows[other] & (1 << vertex)),
                f"{label}: asymmetry",
            )
    return rows


def decode_graph6(text: str) -> tuple[int, ...]:
    require(len(text) >= 2 and ord(text[0]) - 63 == N, "graph6 header")
    bits = []
    for character in text[1:]:
        value = ord(character) - 63
        require(0 <= value < 64, "graph6 character")
        bits.extend((value >> shift) & 1 for shift in range(5, -1, -1))
    rows = [0] * N
    position = 0
    for right in range(1, N):
        for left in range(right):
            require(position < len(bits), "short graph6")
            if bits[position]:
                rows[left] |= 1 << right
                rows[right] |= 1 << left
            position += 1
    require(not any(bits[position:]), "graph6 padding")
    return tuple(rows)


def delete_vertex(rows: Sequence[int], deleted: int) -> tuple[int, ...]:
    surviving = [vertex for vertex in range(len(rows)) if vertex != deleted]
    answer = []
    for old_vertex in surviving:
        answer.append(
            sum(
                1 << new_other
                for new_other, old_other in enumerate(surviving)
                if rows[old_vertex] & (1 << old_other)
            )
        )
    return tuple(answer)


def compress_mask(mask: int, deleted: int, order: int = 19) -> int:
    surviving = [vertex for vertex in range(order) if vertex != deleted]
    return sum(
        1 << new_vertex
        for new_vertex, old_vertex in enumerate(surviving)
        if mask & (1 << old_vertex)
    )


def complement(rows: Sequence[int]) -> tuple[int, ...]:
    full = (1 << len(rows)) - 1
    return tuple(full ^ (1 << vertex) ^ row for vertex, row in enumerate(rows))


def wl_colors(rows: Sequence[int]) -> tuple[int, ...]:
    colors = [row.bit_count() for row in rows]
    while True:
        signatures = [
            (
                colors[vertex],
                tuple(
                    sorted(
                        colors[other]
                        for other in range(len(rows))
                        if rows[vertex] & (1 << other)
                    )
                ),
            )
            for vertex in range(len(rows))
        ]
        labels = {
            signature: label
            for label, signature in enumerate(sorted(set(signatures)))
        }
        updated = [labels[signature] for signature in signatures]
        if updated == colors:
            return tuple(colors)
        colors = updated


def graph_invariant(rows: Sequence[int]) -> tuple[object, ...]:
    colors = wl_colors(rows)
    return (
        tuple(sorted(row.bit_count() for row in rows)),
        tuple(sorted(Counter(colors).values())),
        tuple(sorted((colors[v], rows[v].bit_count()) for v in range(len(rows)))),
    )


def exact_isomorphic(left: Sequence[int], right: Sequence[int]) -> bool:
    require(len(left) == len(right), "isomorphism order")
    if graph_invariant(left) != graph_invariant(right):
        return False
    order = len(left)
    left_colors = wl_colors(left)
    right_colors = wl_colors(right)
    domains = [
        sum(1 << vertex for vertex in range(order) if right_colors[vertex] == left_colors[source])
        for source in range(order)
    ]
    mapping: dict[int, int] = {}
    used = 0

    def search() -> bool:
        nonlocal used
        if len(mapping) == order:
            return True
        best_source = -1
        best_choices = 0
        best_key = None
        for source in range(order):
            if source in mapping:
                continue
            choices = domains[source] & ~used
            for mapped_source, mapped_target in mapping.items():
                required_edge = bool(left[source] & (1 << mapped_source))
                choices &= (
                    right[mapped_target]
                    if required_edge
                    else ~right[mapped_target]
                )
            choices &= (1 << order) - 1
            key = (choices.bit_count(), -left[source].bit_count(), source)
            if best_key is None or key < best_key:
                best_key = key
                best_source = source
                best_choices = choices
        if not best_choices:
            return False
        while best_choices:
            choice = best_choices & -best_choices
            best_choices ^= choice
            target = choice.bit_length() - 1
            mapping[best_source] = target
            used |= choice
            if search():
                return True
            used ^= choice
            del mapping[best_source]
        return False

    return search()


def git_blob_sha256(commit: str, name: str) -> str:
    blob = subprocess.run(
        ["git", "show", f"{commit}:{name}"],
        cwd=ROOT,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    ).stdout
    return bytes_sha256(blob)


def verify_package_boundary(report: dict, enforce: bool) -> dict:
    boundary = report.get("execution", {}).get("package", {})
    if not enforce:
        return {"status": "SKIPPED_TEST_MODE"}
    require(
        boundary.get("status") == "PASS"
        and boundary.get("branch") == "codex/dimension6"
        and boundary.get("package_sources_equal_committed_blobs") is True,
        "package boundary",
    )
    commit = boundary.get("commit")
    sources = boundary.get("package_sources")
    require(isinstance(commit, str) and len(commit) == 40, "package commit")
    require(isinstance(sources, dict) and set(sources) == PACKAGE_FILES, "package files")
    for name, expected in sources.items():
        require(sha256(ROOT / name) == expected, f"working package source: {name}")
        require(git_blob_sha256(commit, name) == expected, f"committed package source: {name}")
    for name, expected in INPUTS.items():
        require(git_blob_sha256(commit, name) == expected, f"committed input: {name}")
    return {
        "status": "PASS",
        "commit": commit,
        "branch": boundary["branch"],
        "source_files": len(sources),
    }


def verify_static_files(configuration: dict, report: dict) -> dict:
    require(report.get("inputs_sha256") == dict(sorted(INPUTS.items())), "report inputs")
    for name, expected in INPUTS.items():
        require(sha256(ROOT / name) == expected, f"working input: {name}")
    sources = configuration.get("sources", {})
    git = configuration.get("git", {})
    require(
        git.get("branch") == "codex/dimension6"
        and git.get("source_commit") == CAMPAIGN_SOURCE_COMMIT
        and git.get("sources") == dict(sorted(CAMPAIGN_SOURCES.items())),
        "raw campaign source boundary",
    )
    for name, expected in CAMPAIGN_SOURCES.items():
        require(sources.get(name) == expected, f"raw source record: {name}")
        require(sha256(ROOT / name) == expected, f"working campaign source: {name}")
        require(
            git_blob_sha256(CAMPAIGN_SOURCE_COMMIT, name) == expected,
            f"campaign commit blob: {name}",
        )
    require(
        sources.get(cdriver6.LIBNAME) == KERNEL_BINARY_SHA256
        and sha256(ROOT / cdriver6.LIBNAME) == KERNEL_BINARY_SHA256,
        "kernel binary hash",
    )
    report_boundary = report.get("campaign_source_boundary", {})
    require(
        report_boundary.get("source_commit") == CAMPAIGN_SOURCE_COMMIT
        and report_boundary.get("source_sha256") == dict(sorted(CAMPAIGN_SOURCES.items()))
        and report_boundary.get("kernel_binary_sha256") == KERNEL_BINARY_SHA256
        and report_boundary.get("input_sha256") == dict(sorted(INPUTS.items())),
        "report campaign sources",
    )
    return {
        "source_commit": CAMPAIGN_SOURCE_COMMIT,
        "sources_checked": len(CAMPAIGN_SOURCES),
        "inputs_checked": len(INPUTS),
        "kernel_binary_sha256": KERNEL_BINARY_SHA256,
    }


def load_residue(name: str, counts: dict, hashes: dict, combined_hash: str) -> list[dict]:
    manifest = load_json(ROOT / name)
    check_name = name.replace(".json", "_verification.json")
    verification = load_json(ROOT / check_name)
    version = "v7" if "v7" in name else "v8"
    require(
        manifest.get("schema") == f"d6-current-certified-residue-{version}"
        and manifest.get("status") == "COMPLETE_MIXED_CERTIFICATE_UNION"
        and manifest.get("class_order") == ["K7", "K6_only"]
        and verification.get("status") == "PASS"
        and verification.get("manifest", {}).get("sha256") == INPUTS[name]
        and all(verification.get("checks", {}).values()),
        f"{version} manifest/verification",
    )
    graphs = []
    ordered_indices = []
    for class_name in ("K7", "K6_only"):
        block = manifest.get("classes", {}).get(class_name, {})
        records = block.get("graphs")
        indices = block.get("indices")
        require(
            isinstance(records, list)
            and isinstance(indices, list)
            and len(records) == len(indices) == counts[class_name]
            and indices == [record.get("index") for record in records]
            and stable_hash(indices) == block.get("indices_sha256") == hashes[class_name]
            and stable_hash(records) == block.get("graphs_sha256"),
            f"{version} {class_name} roots",
        )
        for record in records:
            graphs.append(
                {
                    "class": class_name,
                    "index": int(record["index"]),
                    "adjacency": validate_rows(record.get("adjacency"), 19, "parent"),
                }
            )
            ordered_indices.append(int(record["index"]))
    require(
        stable_hash(ordered_indices)
        == manifest.get("combined", {}).get("ordered_indices_sha256")
        == combined_hash,
        f"{version} combined root",
    )
    return graphs


def known_positive(record: dict) -> bool:
    return bool(record.get("standard18_compatible")) or record.get("class_id") in {
        "u18-761add09c775ad5ba966a0eef97848e076faacc647744c4079b7bee48a65d92f",
        "u18-84a0a88e385b1aa4673cc149e0c2458f37b7cc6aef3def3cf2589f2075be4841",
    }


def parent_key(class_name: str, index: int) -> tuple[int, int]:
    return ({"K7": 0, "K6_only": 1}[class_name], index)


def tagged(parent: tuple[str, int]) -> str:
    return f"{parent[0]}:{parent[1]}"


def reconstruct_production(v7_graphs: list[dict], deletion_manifest: dict, configuration: dict) -> dict:
    records = deletion_manifest.get("unique_deletions")
    require(isinstance(records, list) and len(records) == 11_975, "deletion records")
    parents = {(graph["class"], graph["index"]) for graph in v7_graphs}
    occurrences: dict[tuple[str, int], list[tuple[int, int, dict]]] = defaultdict(list)
    for class_index, record in enumerate(records):
        for occurrence in record.get("occurrences", ()):
            parent = (occurrence.get("parent_class"), occurrence.get("parent_index"))
            if parent in parents:
                occurrences[parent].append((class_index, int(record["edges"]), occurrence))
    require(set(occurrences) == parents, "parent occurrence coverage")
    coverage: dict[int, set[tuple[str, int]]] = defaultdict(set)
    dense_by_parent: dict[str, list[int]] = {}
    for parent in sorted(parents, key=lambda p: parent_key(*p)):
        options = occurrences[parent]
        require(
            len(options) == 19
            and {int(option[2]["deleted_vertex"]) for option in options} == set(range(19)),
            f"rooted deletion partition: {parent}",
        )
        maximum = max(option[1] for option in options)
        dense = sorted({option[0] for option in options if option[1] == maximum})
        dense_by_parent[tagged(parent)] = dense
        for class_index in dense:
            coverage[class_index].add(parent)
    uncovered = set(parents)
    selected = []
    steps = []
    while uncovered:
        choices = [index for index, covered in coverage.items() if covered & uncovered]
        require(bool(choices), "greedy cover stalled")
        chosen = max(
            choices,
            key=lambda index: (len(coverage[index] & uncovered), records[index]["class_id"]),
        )
        new = sorted(coverage[chosen] & uncovered, key=lambda p: parent_key(*p))
        selected.append(chosen)
        steps.append((chosen, [tagged(parent) for parent in new]))
        uncovered.difference_update(new)
    require(len(selected) == 179 and stable_hash(selected) == EXPECTED_GREEDY_HASH, "greedy cover")
    selected_coverers: dict[str, list[int]] = defaultdict(list)
    for index in selected:
        for parent in coverage[index]:
            selected_coverers[tagged(parent)].append(index)
    blocked = sorted(
        [
            parent
            for parent, coverers in selected_coverers.items()
            if coverers and all(known_positive(records[index]) for index in coverers)
        ],
        key=lambda item: parent_key(item.split(":")[0], int(item.split(":")[1])),
    )
    require(blocked == ["K6_only:364827", "K6_only:3335955"], "blocked parents")
    backups = []
    for parent in blocked:
        candidates = []
        for index in dense_by_parent[parent]:
            if known_positive(records[index]):
                continue
            active = coverage[index]
            candidates.append((len(active), records[index]["class_id"], index))
        require(bool(candidates), "missing backup")
        backups.append(max(candidates)[2])
    require(backups == [126, 5_673], "backup classes")
    production = sorted(selected + backups)
    require(
        len(production) == 181
        and len(set(production)) == 181
        and stable_hash(production) == EXPECTED_PRODUCTION_HASH
        and production[TARGET_ORDINAL] == TARGET_CLASS_INDEX,
        "production selection",
    )
    raw_selection = configuration.get("selection", {})
    require(
        raw_selection.get("run_selection", {}).get("class_indices") == production
        and raw_selection.get("run_selection", {}).get("class_indices_sha256")
        == EXPECTED_PRODUCTION_HASH
        and raw_selection.get("actionable_augmentation", {}).get(
            "production_class_indices"
        )
        == production
        and raw_selection.get("full_dense_cover", {}).get(
            "selected_class_indices_greedy_order"
        )
        == selected,
        "raw selection agrees with reconstruction",
    )
    target_step = next(
        position for position, (index, _new) in enumerate(steps) if index == TARGET_CLASS_INDEX
    )
    require(steps[target_step] == (TARGET_CLASS_INDEX, [f"K7:{TARGET_PARENT_INDEX}"]), "target greedy step")
    return {
        "v7_parents": len(parents),
        "eligible_dense_classes": len(coverage),
        "selected_classes": len(selected),
        "greedy_indices_sha256": stable_hash(selected),
        "production_classes": len(production),
        "production_indices_sha256": stable_hash(production),
        "target_production_ordinal": TARGET_ORDINAL,
        "target_greedy_step": target_step,
        "target_newly_covered_parent": f"K7:{TARGET_PARENT_INDEX}",
    }


def reconstruct_occurrences(v7_graphs: list[dict], target: Sequence[int], record: dict) -> dict:
    candidates = []
    occurrences = []
    target_edges = sum(row.bit_count() for row in target) // 2
    target_degrees = sorted(row.bit_count() for row in target)
    for parent in v7_graphs:
        parent_rows = parent["adjacency"]
        parent_edges = sum(row.bit_count() for row in parent_rows) // 2
        for deleted in range(19):
            edge_count = parent_edges - parent_rows[deleted].bit_count()
            if edge_count != target_edges:
                continue
            deletion = delete_vertex(parent_rows, deleted)
            if sorted(row.bit_count() for row in deletion) != target_degrees:
                continue
            candidates.append((parent, deleted, deletion))
            if exact_isomorphic(deletion, target):
                occurrences.append(
                    {
                        "attachment_degree": parent_rows[deleted].bit_count(),
                        "attachment_mask": compress_mask(parent_rows[deleted], deleted),
                        "deleted_vertex": deleted,
                        "parent_class": parent["class"],
                        "parent_index": parent["index"],
                    }
                )
    expected = [
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
    ]
    require(occurrences == expected, "induced occurrence reconstruction")
    require(
        record.get("occurrence_count") == 2
        and record.get("parent_count") == 1
        and record.get("parent_classes") == ["K7"]
        and record.get("occurrences") == occurrences,
        "deletion manifest occurrence record",
    )
    return {
        "edge_and_degree_candidates": len(candidates),
        "exact_isomorphism_occurrences": occurrences,
        "occurrences_sha256": stable_hash(occurrences),
        "unique_parent_indices": [TARGET_PARENT_INDEX],
    }


def clique_list(adjacency: Sequence[int], size: int) -> list[tuple[int, ...]]:
    answer = []

    def extend(clique: list[int], candidates: int, start: int) -> None:
        if len(clique) == size:
            answer.append(tuple(clique))
            return
        for vertex in range(start, len(adjacency)):
            if candidates & (1 << vertex):
                extend(clique + [vertex], candidates & adjacency[vertex], vertex + 1)

    extend([], (1 << len(adjacency)) - 1, 0)
    return answer


def closure_size(adjacency: Sequence[int], vertex: int, placed: set[int]) -> int:
    closure = set(placed)
    closure.add(vertex)
    changed = True
    while changed:
        changed = False
        for candidate in range(len(adjacency)):
            if candidate not in closure and sum(
                bool(adjacency[candidate] & (1 << other)) for other in closure
            ) >= D:
                closure.add(candidate)
                changed = True
    return len(closure)


def layered_order(
    adjacency: Sequence[int],
    seed: tuple[int, ...],
    picker: Callable[[list[int], set[int]], int],
    force_circle_first: bool,
) -> list[int] | None:
    placed = set(seed)
    answer = []
    forced = not force_circle_first
    while len(placed) < len(adjacency):
        if not forced:
            circle = [
                vertex
                for vertex in range(len(adjacency))
                if vertex not in placed
                and sum(
                    bool(adjacency[vertex] & (1 << other)) for other in placed
                )
                == D - 1
            ]
            if circle:
                vertex = picker(circle, placed)
                answer.append(vertex)
                placed.add(vertex)
                forced = True
                continue
        progressed = False
        while True:
            candidates = [
                (
                    sum(
                        bool(adjacency[vertex] & (1 << other))
                        for other in placed
                    ),
                    vertex,
                )
                for vertex in range(len(adjacency))
                if vertex not in placed
            ]
            candidates = [item for item in candidates if item[0] >= D]
            if not candidates:
                break
            _degree, vertex = max(candidates)
            answer.append(vertex)
            placed.add(vertex)
            progressed = True
            if not forced:
                break
        if len(placed) == len(adjacency):
            break
        if not forced and progressed:
            continue
        circle = [
            vertex
            for vertex in range(len(adjacency))
            if vertex not in placed
            and sum(bool(adjacency[vertex] & (1 << other)) for other in placed)
            == D - 1
        ]
        if not circle:
            return None
        vertex = picker(circle, placed)
        answer.append(vertex)
        placed.add(vertex)
        forced = True
    return answer


def circle_positions(
    adjacency: Sequence[int], seed: Sequence[int], order: Sequence[int]
) -> list[int]:
    placed = set(seed)
    answer = []
    for position, vertex in enumerate(order):
        if sum(bool(adjacency[vertex] & (1 << other)) for other in placed) == D - 1:
            answer.append(position)
        placed.add(vertex)
    return answer


def generate_orders(adjacency: Sequence[int], maximum: int) -> list[tuple[tuple[int, ...], list[int]]]:
    seeds = clique_list(adjacency, 7) or clique_list(adjacency, 6)

    def largest_closure(candidates: list[int], placed: set[int]) -> int:
        return max(candidates, key=lambda vertex: closure_size(adjacency, vertex, placed))

    def smallest_closure(candidates: list[int], placed: set[int]) -> int:
        return min(candidates, key=lambda vertex: closure_size(adjacency, vertex, placed))

    pickers: list[Callable[[list[int], set[int]], int]] = [
        largest_closure,
        smallest_closure,
        lambda candidates, _placed: max(candidates),
        lambda candidates, _placed: min(candidates),
    ]
    generated = []
    seen = set()
    for seed in seeds:
        for picker in pickers:
            for force in (False, True):
                order = layered_order(adjacency, seed, picker, force)
                key = None if order is None else (seed, tuple(order))
                if order is not None and key not in seen:
                    seen.add(key)
                    generated.append((seed, order))

    def score(item: tuple[tuple[int, ...], list[int]]) -> tuple[int, list[int]]:
        positions = circle_positions(adjacency, *item)
        return len(positions), [-position for position in positions]

    generated.sort(key=score)
    if len(generated) <= maximum:
        return generated
    zero = [item for item in generated if score(item)[0] == 0]
    circles = [item for item in generated if score(item)[0] > 0]
    circles.sort(key=lambda item: circle_positions(adjacency, *item)[0])
    keep = zero[: maximum - min(len(circles), maximum // 2)] + circles[: maximum // 2]
    return keep[:maximum] if keep else generated[:maximum]


def expected_interval(part: int, slices: int) -> tuple[float, float]:
    return 2.0 * math.pi * part / slices, 2.0 * math.pi * (part + 1) / slices


def validate_slice_records(records: object, required: str | None = None) -> list[dict]:
    require(isinstance(records, list) and len(records) == 24, "slice record count")
    for part, record in enumerate(records):
        require(isinstance(record, dict) and record.get("part") == part, "slice part")
        lo, hi = expected_interval(part, 24)
        require(record.get("lo") == lo and record.get("hi") == hi, "slice tiling")
        require(record.get("status") in {"KILLED", "SURVIVORS", "ABORT"}, "slice status")
        if required is not None:
            require(record.get("status") == required, "required slice status")
        require(
            type(record.get("nodes")) is int
            and record["nodes"] >= 0
            and type(record.get("unresolved_cells")) is int
            and record["unresolved_cells"] >= 0,
            "slice counters",
        )
    require(records[0]["lo"] == 0.0 and records[-1]["hi"] == 2.0 * math.pi, "full circle endpoints")
    require(
        all(records[position]["hi"] == records[position + 1]["lo"] for position in range(23)),
        "slice gaps/overlaps",
    )
    return records


def replay_target(adjacency: Sequence[int], raw_result: dict, replay: bool) -> dict:
    orders = generate_orders(adjacency, 4)
    expected_orders = [
        {
            "seed": [1, 7, 8, 13, 14, 15, 16],
            "order": [5, 12, 17, 6, 11, 10, 9, 4, 3, 2, 0],
            "circle_positions": [0],
        },
        {
            "seed": [1, 7, 8, 13, 14, 15, 16],
            "order": [12, 6, 17, 11, 10, 9, 5, 4, 3, 2, 0],
            "circle_positions": [0],
        },
    ]
    observed_orders = [
        {
            "seed": list(seed),
            "order": list(order),
            "circle_positions": circle_positions(adjacency, seed, order),
        }
        for seed, order in orders
    ]
    require(observed_orders == expected_orders, "independent target order generation")
    require(
        raw_result.get("status") == "KILLED"
        and raw_result.get("ordinal") == TARGET_ORDINAL
        and raw_result.get("index") == TARGET_CLASS_INDEX
        and raw_result.get("population") == "K7"
        and raw_result.get("orders_available") == 2
        and raw_result.get("orders_tried") == 2
        and raw_result.get("winning_order") == 1
        and raw_result.get("winning_circles") == 1
        and raw_result.get("winning_seed") == expected_orders[1]["seed"]
        and raw_result.get("winning_placement_order") == expected_orders[1]["order"],
        "target result identity/winner",
    )
    validate_slice_records(raw_result.get("winning_slice_records"), "KILLED")
    require(
        raw_result.get("error_type") is None
        and raw_result.get("error_message") is None
        and raw_result.get("error_traceback") is None
        and float(raw_result.get("elapsed_seconds", -1)) >= 0,
        "target error/time fields",
    )
    if not replay:
        return {
            "status": "SKIPPED",
            "orders_generated": 2,
            "winning_slices_structurally_checked": 24,
        }
    counts: Counter[str] = Counter()
    total_nodes = 0
    unresolved_cells = 0
    attempts = []
    winner = None
    for position, (seed, order) in enumerate(orders):
        slice_records = []
        for part in range(24):
            lo, hi = expected_interval(part, 24)
            status, nodes, cells = cdriver6.decide6(
                adjacency,
                N,
                seed=seed,
                order=order,
                th0=(lo, hi),
                max_nodes=100_000,
            )
            counts[status] += 1
            total_nodes += nodes
            unresolved_cells += cells
            slice_records.append(
                {
                    "part": part,
                    "lo": lo,
                    "hi": hi,
                    "status": status,
                    "nodes": nodes,
                    "unresolved_cells": cells,
                }
            )
            if status != "KILLED":
                break
        attempts.append(
            {
                "order_position": position,
                "seed": list(seed),
                "placement_order": list(order),
                "circle_positions": circle_positions(adjacency, seed, order),
                "slice_records": slice_records,
            }
        )
        if len(slice_records) == 24 and all(record["status"] == "KILLED" for record in slice_records):
            winner = position
            require(slice_records == raw_result["winning_slice_records"], "winning slice replay")
            break
    require(winner == 1, "replayed winning order")
    require(
        counts
        == Counter(
            {
                "KILLED": raw_result["kernel_killed"],
                "SURVIVORS": raw_result["kernel_survivors"],
                "ABORT": raw_result["kernel_aborts"],
            }
        )
        and sum(counts.values()) == raw_result["kernel_calls"]
        and total_nodes == raw_result["total_nodes"]
        and unresolved_cells == raw_result["unresolved_cells"],
        "complete target replay counters",
    )
    return {
        "status": "PASS",
        "orders_generated": len(orders),
        "orders_tried": len(attempts),
        "kernel_calls": sum(counts.values()),
        "kernel_status_counts": dict(counts),
        "total_nodes": total_nodes,
        "unresolved_cells": unresolved_cells,
        "attempts": attempts,
        "attempts_sha256": stable_hash(attempts),
        "winning_slices_replayed": 24,
        "exact_status_node_cell_match": True,
    }


def control_graph(negative: bool) -> tuple[int, ...]:
    rows = [0] * N
    clique_order = 8 if negative else 7
    for vertex in range(clique_order):
        for other in range(vertex):
            rows[vertex] |= 1 << other
            rows[other] |= 1 << vertex
    for vertex in range(clique_order, N):
        for other in range(5):
            rows[vertex] |= 1 << other
            rows[other] |= 1 << vertex
    return validate_rows(rows, N, "synthetic control")


def replay_slice_control(raw: dict, negative: bool, replay: bool) -> int:
    rows = control_graph(negative)
    orders = generate_orders(rows, 12)
    require(bool(orders), "control orders")
    seed, order = orders[0]
    require(raw.get("seed") == list(seed) and raw.get("placement_order") == list(order), "control order")
    records = validate_slice_records(raw.get("records"), "KILLED" if negative else None)
    if not negative:
        require(all(record["status"] != "KILLED" for record in records), "positive slice control killed")
    require(raw.get("status_counts") == dict(Counter(record["status"] for record in records)), "control counts")
    if replay:
        for record in records:
            observed = cdriver6.decide6(
                rows,
                N,
                seed=seed,
                order=order,
                th0=(record["lo"], record["hi"]),
                max_nodes=1_000,
            )
            require(
                observed
                == (record["status"], record["nodes"], record["unresolved_cells"]),
                "slice control replay",
            )
    return len(records) if replay else 0


def independent_analyze_control(adjacency: Sequence[int], index: int, ordinal: int, replay: bool) -> dict | None:
    if not replay:
        return None
    orders = generate_orders(adjacency, 4)
    counts: Counter[str] = Counter()
    total_nodes = 0
    unresolved = 0
    for position, (seed, order) in enumerate(orders):
        killed = True
        for part in range(24):
            lo, hi = expected_interval(part, 24)
            status, nodes, cells = cdriver6.decide6(
                adjacency,
                N,
                seed=seed,
                order=order,
                th0=(lo, hi),
                max_nodes=1_000,
            )
            counts[status] += 1
            total_nodes += nodes
            unresolved += cells
            if status != "KILLED":
                killed = False
                break
        if killed:
            status = "KILLED"
            return {
                "ordinal": ordinal,
                "index": index,
                "population": "POSITIVE_CONTROL",
                "status": status,
                "winning_order": position,
                "winning_circles": len(circle_positions(adjacency, seed, order)),
                "orders_available": len(orders),
                "orders_tried": position + 1,
                "kernel_calls": sum(counts.values()),
                "kernel_killed": counts["KILLED"],
                "kernel_survivors": counts["SURVIVORS"],
                "kernel_aborts": counts["ABORT"],
                "total_nodes": total_nodes,
                "unresolved_cells": unresolved,
                "error_type": None,
                "error_message": None,
                "error_traceback": None,
                "winning_seed": list(seed),
                "winning_placement_order": list(order),
                "winning_slice_records": None,
            }
    status = "ABORT" if counts["ABORT"] else "UNRESOLVED"
    return {
        "ordinal": ordinal,
        "index": index,
        "population": "POSITIVE_CONTROL",
        "status": status,
        "winning_order": -1,
        "winning_circles": -1,
        "orders_available": len(orders),
        "orders_tried": len(orders),
        "kernel_calls": sum(counts.values()),
        "kernel_killed": counts["KILLED"],
        "kernel_survivors": counts["SURVIVORS"],
        "kernel_aborts": counts["ABORT"],
        "total_nodes": total_nodes,
        "unresolved_cells": unresolved,
        "error_type": None,
        "error_message": None,
        "error_traceback": None,
        "winning_seed": None,
        "winning_placement_order": None,
        "winning_slice_records": None,
    }


def verify_controls(configuration: dict, deletion_manifest: dict, replay: bool) -> dict:
    controls = configuration.get("kernel_controls", {})
    require(controls.get("status") == "PASS" and controls.get("graph_order") == N, "control envelope")
    flexible = controls.get("flexible_circle_and_K8", {})
    require(flexible.get("status") == "PASS" and flexible.get("production_slices_exercised") == 24, "slice controls")
    replayed_slices = 0
    replayed_slices += replay_slice_control(flexible["positive"], False, replay)
    replayed_slices += replay_slice_control(flexible["negative"], True, replay)
    standard = load_json(ROOT / "d6_standard18_geometry_report.json")
    standard_check = load_json(ROOT / "d6_standard18_geometry_verification.json")
    exact = load_json(ROOT / "d6_18_exact_reconstruction.json")
    exact_check = load_json(ROOT / "d6_18_exact_reconstruction_verification.json")
    require(
        standard.get("status") == "PASS"
        and standard_check.get("status") == "PASS"
        and exact.get("status") == "PASS"
        and exact_check.get("status") == "PASS"
        and len(exact.get("configurations", ())) == 2,
        "exact positive input boundary",
    )
    control_inputs = [
        (
            "selected_standard18_compatible_deletion",
            441,
            deletion_manifest["unique_deletions"][441]["adjacency"],
        ),
        ("standard18_full_unit_graph", -101, standard["unit_graph"]["adjacency"]),
        (
            "nonstandard_exact18_1",
            -102,
            exact["configurations"][0]["exact_geometry"]["unit_graph_adjacency"],
        ),
        (
            "nonstandard_exact18_2",
            -103,
            exact["configurations"][1]["exact_geometry"]["unit_graph_adjacency"],
        ),
    ]
    raw_exact = controls.get("known_realizable_exact18")
    require(isinstance(raw_exact, list) and len(raw_exact) == 4, "exact control records")
    checked = []
    for ordinal, ((name, index, values), raw) in enumerate(zip(control_inputs, raw_exact, strict=True)):
        adjacency = validate_rows(values, N, name)
        require(raw.get("name") == name and raw.get("status") != "KILLED", f"{name} raw status")
        observed = independent_analyze_control(adjacency, index, ordinal, replay)
        if replay:
            require(observed == raw.get("result"), f"{name} exact replay")
        else:
            require(raw.get("result", {}).get("status") != "KILLED", f"{name} structural")
        checked.append({"name": name, "status": raw["status"]})
    require(
        controls.get("positive_source_hashes")
        == {name: INPUTS[name] for name in sorted(INPUTS) if name not in {
            "d6_current_residue_manifest_v7.json",
            "d6_current_residue_manifest_v7_verification.json",
            "d6_current_residue_manifest_v8.json",
            "d6_current_residue_manifest_v8_verification.json",
        }},
        "positive source hashes",
    )
    return {
        "status": "PASS" if replay else "STRUCTURAL_PASS_REPLAY_SKIPPED",
        "slice_control_records_checked": 48,
        "slice_control_records_replayed": replayed_slices,
        "exact_positive_controls_checked": checked,
        "exact_positive_controls_replayed": 4 if replay else 0,
        "recorded_controls_sha256": stable_hash(controls),
    }


def containment_scan(v8_graphs: list[dict], target: Sequence[int]) -> dict:
    """Independent complement-deletion enumeration, not producer backtracking."""

    target_edges = sum(row.bit_count() for row in target) // 2
    target_complement = complement(target)
    complement_edges = [
        (left, right)
        for right in range(1, N)
        for left in range(right)
        if target_complement[left] & (1 << right)
    ]
    candidates = []
    maximum = 0
    for parent in v8_graphs:
        rows = parent["adjacency"]
        edge_count = sum(row.bit_count() for row in rows) // 2
        for deleted in range(19):
            deletion_edges = edge_count - rows[deleted].bit_count()
            maximum = max(maximum, deletion_edges)
            if deletion_edges >= target_edges:
                candidates.append(
                    (
                        parent,
                        deleted,
                        deletion_edges,
                        complement(delete_vertex(rows, deleted)),
                    )
                )
    require(maximum == 109 and len(candidates) == 1_085, "v8 eligible deletion profile")
    maximum_removed = maximum - target_edges
    variants: dict[tuple[int, tuple[object, ...]], list[tuple[int, ...]]] = defaultdict(list)
    variant_count = 0
    for removed_count in range(maximum_removed + 1):
        for removed in itertools.combinations(complement_edges, removed_count):
            rows = list(target_complement)
            for left, right in removed:
                rows[left] ^= 1 << right
                rows[right] ^= 1 << left
            variant = tuple(rows)
            variants[(removed_count, graph_invariant(variant))].append(variant)
            variant_count += 1
    occurrences = []
    exact_comparisons = 0
    nonempty_buckets = 0
    for parent, deleted, deletion_edges, candidate_complement in candidates:
        key = (deletion_edges - target_edges, graph_invariant(candidate_complement))
        bucket = variants.get(key, ())
        if bucket:
            nonempty_buckets += 1
        matched = False
        for variant in bucket:
            exact_comparisons += 1
            if exact_isomorphic(candidate_complement, variant):
                matched = True
                break
        if matched:
            occurrences.append(
                {
                    "parent_class": parent["class"],
                    "parent_index": parent["index"],
                    "deleted_vertex": deleted,
                    "deletion_edges": deletion_edges,
                }
            )
    expected = [
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
    ]
    require(occurrences == expected, "v8 containment occurrences")
    require(variant_count == 17_344, "containment variant count")
    return {
        "parents_scanned": len(v8_graphs),
        "deletions_scanned": len(v8_graphs) * 19,
        "eligible_deletions": len(candidates),
        "maximum_deletion_edges": maximum,
        "target_edges": target_edges,
        "target_complement_edges": len(complement_edges),
        "deleted_complement_edge_variants": variant_count,
        "invariant_buckets": len(variants),
        "candidate_nonempty_buckets": nonempty_buckets,
        "exact_isomorphism_comparisons": exact_comparisons,
        "occurrences": occurrences,
        "occurrences_sha256": stable_hash(occurrences),
        "rejected_parent_indices": [TARGET_PARENT_INDEX],
        "method": (
            "enumerate every target-complement graph after deleting 0..3 "
            "edges, bucket by independent 1-WL invariants, then exact "
            "isomorphism backtracking"
        ),
    }


def read_evidence(report: dict, evidence_path: Path) -> tuple[dict, dict, dict]:
    evidence_bytes = evidence_path.read_bytes()
    raw_meta = report.get("raw_evidence", {})
    require(
        raw_meta.get("path") == evidence_path.name
        and raw_meta.get("sha256") == bytes_sha256(evidence_bytes)
        and raw_meta.get("campaign_sha256") == RAW_CAMPAIGN_SHA256
        and raw_meta.get("result_sha256") == RAW_RESULT_SHA256
        and raw_meta.get("config_sha256") == CONFIG_SHA256
        and raw_meta.get("evidence_contains_exact_raw_utf8") is True,
        "report/evidence binding",
    )
    evidence = parse_json_raw(gzip.decompress(evidence_bytes), "evidence")
    require(evidence.get("schema") == "d6-interval-18-class2100-raw-evidence-v1", "evidence schema")
    campaign_record = evidence.get("campaign", {})
    result_record = evidence.get("result", {})
    require(
        campaign_record.get("name") == "campaign.json"
        and campaign_record.get("sha256") == RAW_CAMPAIGN_SHA256
        and result_record.get("name") == "result_000032_0002100.json"
        and result_record.get("sha256") == RAW_RESULT_SHA256,
        "evidence member metadata",
    )
    campaign_raw = campaign_record.get("raw_utf8", "").encode("utf-8")
    result_raw = result_record.get("raw_utf8", "").encode("utf-8")
    require(bytes_sha256(campaign_raw) == RAW_CAMPAIGN_SHA256, "embedded campaign raw hash")
    require(bytes_sha256(result_raw) == RAW_RESULT_SHA256, "embedded result raw hash")
    campaign = parse_json_raw(campaign_raw, "embedded campaign")
    checkpoint = parse_json_raw(result_raw, "embedded result")
    return evidence, campaign, checkpoint


def verify_increment(
    report_path: Path,
    evidence_path: Path,
    *,
    expected_report_sha256: str | None = None,
    replay: bool = True,
    enforce_source_boundary: bool = True,
) -> dict:
    report_path = report_path.resolve()
    evidence_path = evidence_path.resolve()
    report_hash = sha256(report_path)
    if expected_report_sha256 is not None:
        require(report_hash == expected_report_sha256, "report hash")
    report = load_json(report_path)
    require(
        report.get("schema") == "d6-interval-18-class2100-increment-v1"
        and report.get("status") == "COMPLETE_INTERVAL_REJECTION"
        and report.get("trust_assumptions") == TRUST_ASSUMPTIONS
        and report.get("search") == {**EXPECTED_SEARCH, "workers": 3},
        "report envelope",
    )
    package = verify_package_boundary(report, enforce_source_boundary)
    _evidence, campaign, checkpoint = read_evidence(report, evidence_path)
    configuration = campaign.get("configuration")
    require(
        campaign.get("schema") == 1
        and campaign.get("config_sha256") == CONFIG_SHA256
        and isinstance(configuration, dict)
        and stable_hash(configuration) == CONFIG_SHA256
        and configuration.get("graph_order") == N
        and configuration.get("search") == EXPECTED_SEARCH
        and configuration.get("trust_assumptions") == TRUST_ASSUMPTIONS,
        "embedded campaign configuration",
    )
    require(
        checkpoint.get("schema") == 1
        and checkpoint.get("config_sha256") == CONFIG_SHA256
        and checkpoint.get("search_slices") == 24,
        "embedded result envelope",
    )
    static = verify_static_files(configuration, report)
    deletion_manifest = load_json(ROOT / "d6_residue_18_deletions_v2.json")
    deletion_check = load_json(ROOT / "d6_residue_18_deletions_v2_verification.json")
    require(
        deletion_manifest.get("schema") == "d6-residue-18-deletion-manifest-v2"
        and deletion_manifest.get("summary", {}).get("parents") == 911
        and deletion_manifest.get("summary", {}).get("unique_deletions") == 11_975
        and deletion_manifest.get("semantics", {}).get("nonedges")
        == "unconstrained and may also have distance one"
        and deletion_check.get("schema") == "d6-residue-18-deletion-verification-v2"
        and deletion_check.get("status") == "PASS"
        and deletion_check.get("manifest", {}).get("sha256")
        == INPUTS["d6_residue_18_deletions_v2.json"]
        and all(deletion_check.get("checks", {}).values()),
        "deletion corpus boundary",
    )
    target_record = deletion_manifest["unique_deletions"][TARGET_CLASS_INDEX]
    target = validate_rows(target_record.get("adjacency"), N, "target class")
    require(
        target_record.get("class_id") == TARGET_CLASS_ID
        and target_record.get("canonical_graph6") == TARGET_GRAPH6
        and hashlib.sha256(TARGET_GRAPH6.encode("ascii")).hexdigest()
        == TARGET_CLASS_ID.removeprefix("u18-")
        and decode_graph6(TARGET_GRAPH6) == target
        and sum(row.bit_count() for row in target) // 2 == target_record.get("edges") == 106
        and target_record.get("clique_number") == 7
        and target_record.get("standard18_compatible") is False,
        "target class reconstruction",
    )
    report_class = report.get("deletion_class", {})
    require(
        report_class.get("class_index") == TARGET_CLASS_INDEX
        and report_class.get("ordinal_in_production_run") == TARGET_ORDINAL
        and report_class.get("class_id") == TARGET_CLASS_ID
        and report_class.get("canonical_graph6") == TARGET_GRAPH6
        and report_class.get("adjacency") == list(target)
        and report_class.get("adjacency_sha256") == stable_hash(list(target)),
        "report target class",
    )
    v7_graphs = load_residue(
        "d6_current_residue_manifest_v7.json",
        EXPECTED_V7_COUNTS,
        EXPECTED_V7_CLASS_HASHES,
        EXPECTED_V7_COMBINED_HASH,
    )
    production = reconstruct_production(v7_graphs, deletion_manifest, configuration)
    occurrences = reconstruct_occurrences(v7_graphs, target, target_record)
    require(
        report_class.get("parent_occurrences") == occurrences["exact_isomorphism_occurrences"]
        and report_class.get("parent_occurrences_sha256") == occurrences["occurrences_sha256"],
        "report occurrence reconstruction",
    )
    raw_result = checkpoint.get("result")
    require(isinstance(raw_result, dict), "raw result object")
    report_certificate = report.get("interval_certificate", {})
    require(
        report_certificate.get("status") == "KILLED"
        and report_certificate.get("raw_result") == raw_result
        and report_certificate.get("raw_result_sha256") == stable_hash(raw_result)
        and report_certificate.get("winning_slice_records")
        == raw_result.get("winning_slice_records")
        and report_certificate.get("winning_slice_records_sha256")
        == stable_hash(raw_result.get("winning_slice_records")),
        "report/raw interval certificate",
    )
    target_replay = replay_target(target, raw_result, replay)
    controls = verify_controls(configuration, deletion_manifest, replay)
    require(
        report.get("controls", {}).get("recorded_status") == "PASS"
        and report.get("controls", {}).get("recorded_controls_sha256")
        == controls["recorded_controls_sha256"]
        and report.get("controls", {}).get("independent_replay_required") is True,
        "report controls",
    )
    v8_graphs = load_residue(
        "d6_current_residue_manifest_v8.json",
        EXPECTED_V8_COUNTS,
        EXPECTED_V8_CLASS_HASHES,
        EXPECTED_V8_COMBINED_HASH,
    )
    containment = containment_scan(v8_graphs, target)
    report_containment = report.get("current_v8_containment", {})
    require(
        report_containment.get("parents_scanned") == containment["parents_scanned"]
        and report_containment.get("deletions_scanned") == containment["deletions_scanned"]
        and report_containment.get("eligible_deletions") == containment["eligible_deletions"]
        and report_containment.get("maximum_deletion_edges")
        == containment["maximum_deletion_edges"]
        and report_containment.get("occurrences") == containment["occurrences"]
        and report_containment.get("occurrences_sha256")
        == containment["occurrences_sha256"]
        and report_containment.get("rejected_parent_indices") == [TARGET_PARENT_INDEX],
        "producer/checker containment agreement",
    )
    conclusion = report.get("conclusion", {})
    require(
        conclusion.get("rejected_parent_indices") == [TARGET_PARENT_INDEX]
        and conclusion.get("rejected_parent_indices_sha256")
        == stable_hash([TARGET_PARENT_INDEX])
        and conclusion.get("new_interval_rejections") == 1,
        "report conclusion",
    )
    semantics = report.get("semantics", {})
    require(
        semantics.get("candidate_nonedges")
        == "unconstrained and may also have distance one"
        and semantics.get("distinct_points_required") is True
        and semantics.get("only_KILLED_is_a_rejection") is True,
        "report semantics",
    )
    return {
        "schema": "d6-interval-18-class2100-increment-verification-v1",
        "status": "PASS",
        "claim": (
            "Deletion class 2100 is certified non-realizable in R6 by a "
            "complete outward-rounded interval replay; parent 3945564 is "
            "therefore rejected by required-edge subgraph monotonicity."
        ),
        "report": {"path": report_path.name, "sha256": report_hash},
        "evidence": {"path": evidence_path.name, "sha256": sha256(evidence_path)},
        "raw_checkpoint": {
            "campaign_sha256": RAW_CAMPAIGN_SHA256,
            "result_sha256": RAW_RESULT_SHA256,
            "configuration_sha256": CONFIG_SHA256,
        },
        "source_boundary": {"package": package, "campaign": static},
        "selection_reconstruction": production,
        "class_reconstruction": {
            "class_index": TARGET_CLASS_INDEX,
            "class_id": TARGET_CLASS_ID,
            "canonical_graph6_decoded": True,
            "edges": 106,
            "occurrences": occurrences,
        },
        "interval_replay": target_replay,
        "controls": controls,
        "current_v8_containment": containment,
        "conclusion": {
            "certified_deletion_classes": [TARGET_CLASS_INDEX],
            "rejected_parent_indices": [TARGET_PARENT_INDEX],
            "rejected_parent_indices_sha256": stable_hash([TARGET_PARENT_INDEX]),
            "interval_dependent_rejections": 1,
        },
        "trust_assumptions": TRUST_ASSUMPTIONS,
        "execution": {
            "verified_utc": datetime.now(UTC).isoformat(),
            "python": platform.python_version(),
            "platform": platform.platform(),
            "verifier_source_sha256": sha256(Path(__file__).resolve()),
            "replay_requested": replay,
        },
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
    parser.add_argument("--evidence", type=Path, required=True)
    parser.add_argument("--expected-report-sha256")
    parser.add_argument("--no-replay", action="store_true")
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "d6_interval_18_class2100_increment_verification.json",
    )
    args = parser.parse_args()
    result = verify_increment(
        args.report,
        args.evidence,
        expected_report_sha256=args.expected_report_sha256,
        replay=not args.no_replay,
    )
    atomic_json(args.output.resolve(), result)
    print(
        json.dumps(
            {
                "status": result["status"],
                "rejected_parent_indices": [TARGET_PARENT_INDEX],
                "replay": result["interval_replay"]["status"],
                "output": str(args.output),
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
