#!/usr/bin/env python3
"""Independent verifier for the complete pattern-954 containment report.

This file deliberately does not import the containment runner, its C kernel,
the extractor, or the earlier certificate verifier.  It independently:

* reconstructs the 13-vertex required-edge pattern from the n=14 corpus;
* rebuilds the exact 12,839-target selection from the parent selection and
  degree-one dual decisions, and checks every target against the rank input;
* audits the report, campaign, all 51 immutable chunks, compiler metadata,
  and pinned source/binary hashes;
* checks injectivity and all 54 required pattern edges for every one of the
  3,403 reported HIT mappings; and
* records the ordered 9,436-index complement as the new exact residue.

It does not rerun the C search and therefore does not independently establish
the completeness of a NO_HIT result.  The exact mathematical reduction uses
only the independently checked HIT mappings; its residue is their complement.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import shlex
import subprocess
from collections import Counter
from datetime import datetime
from itertools import combinations
from pathlib import Path
from typing import Sequence


ROOT = Path(__file__).resolve().parent
REPORT_NAME = "d6_n14_pattern_954_containment_report.json"
REPORT_SHA256 = "a2d4fb5a07dcea580315730b4155fc0c08ecd0439a9111e733d29e2da0b5be62"
CONFIG_SHA256 = "42410c706b670875f821a40dc3bb42bb19e8c138c7f3387e3faa4a1dd5ee1913"
CAMPAIGN_SHA256 = "afdeb5bea296edb30d9ff65cbe4838f654dd92fa446b9cdf9bda28b06d8c8710"
CHUNK_MANIFEST_SHA256 = "a9c2d3b8c3737b5bfcdfcd4f3d16041666342fdc3959d3f4ad5f73aae4574ea6"
CHUNK_COUNT = 51
CHUNK_BYTES = 2_229_410
TARGET_COUNT = 12_839
HIT_COUNT = 3_403
RESIDUE_COUNT = 9_436
HIT_INDICES_SHA256 = "24e88173a5bef7cdc191cb610b4048e662964066312694a94b180cc5a21edc7e"
RESIDUE_INDICES_SHA256 = "b740e9c14cab4c204abe4e5e8b9b07c1cbc7a5e037f7feab685df0573b7250db"
TARGET_INDICES_SHA256 = "320a68221ef597a1252d0c16cba6b30b4bef2bb6c13b73ce529e860ed431d497"
TARGET_TSV_SHA256 = "9355854172c324f9d93cc4085020974fb9622a010b2e8fe19f5a954d17a7277d"
PARENT_SELECTION_SHA256 = "814c80651746ce05048f72f0cfa7e49a921554abaf167bd5c8cf4063a34d35c0"
DUAL_REPORT_SHA256 = "c2de7e06bbe4f3d163a4f747d42721f7f60da985b867ce74669664431738e345"
DUAL_DECISIONS_SHA256 = "724a97928422fc8705946ed64e3ce9afd37579d3229e5b3ebb7c4f1144c002ec"
RANK_INPUT_SHA256 = "7a0a350142930e4e13830ea44c4fd217f25aeb07bbd6286e51e316f72712c4f9"
N14_CORPUS_SHA256 = "0e3d74c081b272731848d655da0c68cfba09435e39a2fd2107c32c2bd3b378f0"
BINARY_SHA256 = "019cc3a098b1fc66e6667dc0b4b5d01e1e55a2bf3b4400e6a06cdf232204d028"
EXPECTED_CORE = (
    2842, 6777, 3496, 8039, 4931, 7374, 8122,
    5476, 7901, 7515, 7148, 5999, 4090,
)
ORIGINAL_CORE_VERTICES = (0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13)
REMOVED_SOURCE_EDGES = {frozenset((0, 3)), frozenset((5, 8))}
EXPECTED_BUNDLE_HASHES = {
    "d6_n14_pattern_954_extract.py": "972d522d1efd3f38cb7f2d95cdb3e0c2ffd4008268071ab0b0dab18343d1d320",
    "d6_n14_pattern_954_input.json": "e302b818747a08e07fbd137211c05886a7e290cd4b561c1e977d4e0830b953a0",
    "d6_n14_pattern_954_verify.py": "a930033be5cc93c4a93badb3f06e2c34b675105c99a55b195ed7fee1f8cb8309",
    "d6_n14_pattern_954_test.py": "bb8eb3340f8519027b12ba9b43dd02b74402f3f3aac27ddac54678b778bb2e95",
    "d6_n14_pattern_954_certificate.md": "cb3e63a9e0daae63b52037c38e5babd693c736a87b636578061ac4dadcae86e1",
    "d6_n14_pattern_954_containment.c": "038488956021093cfed039973b1f6ad3e4e2e73c8beb7bc0b94c9619ad8021f6",
    "d6_n14_pattern_954_containment.py": "994cd406de16ca4791789900991971c3fcbf25b06f43add73defb972b72da24c",
    "d6_n14_pattern_954_containment_test.py": "acb0e8756ba8b4de616160bff04341cbf690a320e85ecdfe79e97443ebcee474",
}
EXPECTED_COMPILER = {
    "requested": "cc",
    "path": "/usr/bin/cc",
    "version": "Apple clang version 16.0.0 (clang-1600.0.26.6)",
}
EXPECTED_FLAGS = [
    "-O3", "-std=c11", "-Wall", "-Wextra", "-pedantic", "-shared",
    "-Wl,-no_uuid", "-Wl,-install_name,d6_n14_pattern_954_containment.dylib",
]
EXPECTED_CONFIG = {
    "schema": 1,
    "kind": "d6_n14_pattern_954_containment_config",
    "required_edge_only": True,
    "input_json_sha256": EXPECTED_BUNDLE_HASHES["d6_n14_pattern_954_input.json"],
    "c_source_sha256": EXPECTED_BUNDLE_HASHES["d6_n14_pattern_954_containment.c"],
    "binary_sha256": BINARY_SHA256,
    "target_tsv_sha256": TARGET_TSV_SHA256,
    "full_target_indices_sha256": TARGET_INDICES_SHA256,
    "full_target_population": TARGET_COUNT,
    "selected_indices_sha256": TARGET_INDICES_SHA256,
    "selected_population": TARGET_COUNT,
    "node_limit": 50_000,
    "chunk_size": 256,
    "kernel_status_semantics": {
        "HIT": "mapping independently edge-checked; exact rejection",
        "NO_HIT": "exhaustive kernel search ended within node limit",
        "TIMEOUT": "unresolved; never a rejection",
        "INFRA_ERROR": "unresolved; never a rejection",
    },
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_hash(value: object) -> str:
    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


def atomic_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    with temporary.open("w", encoding="ascii", newline="") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def validate_graph(adjacency: Sequence[int], order: int) -> tuple[int, ...]:
    require(len(adjacency) == order, f"graph has order {len(adjacency)}, expected {order}")
    answer = tuple(adjacency)
    full = (1 << order) - 1
    for vertex, row in enumerate(answer):
        require(isinstance(row, int) and not isinstance(row, bool), "noninteger row")
        require(not (row & ~full), "out-of-range adjacency bit")
        require(not (row & (1 << vertex)), "adjacency loop")
        for other in range(vertex):
            require(
                bool(row & (1 << other)) == bool(answer[other] & (1 << vertex)),
                "asymmetric adjacency",
            )
    return answer


def reconstruct_core(corpus: Path) -> tuple[int, ...]:
    require(sha256(corpus) == N14_CORPUS_SHA256, "n=14 corpus hash mismatch")
    lines = corpus.read_text(encoding="ascii").splitlines()
    require(len(lines) == 1_052, "n=14 corpus population mismatch")
    fields = tuple(map(int, lines[954].split()))
    require(len(fields) == 15 and fields[0] == 14, "source pattern row malformed")
    source = validate_graph(fields[1:], 14)
    positions = {vertex: offset for offset, vertex in enumerate(ORIGINAL_CORE_VERTICES)}
    core = [0] * 13
    for first, second in combinations(ORIGINAL_CORE_VERTICES, 2):
        if not (source[first] & (1 << second)):
            continue
        if frozenset((first, second)) in REMOVED_SOURCE_EDGES:
            continue
        left, right = positions[first], positions[second]
        core[left] |= 1 << right
        core[right] |= 1 << left
    answer = validate_graph(core, 13)
    require(answer == EXPECTED_CORE, "independently reconstructed core changed")
    require(sum(row.bit_count() for row in answer) // 2 == 54, "core edge count changed")
    return answer


def rebuild_targets() -> tuple[list[int], dict[int, tuple[int, ...]], dict]:
    parent_path = ROOT / "d6_k7_positive_dual_selection.json"
    dual_report_path = ROOT / "d6_k7_positive_dual_full_report.json"
    decisions_path = ROOT / "d6_k7_positive_dual_full_decisions.tsv.gz"
    rank_path = ROOT / ".runs/d6_k7_rank_survivors.json"
    tsv_path = ROOT / ".runs/d6_n14_pattern_954_targets_12839.tsv"
    expected = {
        parent_path: PARENT_SELECTION_SHA256,
        dual_report_path: DUAL_REPORT_SHA256,
        decisions_path: DUAL_DECISIONS_SHA256,
        rank_path: RANK_INPUT_SHA256,
        tsv_path: TARGET_TSV_SHA256,
    }
    for path, digest in expected.items():
        require(sha256(path) == digest, f"hash mismatch: {path.name}")

    parent_payload = json.loads(parent_path.read_text(encoding="utf-8"))
    parent = parent_payload["selected_indices"]
    require(parent_payload.get("selected") == 12_941, "parent selection count mismatch")
    require(len(parent) == 12_941 and len(set(parent)) == len(parent),
            "parent selection malformed")
    dual_report = json.loads(dual_report_path.read_text(encoding="utf-8"))
    summary = dual_report["summary"]
    require(summary["graphs"] == summary["complete"] == 12_941,
            "dual report coverage mismatch")
    require(summary["marginal_dual_rejected"] == 102, "dual rejection count mismatch")
    require(summary["survivors"] == TARGET_COUNT and summary["infra_errors"] == 0,
            "dual residue accounting mismatch")
    with gzip.open(decisions_path, "rt", encoding="ascii", newline="") as stream:
        decisions = list(csv.DictReader(stream, delimiter="\t"))
    require([int(row["index"]) for row in decisions] == parent,
            "dual decision order/coverage mismatch")
    require(all(row["status"] in ("REJECTED", "SURVIVOR") for row in decisions),
            "dual decision has invalid status")
    rejected = [int(row["index"]) for row in decisions if row["status"] == "REJECTED"]
    selected = [int(row["index"]) for row in decisions if row["status"] == "SURVIVOR"]
    require(rejected == summary["marginal_dual_rejected_indices"],
            "dual report/archive rejection order mismatch")
    require(len(selected) == TARGET_COUNT, "rebuilt target selection count mismatch")
    require(stable_hash(selected) == TARGET_INDICES_SHA256,
            "rebuilt target index hash mismatch")

    rank_payload = json.loads(rank_path.read_text(encoding="utf-8"))
    rank_rows = rank_payload["graphs"]
    require(len(rank_rows) == 17_764, "rank input population mismatch")
    rank_by_index = {
        int(row["index"]): validate_graph(tuple(map(int, row["adjacency"])), 19)
        for row in rank_rows
    }
    require(len(rank_by_index) == len(rank_rows), "rank input repeats an index")
    observed: list[int] = []
    target_by_index: dict[int, tuple[int, ...]] = {}
    with tsv_path.open("r", encoding="ascii") as stream:
        for line_number, line in enumerate(stream, 1):
            fields = tuple(map(int, line.split()))
            require(len(fields) == 20, f"malformed target TSV row {line_number}")
            index = fields[0]
            adjacency = validate_graph(fields[1:], 19)
            require(index not in target_by_index, "target TSV repeats an index")
            require(rank_by_index.get(index) == adjacency,
                    "target TSV adjacency differs from rank input")
            observed.append(index)
            target_by_index[index] = adjacency
    require(observed == selected, "target TSV differs from rebuilt selection")
    return selected, target_by_index, {
        "parent_population": len(parent),
        "dual_rejections": len(rejected),
        "population": len(selected),
        "indices_sha256": stable_hash(selected),
        "tsv_sha256": sha256(tsv_path),
    }


def verify_mapping(
    pattern: Sequence[int], target: Sequence[int], mapping: object
) -> int:
    require(isinstance(mapping, list) and len(mapping) == len(pattern),
            "HIT mapping has wrong type/length")
    require(all(isinstance(value, int) and not isinstance(value, bool)
                for value in mapping), "HIT mapping contains noninteger")
    require(all(0 <= value < len(target) for value in mapping),
            "HIT mapping has out-of-range vertex")
    require(len(set(mapping)) == len(mapping), "HIT mapping is not injective")
    checked = 0
    for first in range(len(pattern)):
        for second in range(first):
            if pattern[first] & (1 << second):
                require(target[mapping[first]] & (1 << mapping[second]),
                        f"HIT mapping loses required edge ({first},{second})")
                checked += 1
    require(checked == 54, "wrong number of pattern edges checked")
    return checked


def validate_result_record(
    result: dict,
    expected_index: int,
    pattern: Sequence[int],
    target: Sequence[int],
) -> int:
    require(isinstance(result, dict), "chunk result is not an object")
    require(result.get("index") == expected_index, "chunk result index/order mismatch")
    status = result.get("status")
    require(status in ("HIT", "NO_HIT", "TIMEOUT", "INFRA_ERROR"),
            "chunk result has invalid status")
    nodes = result.get("nodes")
    require(isinstance(nodes, int) and not isinstance(nodes, bool) and nodes >= 0,
            "chunk result has invalid node count")
    wall = result.get("wall_seconds")
    require(isinstance(wall, (int, float)) and not isinstance(wall, bool)
            and math.isfinite(wall) and wall >= 0, "invalid target wall time")
    if status in ("HIT", "NO_HIT"):
        require(nodes <= 50_000, "terminal result exceeds node limit")
    if status == "TIMEOUT":
        require(nodes > 50_000, "TIMEOUT did not exhaust node limit")
    if status == "HIT":
        require(set(result) == {"index", "mapping", "nodes", "status", "wall_seconds"},
                "HIT result fields changed")
        return verify_mapping(pattern, target, result["mapping"])
    if status in ("NO_HIT", "TIMEOUT"):
        require(set(result) == {"index", "nodes", "status", "wall_seconds"},
                f"{status} result fields changed")
    if status == "INFRA_ERROR":
        require("error" in result, "INFRA_ERROR has no diagnostic")
    return 0


def chunk_manifest(checkpoint_dir: Path) -> list[dict]:
    files = sorted(checkpoint_dir.glob("chunk_*.json"))
    return [
        {"name": path.name, "sha256": sha256(path), "bytes": path.stat().st_size}
        for path in files
    ]


def audit_chunks(
    checkpoint_dir: Path,
    selected: Sequence[int],
    targets: dict[int, tuple[int, ...]],
    pattern: Sequence[int],
) -> tuple[list[dict], list[dict], int]:
    manifest = chunk_manifest(checkpoint_dir)
    require(len(manifest) == CHUNK_COUNT, "wrong number of chunk artifacts")
    require(sum(row["bytes"] for row in manifest) == CHUNK_BYTES,
            "chunk byte total mismatch")
    require(stable_hash(manifest) == CHUNK_MANIFEST_SHA256,
            "immutable chunk manifest hash mismatch")
    expected_names = []
    results: list[dict] = []
    edges_checked = 0
    for start in range(0, TARGET_COUNT, 256):
        stop = min(start + 256, TARGET_COUNT)
        name = f"chunk_{start:06d}_{stop:06d}.json"
        expected_names.append(name)
        path = checkpoint_dir / name
        payload = json.loads(path.read_text(encoding="utf-8"))
        require(payload.get("schema") == 1 and payload.get("kind") == "pattern954_chunk",
                "chunk schema/kind mismatch")
        require(payload.get("config_sha256") == CONFIG_SHA256,
                "chunk configuration hash mismatch")
        require(payload.get("start") == start and payload.get("stop") == stop,
                "chunk range mismatch")
        datetime.fromisoformat(payload["created_utc"])
        rows = payload.get("results")
        require(isinstance(rows, list) and len(rows) == stop - start,
                "chunk result population mismatch")
        for result, expected_index in zip(rows, selected[start:stop]):
            edges_checked += validate_result_record(
                result, expected_index, pattern, targets[expected_index]
            )
        results.extend(rows)
    require([row["name"] for row in manifest] == expected_names,
            "chunk filenames/ranges are not an exact partition")
    require(len(results) == TARGET_COUNT, "chunks do not cover full selection")
    return manifest, results, edges_checked


def verify_compilation(report: dict) -> dict:
    compilation = report["compilation"]
    require(compilation == json.loads(
        (Path(report["checkpoint_directory"]) / "campaign.json").read_text(
            encoding="utf-8"
        )
    )["compilation"], "campaign/report compilation metadata differ")
    require(compilation["compiler"] == EXPECTED_COMPILER, "compiler identity mismatch")
    completed = subprocess.run(
        [EXPECTED_COMPILER["path"], "--version"],
        check=True, capture_output=True, text=True,
    )
    first_line = (completed.stdout or completed.stderr).splitlines()[0]
    require(first_line == EXPECTED_COMPILER["version"], "installed compiler changed")
    require(compilation["flags"] == EXPECTED_FLAGS, "compile flags mismatch")
    require(compilation["stdout"] == compilation["stderr"] == "",
            "compile diagnostics were not clean")
    command = shlex.split(compilation["command"])
    require(command[0] == EXPECTED_COMPILER["path"], "compile command path mismatch")
    require(command[1:1 + len(EXPECTED_FLAGS)] == EXPECTED_FLAGS,
            "compile command flags mismatch")
    require(command[-1] == str((ROOT / "d6_n14_pattern_954_containment.c").resolve()),
            "compile command source mismatch")
    require(command[-3] == "-o", "compile command output syntax mismatch")
    require(Path(command[-2]).name.startswith("d6_n14_pattern_954_containment.dylib.tmp."),
            "compile temporary output name mismatch")
    post = [shlex.split(value) for value in compilation["postprocess_commands"]]
    require(len(post) == 2 and post[0][:2] == ["codesign", "--remove-signature"],
            "signature-removal command mismatch")
    require(post[1][:4] == ["codesign", "--force", "--sign", "-"],
            "signing command mismatch")
    require(post[0][-1] == command[-2] == post[1][-1],
            "compile/sign temporary targets differ")
    binary = ROOT / ".runs/d6_n14_pattern_954_containment.dylib"
    require(sha256(binary) == BINARY_SHA256, "compiled binary hash mismatch")
    require(binary.stat().st_size == compilation["binary_bytes"] == 51_648,
            "compiled binary size mismatch")
    require(compilation["binary_sha256"] == BINARY_SHA256,
            "reported compiled binary hash mismatch")
    signature = subprocess.run(
        ["codesign", "--verify", "--strict", str(binary)],
        capture_output=True, text=True,
    )
    require(signature.returncode == 0, "compiled binary signature verification failed")
    return {
        "compiler": EXPECTED_COMPILER,
        "flags": EXPECTED_FLAGS,
        "binary_sha256": sha256(binary),
        "binary_bytes": binary.stat().st_size,
        "code_signature_valid": True,
    }


def verify(report_path: Path) -> dict:
    require(sha256(report_path) == REPORT_SHA256, "containment report hash mismatch")
    for filename, digest in EXPECTED_BUNDLE_HASHES.items():
        require(sha256(ROOT / filename) == digest, f"bundle hash mismatch: {filename}")
    input_payload = json.loads(
        (ROOT / "d6_n14_pattern_954_input.json").read_text(encoding="utf-8")
    )
    pattern = reconstruct_core(ROOT / "aeq_d6_n14.txt")
    require(tuple(input_payload["obstruction_core"]["adjacency"]) == pattern,
            "input JSON core differs from corpus reconstruction")
    require(input_payload["obstruction_core"]["nonedges_used_as_distance_constraints"] is False,
            "input JSON does not preserve required-edge-only semantics")
    selected, targets, selection_summary = rebuild_targets()
    require(input_payload["current_K7_targets"]["population"] == TARGET_COUNT,
            "input JSON target population mismatch")
    require(input_payload["current_K7_targets"]["indices_sha256"] == TARGET_INDICES_SHA256,
            "input JSON target index hash mismatch")

    report = json.loads(report_path.read_text(encoding="utf-8"))
    require(report.get("schema") == 1 and
            report.get("kind") == "d6_n14_pattern_954_containment_report",
            "report schema/kind mismatch")
    require(report.get("status") == "COMPLETE", "report is not COMPLETE")
    require(report.get("required_edge_only") is True, "report semantics changed")
    require(report.get("processed") == TARGET_COUNT, "report processed count mismatch")
    require(report.get("config") == EXPECTED_CONFIG, "report configuration changed")
    require(stable_hash(report["config"]) == report.get("config_sha256") == CONFIG_SHA256,
            "report configuration hash mismatch")
    require(report.get("counts") == {
        "HIT": HIT_COUNT, "NO_HIT": RESIDUE_COUNT,
        "TIMEOUT": 0, "INFRA_ERROR": 0,
    }, "report status accounting mismatch")
    require(report.get("unresolved_indices") == [], "report has unresolved indices")
    require(stable_hash([]) == report.get("unresolved_indices_sha256"),
            "empty unresolved-list hash mismatch")
    datetime.fromisoformat(report["finished_utc"])

    checkpoint_dir = Path(report["checkpoint_directory"])
    require(checkpoint_dir.name == CONFIG_SHA256, "checkpoint directory/config mismatch")
    campaign_path = checkpoint_dir / "campaign.json"
    require(sha256(campaign_path) == report.get("campaign_sha256") == CAMPAIGN_SHA256,
            "campaign artifact hash mismatch")
    campaign = json.loads(campaign_path.read_text(encoding="utf-8"))
    require(campaign.get("schema") == 1 and
            campaign.get("kind") == "d6_n14_pattern_954_containment_campaign",
            "campaign schema/kind mismatch")
    require(campaign.get("config") == EXPECTED_CONFIG and
            campaign.get("config_sha256") == CONFIG_SHA256,
            "campaign configuration mismatch")
    require(campaign.get("command") == report.get("command"),
            "campaign/report commands differ")
    require(campaign.get("controls") == report.get("controls"),
            "campaign/report controls differ")
    require(campaign.get("workers") == report.get("workers") == 11,
            "campaign worker count mismatch")
    require(campaign["machine"] == {
        "logical_cpus": 12, "machine": "arm64",
        "platform": "macOS-14.5-arm64-arm-64bit", "python": "3.11.15",
    }, "campaign machine identity mismatch")
    require(datetime.fromisoformat(campaign["started_utc"])
            <= datetime.fromisoformat(report["finished_utc"]),
            "campaign timestamps are reversed")
    compilation_summary = verify_compilation(report)

    manifest, results, edge_checks = audit_chunks(
        checkpoint_dir, selected, targets, pattern
    )
    statuses = Counter(result["status"] for result in results)
    require(dict(statuses) == {"NO_HIT": RESIDUE_COUNT, "HIT": HIT_COUNT},
            "chunk status accounting mismatch")
    require(sum(result["nodes"] for result in results) == report["total_kernel_nodes"],
            "chunk/report node totals differ")
    require(math.isclose(
        sum(result["wall_seconds"] for result in results),
        report["sum_target_wall_seconds"], rel_tol=0, abs_tol=1e-12,
    ), "chunk/report target wall totals differ")
    hit_records = [result for result in results if result["status"] == "HIT"]
    hit_indices = [result["index"] for result in hit_records]
    residue_indices = [result["index"] for result in results
                       if result["status"] == "NO_HIT"]
    require(hit_indices == report["certified_hit_indices"],
            "chunk/report HIT index order differs")
    require(hit_records == report["certified_hits"],
            "chunk/report HIT records differ")
    require(stable_hash(hit_indices) == report["certified_hit_indices_sha256"]
            == HIT_INDICES_SHA256, "HIT index hash mismatch")
    require(stable_hash(residue_indices) == RESIDUE_INDICES_SHA256,
            "NO_HIT/residue index hash mismatch")
    require([index for index in selected if index not in set(hit_indices)] == residue_indices,
            "residue is not the ordered complement of exact HITs")
    require(edge_checks == HIT_COUNT * 54,
            "not every edge of every HIT mapping was checked")

    allowed_files = {"campaign.json", *(row["name"] for row in manifest)}
    actual_files = {path.name for path in checkpoint_dir.iterdir() if path.is_file()}
    require(actual_files == allowed_files, "checkpoint directory has missing/extra files")
    source_hashes = dict(EXPECTED_BUNDLE_HASHES)
    source_hashes.update({
        "aeq_d6_n14.txt": N14_CORPUS_SHA256,
        "d6_k7_positive_dual_selection.json": PARENT_SELECTION_SHA256,
        "d6_k7_positive_dual_full_report.json": DUAL_REPORT_SHA256,
        "d6_k7_positive_dual_full_decisions.tsv.gz": DUAL_DECISIONS_SHA256,
        ".runs/d6_k7_rank_survivors.json": RANK_INPUT_SHA256,
        ".runs/d6_n14_pattern_954_targets_12839.tsv": TARGET_TSV_SHA256,
        ".runs/d6_n14_pattern_954_containment.dylib": BINARY_SHA256,
        REPORT_NAME: REPORT_SHA256,
        str(campaign_path.relative_to(ROOT)): CAMPAIGN_SHA256,
    })
    return {
        "schema": 1,
        "kind": "d6_n14_pattern_954_containment_independent_verification",
        "status": "PASS",
        "claim": (
            "All 3,403 HIT mappings are injective non-induced embeddings of "
            "the 54-edge obstruction. Their ordered complement is the exact "
            "9,436-graph current residue. NO_HIT search completeness was "
            "audited from artifacts but not independently recomputed."
        ),
        "required_edge_only": True,
        "source_hashes": source_hashes,
        "pattern": {"order": 13, "edges": 54, "adjacency": list(pattern)},
        "selection": selection_summary,
        "report": {
            "sha256": REPORT_SHA256,
            "config_sha256": CONFIG_SHA256,
            "campaign_sha256": CAMPAIGN_SHA256,
            "processed": TARGET_COUNT,
            "counts": dict(report["counts"]),
            "total_kernel_nodes": report["total_kernel_nodes"],
        },
        "mapping_verification": {
            "mappings_checked": len(hit_indices),
            "required_edges_checked_per_mapping": 54,
            "required_edge_incidence_checks": edge_checks,
            "all_injective": True,
        },
        "certified_hit_indices": hit_indices,
        "certified_hit_indices_sha256": stable_hash(hit_indices),
        "residue_indices": residue_indices,
        "residue_indices_sha256": stable_hash(residue_indices),
        "chunk_artifacts": {
            "count": len(manifest),
            "bytes": sum(row["bytes"] for row in manifest),
            "manifest_sha256": stable_hash(manifest),
            "manifest": manifest,
            "ordered_coverage": True,
        },
        "compilation": compilation_summary,
        "trust_scope": {
            "exactly_reverified": "all positive HIT mappings and residue complement",
            "artifact_audited_only": "C-kernel exhaustiveness of 9,436 NO_HIT searches",
            "floating_point_used_for_mathematical_claim": False,
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, default=ROOT / REPORT_NAME)
    parser.add_argument(
        "--output", type=Path,
        default=ROOT / "d6_n14_pattern_954_containment_verification.json",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    result = verify(args.report)
    result["verifier"] = {
        "file": Path(__file__).name,
        "sha256": sha256(Path(__file__)),
        "command": "python3 d6_n14_pattern_954_containment_report_verify.py",
    }
    atomic_json(args.output, result)
    print(json.dumps({
        "status": result["status"],
        "mappings_checked": result["mapping_verification"]["mappings_checked"],
        "edge_checks": result["mapping_verification"]["required_edge_incidence_checks"],
        "residue": len(result["residue_indices"]),
        "residue_indices_sha256": result["residue_indices_sha256"],
        "chunk_manifest_sha256": result["chunk_artifacts"]["manifest_sha256"],
        "output": str(args.output),
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
