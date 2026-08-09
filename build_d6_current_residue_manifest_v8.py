#!/usr/bin/env python3
"""Build the source-bound dimension-six current residue manifest v8.

The v8 boundary is a deliberately small exact delta.  It consumes the frozen,
independently verified v7 manifest and the independently verified K6
saturated-singleton-basis result, removes exactly its two rejected K6-only
indices, and carries the K7 class through unchanged.  No production search
kernel is imported.
"""

from __future__ import annotations

import argparse
import copy
import gzip
import hashlib
import json
import os
import platform
import subprocess
import sys
import time
from datetime import UTC, datetime
from pathlib import Path
from typing import Sequence


ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / "d6_current_residue_manifest_v8.json"
PACKAGE_FILES = (
    "build_d6_current_residue_manifest_v8.py",
    "verify_d6_current_residue_manifest_v8.py",
    "test_d6_current_residue_manifest_v8.py",
    "d6_current_residue_manifest_v8.md",
)
PRIMARY_ARTIFACTS = {
    "d6_current_residue_manifest_v7.json": (
        "1ab0948780d73cdd6aca2925107fc240dcae5df51d7cf5393accd8de3ae79f95"
    ),
    "d6_current_residue_manifest_v7_verification.json": (
        "ff02a9ca5bfe769ce77bbf67dedbc8e928561fb4f291a78d2ea577e06e784255"
    ),
    "d6_k6_saturated_singleton_basis_report.json": (
        "5f25193311ad769a402ec5245a6a690619cc41225e63a56546fe50a30d7a93e5"
    ),
    "d6_k6_saturated_singleton_basis_certificates.json.gz": (
        "0d3706ccd58cced274cdcd1cdf427ba72c529980caf641c6c37bfaa52ee139c1"
    ),
    "d6_k6_saturated_singleton_basis_verification.json": (
        "9524b2f92baf115406e26d02f62598d6726e09d7fb8622ee750aab2ef7b36eb4"
    ),
}
EXPECTED_REJECTED = [3_138_618, 3_673_988]
EXPECTED_REJECTED_SHA256 = (
    "45ce66ccf4db5fe1a565a576cd16fc719e51b20c6ea9fee858dc3682463e9e84"
)
EXPECTED_V7_K7_SHA256 = (
    "e545ea5a8fc9691ed98e885919e149f84dc7e8b4796d0fe9e6cb3a519df6fd8c"
)
EXPECTED_V7_K6_SHA256 = (
    "a23da7bd9ae21eb7c078a08a7bce92f22d397d5c71f1e4f51d9cf8e59bbfa907"
)
EXPECTED_K6_SHA256 = (
    "f274aab3f2d408fce34c7cd7eb4518621fa9d634a8943e3ee1a38d76df9f3d31"
)
EXPECTED_K6_GRAPHS_SHA256 = (
    "92894350a8a3f76e791cd419049e34a52903c0ff6a4a45ed7a262ab801ca630e"
)
EXPECTED_COMBINED_SHA256 = (
    "2391a93a3629517363106603bdad00be9b6f960966d9089ef37d93be2988c213"
)
EXPECTED_SORTED_SHA256 = (
    "2ac404f7260a130549f572ba5c34812b95b9df96bab7ab546b2f3f4c6656ce65"
)
EXPECTED_TAGGED_INDICES_SHA256 = (
    "eb39e7f43e1bc4bdfbada948af846c6ba33a087e0f1df9fc139f1226699c37bc"
)
EXPECTED_TAGGED_GRAPHS_SHA256 = (
    "f54699e5cb758864cda879c347ba78b90f4b14390e24d8328858d4b13927d278"
)
EXPECTED_POST_V6_K6_REJECTED_SHA256 = (
    "621e3d8f9859cbab831f307ef20d57a0a1ffeba988b5c3178a0352f26ea7fa0b"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


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
    return hashlib.sha256(raw).hexdigest()


def load(name: str) -> dict:
    value = json.loads((ROOT / name).read_text(encoding="utf-8"))
    require(isinstance(value, dict), f"{name} is not an object")
    return value


def verify_named_hashes(mapping: object, label: str) -> dict[str, str]:
    require(isinstance(mapping, dict) and mapping, f"empty {label}")
    normalized: dict[str, str] = {}
    for name, expected in mapping.items():
        require(
            isinstance(name, str)
            and Path(name).name == name
            and isinstance(expected, str)
            and len(expected) == 64,
            f"malformed {label}",
        )
        require(sha256(ROOT / name) == expected, f"{label} hash: {name}")
        normalized[name] = expected
    return dict(sorted(normalized.items()))


def git_blob_sha256(commit: str, name: str) -> str:
    blob = subprocess.run(
        ["git", "show", f"{commit}:{name}"],
        cwd=ROOT,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    ).stdout
    return hashlib.sha256(blob).hexdigest()


def committed_package_boundary() -> dict:
    """Bind all four v8 sources to the launch commit."""

    def git(*arguments: str) -> str:
        return subprocess.run(
            ["git", *arguments],
            cwd=ROOT,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        ).stdout.strip()

    commit = git("rev-parse", "HEAD")
    branch = git("branch", "--show-current")
    require(branch == "codex/dimension6", "official v8 launch branch")
    sources = {name: sha256(ROOT / name) for name in PACKAGE_FILES}
    for name, expected in sources.items():
        require(
            git_blob_sha256(commit, name) == expected,
            f"uncommitted v8 package source: {name}",
        )
    porcelain = git("status", "--porcelain=v1", "--untracked-files=all")
    return {
        "commit": commit,
        "branch": branch,
        "package_sources": dict(sorted(sources.items())),
        "porcelain_lines": porcelain.splitlines(),
        "porcelain_sha256": hashlib.sha256(
            porcelain.encode("utf-8")
        ).hexdigest(),
        "package_sources_equal_committed_blobs": True,
    }


def validate_rows(value: object, label: str) -> tuple[int, ...]:
    require(
        isinstance(value, list)
        and len(value) == 19
        and all(isinstance(row, int) and row >= 0 for row in value),
        f"malformed adjacency: {label}",
    )
    rows = tuple(value)
    full = (1 << len(rows)) - 1
    for vertex, row in enumerate(rows):
        require(not row & ~full and not row & (1 << vertex), f"loop/range: {label}")
        for other in range(vertex):
            require(
                bool(row & (1 << other))
                == bool(rows[other] & (1 << vertex)),
                f"asymmetry: {label}",
            )
    return rows


def alpha_at_most_two(rows: Sequence[int]) -> bool:
    full = (1 << len(rows)) - 1
    for right in range(len(rows)):
        non_right = full & ~(rows[right] | (1 << right))
        for left in range(right):
            if rows[right] & (1 << left):
                continue
            common = non_right & ~(rows[left] | (1 << left))
            if common & ~((1 << right) - 1):
                return False
    return True


def contains_clique(rows: Sequence[int], size: int) -> bool:
    def search(candidates: int, needed: int) -> bool:
        if needed == 0:
            return True
        if candidates.bit_count() < needed:
            return False
        while candidates:
            bit = candidates & -candidates
            candidates ^= bit
            vertex = bit.bit_length() - 1
            if search(candidates & rows[vertex], needed - 1):
                return True
            if candidates.bit_count() < needed:
                break
        return False

    return search((1 << len(rows)) - 1, size)


def validate_embedded_class(records: object, name: str) -> list[int]:
    require(isinstance(records, list), f"{name} graphs are not a list")
    indices = []
    for record in records:
        require(isinstance(record, dict), f"{name} graph is not an object")
        index = record.get("index")
        require(isinstance(index, int), f"{name} graph index")
        rows = validate_rows(record.get("adjacency"), f"{name}/{index}")
        require(alpha_at_most_two(rows), f"{name}/{index}: alpha exceeds two")
        require(contains_clique(rows, 6), f"{name}/{index}: no K6")
        require(
            contains_clique(rows, 7) == (name == "K7"),
            f"{name}/{index}: wrong clique class",
        )
        indices.append(index)
    require(len(indices) == len(set(indices)), f"duplicate {name} index")
    return indices


def validate_v7() -> tuple[dict, dict, list[int], list[int], dict]:
    v7 = load("d6_current_residue_manifest_v7.json")
    check = load("d6_current_residue_manifest_v7_verification.json")
    require(
        v7.get("schema") == "d6-current-certified-residue-v7"
        and v7.get("status") == "COMPLETE_MIXED_CERTIFICATE_UNION"
        and v7.get("class_order") == ["K7", "K6_only"]
        and check.get("schema")
        == "d6-current-certified-residue-v7-verification-v1"
        and check.get("status") == "PASS"
        and check.get("manifest", {}).get("sha256")
        == PRIMARY_ARTIFACTS["d6_current_residue_manifest_v7.json"]
        and all(check.get("checks", {}).values()),
        "frozen v7 manifest/verification boundary",
    )
    k7_graphs = v7.get("classes", {}).get("K7", {}).get("graphs")
    k6_graphs = v7.get("classes", {}).get("K6_only", {}).get("graphs")
    k7 = validate_embedded_class(k7_graphs, "K7")
    k6 = validate_embedded_class(k6_graphs, "K6_only")
    require(
        len(k7) == 12
        and stable_hash(k7) == EXPECTED_V7_K7_SHA256
        and v7["classes"]["K7"].get("indices") == k7
        and v7["classes"]["K7"].get("indices_sha256")
        == EXPECTED_V7_K7_SHA256
        and len(k6) == 251
        and stable_hash(k6) == EXPECTED_V7_K6_SHA256
        and v7["classes"]["K6_only"].get("indices") == k6
        and v7["classes"]["K6_only"].get("indices_sha256")
        == EXPECTED_V7_K6_SHA256
        and not set(k7) & set(k6)
        and check.get("ordered_residue_indices_sha256")
        == stable_hash(k7 + k6)
        and check.get("counts")
        == {
            "K7": 12,
            "K6_only": 251,
            "combined": 263,
            "exact_rejections_from_v5": 643,
            "interval_incremental_rejections": 5,
        },
        "frozen v7 class/accounting boundary",
    )
    accounting = v7.get("certificate_accounting", {})
    require(
        accounting.get("v5_graphs") == 911
        and accounting.get("exact", {}).get("total_rejections_from_v5") == 643
        and accounting.get("interval", {}).get(
            "total_incremental_rejections_from_v5"
        )
        == 5
        and accounting.get("union_rejections_from_v5") == 648
        and accounting.get("current_residue") == 263,
        "v7 certificate accounting",
    )
    package = v7.get("source_boundary", {}).get("package", {})
    package_sources = package.get("package_sources", {})
    require(
        package.get("package_sources_equal_committed_blobs") is True
        and package.get("branch") == "codex/dimension6"
        and set(package_sources) == {
            "build_d6_current_residue_manifest_v7.py",
            "verify_d6_current_residue_manifest_v7.py",
            "test_d6_current_residue_manifest_v7.py",
            "d6_current_residue_manifest_v7.md",
        }
        and check.get("execution", {}).get("verifier_sha256")
        == package_sources["verify_d6_current_residue_manifest_v7.py"],
        "v7 source package boundary",
    )
    for name, expected in package_sources.items():
        require(sha256(ROOT / name) == expected, f"v7 working source: {name}")
        require(
            git_blob_sha256(package["commit"], name) == expected,
            f"v7 committed source: {name}",
        )
    source_root = {
        "manifest_source_boundary": v7["source_boundary"],
        "builder": {
            "build_d6_current_residue_manifest_v7.py": package_sources[
                "build_d6_current_residue_manifest_v7.py"
            ]
        },
        "verifier": {
            "verify_d6_current_residue_manifest_v7.py": check["execution"][
                "verifier_sha256"
            ]
        },
    }
    return v7, check, k7, k6, source_root


def validate_singleton_layer(v7_k6: list[int]) -> tuple[list[int], list[int], dict]:
    report = load("d6_k6_saturated_singleton_basis_report.json")
    verification = load("d6_k6_saturated_singleton_basis_verification.json")
    rejected = list(map(int, report.get("rejected_indices", ())))
    residue = list(map(int, report.get("ordered_residue_indices", ())))
    require(
        report.get("schema") == "d6-k6-saturated-singleton-basis-v1"
        and report.get("status") == "COMPLETE"
        and report.get("ordered_input_indices") == v7_k6
        and report.get("ordered_input_indices_sha256")
        == EXPECTED_V7_K6_SHA256
        and report.get("input_graphs") == 251
        and report.get("graphs_rejected") == 2
        and report.get("graphs_surviving") == 249
        and rejected == EXPECTED_REJECTED
        and report.get("rejected_indices_sha256")
        == EXPECTED_REJECTED_SHA256
        and residue
        == [index for index in v7_k6 if index not in set(EXPECTED_REJECTED)]
        and report.get("ordered_residue_indices_sha256") == EXPECTED_K6_SHA256,
        "saturated-singleton exact partition",
    )
    graph_results = report.get("graph_results", ())
    require(
        [int(row["index"]) for row in graph_results] == v7_k6
        and [int(row["index"]) for row in graph_results if row.get("rejected")]
        == EXPECTED_REJECTED
        and report.get("totals", {}).get(
            "saturated_singleton_branches_checked"
        )
        == 291
        and report.get("totals", {}).get(
            "saturated_singleton_branches_failed"
        )
        == 291,
        "saturated-singleton graph coverage",
    )
    require(
        verification.get("schema")
        == "d6-k6-saturated-singleton-basis-verification-v1"
        and verification.get("status") == "PASS"
        and verification.get("report", {}).get("sha256")
        == PRIMARY_ARTIFACTS["d6_k6_saturated_singleton_basis_report.json"]
        and verification.get("certificates", {}).get("sha256")
        == PRIMARY_ARTIFACTS[
            "d6_k6_saturated_singleton_basis_certificates.json.gz"
        ]
        and verification.get("graphs_recomputed") == 251
        and verification.get("graphs_rejected") == 2
        and verification.get("graphs_surviving") == 249
        and verification.get("rejected_indices") == EXPECTED_REJECTED
        and verification.get("rejected_indices_sha256")
        == EXPECTED_REJECTED_SHA256
        and verification.get("ordered_residue_indices_sha256")
        == EXPECTED_K6_SHA256
        and all(verification.get("checks", {}).values()),
        "saturated-singleton independent verification",
    )
    require(
        report.get("semantics", {}).get("candidate_nonedges")
        == "unconstrained_and_may_be_unit"
        and report.get("semantics", {}).get("allowed_supports")
        == "upper_bounds_actual_subsets_exhausted"
        and report.get("semantics", {}).get("actual_support_quantifier")
        == "basis_rule_prunes_only_one_complete_support_leaf"
        and report.get("positive_18_control")
        == {"passed": True, "K6_seeds": 32}
        and verification.get("positive_18_control")
        == {"passed": True, "K6_seeds": 32},
        "saturated-singleton semantics/control",
    )
    archive_entry = report.get("certificate_archive", {})
    require(
        archive_entry.get("path")
        == "d6_k6_saturated_singleton_basis_certificates.json.gz"
        and archive_entry.get("sha256")
        == PRIMARY_ARTIFACTS[
            "d6_k6_saturated_singleton_basis_certificates.json.gz"
        ]
        and archive_entry.get("rejected_graphs") == 2,
        "saturated-singleton certificate root",
    )
    with gzip.open(ROOT / archive_entry["path"], "rt", encoding="utf-8") as stream:
        archive = json.load(stream)
    require(
        archive.get("schema")
        == "d6-k6-saturated-singleton-basis-certificates-v1"
        and archive.get("ordered_input_indices_sha256")
        == EXPECTED_V7_K6_SHA256
        and [int(row["index"]) for row in archive.get("rejected_graphs", ())]
        == EXPECTED_REJECTED,
        "saturated-singleton certificate coverage",
    )

    dependencies = verify_named_hashes(
        report.get("dependencies"), "saturated-singleton dependencies"
    )
    production_sources = verify_named_hashes(
        report.get("source_sha256"), "saturated-singleton production sources"
    )
    verifier_sources = verify_named_hashes(
        verification.get("sources"), "saturated-singleton verifier sources"
    )
    git = report.get("execution", {}).get("git", {})
    commit = git.get("commit")
    require(
        isinstance(commit, str)
        and len(commit) == 40
        and git.get("branch") == "codex/dimension6"
        and git.get("tracked_clean") is True
        and git.get("committed_source_sha256") == production_sources
        and verification.get("provenance", {}).get("commit") == commit
        and verification.get("provenance", {}).get(
            "tracked_sources_clean_at_launch"
        )
        is True,
        "saturated-singleton source launch",
    )
    for name, expected in production_sources.items():
        require(
            git_blob_sha256(commit, name) == expected,
            f"saturated-singleton committed source: {name}",
        )
    require(
        archive.get("source_sha256") == production_sources
        and archive.get("production_source_sha256")
        == report.get("production_source_sha256"),
        "saturated-singleton archive source boundary",
    )
    return rejected, residue, {
        "launch_commit": commit,
        "dependencies": dependencies,
        "production_sources": production_sources,
        "verifier_sources": verifier_sources,
        "certificate_archive": copy.deepcopy(archive_entry),
    }


def build_manifest(package_boundary: dict | None = None) -> dict:
    verify_named_hashes(PRIMARY_ARTIFACTS, "primary artifact")
    v7, v7_check, k7_indices, v7_k6, inherited_sources = validate_v7()
    rejected, k6_indices, singleton_sources = validate_singleton_layer(v7_k6)
    require(
        rejected == EXPECTED_REJECTED
        and all(index in set(v7_k6) for index in rejected)
        and not set(rejected) & set(k7_indices),
        "v8 rejection class/disjointness",
    )

    k7_class = copy.deepcopy(v7["classes"]["K7"])
    require(
        k7_class.get("count") == 12
        and k7_class.get("indices") == k7_indices
        and k7_class.get("indices_sha256") == EXPECTED_V7_K7_SHA256,
        "unchanged K7 class",
    )
    v7_k6_class = v7["classes"]["K6_only"]
    k6_graphs = [
        copy.deepcopy(record)
        for record in v7_k6_class["graphs"]
        if int(record["index"]) not in set(rejected)
    ]
    require(
        [int(record["index"]) for record in k6_graphs] == k6_indices
        and stable_hash(k6_indices) == EXPECTED_K6_SHA256
        and stable_hash(k6_graphs) == EXPECTED_K6_GRAPHS_SHA256,
        "v8 K6 embedded residue",
    )
    for record in k6_graphs:
        rows = validate_rows(record.get("adjacency"), f"K6_only/{record.get('index')}")
        require(
            alpha_at_most_two(rows)
            and contains_clique(rows, 6)
            and not contains_clique(rows, 7),
            f"K6_only/{record.get('index')}: graph invariant",
        )

    new_layer = {
        "name": "saturated_singleton_basis",
        "input_count": 251,
        "input_indices_sha256": EXPECTED_V7_K6_SHA256,
        "rejection_count": 2,
        "rejected_indices": rejected,
        "rejected_indices_sha256": EXPECTED_REJECTED_SHA256,
        "residue_count": 249,
        "residue_indices_sha256": EXPECTED_K6_SHA256,
    }
    post_v6_k6_rejected = [
        *map(int, v7_k6_class["post_v6_exact_rejected_indices"]),
        *rejected,
    ]
    require(
        len(post_v6_k6_rejected) == 376
        and len(set(post_v6_k6_rejected)) == 376
        and stable_hash(post_v6_k6_rejected)
        == EXPECTED_POST_V6_K6_REJECTED_SHA256,
        "cumulative K6 rejection accounting",
    )
    k6_class = copy.deepcopy(v7_k6_class)
    k6_class.update({
        "v7_count": 251,
        "exact_layers": [*copy.deepcopy(v7_k6_class["exact_layers"]), new_layer],
        "post_v7_exact_rejection_count": 2,
        "post_v7_exact_rejected_indices": rejected,
        "post_v7_exact_rejected_indices_sha256": EXPECTED_REJECTED_SHA256,
        "post_v6_exact_rejection_count": 376,
        "post_v6_exact_rejected_indices": post_v6_k6_rejected,
        "post_v6_exact_rejected_indices_sha256": (
            EXPECTED_POST_V6_K6_REJECTED_SHA256
        ),
        "count": 249,
        "indices": k6_indices,
        "indices_sha256": EXPECTED_K6_SHA256,
        "graphs": k6_graphs,
        "graphs_sha256": EXPECTED_K6_GRAPHS_SHA256,
    })

    combined_indices = k7_indices + k6_indices
    require(
        len(combined_indices) == 261
        and stable_hash(combined_indices) == EXPECTED_COMBINED_SHA256
        and stable_hash(sorted(combined_indices)) == EXPECTED_SORTED_SHA256
        and not set(k7_indices) & set(k6_indices),
        "combined v8 residue",
    )
    tagged_indices = [
        {"class": class_name, "index": int(record["index"])}
        for class_name, records in (("K7", k7_class["graphs"]), ("K6_only", k6_graphs))
        for record in records
    ]
    tagged_graphs = [
        {"class": class_name, **record}
        for class_name, records in (("K7", k7_class["graphs"]), ("K6_only", k6_graphs))
        for record in records
    ]
    require(
        stable_hash(tagged_indices) == EXPECTED_TAGGED_INDICES_SHA256
        and stable_hash(tagged_graphs) == EXPECTED_TAGGED_GRAPHS_SHA256,
        "tagged v8 roots",
    )

    if package_boundary is None:
        package_boundary = committed_package_boundary()
    require(
        set(package_boundary.get("package_sources", {})) == set(PACKAGE_FILES),
        "v8 package source set",
    )
    accounting = copy.deepcopy(v7["certificate_accounting"])
    accounting["exact"]["post_v6_rejections"] = 381
    accounting["exact"]["total_rejections_from_v5"] = 645
    accounting["exact"]["layers"]["K6_saturated_singleton_basis"] = 2
    accounting["union_rejections_from_v5"] = 650
    accounting["current_residue"] = 261
    controls = copy.deepcopy(v7["controls"])
    singleton_report = load("d6_k6_saturated_singleton_basis_report.json")
    controls["known_realizable_18"]["exact_K6_saturated_singleton_basis"] = (
        copy.deepcopy(singleton_report["positive_18_control"])
    )
    controls["exact_K6_saturated_singleton_synthetic"] = copy.deepcopy(
        singleton_report["synthetic_controls"]
    )
    source_roots = {
        "inherited_v7": inherited_sources,
        "exact_K6_saturated_singleton_basis": singleton_sources,
    }
    return {
        "schema": "d6-current-certified-residue-v8",
        "status": "COMPLETE_MIXED_CERTIFICATE_UNION",
        "description": (
            "The frozen independently verified v7 residue after the exact "
            "K6 saturated-singleton-basis increment. Survival is not "
            "realizability."
        ),
        "class_order": ["K7", "K6_only"],
        "classes": {"K7": k7_class, "K6_only": k6_class},
        "combined": {
            "v6_count": 644,
            "v7_count": 263,
            "post_v7_exact_rejection_count": 2,
            "post_v6_exact_rejection_count": 381,
            "post_v6_interval_incremental_rejection_count": 2,
            "post_v6_union_rejection_count": 383,
            "count": 261,
            "class_counts": {"K7": 12, "K6_only": 249},
            "ordered_indices_sha256": EXPECTED_COMBINED_SHA256,
            "sorted_indices_sha256": EXPECTED_SORTED_SHA256,
            "tagged_indices_sha256": EXPECTED_TAGGED_INDICES_SHA256,
            "tagged_graphs_sha256": EXPECTED_TAGGED_GRAPHS_SHA256,
            "cross_class_overlap": 0,
        },
        "certificate_accounting": accounting,
        "semantics": copy.deepcopy(v7["semantics"]),
        "controls": controls,
        "nonclaims": copy.deepcopy(v7["nonclaims"]),
        "execution": {
            "command": " ".join(map(str, sys.argv)),
            "generated_utc": datetime.now(UTC).isoformat(),
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python_version": sys.version,
            "logical_cpus": os.cpu_count(),
        },
        "source_boundary": {
            "primary_artifacts": dict(sorted(PRIMARY_ARTIFACTS.items())),
            "layer_source_roots": source_roots,
            "package": package_boundary,
        },
    }


def atomic_json(path: Path, value: object) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    manifest = build_manifest()
    atomic_json(args.output, manifest)
    print(json.dumps({
        "status": manifest["status"],
        "manifest": str(args.output),
        "sha256": sha256(args.output),
        "K7": manifest["classes"]["K7"]["count"],
        "K6_only": manifest["classes"]["K6_only"]["count"],
        "combined": manifest["combined"]["count"],
        "source_commit": manifest["source_boundary"]["package"]["commit"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
