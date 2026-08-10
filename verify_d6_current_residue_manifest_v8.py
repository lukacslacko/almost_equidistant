#!/usr/bin/env python3
"""Import-independent structural checker for the d=6 residue manifest v8."""

from __future__ import annotations

import argparse
import ast
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
PACKAGE_FILES = {
    "build_d6_current_residue_manifest_v8.py",
    "verify_d6_current_residue_manifest_v8.py",
    "test_d6_current_residue_manifest_v8.py",
    "d6_current_residue_manifest_v8.md",
}
PRIMARY = {
    "d6_current_residue_manifest_v7.json": (
        "1ab0948780d73cdd6aca2925107fc240dcae5df51d7cf5393accd8de3ae79f95"
    ),
    "d6_current_residue_manifest_v7_verification.json": (
        "ff02a9ca5bfe769ce77bbf67dedbc8e928561fb4f291a78d2ea577e06e784255"
    ),
    "d6_k6_saturated_singleton_basis_report.json": (
        "5f25193311ad769a402ec5245a6a690619cc41225e63a56546fe50a30d7a93e5"
    ),
    "d6_k6_saturated_singleton_basis_verification.json": (
        "9524b2f92baf115406e26d02f62598d6726e09d7fb8622ee750aab2ef7b36eb4"
    ),
    "d6_k6_saturated_singleton_basis_certificates.json.gz": (
        "0d3706ccd58cced274cdcd1cdf427ba72c529980caf641c6c37bfaa52ee139c1"
    ),
}
EXPECTED_REJECTED = [3138618, 3673988]
EXPECTED_REJECTED_SHA256 = (
    "45ce66ccf4db5fe1a565a576cd16fc719e51b20c6ea9fee858dc3682463e9e84"
)
EXPECTED_K7_SHA256 = (
    "e545ea5a8fc9691ed98e885919e149f84dc7e8b4796d0fe9e6cb3a519df6fd8c"
)
EXPECTED_V7_K6_SHA256 = (
    "a23da7bd9ae21eb7c078a08a7bce92f22d397d5c71f1e4f51d9cf8e59bbfa907"
)
EXPECTED_K6_SHA256 = (
    "f274aab3f2d408fce34c7cd7eb4518621fa9d634a8943e3ee1a38d76df9f3d31"
)
EXPECTED_K7_GRAPHS_SHA256 = (
    "f0ad41194d3e16e58d3210ff9f60507bd5e4332671f25c02e7bd4924b1290680"
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
EXPECTED_CUMULATIVE_K6_REJECTIONS_SHA256 = (
    "621e3d8f9859cbab831f307ef20d57a0a1ffeba988b5c3178a0352f26ea7fa0b"
)
EXPECTED_UNCOMPRESSED_CERTIFICATES_SHA256 = (
    "36d3bf52be35f5b3eaf08106c4895b7fef69176b082285d445405f4382dbc49c"
)
INTERVAL_TRUST = {
    "basic_arithmetic": (
        "IEEE-754 binary64 basic operations and sqrt are correctly rounded; "
        "each interval endpoint is expanded with nextafter."
    ),
    "candidate_nonedges": "unconstrained and allowed to be unit distance",
    "transcendentals": (
        "macOS libm cos endpoint values are assumed within 8 ulps; the kernel "
        "pads both directions by 8 nextafter steps and includes all interior extrema."
    ),
}


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


def load(name_or_path: str | Path) -> dict:
    path = Path(name_or_path)
    if not path.is_absolute():
        path = ROOT / path
    value = json.loads(path.read_text(encoding="utf-8"))
    require(isinstance(value, dict), f"{path.name}: object")
    return value


def named_hashes(mapping: object, label: str) -> dict[str, str]:
    require(isinstance(mapping, dict) and mapping, f"{label}: map")
    answer = {}
    for name, expected in mapping.items():
        require(
            isinstance(name, str)
            and Path(name).name == name
            and isinstance(expected, str)
            and len(expected) == 64,
            f"{label}: malformed",
        )
        require(sha256(ROOT / name) == expected, f"{label}: {name}")
        answer[name] = expected
    return dict(sorted(answer.items()))


def git_blob_sha256(commit: str, name: str) -> str:
    blob = subprocess.run(
        ["git", "show", f"{commit}:{name}"],
        cwd=ROOT,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    ).stdout
    return hashlib.sha256(blob).hexdigest()


def validate_rows(value: object, label: str) -> tuple[int, ...]:
    require(
        isinstance(value, list)
        and len(value) == 19
        and all(type(row) is int for row in value),
        f"{label}: adjacency",
    )
    rows = tuple(value)
    for first, row in enumerate(rows):
        require(0 <= row < (1 << 19) and not row & (1 << first), f"{label}: row")
        for second in range(first):
            require(
                bool(row & (1 << second))
                == bool(rows[second] & (1 << first)),
                f"{label}: symmetry",
            )
    return rows


def alpha_at_most_two(rows: Sequence[int]) -> bool:
    for first in range(19):
        for second in range(first):
            if rows[first] & (1 << second):
                continue
            for third in range(second):
                if not rows[first] & (1 << third) and not rows[second] & (
                    1 << third
                ):
                    return False
    return True


def contains_clique(rows: Sequence[int], size: int) -> bool:
    def visit(candidates: int, need: int) -> bool:
        if not need:
            return True
        while candidates.bit_count() >= need:
            bit = candidates & -candidates
            candidates ^= bit
            if visit(candidates & rows[bit.bit_length() - 1], need - 1):
                return True
        return False

    return visit((1 << len(rows)) - 1, size)


def package_boundary(manifest: dict) -> dict:
    boundary = manifest.get("source_boundary", {}).get("package", {})
    sources = boundary.get("package_sources")
    require(isinstance(sources, dict) and set(sources) == PACKAGE_FILES, "package sources")
    require(
        boundary.get("branch") == "codex/dimension6"
        and boundary.get("package_sources_equal_committed_blobs") is True,
        "package launch",
    )
    commit = boundary.get("commit")
    require(isinstance(commit, str) and len(commit) == 40, "package commit")
    for name, expected in sources.items():
        require(sha256(ROOT / name) == expected, f"working package source: {name}")
        require(git_blob_sha256(commit, name) == expected, f"committed package source: {name}")
    porcelain = "\n".join(boundary.get("porcelain_lines", ()))
    require(
        hashlib.sha256(porcelain.encode("utf-8")).hexdigest()
        == boundary.get("porcelain_sha256"),
        "package porcelain",
    )
    return {"commit": commit, "branch": boundary["branch"], "source_files": 4}


def frozen_v7_boundary(v7: dict, verification: dict) -> dict:
    require(
        v7.get("schema") == "d6-current-certified-residue-v7"
        and v7.get("status") == "COMPLETE_MIXED_CERTIFICATE_UNION"
        and v7.get("class_order") == ["K7", "K6_only"]
        and verification.get("schema")
        == "d6-current-certified-residue-v7-verification-v1"
        and verification.get("status") == "PASS"
        and verification.get("manifest", {}).get("sha256")
        == PRIMARY["d6_current_residue_manifest_v7.json"]
        and all(verification.get("checks", {}).values()),
        "frozen v7 report/check",
    )
    k7 = v7["classes"]["K7"]
    k6 = v7["classes"]["K6_only"]
    require(
        k7.get("count") == 12
        and k7.get("indices_sha256") == EXPECTED_K7_SHA256
        and k7.get("graphs_sha256") == EXPECTED_K7_GRAPHS_SHA256
        and stable_hash(k7.get("indices")) == EXPECTED_K7_SHA256
        and stable_hash(k7.get("graphs")) == EXPECTED_K7_GRAPHS_SHA256
        and k6.get("count") == 251
        and k6.get("indices_sha256") == EXPECTED_V7_K6_SHA256
        and stable_hash(k6.get("indices")) == EXPECTED_V7_K6_SHA256
        and v7.get("combined", {}).get("count") == 263
        and verification.get("counts")
        == {
            "K7": 12,
            "K6_only": 251,
            "combined": 263,
            "exact_rejections_from_v5": 643,
            "interval_incremental_rejections": 5,
        }
        and verification.get("class_indices_sha256")
        == {"K7": EXPECTED_K7_SHA256, "K6_only": EXPECTED_V7_K6_SHA256},
        "frozen v7 counts/hashes",
    )
    require(
        v7["certificate_accounting"]["exact"]["total_rejections_from_v5"]
        == 643
        and v7["certificate_accounting"]["interval"]
        ["total_incremental_rejections_from_v5"]
        == 5
        and v7["certificate_accounting"]["interval"]["trust_assumptions"]
        == INTERVAL_TRUST,
        "frozen v7 trust tiers",
    )
    package = v7["source_boundary"]["package"]
    package_sources = named_hashes(package["package_sources"], "v7 package sources")
    commit = package.get("commit")
    require(
        isinstance(commit, str)
        and len(commit) == 40
        and package.get("package_sources_equal_committed_blobs") is True
        and verification.get("source_boundary", {}).get("commit") == commit,
        "v7 package boundary",
    )
    for name, expected in package_sources.items():
        require(git_blob_sha256(commit, name) == expected, f"v7 package blob: {name}")
    verifier_hash = verification.get("execution", {}).get("verifier_sha256")
    require(
        verifier_hash == package_sources["verify_d6_current_residue_manifest_v7.py"],
        "v7 verifier root",
    )
    return {
        "manifest_source_boundary": v7["source_boundary"],
        "builder": {
            "build_d6_current_residue_manifest_v7.py": package_sources[
                "build_d6_current_residue_manifest_v7.py"
            ]
        },
        "verifier": {
            "verify_d6_current_residue_manifest_v7.py": verifier_hash
        },
    }


def singleton_boundary(report: dict, verification: dict) -> tuple[list[int], list[int], dict]:
    inputs = list(map(int, report.get("ordered_input_indices", ())))
    rejected = list(map(int, report.get("rejected_indices", ())))
    residue = list(map(int, report.get("ordered_residue_indices", ())))
    require(
        report.get("schema") == "d6-k6-saturated-singleton-basis-v1"
        and report.get("status") == "COMPLETE"
        and report.get("input_graphs") == 251
        and report.get("graphs_rejected") == 2
        and report.get("graphs_surviving") == 249
        and verification.get("schema")
        == "d6-k6-saturated-singleton-basis-verification-v1"
        and verification.get("status") == "PASS"
        and verification.get("report", {}).get("sha256")
        == PRIMARY["d6_k6_saturated_singleton_basis_report.json"]
        and verification.get("graphs_recomputed") == 251
        and verification.get("graphs_rejected") == 2
        and verification.get("graphs_surviving") == 249
        and all(verification.get("checks", {}).values()),
        "singleton report/check",
    )
    rejected_set = set(rejected)
    require(
        len(inputs) == len(set(inputs)) == 251
        and rejected == EXPECTED_REJECTED
        and residue == [index for index in inputs if index not in rejected_set]
        and len(residue) == len(set(residue)) == 249
        and not rejected_set & set(residue)
        and stable_hash(inputs) == EXPECTED_V7_K6_SHA256
        and report.get("ordered_input_indices_sha256") == EXPECTED_V7_K6_SHA256
        and stable_hash(rejected) == EXPECTED_REJECTED_SHA256
        and report.get("rejected_indices_sha256") == EXPECTED_REJECTED_SHA256
        and stable_hash(residue) == EXPECTED_K6_SHA256
        and report.get("ordered_residue_indices_sha256") == EXPECTED_K6_SHA256
        and verification.get("rejected_indices") == rejected
        and verification.get("rejected_indices_sha256") == EXPECTED_REJECTED_SHA256
        and verification.get("ordered_residue_indices_sha256")
        == EXPECTED_K6_SHA256,
        "singleton exact partition",
    )
    graph_results = report.get("graph_results")
    require(
        isinstance(graph_results, list)
        and [int(row["index"]) for row in graph_results] == inputs
        and [int(row["index"]) for row in graph_results if row.get("rejected")]
        == rejected,
        "singleton graph decisions",
    )
    totals = report.get("totals", {})
    require(
        totals.get("saturated_singleton_branches_checked") == 291
        and totals.get("saturated_singleton_branches_failed") == 291
        and report.get("positive_18_control") == {"K6_seeds": 32, "passed": True}
        and verification.get("positive_18_control")
        == {"K6_seeds": 32, "passed": True}
        and report.get("synthetic_controls", {}).get(
            "odd_opposite_sign_cycle_fails"
        )
        is True
        and all(verification.get("synthetic_controls", {}).values()),
        "singleton exact controls",
    )
    semantics = report.get("semantics", {})
    require(
        semantics.get("actual_support_quantifier")
        == "basis_rule_prunes_only_one_complete_support_leaf"
        and semantics.get("allowed_supports")
        == "upper_bounds_actual_subsets_exhausted"
        and semantics.get("candidate_nonedges")
        == "unconstrained_and_may_be_unit"
        and semantics.get("arithmetic")
        == "exact_integer_bitmasks_and_exhaustive_two_sign_CSP",
        "singleton semantics",
    )

    archive = report.get("certificate_archive", {})
    require(
        archive.get("path")
        == "d6_k6_saturated_singleton_basis_certificates.json.gz"
        and archive.get("sha256") == PRIMARY[archive["path"]]
        and archive.get("rejected_graphs") == 2
        and archive.get("uncompressed_bytes") == 71169
        and archive.get("uncompressed_sha256")
        == EXPECTED_UNCOMPRESSED_CERTIFICATES_SHA256
        and verification.get("certificates", {}).get("sha256")
        == PRIMARY[archive["path"]]
        and verification.get("certificates", {}).get("uncompressed_sha256")
        == EXPECTED_UNCOMPRESSED_CERTIFICATES_SHA256,
        "singleton certificate roots",
    )
    with gzip.open(ROOT / archive["path"], "rb") as stream:
        certificate_bytes = stream.read()
    require(
        len(certificate_bytes) == archive["uncompressed_bytes"]
        and hashlib.sha256(certificate_bytes).hexdigest()
        == EXPECTED_UNCOMPRESSED_CERTIFICATES_SHA256,
        "singleton decompressed archive root",
    )
    certificates = json.loads(certificate_bytes)
    require(
        certificates.get("schema")
        == "d6-k6-saturated-singleton-basis-certificates-v1"
        and certificates.get("ordered_input_indices_sha256")
        == EXPECTED_V7_K6_SHA256
        and certificates.get("production_source_sha256")
        == report.get("production_source_sha256")
        and certificates.get("source_sha256") == report.get("source_sha256")
        and [int(row["index"]) for row in certificates.get("rejected_graphs", ())]
        == rejected,
        "singleton certificate archive structure",
    )

    dependencies = named_hashes(report.get("dependencies"), "singleton dependencies")
    production_sources = named_hashes(
        report.get("source_sha256"), "singleton production sources"
    )
    verifier_sources = named_hashes(
        verification.get("sources"), "singleton verifier sources"
    )
    execution_git = report.get("execution", {}).get("git", {})
    commit = execution_git.get("commit")
    require(
        isinstance(commit, str)
        and len(commit) == 40
        and execution_git.get("branch") == "codex/dimension6"
        and execution_git.get("tracked_clean") is True
        and execution_git.get("committed_source_sha256") == production_sources
        and verification.get("provenance", {}).get("commit") == commit
        and verification.get("provenance", {}).get(
            "tracked_sources_clean_at_launch"
        )
        is True,
        "singleton source launch",
    )
    for name, expected in production_sources.items():
        require(git_blob_sha256(commit, name) == expected, f"singleton blob: {name}")
    require(
        report.get("production_source_sha256")
        == production_sources["d6_k6_saturated_singleton_basis.py"]
        and verifier_sources["d6_k6_saturated_singleton_basis.py"]
        == production_sources["d6_k6_saturated_singleton_basis.py"]
        and verifier_sources["verify_d6_k6_saturated_singleton_basis.py"]
        == production_sources["verify_d6_k6_saturated_singleton_basis.py"],
        "singleton producer/verifier roots",
    )
    source_root = {
        "launch_commit": commit,
        "dependencies": dependencies,
        "production_sources": production_sources,
        "verifier_sources": verifier_sources,
        "certificate_archive": archive,
    }
    return rejected, residue, source_root


def verify(
    manifest_path: Path,
    expected_manifest_sha256: str | None,
    *,
    enforce_source_boundary: bool = True,
) -> dict:
    manifest_path = manifest_path.resolve()
    manifest_hash = sha256(manifest_path)
    if expected_manifest_sha256 is not None:
        require(manifest_hash == expected_manifest_sha256, "manifest hash")
    named_hashes(PRIMARY, "primary")
    checks: dict[str, bool] = {}

    v7 = load("d6_current_residue_manifest_v7.json")
    v7_check = load("d6_current_residue_manifest_v7_verification.json")
    inherited_v7_sources = frozen_v7_boundary(v7, v7_check)
    checks["frozen_v7_boundary_sources"] = True

    report = load("d6_k6_saturated_singleton_basis_report.json")
    report_check = load("d6_k6_saturated_singleton_basis_verification.json")
    new_rejected, expected_k6_indices, singleton_sources = singleton_boundary(
        report, report_check
    )
    checks["singleton_partition_sources_controls"] = True

    v7_k7 = v7["classes"]["K7"]
    v7_k6 = v7["classes"]["K6_only"]
    v7_k7_indices = list(map(int, v7_k7["indices"]))
    v7_k6_indices = list(map(int, v7_k6["indices"]))
    require(
        report["ordered_input_indices"] == v7_k6_indices
        and not set(new_rejected) & set(v7_k7_indices)
        and set(new_rejected) <= set(v7_k6_indices)
        and not set(v7_k7_indices) & set(v7_k6_indices)
        and not set(v7_k7_indices) & set(expected_k6_indices),
        "v7/singleton class and layer disjointness",
    )
    checks["layer_input_and_class_disjointness"] = True

    manifest = load(manifest_path)
    require(
        manifest.get("schema") == "d6-current-certified-residue-v8"
        and manifest.get("status") == "COMPLETE_MIXED_CERTIFICATE_UNION"
        and manifest.get("class_order") == ["K7", "K6_only"],
        "v8 schema",
    )
    k7 = manifest["classes"]["K7"]
    k6 = manifest["classes"]["K6_only"]
    require(k7 == v7_k7, "K7 class changed from frozen v7")
    checks["K7_byte_for_value_unchanged"] = True

    expected_k6_graphs = [
        record
        for record in v7_k6["graphs"]
        if int(record["index"]) not in set(new_rejected)
    ]
    require(
        [int(record["index"]) for record in expected_k6_graphs]
        == expected_k6_indices
        and k6.get("graphs") == expected_k6_graphs,
        "K6 embedded graph filter",
    )
    for class_name, records in (
        ("K7", k7["graphs"]),
        ("K6_only", k6["graphs"]),
    ):
        for record in records:
            rows = validate_rows(
                record.get("adjacency"), f"{class_name}/{record.get('index')}"
            )
            require(
                alpha_at_most_two(rows) and contains_clique(rows, 6),
                "graph invariant",
            )
            require(
                contains_clique(rows, 7) == (class_name == "K7"),
                "graph class",
            )
    checks["embedded_graphs_and_classes"] = True

    new_layer = {
        "name": "saturated_singleton_basis",
        "input_count": 251,
        "input_indices_sha256": EXPECTED_V7_K6_SHA256,
        "rejection_count": 2,
        "rejected_indices": EXPECTED_REJECTED,
        "rejected_indices_sha256": EXPECTED_REJECTED_SHA256,
        "residue_count": 249,
        "residue_indices_sha256": EXPECTED_K6_SHA256,
    }
    cumulative_k6_rejected = list(v7_k6["post_v6_exact_rejected_indices"]) + new_rejected
    require(
        len(cumulative_k6_rejected) == len(set(cumulative_k6_rejected)) == 376
        and stable_hash(cumulative_k6_rejected)
        == EXPECTED_CUMULATIVE_K6_REJECTIONS_SHA256
        and k6.get("v6_count") == 625
        and k6.get("v7_count") == 251
        and k6.get("exact_layers") == v7_k6["exact_layers"] + [new_layer]
        and k6.get("post_v6_exact_rejection_count") == 376
        and k6.get("post_v6_exact_rejected_indices")
        == cumulative_k6_rejected
        and k6.get("post_v6_exact_rejected_indices_sha256")
        == EXPECTED_CUMULATIVE_K6_REJECTIONS_SHA256
        and k6.get("post_v7_exact_rejection_count") == 2
        and k6.get("post_v7_exact_rejected_indices") == EXPECTED_REJECTED
        and k6.get("post_v7_exact_rejected_indices_sha256")
        == EXPECTED_REJECTED_SHA256
        and k6.get("count") == 249
        and k6.get("indices") == expected_k6_indices
        and k6.get("indices_sha256") == EXPECTED_K6_SHA256
        and k6.get("graphs_sha256") == EXPECTED_K6_GRAPHS_SHA256
        and stable_hash(k6["graphs"]) == EXPECTED_K6_GRAPHS_SHA256,
        "K6 exact-layer accounting",
    )
    checks["K6_delta_partition_and_hashes"] = True

    combined_indices = v7_k7_indices + expected_k6_indices
    tagged_indices = [
        {"class": class_name, "index": int(record["index"])}
        for class_name, records in (("K7", k7["graphs"]), ("K6_only", k6["graphs"]))
        for record in records
    ]
    tagged_graphs = [
        {"class": class_name, **record}
        for class_name, records in (("K7", k7["graphs"]), ("K6_only", k6["graphs"]))
        for record in records
    ]
    combined = manifest["combined"]
    require(
        stable_hash(combined_indices) == EXPECTED_COMBINED_SHA256
        and stable_hash(sorted(combined_indices)) == EXPECTED_SORTED_SHA256
        and stable_hash(tagged_indices) == EXPECTED_TAGGED_INDICES_SHA256
        and stable_hash(tagged_graphs) == EXPECTED_TAGGED_GRAPHS_SHA256
        and combined.get("v6_count") == 644
        and combined.get("v7_count") == 263
        and combined.get("post_v6_exact_rejection_count") == 381
        and combined.get("post_v6_interval_incremental_rejection_count") == 2
        and combined.get("post_v6_union_rejection_count") == 383
        and combined.get("post_v7_exact_rejection_count") == 2
        and combined.get("count") == 261
        and combined.get("class_counts") == {"K7": 12, "K6_only": 249}
        and combined.get("ordered_indices_sha256") == EXPECTED_COMBINED_SHA256
        and combined.get("sorted_indices_sha256") == EXPECTED_SORTED_SHA256
        and combined.get("tagged_indices_sha256")
        == EXPECTED_TAGGED_INDICES_SHA256
        and combined.get("tagged_graphs_sha256")
        == EXPECTED_TAGGED_GRAPHS_SHA256
        and combined.get("cross_class_overlap") == 0,
        "combined accounting/hashes",
    )
    checks["combined_counts_disjointness_and_hashes"] = True

    expected_accounting = copy.deepcopy(v7["certificate_accounting"])
    expected_accounting["current_residue"] = 261
    expected_accounting["exact"]["layers"]["K6_saturated_singleton_basis"] = 2
    expected_accounting["exact"]["post_v6_rejections"] = 381
    expected_accounting["exact"]["total_rejections_from_v5"] = 645
    expected_accounting["union_rejections_from_v5"] = 650
    require(
        manifest.get("certificate_accounting") == expected_accounting
        and expected_accounting["interval"]["total_incremental_rejections_from_v5"]
        == 5
        and expected_accounting["interval"]["trust_assumptions"] == INTERVAL_TRUST,
        "certificate accounting/trust tiers",
    )
    expected_controls = copy.deepcopy(v7["controls"])
    expected_controls["known_realizable_18"][
        "exact_K6_saturated_singleton_basis"
    ] = report["positive_18_control"]
    expected_controls["exact_K6_saturated_singleton_synthetic"] = report[
        "synthetic_controls"
    ]
    require(manifest.get("controls") == expected_controls, "manifest controls")
    semantics = manifest.get("semantics", {})
    require(
        semantics == v7["semantics"]
        and "unconstrained" in semantics.get("nonedges", "")
        and "not a realizability" in semantics.get("residue", ""),
        "manifest semantics",
    )
    checks["trust_tiers_semantics_and_controls"] = True

    expected_source_roots = {
        "inherited_v7": inherited_v7_sources,
        "exact_K6_saturated_singleton_basis": singleton_sources,
    }
    require(
        manifest["source_boundary"].get("primary_artifacts")
        == dict(sorted(PRIMARY.items()))
        and manifest["source_boundary"].get("layer_source_roots")
        == expected_source_roots,
        "manifest evidence/source roots",
    )
    provenance = (
        package_boundary(manifest)
        if enforce_source_boundary
        else {"status": "SKIPPED_FOR_TEST"}
    )
    checks["source_and_artifact_roots"] = True

    execution = manifest.get("execution", {})
    require(
        isinstance(execution.get("command"), str)
        and bool(execution["command"])
        and isinstance(execution.get("generated_utc"), str)
        and "+00:00" in execution["generated_utc"]
        and isinstance(execution.get("platform"), str)
        and bool(execution["platform"])
        and isinstance(execution.get("machine"), str)
        and bool(execution["machine"])
        and isinstance(execution.get("python_version"), str)
        and bool(execution["python_version"])
        and isinstance(execution.get("logical_cpus"), int)
        and execution["logical_cpus"] > 0,
        "manifest execution provenance",
    )
    checks["execution_provenance"] = True

    syntax = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    imports = []
    for node in ast.walk(syntax):
        if isinstance(node, ast.Import):
            imports.extend(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imports.append(node.module)
    require(
        "build_d6_current_residue_manifest_v8" not in imports,
        "checker imports v8 builder",
    )
    require(
        not any(
            name.startswith(
                ("d6_k6_", "d6_k7_", "verify_d6_k6_", "verify_d6_k7_")
            )
            for name in imports
        ),
        "checker imports production module",
    )
    checks["import_independence"] = True
    require(all(checks.values()), "false v8 check")
    return {
        "schema": "d6-current-certified-residue-v8-verification-v1",
        "status": "PASS",
        "manifest": {"path": manifest_path.name, "sha256": manifest_hash},
        "counts": {
            "K7": 12,
            "K6_only": 249,
            "combined": 261,
            "exact_rejections_from_v5": 645,
            "interval_incremental_rejections": 5,
        },
        "class_indices_sha256": {
            "K7": EXPECTED_K7_SHA256,
            "K6_only": EXPECTED_K6_SHA256,
        },
        "ordered_residue_indices_sha256": EXPECTED_COMBINED_SHA256,
        "new_exact_layer": {
            "name": "saturated_singleton_basis",
            "rejected_indices": EXPECTED_REJECTED,
            "rejected_indices_sha256": EXPECTED_REJECTED_SHA256,
        },
        "trust_tiers": {
            "exact_floating_point": False,
            "interval_incremental_rejections": 5,
            "interval_trust_assumptions": INTERVAL_TRUST,
        },
        "source_boundary": provenance,
        "checks": checks,
        "execution": {
            "command": " ".join([sys.executable, *sys.argv]),
            "finished_utc": datetime.now(UTC).isoformat(),
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python_version": sys.version,
            "logical_cpus": os.cpu_count(),
            "verifier_sha256": sha256(Path(__file__).resolve()),
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
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        type=Path,
        default=ROOT / "d6_current_residue_manifest_v8.json",
    )
    parser.add_argument("--expected-manifest-sha256", required=True)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "d6_current_residue_manifest_v8_verification.json",
    )
    args = parser.parse_args()
    result = verify(args.manifest, args.expected_manifest_sha256)
    atomic_json(args.output.resolve(), result)
    print(
        json.dumps(
            {
                "output": str(args.output.resolve()),
                "sha256": sha256(args.output.resolve()),
                "status": result["status"],
                **result["counts"],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
