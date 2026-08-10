#!/usr/bin/env python3
"""Import-independent structural checker for the d=6 residue manifest v7."""

from __future__ import annotations

import argparse
import ast
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
    "build_d6_current_residue_manifest_v7.py",
    "verify_d6_current_residue_manifest_v7.py",
    "test_d6_current_residue_manifest_v7.py",
    "d6_current_residue_manifest_v7.md",
}
PRIMARY = {
    "d6_current_residue_manifest_v6.json": "c513bcc40c037af23340dce3cc91d593a101dd1d073f93feab82abd556efadc9",
    "d6_current_residue_manifest_v6_verification.json": "3d11cbd20cb59b42a8c8fd5252fb3ef3dea649576982ea8017622fad3279a7bc",
    "d6_interval_v6_k7_cap500000_report.json": "0b04a7bbdcbb5bba264099e3cb9f2c9e39d4bc50f7bf7c211ba5ce81e4660704",
    "d6_interval_v6_k7_cap500000_verification.json": "105550c0ba2b181b1cfec64afd6fefa1324d4cfcde41aee385c99a769d8f1fb7",
    "d6_interval_v6_k7_cap500000_decisions.tsv": "61eb84688c2a24ae8ba8f727c5e2f2250693f15ab74452c21e47beb6c126a869",
    "d6_interval_v6_k7_cap500000_checkpoints.tar.gz": "e4dab680fbc1315be24271dff95769f9d6b7c6f172544fe232b1345ba07589c0",
    "d6_k7_one_two_star_increment_report.json": "12b3d18b1ea81961f58d831d4c7c322fbceb2300e7533b64161e8011fbe9d1ec",
    "d6_k7_one_two_star_increment_verification.json": "c6f40578685036cb1156cb0ea8bf06d6004a70d5f376aff5ac1e16a1c0ca09eb",
    "d6_interval_k7_3936177_cap2000000_report.json": "4135ffbb2a54ab61b28061f0a8d879240b20b1e363c7751337c2bde0559de70e",
    "d6_interval_k7_3936177_cap2000000_verification.json": "b2d15565fb246f152e5d12bc5cfd188bbe187bd2a554ee2f67de821f57b2553c",
    "d6_interval_k7_3936177_cap2000000_decisions.tsv": "960009a2343e719d7d0db1725e5bce1aae858f6c5c908efd3b14aacedd657866",
    "d6_interval_k7_3936177_cap2000000_checkpoints.tar.gz": "31de6fb88ed8f1477fc7b4b57fcb4f3e039ca518c88868ff10ae1db2c1ecda76",
    "d6_k7_schur_3949382_increment_report.json": "30efa40d1db8a3aef03af9e56aa915a04594aabe8b72e5a1c72489e2a8c952b8",
    "d6_k7_schur_3949382_increment_verification.json": "18c9980baa74596523eba9b447af1ae5830a9558d80981e1fa3ff6c984bcb5bc",
    "d6_k7_star_2593240_increment_report.json": "dab2e54c0d09c687cbe9c8488f0580e7777a6ef1fe02eeda675a171c746cabc5",
    "d6_k7_star_2593240_increment_verification.json": "384f821307f1860bed2d4227382c05d9fe4fd14837ffeec6a1d16c5d98b84750",
    "d6_k6_tight_same_z0_report.json": "2f2fb63e91518717f3a340c4e9401f720a3393ee6816535eadff1b89d681be04",
    "d6_k6_tight_same_z0_verification.json": "d87ecd688301e2e811e8eb3f3fe52111a050765643de6a5dbb1f1c06582b311e",
    "d6_k6_tight_same_z0_certificates.json": "ad408071d49d95724c6b9a45afca947b937ea1ce8b79fb2b67cd810a91756a11",
    "d6_k6_opposite_ray_rank_conjunction_report.json": "3415f48f498758a5d0daac9403012e35186eb7b620fe9c0324e4674f7a11fcc8",
    "d6_k6_opposite_ray_rank_conjunction_verification.json": "e9e64d93faad85af321011b7f48842d13c416f3e61ec802455972d13d48b334e",
    "d6_k6_opposite_ray_rank_conjunction_certificates.json.gz": "99a86290da7221f0d2a07ee15ed93432b19dcd58a9ba2bcebaa2956167ac1e37",
}
EXPECTED_K7 = [
    316173,
    2581209,
    3648882,
    3729907,
    3935560,
    3936310,
    3936435,
    3945490,
    3945555,
    3945557,
    3945564,
    3947605,
]
EXPECTED_K7_SHA256 = "e545ea5a8fc9691ed98e885919e149f84dc7e8b4796d0fe9e6cb3a519df6fd8c"
EXPECTED_K6_SHA256 = "a23da7bd9ae21eb7c078a08a7bce92f22d397d5c71f1e4f51d9cf8e59bbfa907"
EXPECTED_COMBINED_SHA256 = "a50baf54bfff642e7c1711fe7041b5184df12f90de2cf52557cd086f10875127"
EXPECTED_SORTED_SHA256 = "ad495b15b84913cb49c8ef1633d27baaadec37e23176c8d89a6dce06396f39cd"
EXPECTED_V6_VERIFIER_SHA256 = (
    "664aa48f3bfed1acc983e2b4ccd1b5546e86ad8e5a57f89e73d20f20bd2c63fa"
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
                if not rows[first] & (1 << third) and not rows[second] & (1 << third):
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


def exact_partition(
    report: dict,
    verification: dict,
    *,
    schema: str,
    verification_schema: str,
    report_name: str,
    input_count: int,
    rejection_count: int,
    residue_count: int,
) -> tuple[list[int], list[int], list[int]]:
    inputs = list(map(int, report.get("ordered_input_indices", ())))
    rejected = list(map(int, report.get("rejected_indices", ())))
    residue = list(map(int, report.get("ordered_residue_indices", ())))
    require(
        report.get("schema") == schema
        and report.get("status") == "COMPLETE"
        and report.get("input_graphs") == input_count
        and report.get("graphs_rejected") == rejection_count
        and report.get("graphs_surviving") == residue_count
        and verification.get("schema") == verification_schema
        and verification.get("status") == "PASS"
        and verification.get("report", {}).get("sha256") == PRIMARY[report_name]
        and verification.get("graphs_recomputed") == input_count
        and verification.get("graphs_rejected") == rejection_count
        and verification.get("graphs_surviving") == residue_count,
        f"{report_name}: report/check",
    )
    require(
        len(inputs) == input_count
        and len(set(inputs)) == input_count
        and len(rejected) == rejection_count
        and len(set(rejected)) == rejection_count
        and residue == [index for index in inputs if index not in set(rejected)]
        and stable_hash(inputs) == report.get("ordered_input_indices_sha256")
        and stable_hash(rejected) == report.get("rejected_indices_sha256")
        and stable_hash(residue) == report.get("ordered_residue_indices_sha256")
        and verification.get("rejected_indices") == rejected
        and verification.get("rejected_indices_sha256") == stable_hash(rejected)
        and verification.get("ordered_residue_indices_sha256") == stable_hash(residue),
        f"{report_name}: partition",
    )
    archive = report["certificate_archive"]
    require(
        archive.get("sha256") == PRIMARY[archive["path"]]
        and archive.get("rejected_graphs") == rejection_count,
        f"{report_name}: archive",
    )
    require(
        "exact" in report.get("semantics", {}).get("arithmetic", "")
        and "unconstrained" in report.get("semantics", {}).get("candidate_nonedges", "")
        and report.get("positive_18_control", {}).get("passed") is True
        and verification.get("positive_18_control", {}).get("passed") is True,
        f"{report_name}: semantics/control",
    )
    return inputs, rejected, residue


def interval_layer(
    report: dict,
    verification: dict,
    *,
    report_name: str,
    input_indices: list[int],
    killed: list[int],
    cap: int,
    focused: bool,
) -> dict:
    require(
        report.get("schema") == 2
        and report.get("total_graphs") == len(input_indices)
        and report.get("completed_graphs") == len(input_indices)
        and report.get("certified_killed") == len(killed)
        and report.get("unresolved_total") == len(input_indices) - len(killed)
        and verification.get("status") == "PASS"
        and verification.get("report", {}).get("sha256") == PRIMARY[report_name]
        and verification.get("certified_killed") == len(killed)
        and verification.get("search", {}).get("cap") == cap
        and verification.get("trust_assumptions") == INTERVAL_TRUST,
        f"{report_name}: boundary",
    )
    results = report.get("results", [])
    require([int(row["index"]) for row in results] == input_indices, "interval order")
    require(
        [int(row["index"]) for row in results if row.get("status") == "KILLED"] == killed,
        "interval KILLED",
    )
    require(
        all(row.get("status") in {"KILLED", "ABORT", "UNRESOLVED", "INFRA_ERROR"} for row in results),
        "interval statuses",
    )
    selection = verification["selection"]
    require(
        selection.get("graphs") == len(input_indices)
        and selection.get("indices_sha256") == stable_hash(input_indices)
        and selection.get("v6_manifest_sha256") == PRIMARY["d6_current_residue_manifest_v6.json"]
        and selection.get("v6_verification_sha256") == PRIMARY["d6_current_residue_manifest_v6_verification.json"],
        "interval selection",
    )
    if focused:
        require(selection.get("index") == killed[0], "focused target")
        decisions = "d6_interval_k7_3936177_cap2000000_decisions.tsv"
    else:
        require(verification.get("killed_indices") == killed, "cap killed list")
        decisions = "d6_interval_v6_k7_cap500000_decisions.tsv"
    require(
        report.get("decisions", {}).get("sha256") == PRIMARY[decisions]
        and verification.get("checkpoint", {}).get("decision_TSV_sha256") == PRIMARY[decisions]
        and verification.get("kernel_controls", {}).get("status") == "PASS"
        and verification.get("kernel_controls", {}).get("records_replayed") == 48
        and verification.get("witness_replay", {}).get("status") == "PASS"
        and verification.get("witness_replay", {}).get("slices") == 24,
        "interval controls/roots",
    )
    sources = report["configuration"]["sources"]
    commit = report["configuration"]["git"]["commit"]
    require(verification["provenance"]["commit"] == commit, "interval commit")
    paths = {
        "cdriver6.py": "cdriver6.py",
        "ckernel6.c": "ckernel6.c",
        "interval_runner": "run_d6_interval_residue.py",
        "ival.py": "ival.py",
        "v6_wrapper" if focused else "wrapper": "run_d6_interval_v6.py",
    }
    if focused:
        paths["focused_wrapper"] = "run_d6_interval_k7_3936177.py"
    for key, path in paths.items():
        require(git_blob_sha256(commit, path) == sources[key], f"interval blob {key}")
    require(sha256(ROOT / "ckernel6.dylib") == sources["ckernel6.dylib"], "kernel dylib")
    verifier_name = "verify_d6_interval_k7_3936177.py" if focused else "verify_d6_interval_v6.py"
    verifier_hash = verification["verifier_source_sha256"]
    common_hash = verification["common_helper_sha256"]
    require(
        sha256(ROOT / verifier_name) == verifier_hash
        and git_blob_sha256(commit, verifier_name) == verifier_hash
        and sha256(ROOT / "verify_d6_interval_benchmarks.py") == common_hash
        and git_blob_sha256(commit, "verify_d6_interval_benchmarks.py")
        == common_hash,
        "interval verifier/helper",
    )
    return {
        "launch_commit": commit,
        "configuration_sources": dict(sorted(sources.items())),
        "verifier": {
            verifier_name: verifier_hash,
            "verify_d6_interval_benchmarks.py": common_hash,
        },
    }


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

    v6 = load("d6_current_residue_manifest_v6.json")
    v6_check = load("d6_current_residue_manifest_v6_verification.json")
    require(
        v6.get("schema") == "d6-current-certified-residue-v6"
        and v6.get("status") == "COMPLETE_MIXED_CERTIFICATE_UNION"
        and v6.get("combined", {}).get("count") == 644
        and v6_check.get("status") == "PASS"
        and v6_check.get("manifest", {}).get("sha256") == PRIMARY["d6_current_residue_manifest_v6.json"]
        and all(v6_check.get("checks", {}).values()),
        "v6 boundary",
    )
    v6_k7 = list(map(int, v6["classes"]["K7"]["indices"]))
    v6_k6 = list(map(int, v6["classes"]["K6_only"]["indices"]))
    v6_k7_graphs = v6["classes"]["K7"]["graphs"]
    v6_k6_graphs = v6["classes"]["K6_only"]["graphs"]
    require(len(v6_k7) == 19 and len(v6_k6) == 625, "v6 class counts")
    require(v6["certificate_accounting"]["interval"]["trust_assumptions"] == INTERVAL_TRUST, "v6 trust")
    named_hashes(v6["source_boundary"]["transitive_files"], "v6 transitive")
    require(sha256(ROOT / "build_d6_current_residue_manifest_v6.py") == v6_check["builder_sha256"], "v6 builder")
    require(
        sha256(ROOT / "verify_d6_current_residue_manifest_v6.py")
        == EXPECTED_V6_VERIFIER_SHA256,
        "v6 verifier",
    )
    checks["v6_boundary_and_roots"] = True

    cap = load("d6_interval_v6_k7_cap500000_report.json")
    cap_check = load("d6_interval_v6_k7_cap500000_verification.json")
    cap_sources = interval_layer(
        cap,
        cap_check,
        report_name="d6_interval_v6_k7_cap500000_report.json",
        input_indices=v6_k7,
        killed=[3595554],
        cap=500_000,
        focused=False,
    )
    focused = load("d6_interval_k7_3936177_cap2000000_report.json")
    focused_check = load("d6_interval_k7_3936177_cap2000000_verification.json")
    focused_sources = interval_layer(
        focused,
        focused_check,
        report_name="d6_interval_k7_3936177_cap2000000_report.json",
        input_indices=[3936177],
        killed=[3936177],
        cap=2_000_000,
        focused=True,
    )
    checks["interval_layers_sources_controls"] = True

    star = load("d6_k7_one_two_star_increment_report.json")
    star_check = load("d6_k7_one_two_star_increment_verification.json")
    star_bad = list(map(int, star["summary"]["rejected_indices"]))
    star_left = list(map(int, star["summary"]["ordered_survivor_indices"]))
    require(
        star.get("kind") == "d6_k7_one_free_neighbour_two_free_center_increment"
        and star.get("status") == "COMPLETE"
        and star.get("input", {}).get("ordered_indices") == v6_k7
        and star.get("input", {}).get("embedded_graphs_sha256") == stable_hash(v6_k7_graphs)
        and star_bad == [2592657, 3785980, 3888410]
        and star_left == [index for index in v6_k7 if index not in set(star_bad)]
        and star_check.get("status") == "PASS"
        and star_check.get("report", {}).get("sha256") == PRIMARY["d6_k7_one_two_star_increment_report.json"]
        and star_check.get("rejected_indices") == star_bad
        and star_check.get("checked", {}).get("positive_18_control") is True
        and star_check.get("checked", {}).get("optional_zero_control") is True
        and star.get("semantics", {}).get("floating_point_enters_rejection") is False
        and star.get("semantics", {}).get("candidate_nonedges_optional") is True,
        "star exact layer",
    )
    star_sources_map = named_hashes(star["source_sha256"], "star source")
    star_commit = star["execution"]["git"]["commit"]
    for name, expected in star_sources_map.items():
        require(git_blob_sha256(star_commit, name) == expected, f"star blob: {name}")
    star_verifier = star_check["verifier_source_sha256"]
    require(
        sha256(ROOT / "verify_d6_k7_one_two_star_increment.py") == star_verifier
        and git_blob_sha256(star_commit, "verify_d6_k7_one_two_star_increment.py") == star_verifier,
        "star verifier",
    )
    star_sources = {
        "launch_commit": star_commit,
        "producer_and_kernels": star_sources_map,
        "verifier": {"verify_d6_k7_one_two_star_increment.py": star_verifier},
        "upstream": dict(sorted(star["upstream_artifact_sha256"].items())),
    }
    checks["star_exact_partition_sources_controls"] = True

    schur = load("d6_k7_schur_3949382_increment_report.json")
    schur_check = load("d6_k7_schur_3949382_increment_verification.json")
    schur_bad = list(map(int, schur["summary"]["rejected_indices"]))
    require(
        schur.get("kind") == "d6_k7_schur_3949382_exact_increment"
        and schur.get("status") == "COMPLETE_EXACT_REJECTION"
        and schur.get("input", {}).get("ordered_indices") == star_left
        and schur_bad == [3949382]
        and schur.get("certificate_sha256") == stable_hash(schur.get("certificate"))
        and schur_check.get("status") == "PASS"
        and schur_check.get("report", {}).get("sha256") == PRIMARY["d6_k7_schur_3949382_increment_report.json"]
        and schur_check.get("checked", {}).get("unresolved_optional_entries") == 0
        and schur.get("semantics", {}).get("floating_point_enters_rejection") is False
        and schur.get("semantics", {}).get("candidate_nonedges_optional") is True
        and schur.get("controls", {}).get("nonnegative_boundary_solution", {}).get("rejected_by_strict_domain") is True,
        "Schur exact layer",
    )
    schur_sources_map = named_hashes(schur["source_sha256"], "Schur source")
    schur_commit = schur["execution"]["git"]["commit"]
    for name, expected in schur_sources_map.items():
        require(git_blob_sha256(schur_commit, name) == expected, f"Schur blob: {name}")
    schur_sources = {
        "launch_commit": schur_commit,
        "package_sources": schur_sources_map,
        "upstream": dict(sorted(schur["upstream_sha256"].items())),
    }
    checks["Schur_exact_source_bound"] = True

    tight = load("d6_k6_tight_same_z0_report.json")
    tight_check = load("d6_k6_tight_same_z0_verification.json")
    opposite = load("d6_k6_opposite_ray_rank_conjunction_report.json")
    opposite_check = load("d6_k6_opposite_ray_rank_conjunction_verification.json")
    tight_input, tight_bad, tight_left = exact_partition(
        tight,
        tight_check,
        schema="d6-k6-tight-same-z0-v1",
        verification_schema="d6-k6-tight-same-z0-verification-v1",
        report_name="d6_k6_tight_same_z0_report.json",
        input_count=625,
        rejection_count=2,
        residue_count=623,
    )
    opposite_input, opposite_bad, k6_indices = exact_partition(
        opposite,
        opposite_check,
        schema="d6-k6-opposite-ray-rank-conjunction-v1",
        verification_schema="d6-k6-opposite-ray-rank-conjunction-verification-v1",
        report_name="d6_k6_opposite_ray_rank_conjunction_report.json",
        input_count=623,
        rejection_count=372,
        residue_count=251,
    )
    require(
        tight_input == v6_k6
        and tight_left == opposite_input
        and not set(tight_bad) & set(opposite_bad)
        and stable_hash(k6_indices) == EXPECTED_K6_SHA256,
        "K6 exact chain",
    )
    tight_dependencies = named_hashes(tight["dependencies"], "tight dependencies")
    tight_sources_map = named_hashes(tight_check["sources"], "tight sources")
    require(sha256(ROOT / "d6_k6_tight_same_z0.py") == tight["production_source_sha256"], "tight producer")
    opposite_dependencies = named_hashes(opposite["dependencies"], "opposite dependencies")
    opposite_sources_map = named_hashes(opposite["source_sha256"], "opposite sources")
    opposite_check_sources = named_hashes(opposite_check["sources"], "opposite check sources")
    opposite_commit = opposite["execution"]["git"]["commit"]
    for name, expected in opposite_sources_map.items():
        require(git_blob_sha256(opposite_commit, name) == expected, f"opposite blob: {name}")
    require(
        opposite["execution"]["git"].get("tracked_clean") is True
        and opposite_check["provenance"].get("tracked_sources_clean_at_launch") is True,
        "opposite provenance",
    )
    k6_sources = {
        "tight_same_Z0": {
            "dependencies": tight_dependencies,
            "producer_and_verifier_sources": tight_sources_map,
        },
        "opposite_ray_conjunction": {
            "launch_commit": opposite_commit,
            "dependencies": opposite_dependencies,
            "production_sources": opposite_sources_map,
            "verifier_sources": opposite_check_sources,
        },
    }
    checks["K6_exact_chain_sources_controls"] = True

    star_2593240 = load("d6_k7_star_2593240_increment_report.json")
    star_2593240_check = load(
        "d6_k7_star_2593240_increment_verification.json"
    )
    schur_left = list(map(int, schur["summary"]["ordered_survivor_indices"]))
    new_2593240_bad = list(
        map(int, star_2593240["summary"]["rejected_indices"])
    )
    new_2593240_left = list(
        map(int, star_2593240["summary"]["ordered_survivor_indices"])
    )
    require(
        star_2593240.get("kind")
        == "d6_k7_rankone_star_2593240_exact_increment"
        and star_2593240.get("status") == "COMPLETE_EXACT_REJECTION"
        and star_2593240.get("input", {}).get("ordered_indices") == schur_left
        and star_2593240.get("input", {}).get("ordered_indices_sha256")
        == stable_hash(schur_left)
        and new_2593240_bad == [2593240]
        and new_2593240_left
        == [index for index in schur_left if index != 2593240]
        and star_2593240.get("certificate_sha256")
        == stable_hash(star_2593240.get("certificate"))
        and star_2593240_check.get("kind")
        == "d6_k7_star_2593240_increment_verification"
        and star_2593240_check.get("status") == "PASS"
        and star_2593240_check.get("report", {}).get("sha256")
        == PRIMARY["d6_k7_star_2593240_increment_report.json"]
        and star_2593240_check.get("conclusion", {}).get("rejected_indices")
        == new_2593240_bad
        and star_2593240_check.get("checked", {}).get("required_k7_seeds") == 2
        and star_2593240_check.get("checked", {}).get("new_families") == 2
        and star_2593240_check.get("checked", {}).get(
            "unresolved_optional_entries"
        )
        == 0,
        "2593240 exact layer",
    )
    require(
        star_2593240.get("semantics", {}).get("candidate_nonedges_optional")
        is True
        and star_2593240.get("semantics", {}).get(
            "floating_point_enters_rejection"
        )
        is False
        and star_2593240.get("semantics", {}).get("both_required_K7_seeds_checked")
        is True
        and star_2593240_check.get("semantics")
        == star_2593240.get("semantics")
        and star_2593240.get("controls", {})
        .get("positive_tetrad_star_control", {})
        .get("all_thirty_tetrads_zero")
        is True
        and star_2593240.get("controls", {})
        .get("positive_tetrad_star_control", {})
        .get("rank_one_completion_exists")
        is False,
        "2593240 semantics/control",
    )
    new_2593240_source_map = named_hashes(
        star_2593240["source_sha256"], "2593240 sources"
    )
    new_2593240_commit = star_2593240["execution"]["git"]["commit"]
    require(
        new_2593240_commit
        == "64cde398d8d116632551b6f40edfc5de1b84405b"
        and star_2593240["execution"]["git"].get(
            "proof_and_checker_sources_equal_committed_blobs"
        )
        is True
        and star_2593240_check.get("source_boundary", {}).get("commit")
        == new_2593240_commit
        and star_2593240_check.get("source_boundary", {}).get(
            "proof_and_checker_sources_equal_committed_blobs"
        )
        is True,
        "2593240 source boundary",
    )
    for name, expected in new_2593240_source_map.items():
        require(
            git_blob_sha256(new_2593240_commit, name) == expected,
            f"2593240 blob: {name}",
        )
    new_2593240_sources = {
        "launch_commit": new_2593240_commit,
        "package_sources": new_2593240_source_map,
        "upstream": dict(sorted(star_2593240["upstream_sha256"].items())),
    }
    checks["star_2593240_exact_source_bound"] = True

    exact_k7_set = set(star_bad + schur_bad + new_2593240_bad)
    interval_k7_set = {3595554, 3936177}
    require(not exact_k7_set & interval_k7_set, "K7 trust-tier overlap")
    k7_indices = [index for index in v6_k7 if index not in exact_k7_set | interval_k7_set]
    require(k7_indices == EXPECTED_K7 and stable_hash(k7_indices) == EXPECTED_K7_SHA256, "K7 final")
    checks["K7_layer_union_disjointness"] = True

    manifest = load(manifest_path)
    require(
        manifest.get("schema") == "d6-current-certified-residue-v7"
        and manifest.get("status") == "COMPLETE_MIXED_CERTIFICATE_UNION"
        and manifest.get("class_order") == ["K7", "K6_only"],
        "v7 schema",
    )
    embedded_k7 = manifest["classes"]["K7"]["graphs"]
    embedded_k6 = manifest["classes"]["K6_only"]["graphs"]
    expected_k7_graphs = [record for record in v6_k7_graphs if int(record["index"]) in set(k7_indices)]
    expected_k6_graphs = [record for record in v6_k6_graphs if int(record["index"]) in set(k6_indices)]
    require(embedded_k7 == expected_k7_graphs and embedded_k6 == expected_k6_graphs, "embedded graphs")
    for class_name, records in (("K7", embedded_k7), ("K6_only", embedded_k6)):
        for record in records:
            rows = validate_rows(record.get("adjacency"), f"{class_name}/{record.get('index')}")
            require(alpha_at_most_two(rows) and contains_clique(rows, 6), "graph invariant")
            require(contains_clique(rows, 7) == (class_name == "K7"), "graph class")
    checks["embedded_graphs_and_classes"] = True

    k7 = manifest["classes"]["K7"]
    k6 = manifest["classes"]["K6_only"]
    k6_bad = [index for index in v6_k6 if index in set(tight_bad + opposite_bad)]
    all_k7_bad = [index for index in v6_k7 if index in exact_k7_set | interval_k7_set]
    exact_k7 = [index for index in v6_k7 if index in exact_k7_set]
    interval_k7 = [index for index in v6_k7 if index in interval_k7_set]
    expected_k7_layers = [
        {
            "name": "interval_cap500000",
            "trust_tier": "interval",
            "input_count": 19,
            "input_indices_sha256": stable_hash(v6_k7),
            "rejection_count": 1,
            "rejected_indices": [3595554],
            "rejected_indices_sha256": stable_hash([3595554]),
        },
        {
            "name": "one_two_star",
            "trust_tier": "exact",
            "input_count": 19,
            "input_indices_sha256": stable_hash(v6_k7),
            "rejection_count": 3,
            "rejected_indices": star_bad,
            "rejected_indices_sha256": stable_hash(star_bad),
            "residue_count": len(star_left),
            "residue_indices_sha256": stable_hash(star_left),
        },
        {
            "name": "focused_interval_cap2000000",
            "trust_tier": "interval",
            "input_count": 1,
            "input_indices_sha256": stable_hash([3936177]),
            "rejection_count": 1,
            "rejected_indices": [3936177],
            "rejected_indices_sha256": stable_hash([3936177]),
        },
        {
            "name": "saturating_clique_Schur",
            "trust_tier": "exact",
            "input_count": len(star_left),
            "input_indices_sha256": stable_hash(star_left),
            "rejection_count": 1,
            "rejected_indices": schur_bad,
            "rejected_indices_sha256": stable_hash(schur_bad),
            "residue_count": len(schur_left),
            "residue_indices_sha256": stable_hash(schur_left),
        },
        {
            "name": "rankone_star_2593240",
            "trust_tier": "exact",
            "input_count": len(schur_left),
            "input_indices_sha256": stable_hash(schur_left),
            "rejection_count": 1,
            "rejected_indices": new_2593240_bad,
            "rejected_indices_sha256": stable_hash(new_2593240_bad),
            "residue_count": len(new_2593240_left),
            "residue_indices_sha256": stable_hash(new_2593240_left),
        },
    ]
    expected_k6_layers = [
        {
            "name": "tight_same_Z0",
            "input_count": 625,
            "input_indices_sha256": stable_hash(v6_k6),
            "rejection_count": 2,
            "rejected_indices_sha256": stable_hash(tight_bad),
            "residue_count": 623,
            "residue_indices_sha256": stable_hash(tight_left),
        },
        {
            "name": "opposite_ray_rank_conjunction",
            "input_count": 623,
            "input_indices_sha256": stable_hash(opposite_input),
            "rejection_count": 372,
            "rejected_indices_sha256": stable_hash(opposite_bad),
            "residue_count": 251,
            "residue_indices_sha256": EXPECTED_K6_SHA256,
        },
    ]
    require(
        k7.get("v6_count") == 19
        and k7.get("layers") == expected_k7_layers
        and k7.get("post_v6_exact_rejection_count") == 5
        and k7.get("post_v6_exact_rejected_indices") == exact_k7
        and k7.get("post_v6_exact_rejected_indices_sha256")
        == stable_hash(exact_k7)
        and k7.get("post_v6_interval_rejection_count") == 2
        and k7.get("post_v6_interval_rejected_indices") == interval_k7
        and k7.get("post_v6_interval_rejected_indices_sha256")
        == stable_hash(interval_k7)
        and k7.get("post_v6_union_rejection_count") == 7
        and k7.get("post_v6_union_rejected_indices") == all_k7_bad
        and k7.get("post_v6_union_rejected_indices_sha256")
        == stable_hash(all_k7_bad)
        and k7.get("count") == 12
        and k7.get("indices") == k7_indices
        and k7.get("indices_sha256") == EXPECTED_K7_SHA256
        and k7.get("graphs_sha256") == stable_hash(embedded_k7)
        and k6.get("v6_count") == 625
        and k6.get("exact_layers") == expected_k6_layers
        and k6.get("post_v6_exact_rejected_indices") == k6_bad
        and k6.get("post_v6_exact_rejection_count") == 374
        and k6.get("post_v6_exact_rejected_indices_sha256")
        == stable_hash(k6_bad)
        and k6.get("count") == 251
        and k6.get("indices") == k6_indices
        and k6.get("indices_sha256") == EXPECTED_K6_SHA256
        and k6.get("graphs_sha256") == stable_hash(embedded_k6),
        "class accounting",
    )
    combined_indices = k7_indices + k6_indices
    tagged_indices = [
        {"class": class_name, "index": int(record["index"])}
        for class_name, records in (("K7", embedded_k7), ("K6_only", embedded_k6))
        for record in records
    ]
    tagged_graphs = [
        {"class": class_name, **record}
        for class_name, records in (("K7", embedded_k7), ("K6_only", embedded_k6))
        for record in records
    ]
    combined = manifest["combined"]
    require(
        combined.get("v6_count") == 644
        and combined.get("post_v6_exact_rejection_count") == 379
        and combined.get("post_v6_interval_incremental_rejection_count") == 2
        and combined.get("post_v6_union_rejection_count") == 381
        and combined.get("count") == 263
        and combined.get("class_counts") == {"K7": 12, "K6_only": 251}
        and combined.get("ordered_indices_sha256") == EXPECTED_COMBINED_SHA256
        and combined.get("sorted_indices_sha256") == EXPECTED_SORTED_SHA256
        and combined.get("tagged_indices_sha256") == stable_hash(tagged_indices)
        and combined.get("tagged_graphs_sha256") == stable_hash(tagged_graphs)
        and combined.get("cross_class_overlap") == 0,
        "combined accounting",
    )
    checks["counts_and_hashes"] = True

    accounting = manifest["certificate_accounting"]
    require(
        accounting["exact"].get("post_v6_rejections") == 379
        and accounting["exact"].get("total_rejections_from_v5") == 643
        and accounting["exact"].get("floating_point_enters_rejection") is False
        and accounting["interval"].get("post_v6_incremental_rejections") == 2
        and accounting["interval"].get("total_incremental_rejections_from_v5") == 5
        and accounting["interval"].get("only_KILLED_used_for_rejection") is True
        and accounting["interval"].get("ABORT_UNRESOLVED_INFRA_ERROR_used_for_rejection") is False
        and accounting["interval"].get("trust_assumptions") == INTERVAL_TRUST
        and accounting.get("union_rejections_from_v5") == 648
        and accounting.get("current_residue") == 263,
        "trust-tier accounting",
    )
    semantics = manifest["semantics"]
    require(
        "unconstrained" in semantics.get("nonedges", "")
        and "KILLED" in semantics.get("interval_rejection", "")
        and "not a realizability" in semantics.get("residue", ""),
        "manifest semantics",
    )
    controls = manifest["controls"]
    require(
        controls["known_realizable_18"]["inherited_v6"]["passed"] is True
        and controls["known_realizable_18"]["exact_star"]["passed"] is True
        and controls["known_realizable_18"]["exact_K6_tight_same_Z0"]["passed"]
        is True
        and controls["known_realizable_18"][
            "exact_K6_opposite_ray_conjunction"
        ]["passed"]
        is True
        and controls["optional_zero_star"]["passed"] is True
        and controls["interval_cap500000_kernel"]["status"] == "PASS"
        and controls["interval_focused_3936177_kernel"]["status"] == "PASS"
        and controls["exact_Schur_nonnegative_boundary"]
        == {"checked": True, "strict_domain_required": True}
        and controls["exact_2593240_tetrad_star"]
        == {
            "two_required_K7_seeds_checked": True,
            "all_thirty_tetrads_zero_at_control": True,
            "rank_one_completion_exists": False,
        },
        "manifest controls",
    )
    checks["trust_semantics_controls"] = True

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

    expected_source_roots = {
        "inherited_v6": {
            "manifest_source_boundary": v6["source_boundary"],
            "builder": {
                "build_d6_current_residue_manifest_v6.py": v6_check[
                    "builder_sha256"
                ]
            },
            "verifier": {
                "verify_d6_current_residue_manifest_v6.py": (
                    EXPECTED_V6_VERIFIER_SHA256
                )
            },
        },
        "interval_cap500000": cap_sources,
        "exact_star": star_sources,
        "interval_focused_3936177": focused_sources,
        "exact_schur_3949382": schur_sources,
        "exact_2593240": new_2593240_sources,
        "exact_K6_chain": k6_sources,
    }
    require(
        manifest["source_boundary"].get("primary_artifacts") == dict(sorted(PRIMARY.items()))
        and manifest["source_boundary"].get("layer_source_roots") == expected_source_roots,
        "manifest evidence/source roots",
    )
    provenance = package_boundary(manifest) if enforce_source_boundary else {"status": "SKIPPED_FOR_TEST"}
    checks["source_and_artifact_roots"] = True

    syntax = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    imports = []
    for node in ast.walk(syntax):
        if isinstance(node, ast.Import):
            imports.extend(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imports.append(node.module)
    require("build_d6_current_residue_manifest_v7" not in imports, "checker imports builder")
    require(
        not any(name.startswith(("d6_k6_", "d6_k7_", "verify_d6_k6_", "verify_d6_k7_")) for name in imports),
        "checker imports production module",
    )
    checks["import_independence"] = True
    require(all(checks.values()), "false v7 check")
    return {
        "schema": "d6-current-certified-residue-v7-verification-v1",
        "status": "PASS",
        "manifest": {"path": manifest_path.name, "sha256": manifest_hash},
        "counts": {
            "K7": 12,
            "K6_only": 251,
            "combined": 263,
            "exact_rejections_from_v5": 643,
            "interval_incremental_rejections": 5,
        },
        "class_indices_sha256": {"K7": EXPECTED_K7_SHA256, "K6_only": EXPECTED_K6_SHA256},
        "ordered_residue_indices_sha256": EXPECTED_COMBINED_SHA256,
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
    parser.add_argument("--manifest", type=Path, default=ROOT / "d6_current_residue_manifest_v7.json")
    parser.add_argument("--expected-manifest-sha256", required=True)
    parser.add_argument("--output", type=Path, default=ROOT / "d6_current_residue_manifest_v7_verification.json")
    args = parser.parse_args()
    result = verify(args.manifest, args.expected_manifest_sha256)
    atomic_json(args.output.resolve(), result)
    print(json.dumps({"output": str(args.output.resolve()), "sha256": sha256(args.output.resolve()), "status": result["status"], **result["counts"]}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
