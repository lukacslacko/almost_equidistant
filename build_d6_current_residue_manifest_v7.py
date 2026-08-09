#!/usr/bin/env python3
"""Build the source-bound dimension-six current residue manifest v7.

The v7 boundary unions every independently checked certificate added after
v6.  It deliberately keeps exact and interval-dependent rejection accounting
separate.  No production kernel is imported: this file consumes immutable
reports, checks their roots and partitions, and embeds the surviving graphs.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import subprocess
import sys
import time
from datetime import UTC, datetime
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / "d6_current_residue_manifest_v7.json"
PACKAGE_FILES = (
    "build_d6_current_residue_manifest_v7.py",
    "verify_d6_current_residue_manifest_v7.py",
    "test_d6_current_residue_manifest_v7.py",
    "d6_current_residue_manifest_v7.md",
)
PRIMARY_ARTIFACTS = {
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


def load(name: str) -> dict:
    value = json.loads((ROOT / name).read_text(encoding="utf-8"))
    require(isinstance(value, dict), f"{name} is not an object")
    return value


def verify_named_hashes(mapping: dict, label: str) -> dict[str, str]:
    require(isinstance(mapping, dict) and mapping, f"empty {label}")
    normalized = {}
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
    require(branch == "codex/dimension6", "official v7 launch branch")
    sources = {name: sha256(ROOT / name) for name in PACKAGE_FILES}
    for name, expected in sources.items():
        require(git_blob_sha256(commit, name) == expected, f"uncommitted package source: {name}")
    porcelain = git("status", "--porcelain=v1", "--untracked-files=all")
    return {
        "commit": commit,
        "branch": branch,
        "package_sources": dict(sorted(sources.items())),
        "porcelain_lines": porcelain.splitlines(),
        "porcelain_sha256": hashlib.sha256(porcelain.encode("utf-8")).hexdigest(),
        "package_sources_equal_committed_blobs": True,
    }


def validate_interval_sources(report: dict, verification: dict, *, focused: bool) -> dict:
    sources = report["configuration"]["sources"]
    commit = report["configuration"]["git"]["commit"]
    require(verification["provenance"]["commit"] == commit, "interval source commit")
    paths = {
        "cdriver6.py": "cdriver6.py",
        "ckernel6.c": "ckernel6.c",
        "interval_runner": "run_d6_interval_residue.py",
        "ival.py": "ival.py",
        "v6_wrapper" if focused else "wrapper": "run_d6_interval_v6.py",
    }
    if focused:
        paths["focused_wrapper"] = "run_d6_interval_k7_3936177.py"
    for key, name in paths.items():
        require(git_blob_sha256(commit, name) == sources[key], f"interval source blob: {key}")
    require(sha256(ROOT / "ckernel6.dylib") == sources["ckernel6.dylib"], "kernel binary")
    verifier_name = (
        "verify_d6_interval_k7_3936177.py" if focused else "verify_d6_interval_v6.py"
    )
    require(
        sha256(ROOT / verifier_name) == verification["verifier_source_sha256"]
        and git_blob_sha256(commit, verifier_name)
        == verification["verifier_source_sha256"],
        "interval verifier source",
    )
    common_hash = verification["common_helper_sha256"]
    require(
        sha256(ROOT / "verify_d6_interval_benchmarks.py") == common_hash
        and git_blob_sha256(commit, "verify_d6_interval_benchmarks.py")
        == common_hash,
        "interval common helper root",
    )
    return {
        "launch_commit": commit,
        "configuration_sources": dict(sorted(sources.items())),
        "verifier": {
            verifier_name: verification["verifier_source_sha256"],
            "verify_d6_interval_benchmarks.py": common_hash,
        },
    }


def validate_interval_layer(
    report: dict,
    verification: dict,
    *,
    report_name: str,
    input_indices: list[int],
    killed_indices: list[int],
    cap: int,
    focused: bool,
) -> dict:
    require(
        report.get("schema") == 2
        and report.get("total_graphs") == len(input_indices)
        and report.get("completed_graphs") == len(input_indices)
        and report.get("certified_killed") == len(killed_indices)
        and report.get("unresolved_total")
        == len(input_indices) - len(killed_indices)
        and report.get("status_counts", {}).get("KILLED") == len(killed_indices)
        and report.get("status_counts", {}).get("INFRA_ERROR") == 0
        and verification.get("status") == "PASS"
        and verification.get("report", {}).get("sha256")
        == PRIMARY_ARTIFACTS[report_name]
        and verification.get("certified_killed") == len(killed_indices)
        and verification.get("search", {}).get("cap") == cap
        and verification.get("trust_assumptions") == INTERVAL_TRUST,
        "interval report/check boundary",
    )
    results = report.get("results", [])
    require([int(row["index"]) for row in results] == input_indices, "interval result order")
    observed_killed = [int(row["index"]) for row in results if row.get("status") == "KILLED"]
    require(observed_killed == killed_indices, "interval KILLED list")
    require(
        all(row.get("status") in {"KILLED", "ABORT", "UNRESOLVED", "INFRA_ERROR"} for row in results),
        "interval status vocabulary",
    )
    selection = verification["selection"]
    require(
        selection.get("graphs") == len(input_indices)
        and selection.get("indices_sha256") == stable_hash(input_indices)
        and selection.get("v6_manifest_sha256")
        == PRIMARY_ARTIFACTS["d6_current_residue_manifest_v6.json"]
        and selection.get("v6_verification_sha256")
        == PRIMARY_ARTIFACTS["d6_current_residue_manifest_v6_verification.json"],
        "interval selection",
    )
    if focused:
        require(selection.get("index") == killed_indices[0], "focused index")
    else:
        require(verification.get("killed_indices") == killed_indices, "cap KILLED list")
    decisions_name = (
        "d6_interval_k7_3936177_cap2000000_decisions.tsv"
        if focused
        else "d6_interval_v6_k7_cap500000_decisions.tsv"
    )
    require(
        report.get("decisions", {}).get("sha256") == PRIMARY_ARTIFACTS[decisions_name]
        and verification.get("checkpoint", {}).get("decision_TSV_sha256")
        == PRIMARY_ARTIFACTS[decisions_name],
        "interval decisions root",
    )
    require(
        verification.get("kernel_controls", {}).get("status") == "PASS"
        and verification.get("kernel_controls", {}).get("records_replayed") == 48
        and verification.get("witness_replay", {}).get("status") == "PASS"
        and verification.get("witness_replay", {}).get("slices") == 24,
        "interval controls/replay",
    )
    return validate_interval_sources(report, verification, focused=focused)


def validate_current_file_map(mapping: dict, label: str) -> dict[str, str]:
    return verify_named_hashes(mapping, label)


def validate_star(v6_k7: list[int], v6_k7_graphs: list[dict]) -> tuple[list[int], list[int], dict]:
    report = load("d6_k7_one_two_star_increment_report.json")
    verification = load("d6_k7_one_two_star_increment_verification.json")
    summary = report.get("summary", {})
    rejected = list(map(int, summary.get("rejected_indices", ())))
    survivors = list(map(int, summary.get("ordered_survivor_indices", ())))
    require(
        report.get("kind") == "d6_k7_one_free_neighbour_two_free_center_increment"
        and report.get("status") == "COMPLETE"
        and report.get("input", {}).get("ordered_indices") == v6_k7
        and report.get("input", {}).get("ordered_indices_sha256") == stable_hash(v6_k7)
        and report.get("input", {}).get("embedded_graphs_sha256") == stable_hash(v6_k7_graphs)
        and rejected == [2592657, 3785980, 3888410]
        and survivors == [index for index in v6_k7 if index not in set(rejected)]
        and summary.get("rejected_indices_sha256") == stable_hash(rejected)
        and summary.get("ordered_survivor_indices_sha256") == stable_hash(survivors)
        and verification.get("kind") == "d6_k7_one_two_star_increment_verification"
        and verification.get("status") == "PASS"
        and verification.get("report", {}).get("sha256")
        == PRIMARY_ARTIFACTS["d6_k7_one_two_star_increment_report.json"]
        and verification.get("rejected_indices") == rejected
        and verification.get("ordered_survivor_indices") == survivors,
        "exact star partition",
    )
    semantics = report.get("semantics", {})
    require(
        semantics.get("candidate_nonedges_optional") is True
        and semantics.get("floating_point_enters_rejection") is False
        and verification.get("semantics") == semantics
        and report.get("positive_controls", {}).get("known_realizable_18", {}).get("status")
        == "PASS_NOT_APPLICABLE_NO_K7"
        and report.get("positive_controls", {}).get("optional_zero_line_hyperbola", {}).get("status")
        == "PASS_FEASIBLE"
        and verification.get("checked", {}).get("positive_18_control") is True
        and verification.get("checked", {}).get("optional_zero_control") is True,
        "star semantics/controls",
    )
    sources = validate_current_file_map(report["source_sha256"], "star source")
    commit = report["execution"]["git"]["commit"]
    for name, expected in sources.items():
        require(git_blob_sha256(commit, name) == expected, f"star committed source: {name}")
    verifier_hash = verification["verifier_source_sha256"]
    require(
        sha256(ROOT / "verify_d6_k7_one_two_star_increment.py") == verifier_hash
        and git_blob_sha256(commit, "verify_d6_k7_one_two_star_increment.py")
        == verifier_hash,
        "star verifier source",
    )
    return rejected, survivors, {
        "launch_commit": commit,
        "producer_and_kernels": sources,
        "verifier": {"verify_d6_k7_one_two_star_increment.py": verifier_hash},
        "upstream": dict(sorted(report["upstream_artifact_sha256"].items())),
    }


def validate_schur(star_survivors: list[int]) -> tuple[list[int], list[int], dict]:
    report = load("d6_k7_schur_3949382_increment_report.json")
    verification = load("d6_k7_schur_3949382_increment_verification.json")
    rejected = list(map(int, report["summary"]["rejected_indices"]))
    survivors = list(map(int, report["summary"]["ordered_survivor_indices"]))
    require(
        report.get("kind") == "d6_k7_schur_3949382_exact_increment"
        and report.get("status") == "COMPLETE_EXACT_REJECTION"
        and report.get("input", {}).get("ordered_indices") == star_survivors
        and report.get("input", {}).get("ordered_indices_sha256") == stable_hash(star_survivors)
        and rejected == [3949382]
        and survivors == [index for index in star_survivors if index != 3949382]
        and report.get("certificate_sha256") == stable_hash(report.get("certificate"))
        and verification.get("kind") == "d6_k7_schur_3949382_increment_verification"
        and verification.get("status") == "PASS"
        and verification.get("report", {}).get("sha256")
        == PRIMARY_ARTIFACTS["d6_k7_schur_3949382_increment_report.json"]
        and verification.get("conclusion", {}).get("rejected_indices") == rejected
        and verification.get("checked", {}).get("unresolved_optional_entries") == 0,
        "Schur exact increment",
    )
    semantics = report["semantics"]
    require(
        semantics.get("candidate_nonedges_optional") is True
        and semantics.get("floating_point_enters_rejection") is False
        and verification.get("semantics") == semantics
        and report.get("controls", {}).get("nonnegative_boundary_solution", {}).get("rejected_by_strict_domain") is True
        and report.get("controls", {}).get("upstream_independent_replay", {}).get("status") == "PASS",
        "Schur semantics/controls",
    )
    sources = validate_current_file_map(report["source_sha256"], "Schur source")
    commit = report["execution"]["git"]["commit"]
    require(report["execution"]["git"].get("proof_and_checker_sources_equal_committed_blobs") is True, "Schur source claim")
    for name, expected in sources.items():
        require(git_blob_sha256(commit, name) == expected, f"Schur committed source: {name}")
    return rejected, survivors, {
        "launch_commit": commit,
        "package_sources": sources,
        "upstream": dict(sorted(report["upstream_sha256"].items())),
    }


def validate_star_2593240(
    schur_survivors: list[int],
) -> tuple[list[int], list[int], dict]:
    report = load("d6_k7_star_2593240_increment_report.json")
    verification = load("d6_k7_star_2593240_increment_verification.json")
    rejected = list(map(int, report["summary"]["rejected_indices"]))
    survivors = list(map(int, report["summary"]["ordered_survivor_indices"]))
    require(
        report.get("kind") == "d6_k7_rankone_star_2593240_exact_increment"
        and report.get("status") == "COMPLETE_EXACT_REJECTION"
        and report.get("input", {}).get("ordered_indices") == schur_survivors
        and report.get("input", {}).get("ordered_indices_sha256")
        == stable_hash(schur_survivors)
        and rejected == [2593240]
        and survivors == [index for index in schur_survivors if index != 2593240]
        and report.get("certificate_sha256") == stable_hash(report.get("certificate"))
        and verification.get("kind")
        == "d6_k7_star_2593240_increment_verification"
        and verification.get("status") == "PASS"
        and verification.get("report", {}).get("sha256")
        == PRIMARY_ARTIFACTS["d6_k7_star_2593240_increment_report.json"]
        and verification.get("conclusion", {}).get("rejected_indices")
        == rejected
        and verification.get("checked", {}).get("target_index") == 2593240
        and verification.get("checked", {}).get("required_k7_seeds") == 2
        and verification.get("checked", {}).get("new_families") == 2
        and verification.get("checked", {}).get("unresolved_optional_entries")
        == 0,
        "2593240 exact increment",
    )
    semantics = report["semantics"]
    require(
        semantics.get("candidate_nonedges_optional") is True
        and semantics.get("floating_point_enters_rejection") is False
        and semantics.get("both_required_K7_seeds_checked") is True
        and semantics.get("rank_one_schur_support_closure_checked") is True
        and verification.get("semantics") == semantics
        and report.get("controls", {}).get("positive_tetrad_star_control", {}).get(
            "all_thirty_tetrads_zero"
        )
        is True
        and report.get("controls", {}).get("positive_tetrad_star_control", {}).get(
            "rank_one_completion_exists"
        )
        is False
        and report.get("controls", {}).get("upstream_independent_replay", {}).get(
            "preceding_exact_increment_status"
        )
        == "PASS",
        "2593240 semantics/controls",
    )
    sources = validate_current_file_map(report["source_sha256"], "2593240 source")
    commit = report["execution"]["git"]["commit"]
    require(
        commit == "64cde398d8d116632551b6f40edfc5de1b84405b"
        and report["execution"]["git"].get(
            "proof_and_checker_sources_equal_committed_blobs"
        )
        is True
        and verification.get("source_boundary", {}).get("commit") == commit
        and verification.get("source_boundary", {}).get(
            "proof_and_checker_sources_equal_committed_blobs"
        )
        is True,
        "2593240 source boundary",
    )
    for name, expected in sources.items():
        require(
            git_blob_sha256(commit, name) == expected,
            f"2593240 committed source: {name}",
        )
    return rejected, survivors, {
        "launch_commit": commit,
        "package_sources": sources,
        "upstream": dict(sorted(report["upstream_sha256"].items())),
    }


def validate_partition(
    report: dict,
    verification: dict,
    *,
    schema: str,
    verification_schema: str,
    report_name: str,
    input_count: int,
    rejected_count: int,
    survivor_count: int,
) -> tuple[list[int], list[int], list[int]]:
    inputs = list(map(int, report.get("ordered_input_indices", ())))
    rejected = list(map(int, report.get("rejected_indices", ())))
    survivors = list(map(int, report.get("ordered_residue_indices", ())))
    require(
        report.get("schema") == schema
        and report.get("status") == "COMPLETE"
        and report.get("input_graphs") == input_count
        and report.get("graphs_rejected") == rejected_count
        and report.get("graphs_surviving") == survivor_count
        and verification.get("schema") == verification_schema
        and verification.get("status") == "PASS"
        and verification.get("report", {}).get("sha256") == PRIMARY_ARTIFACTS[report_name]
        and verification.get("graphs_recomputed") == input_count
        and verification.get("graphs_rejected") == rejected_count
        and verification.get("graphs_surviving") == survivor_count,
        f"partition boundary: {report_name}",
    )
    require(
        len(inputs) == input_count
        and len(set(inputs)) == input_count
        and len(rejected) == rejected_count
        and len(set(rejected)) == rejected_count
        and survivors == [index for index in inputs if index not in set(rejected)]
        and stable_hash(inputs) == report.get("ordered_input_indices_sha256")
        and stable_hash(rejected) == report.get("rejected_indices_sha256")
        and stable_hash(survivors) == report.get("ordered_residue_indices_sha256")
        and verification.get("rejected_indices") == rejected
        and verification.get("rejected_indices_sha256") == stable_hash(rejected)
        and verification.get("ordered_residue_indices_sha256") == stable_hash(survivors),
        f"partition lists: {report_name}",
    )
    archive = report["certificate_archive"]
    require(
        archive.get("sha256") == PRIMARY_ARTIFACTS[archive["path"]]
        and archive.get("rejected_graphs") == rejected_count,
        f"certificate archive: {report_name}",
    )
    require(
        "exact" in report.get("semantics", {}).get("arithmetic", "")
        and "unconstrained" in report.get("semantics", {}).get("candidate_nonedges", "")
        and report.get("positive_18_control", {}).get("passed") is True
        and verification.get("positive_18_control", {}).get("passed") is True,
        f"exact semantics/control: {report_name}",
    )
    return inputs, rejected, survivors


def validate_k6(v6_k6: list[int]) -> tuple[list[int], list[int], dict]:
    tight = load("d6_k6_tight_same_z0_report.json")
    tight_check = load("d6_k6_tight_same_z0_verification.json")
    opposite = load("d6_k6_opposite_ray_rank_conjunction_report.json")
    opposite_check = load("d6_k6_opposite_ray_rank_conjunction_verification.json")
    tight_input, tight_bad, tight_left = validate_partition(
        tight,
        tight_check,
        schema="d6-k6-tight-same-z0-v1",
        verification_schema="d6-k6-tight-same-z0-verification-v1",
        report_name="d6_k6_tight_same_z0_report.json",
        input_count=625,
        rejected_count=2,
        survivor_count=623,
    )
    opposite_input, opposite_bad, survivors = validate_partition(
        opposite,
        opposite_check,
        schema="d6-k6-opposite-ray-rank-conjunction-v1",
        verification_schema="d6-k6-opposite-ray-rank-conjunction-verification-v1",
        report_name="d6_k6_opposite_ray_rank_conjunction_report.json",
        input_count=623,
        rejected_count=372,
        survivor_count=251,
    )
    require(tight_input == v6_k6 and tight_left == opposite_input, "K6 sequential chain")
    require(not set(tight_bad) & set(opposite_bad), "K6 layer overlap")
    require(stable_hash(survivors) == EXPECTED_K6_SHA256, "K6 final hash")
    tight_dependencies = validate_current_file_map(tight["dependencies"], "tight-Z0 dependencies")
    tight_sources = validate_current_file_map(tight_check["sources"], "tight-Z0 verifier sources")
    require(
        sha256(ROOT / "d6_k6_tight_same_z0.py") == tight["production_source_sha256"],
        "tight-Z0 producer source",
    )
    opposite_dependencies = validate_current_file_map(opposite["dependencies"], "opposite dependencies")
    opposite_sources = validate_current_file_map(opposite["source_sha256"], "opposite source")
    opposite_check_sources = validate_current_file_map(opposite_check["sources"], "opposite verifier sources")
    commit = opposite["execution"]["git"]["commit"]
    require(
        opposite["execution"]["git"].get("tracked_clean") is True
        and opposite_check.get("provenance", {}).get("tracked_sources_clean_at_launch") is True
        and opposite_check.get("provenance", {}).get("commit") == commit,
        "opposite launch provenance",
    )
    for name, expected in opposite_sources.items():
        require(git_blob_sha256(commit, name) == expected, f"opposite committed source: {name}")
    return tight_bad + opposite_bad, survivors, {
        "tight_same_Z0": {
            "dependencies": tight_dependencies,
            "producer_and_verifier_sources": tight_sources,
        },
        "opposite_ray_conjunction": {
            "launch_commit": commit,
            "dependencies": opposite_dependencies,
            "production_sources": opposite_sources,
            "verifier_sources": opposite_check_sources,
        },
    }


def validate_v6() -> tuple[dict, dict, list[int], list[int]]:
    v6 = load("d6_current_residue_manifest_v6.json")
    check = load("d6_current_residue_manifest_v6_verification.json")
    require(
        v6.get("schema") == "d6-current-certified-residue-v6"
        and v6.get("status") == "COMPLETE_MIXED_CERTIFICATE_UNION"
        and v6.get("class_order") == ["K7", "K6_only"]
        and v6.get("combined", {}).get("count") == 644
        and check.get("schema") == "d6-current-certified-residue-v6-verification-v1"
        and check.get("status") == "PASS"
        and check.get("manifest", {}).get("sha256")
        == PRIMARY_ARTIFACTS["d6_current_residue_manifest_v6.json"]
        and check.get("counts")
        == {
            "K7": 19,
            "K6_only": 625,
            "combined": 644,
            "exact_rejections_from_v5": 264,
            "interval_incremental_rejections": 3,
        }
        and all(check.get("checks", {}).values()),
        "v6 verified boundary",
    )
    k7 = list(map(int, v6["classes"]["K7"]["indices"]))
    k6 = list(map(int, v6["classes"]["K6_only"]["indices"]))
    require(
        len(k7) == 19
        and len(k6) == 625
        and stable_hash(k7) == v6["classes"]["K7"]["indices_sha256"]
        and stable_hash(k6) == v6["classes"]["K6_only"]["indices_sha256"]
        and check.get("ordered_residue_indices_sha256") == stable_hash(k7 + k6)
        and v6.get("certificate_accounting", {}).get("interval", {}).get("trust_assumptions")
        == INTERVAL_TRUST,
        "v6 lists/trust",
    )
    verify_named_hashes(v6["source_boundary"]["transitive_files"], "v6 transitive root")
    require(
        sha256(ROOT / "build_d6_current_residue_manifest_v6.py") == check["builder_sha256"],
        "v6 builder source",
    )
    require(
        sha256(ROOT / "verify_d6_current_residue_manifest_v6.py")
        == EXPECTED_V6_VERIFIER_SHA256,
        "v6 verifier source",
    )
    return v6, check, k7, k6


def build_manifest(package_boundary: dict | None = None) -> dict:
    verify_named_hashes(PRIMARY_ARTIFACTS, "primary artifact")
    v6, v6_check, v6_k7, v6_k6 = validate_v6()
    v6_k7_graphs = list(v6["classes"]["K7"]["graphs"])
    v6_k6_graphs = list(v6["classes"]["K6_only"]["graphs"])

    cap_report = load("d6_interval_v6_k7_cap500000_report.json")
    cap_check = load("d6_interval_v6_k7_cap500000_verification.json")
    cap_sources = validate_interval_layer(
        cap_report,
        cap_check,
        report_name="d6_interval_v6_k7_cap500000_report.json",
        input_indices=v6_k7,
        killed_indices=[3595554],
        cap=500_000,
        focused=False,
    )
    star_bad, star_left, star_sources = validate_star(v6_k7, v6_k7_graphs)
    focused_report = load("d6_interval_k7_3936177_cap2000000_report.json")
    focused_check = load("d6_interval_k7_3936177_cap2000000_verification.json")
    focused_sources = validate_interval_layer(
        focused_report,
        focused_check,
        report_name="d6_interval_k7_3936177_cap2000000_report.json",
        input_indices=[3936177],
        killed_indices=[3936177],
        cap=2_000_000,
        focused=True,
    )
    schur_bad, schur_left, schur_sources = validate_schur(star_left)
    star_2593240_bad, _star_2593240_left, star_2593240_sources = (
        validate_star_2593240(schur_left)
    )

    exact_k7_set = set(star_bad + schur_bad + star_2593240_bad)
    interval_k7_set = {3595554, 3936177}
    require(not exact_k7_set & interval_k7_set, "K7 exact/interval overlap")
    all_k7_bad = [index for index in v6_k7 if index in exact_k7_set | interval_k7_set]
    k7_indices = [index for index in v6_k7 if index not in exact_k7_set | interval_k7_set]
    require(k7_indices == EXPECTED_K7 and stable_hash(k7_indices) == EXPECTED_K7_SHA256, "K7 v7 residue")
    k7_graphs = [record for record in v6_k7_graphs if int(record["index"]) in set(k7_indices)]
    require([int(record["index"]) for record in k7_graphs] == k7_indices, "K7 embedded order")

    k6_bad, k6_indices, k6_sources = validate_k6(v6_k6)
    require(len(k6_bad) == 374 and len(set(k6_bad)) == 374, "K6 rejection union")
    k6_graphs = [record for record in v6_k6_graphs if int(record["index"]) in set(k6_indices)]
    require([int(record["index"]) for record in k6_graphs] == k6_indices, "K6 embedded order")

    combined_indices = k7_indices + k6_indices
    require(
        len(combined_indices) == 263
        and stable_hash(combined_indices) == EXPECTED_COMBINED_SHA256
        and stable_hash(sorted(combined_indices)) == EXPECTED_SORTED_SHA256,
        "combined v7 residue",
    )
    if package_boundary is None:
        package_boundary = committed_package_boundary()
    require(
        set(package_boundary.get("package_sources", {})) == set(PACKAGE_FILES),
        "v7 package source set",
    )
    tagged_indices = [
        {"class": class_name, "index": int(record["index"])}
        for class_name, records in (("K7", k7_graphs), ("K6_only", k6_graphs))
        for record in records
    ]
    tagged_graphs = [
        {"class": class_name, **record}
        for class_name, records in (("K7", k7_graphs), ("K6_only", k6_graphs))
        for record in records
    ]
    exact_k7 = [index for index in v6_k7 if index in exact_k7_set]
    interval_k7 = [index for index in v6_k7 if index in interval_k7_set]
    source_roots = {
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
        "exact_2593240": star_2593240_sources,
        "exact_K6_chain": k6_sources,
    }
    return {
        "schema": "d6-current-certified-residue-v7",
        "status": "COMPLETE_MIXED_CERTIFICATE_UNION",
        "description": (
            "The independently verified v6 residue after every subsequent "
            "certified K7 and K6 layer. Survival is not realizability."
        ),
        "class_order": ["K7", "K6_only"],
        "classes": {
            "K7": {
                "v6_count": 19,
                "layers": [
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
                        "rejected_indices": star_2593240_bad,
                        "rejected_indices_sha256": stable_hash(
                            star_2593240_bad
                        ),
                        "residue_count": len(_star_2593240_left),
                        "residue_indices_sha256": stable_hash(
                            _star_2593240_left
                        ),
                    },
                ],
                "post_v6_exact_rejection_count": 5,
                "post_v6_exact_rejected_indices": exact_k7,
                "post_v6_exact_rejected_indices_sha256": stable_hash(exact_k7),
                "post_v6_interval_rejection_count": 2,
                "post_v6_interval_rejected_indices": interval_k7,
                "post_v6_interval_rejected_indices_sha256": stable_hash(interval_k7),
                "post_v6_union_rejection_count": 7,
                "post_v6_union_rejected_indices": all_k7_bad,
                "post_v6_union_rejected_indices_sha256": stable_hash(all_k7_bad),
                "count": 12,
                "indices": k7_indices,
                "indices_sha256": EXPECTED_K7_SHA256,
                "graphs": k7_graphs,
                "graphs_sha256": stable_hash(k7_graphs),
            },
            "K6_only": {
                "v6_count": 625,
                "exact_layers": [
                    {
                        "name": "tight_same_Z0",
                        "input_count": 625,
                        "input_indices_sha256": stable_hash(v6_k6),
                        "rejection_count": 2,
                        "rejected_indices_sha256": stable_hash(k6_bad[:2]),
                        "residue_count": 623,
                        "residue_indices_sha256": stable_hash(
                            [index for index in v6_k6 if index not in set(k6_bad[:2])]
                        ),
                    },
                    {
                        "name": "opposite_ray_rank_conjunction",
                        "input_count": 623,
                        "input_indices_sha256": stable_hash(
                            [index for index in v6_k6 if index not in set(k6_bad[:2])]
                        ),
                        "rejection_count": 372,
                        "rejected_indices_sha256": stable_hash(k6_bad[2:]),
                        "residue_count": 251,
                        "residue_indices_sha256": EXPECTED_K6_SHA256,
                    },
                ],
                "post_v6_exact_rejection_count": 374,
                "post_v6_exact_rejected_indices": [index for index in v6_k6 if index in set(k6_bad)],
                "post_v6_exact_rejected_indices_sha256": stable_hash([index for index in v6_k6 if index in set(k6_bad)]),
                "count": 251,
                "indices": k6_indices,
                "indices_sha256": EXPECTED_K6_SHA256,
                "graphs": k6_graphs,
                "graphs_sha256": stable_hash(k6_graphs),
            },
        },
        "combined": {
            "v6_count": 644,
            "post_v6_exact_rejection_count": 379,
            "post_v6_interval_incremental_rejection_count": 2,
            "post_v6_union_rejection_count": 381,
            "count": 263,
            "class_counts": {"K7": 12, "K6_only": 251},
            "ordered_indices_sha256": EXPECTED_COMBINED_SHA256,
            "sorted_indices_sha256": EXPECTED_SORTED_SHA256,
            "tagged_indices_sha256": stable_hash(tagged_indices),
            "tagged_graphs_sha256": stable_hash(tagged_graphs),
            "cross_class_overlap": 0,
        },
        "certificate_accounting": {
            "inherited_v6": {"exact_rejections_from_v5": 264, "interval_incremental_rejections": 3},
            "exact": {
                "post_v6_rejections": 379,
                "total_rejections_from_v5": 643,
                "floating_point_enters_rejection": False,
                "layers": {"K7_one_two_star": 3, "K7_Schur": 1, "K7_rankone_star_2593240": 1, "K6_tight_same_Z0": 2, "K6_opposite_ray_conjunction": 372},
            },
            "interval": {
                "post_v6_incremental_rejections": 2,
                "total_incremental_rejections_from_v5": 5,
                "only_KILLED_used_for_rejection": True,
                "ABORT_UNRESOLVED_INFRA_ERROR_used_for_rejection": False,
                "trust_assumptions": INTERVAL_TRUST,
                "layers": {"v6_cap100000_inherited": 3, "v6_K7_cap500000": 1, "focused_3936177_cap2000000": 1},
            },
            "union_rejections_from_v5": 648,
            "v5_graphs": 911,
            "current_residue": 263,
        },
        "semantics": {
            "edges": "required unit distances",
            "nonedges": "unconstrained and may also have distance one",
            "points": "distinct",
            "exact_rejection": "exact integer, rational, algebraic, polynomial, and graph-logic certificates",
            "interval_rejection": "verified KILLED only, conditional on the recorded IEEE-754 and macOS libm assumptions",
            "residue": "certificate-union non-rejection only; not a realizability claim",
        },
        "controls": {
            "known_realizable_18": {
                "inherited_v6": {
                    "passed": True,
                    "adjacency_sha256": v6["positive_18_control"]["v5"][
                        "adjacency_sha256"
                    ],
                },
                "exact_star": {"passed": True},
                "exact_K6_tight_same_Z0": load(
                    "d6_k6_tight_same_z0_verification.json"
                )["positive_18_control"],
                "exact_K6_opposite_ray_conjunction": load(
                    "d6_k6_opposite_ray_rank_conjunction_verification.json"
                )["positive_18_control"],
            },
            "optional_zero_star": {"passed": True},
            "interval_cap500000_kernel": cap_check["kernel_controls"],
            "interval_focused_3936177_kernel": focused_check["kernel_controls"],
            "exact_Schur_nonnegative_boundary": {
                "checked": True,
                "strict_domain_required": True,
            },
            "exact_2593240_tetrad_star": {
                "two_required_K7_seeds_checked": True,
                "all_thirty_tetrads_zero_at_control": True,
                "rank_one_completion_exists": False,
            },
        },
        "source_boundary": {
            "primary_artifacts": dict(sorted(PRIMARY_ARTIFACTS.items())),
            "layer_source_roots": source_roots,
            "package": package_boundary,
        },
        "execution": {
            "command": " ".join([sys.executable, *sys.argv]),
            "generated_utc": datetime.now(UTC).isoformat(),
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python_version": sys.version,
            "logical_cpus": os.cpu_count(),
        },
        "nonclaims": [
            "A surviving graph is not asserted realizable.",
            "This residue does not by itself prove f(6)=18.",
            "Interval rejections retain their stated platform arithmetic assumptions.",
        ],
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
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    manifest = build_manifest()
    atomic_json(args.output.resolve(), manifest)
    print(
        json.dumps(
            {
                "output": str(args.output.resolve()),
                "sha256": sha256(args.output.resolve()),
                "status": manifest["status"],
                "K7": manifest["classes"]["K7"]["count"],
                "K6_only": manifest["classes"]["K6_only"]["count"],
                "combined": manifest["combined"]["count"],
                "source_commit": manifest["source_boundary"]["package"]["commit"],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
