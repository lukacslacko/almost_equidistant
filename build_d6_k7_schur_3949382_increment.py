#!/usr/bin/env python3
"""Build the source-bound exact Schur increment for K7 graph 3949382.

The mathematical certificate is constructed by the already frozen exploratory
probe.  This production wrapper pins that probe, the independently verified
upstream boundary, and every file in the production checking package to a Git
commit.  The separate verifier deliberately does not import this builder or
the probe.
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

import probe_d6_k7_schur_3949382 as probe


ROOT = Path(__file__).resolve().parent
TARGET_INDEX = 3_949_382
REPORT = ROOT / "d6_k7_schur_3949382_increment_report.json"

UPSTREAM = {
    "d6_current_residue_manifest_v6.json": (
        "c513bcc40c037af23340dce3cc91d593a101dd1d073f93feab82abd556efadc9"
    ),
    "d6_current_residue_manifest_v6_verification.json": (
        "3d11cbd20cb59b42a8c8fd5252fb3ef3dea649576982ea8017622fad3279a7bc"
    ),
    "d6_k7_one_two_star_increment_report.json": (
        "12b3d18b1ea81961f58d831d4c7c322fbceb2300e7533b64161e8011fbe9d1ec"
    ),
    "d6_k7_one_two_star_increment_verification.json": (
        "c6f40578685036cb1156cb0ea8bf06d6004a70d5f376aff5ac1e16a1c0ca09eb"
    ),
}
FROZEN_EXPLORATORY = {
    "probe_d6_k7_schur_3949382.py": (
        "96750e6591aae18dd25e59d382ed64c316575431942780b002be12da4594cb86"
    ),
    "test_probe_d6_k7_schur_3949382.py": (
        "ab3842c4f55e1fa6216ad10e18a197b99dc4c470e9c0e00cba7c26e6eab911ab"
    ),
    "d6_k7_schur_3949382.md": (
        "660de8bb17947dc7848886b51e11327f0cc4645e55dc993ae5df423a8de0aaf4"
    ),
}
PACKAGE_FILES = (
    "build_d6_k7_schur_3949382_increment.py",
    "verify_d6_k7_schur_3949382_increment.py",
    "test_d6_k7_schur_3949382_increment.py",
    "d6_k7_schur_3949382_increment.md",
)
SOURCE_FILES = tuple(FROZEN_EXPLORATORY) + PACKAGE_FILES
EXPECTED_INPUT_INDICES = (
    316173,
    2581209,
    2593240,
    3595554,
    3648882,
    3729907,
    3935560,
    3936177,
    3936310,
    3936435,
    3945490,
    3945555,
    3945557,
    3945564,
    3947605,
    3949382,
)
EXPECTED_INPUT_SHA256 = (
    "b0235ce7e95918197f0f3b4a57d26c8c83632ff51fd678a989f50197591eee8e"
)
EXPECTED_OUTPUT_SHA256 = (
    "715c52011421fb6d8721d52339ed3437ef1cec344ddf3f03a3b0b618f35b00d7"
)
EXPECTED_PROBE_REPORT_STABLE_SHA256 = (
    "7394a569982ae9453acac094a0de3e0818f8c31fbc01b8e21706d14194444217"
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


def atomic_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def git(*arguments: str) -> str:
    return subprocess.run(
        ["git", *arguments],
        cwd=ROOT,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    ).stdout.strip()


def committed_source_boundary() -> dict:
    """Require every proof/checking source to equal its committed blob."""

    commit = git("rev-parse", "HEAD")
    branch = git("branch", "--show-current")
    require(branch == "codex/dimension6", "official launch must use codex/dimension6")
    source_hashes = {name: sha256(ROOT / name) for name in SOURCE_FILES}
    for name, expected in FROZEN_EXPLORATORY.items():
        require(source_hashes[name] == expected, f"frozen source hash: {name}")
    for name, expected in source_hashes.items():
        blob = subprocess.run(
            ["git", "show", f"{commit}:{name}"],
            cwd=ROOT,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        ).stdout
        require(hashlib.sha256(blob).hexdigest() == expected, f"uncommitted source: {name}")
    porcelain = git("status", "--porcelain=v1", "--untracked-files=all")
    return {
        "commit": commit,
        "branch": branch,
        "source_sha256": source_hashes,
        "porcelain_lines": porcelain.splitlines(),
        "porcelain_sha256": hashlib.sha256(porcelain.encode("utf-8")).hexdigest(),
        "proof_and_checker_sources_equal_committed_blobs": True,
    }


def load_and_validate_upstream() -> tuple[dict, dict]:
    for name, expected in UPSTREAM.items():
        require(sha256(ROOT / name) == expected, f"upstream hash: {name}")
    manifest = json.loads(
        (ROOT / "d6_current_residue_manifest_v6.json").read_text(encoding="utf-8")
    )
    manifest_check = json.loads(
        (ROOT / "d6_current_residue_manifest_v6_verification.json").read_text(
            encoding="utf-8"
        )
    )
    star = json.loads(
        (ROOT / "d6_k7_one_two_star_increment_report.json").read_text(
            encoding="utf-8"
        )
    )
    star_check = json.loads(
        (ROOT / "d6_k7_one_two_star_increment_verification.json").read_text(
            encoding="utf-8"
        )
    )
    require(
        manifest.get("schema") == "d6-current-certified-residue-v6"
        and manifest.get("status") == "COMPLETE_MIXED_CERTIFICATE_UNION",
        "manifest semantic boundary",
    )
    require(
        manifest_check.get("status") == "PASS"
        and manifest_check.get("manifest", {}).get("sha256")
        == UPSTREAM["d6_current_residue_manifest_v6.json"],
        "manifest verification boundary",
    )
    require(
        star.get("kind")
        == "d6_k7_one_free_neighbour_two_free_center_increment"
        and star.get("status") == "COMPLETE"
        and star.get("semantics") == probe.EXPECTED_SEMANTICS,
        "star semantic boundary",
    )
    require(
        star_check.get("kind") == "d6_k7_one_two_star_increment_verification"
        and star_check.get("status") == "PASS"
        and star_check.get("report", {}).get("sha256")
        == UPSTREAM["d6_k7_one_two_star_increment_report.json"]
        and star_check.get("semantics") == probe.EXPECTED_SEMANTICS
        and star_check.get("checked", {}).get("eligible_covers") == 19_932
        and star_check.get("checked", {}).get("current_passing_families") == 88,
        "star independent-verification boundary",
    )
    indices = tuple(map(int, star["summary"]["ordered_survivor_indices"]))
    require(indices == EXPECTED_INPUT_INDICES, "star survivor list")
    require(stable_hash(list(indices)) == EXPECTED_INPUT_SHA256, "star survivor hash")
    require(
        star["summary"].get("ordered_survivor_indices_sha256")
        == EXPECTED_INPUT_SHA256,
        "archived star survivor hash",
    )
    return manifest, star


def validated_exploratory_certificate() -> dict:
    """Rebuild and pin the frozen exploratory certificate exactly."""

    for name, expected in FROZEN_EXPLORATORY.items():
        require(sha256(ROOT / name) == expected, f"frozen exploratory file: {name}")
    certificate = probe.build_report()
    require(
        stable_hash(certificate) == EXPECTED_PROBE_REPORT_STABLE_SHA256,
        "probe report stable hash",
    )
    require(certificate.get("status") == "EXACT_CONTRADICTION", "probe status")
    require(certificate.get("target", {}).get("index") == TARGET_INDEX, "probe target")
    require(certificate.get("quantifier", {}).get("raw_eligible_covers") == 128, "cover count")
    require(
        certificate.get("quantifier", {}).get("current_covers") == [0, 2048, 3072],
        "current cover boundary",
    )
    require(
        certificate.get("quantifier", {}).get("sole_new_branch_passing_families") == 1,
        "sole-family boundary",
    )
    require(
        certificate.get("exact_pattern", {}).get("pair_target_string")
        == "0111011011",
        "exact Schur pattern",
    )
    require(
        certificate.get("semantics", {}).get("candidate_nonedges_optional") is True
        and certificate.get("semantics", {}).get("floating_point_enters_contradiction")
        is False,
        "probe semantics",
    )
    return certificate


def build_report(source_boundary: dict | None = None) -> dict:
    """Construct the production report; tests may inject a synthetic boundary."""

    _manifest, star = load_and_validate_upstream()
    probe_report = validated_exploratory_certificate()
    if source_boundary is None:
        source_boundary = committed_source_boundary()
    source_hashes = source_boundary.get("source_sha256", {})
    require(set(source_hashes) == set(SOURCE_FILES), "production source set")

    input_indices = list(EXPECTED_INPUT_INDICES)
    output_indices = [index for index in input_indices if index != TARGET_INDEX]
    require(len(output_indices) == 15, "output cardinality")
    require(stable_hash(output_indices) == EXPECTED_OUTPUT_SHA256, "output hash")
    target_record = next(
        record for record in star["records"] if int(record["index"]) == TARGET_INDEX
    )
    require(target_record.get("decision") == "SURVIVOR", "target upstream decision")

    certificate = {
        "target": probe_report["target"],
        "quantifier": probe_report["quantifier"],
        "exact_pattern": probe_report["exact_pattern"],
        "algebra": probe_report["algebra"],
        "rank_argument": {
            "normalized_gram_rank_upper_bound": 7,
            "required_clique_order": 7,
            "clique_principal_block_positive_definite": True,
            "schur_complement_for_five_remainder_vertices": "identically zero",
            "clique_inverse_method": "Sherman--Morrison over exact rationals",
            "strict_domain": "a,b,c,d,e,f,g > 0",
        },
    }
    return {
        "schema": 1,
        "kind": "d6_k7_schur_3949382_exact_increment",
        "status": "COMPLETE_EXACT_REJECTION",
        "claim": (
            "Graph 3949382 is not the unit-edge graph of distinct points in R^6. "
            "A required K7 seed has an exhaustive prior support boundary, and its "
            "sole surviving family contradicts exact rank-seven Schur equations "
            "in the strict positive orthant."
        ),
        "semantics": {
            "candidate_nonedges_optional": True,
            "zero_K_entries_require_disjoint_propagated_support_supersets": True,
            "required_edges_give_unit_K_entries": True,
            "propagated_masks_are_support_supersets": True,
            "one_infeasible_required_K7_seed_rejects_graph": True,
            "all_actual_support_branches_for_seed_are_quantified": True,
            "strict_positive_diagonal_parameters_encode_distinct_points": True,
            "floating_point_enters_rejection": False,
        },
        "upstream_sha256": dict(sorted(UPSTREAM.items())),
        "frozen_exploratory_sha256": dict(sorted(FROZEN_EXPLORATORY.items())),
        "source_sha256": dict(sorted(source_hashes.items())),
        "input": {
            "artifact": "d6_k7_one_two_star_increment_report.json",
            "graphs": len(input_indices),
            "ordered_indices": input_indices,
            "ordered_indices_sha256": EXPECTED_INPUT_SHA256,
        },
        "summary": {
            "graphs_rejected": 1,
            "rejected_indices": [TARGET_INDEX],
            "rejected_indices_sha256": stable_hash([TARGET_INDEX]),
            "graphs_surviving": len(output_indices),
            "ordered_survivor_indices": output_indices,
            "ordered_survivor_indices_sha256": EXPECTED_OUTPUT_SHA256,
        },
        "certificate": certificate,
        "certificate_sha256": stable_hash(certificate),
        "controls": {
            "nonnegative_boundary_solution": {
                "point": [1, 3, 1, 1, 2, 0, 0],
                "all_ten_equations_zero": True,
                "rejected_by_strict_domain": True,
            },
            "upstream_independent_replay": {
                "eligible_covers": 19_932,
                "current_passing_families": 88,
                "status": "PASS",
            },
        },
        "execution": {
            "command": " ".join([sys.executable, *sys.argv]),
            "started_utc": datetime.now(UTC).isoformat(),
            "finished_utc": datetime.now(UTC).isoformat(),
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python_version": sys.version,
            "logical_cpus": os.cpu_count(),
            "git": source_boundary,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=REPORT)
    args = parser.parse_args()
    report = build_report()
    atomic_json(args.output.resolve(), report)
    print(
        json.dumps(
            {
                "output": str(args.output.resolve()),
                "sha256": sha256(args.output.resolve()),
                "status": report["status"],
                "rejected_indices": report["summary"]["rejected_indices"],
                "survivors": report["summary"]["graphs_surviving"],
                "source_commit": report["execution"]["git"]["commit"],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
