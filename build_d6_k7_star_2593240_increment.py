#!/usr/bin/env python3
"""Build the source-bound exact rank-one-star increment for graph 2593240.

The mathematical certificate is constructed by the frozen exploratory probe.
This production wrapper pins the probe, the independently checked support
boundary, the preceding graph-3949382 exact increment, and every file in this
checking package to a Git commit.  The separate verifier imports neither this
builder nor the probe.
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

import probe_d6_k7_2593240_star_boundary as probe


ROOT = Path(__file__).resolve().parent
TARGET_INDEX = 2_593_240
REPORT = ROOT / "d6_k7_star_2593240_increment_report.json"

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
    "d6_k7_schur_3949382_increment_report.json": (
        "30efa40d1db8a3aef03af9e56aa915a04594aabe8b72e5a1c72489e2a8c952b8"
    ),
    "d6_k7_schur_3949382_increment_verification.json": (
        "18c9980baa74596523eba9b447af1ae5830a9558d80981e1fa3ff6c984bcb5bc"
    ),
}
FROZEN_EXPLORATORY = {
    "probe_d6_k7_2593240_star_boundary.py": (
        "7d88aafdda609a85255476cb81d72ee0685880c8b17e0823f3ffc4c690ad912f"
    ),
    "test_probe_d6_k7_2593240_star_boundary.py": (
        "30e5949baa4180edc6578d25054fd175c0bd81b248470158bb39e8bcc7b6c2a1"
    ),
    "d6_k7_2593240_star_boundary.md": (
        "cbaa5e0bcb1f313136fff96c17cb33b251465232f062df7efa1fa692c8b93d2e"
    ),
}
PACKAGE_FILES = (
    "build_d6_k7_star_2593240_increment.py",
    "verify_d6_k7_star_2593240_increment.py",
    "test_d6_k7_star_2593240_increment.py",
    "d6_k7_star_2593240_increment.md",
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
)
EXPECTED_INPUT_SHA256 = (
    "715c52011421fb6d8721d52339ed3437ef1cec344ddf3f03a3b0b618f35b00d7"
)
EXPECTED_OUTPUT_SHA256 = (
    "3f038f15c97483aaf529300a86da00459d04e4a8f49823b8849fa7a4a91dd7cb"
)
EXPECTED_PROBE_REPORT_STABLE_SHA256 = (
    "0afebed8f02bb5f7db055380b1deefc78b00eaa14f2c36a5482f268eb603d4a4"
)
REPORT_SEMANTICS = {
    "all_actual_support_branches_for_both_seeds_are_quantified": True,
    "both_required_K7_seeds_checked": True,
    "candidate_nonedges_optional": True,
    "floating_point_enters_rejection": False,
    "one_infeasible_required_K7_seed_rejects_graph": True,
    "propagated_masks_are_support_supersets": True,
    "rank_one_schur_support_closure_checked": True,
    "required_edges_give_unit_K_entries": True,
    "strict_positive_diagonal_parameters_encode_distinct_points": True,
    "zero_K_entries_require_disjoint_propagated_support_supersets": True,
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


def load_and_validate_upstream() -> tuple[dict, dict, dict]:
    """Validate the base support boundary and preceding exact increment."""

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
    preceding = json.loads(
        (ROOT / "d6_k7_schur_3949382_increment_report.json").read_text(
            encoding="utf-8"
        )
    )
    preceding_check = json.loads(
        (ROOT / "d6_k7_schur_3949382_increment_verification.json").read_text(
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
    require(
        preceding.get("kind") == "d6_k7_schur_3949382_exact_increment"
        and preceding.get("status") == "COMPLETE_EXACT_REJECTION"
        and preceding.get("summary", {}).get("ordered_survivor_indices")
        == list(EXPECTED_INPUT_INDICES)
        and preceding.get("summary", {}).get("ordered_survivor_indices_sha256")
        == EXPECTED_INPUT_SHA256,
        "preceding exact increment boundary",
    )
    require(
        preceding_check.get("kind")
        == "d6_k7_schur_3949382_increment_verification"
        and preceding_check.get("status") == "PASS"
        and preceding_check.get("report", {}).get("sha256")
        == UPSTREAM["d6_k7_schur_3949382_increment_report.json"]
        and preceding_check.get("conclusion", {}).get(
            "ordered_survivor_indices_sha256"
        )
        == EXPECTED_INPUT_SHA256,
        "preceding independent verification boundary",
    )
    return manifest, star, preceding


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
    require(
        [row.get("raw_eligible_covers") for row in certificate.get("quantifier", [])]
        == [502, 502],
        "probe raw-cover boundary",
    )
    require(
        [row.get("current_covers") for row in certificate.get("quantifier", [])]
        == [[0, 8], [0, 32]],
        "probe current-cover boundary",
    )
    require(
        [row.get("pair_target_string") for row in certificate.get("exact_patterns", [])]
        == ["011111111110111", "111110111111101"],
        "probe exact Schur patterns",
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

    _manifest, star, preceding = load_and_validate_upstream()
    probe_report = validated_exploratory_certificate()
    if source_boundary is None:
        source_boundary = committed_source_boundary()
    source_hashes = source_boundary.get("source_sha256", {})
    require(set(source_hashes) == set(SOURCE_FILES), "production source set")

    input_indices = list(EXPECTED_INPUT_INDICES)
    output_indices = [index for index in input_indices if index != TARGET_INDEX]
    require(len(output_indices) == 14, "output cardinality")
    require(stable_hash(output_indices) == EXPECTED_OUTPUT_SHA256, "output hash")
    require(
        preceding["summary"]["ordered_survivor_indices"] == input_indices,
        "preceding input order",
    )
    target_record = next(
        record for record in star["records"] if int(record["index"]) == TARGET_INDEX
    )
    require(target_record.get("decision") == "SURVIVOR", "target upstream decision")

    certificate = {
        "target": probe_report["target"],
        "quantifier": probe_report["quantifier"],
        "exact_patterns": probe_report["exact_patterns"],
        "system_isomorphism": probe_report["system_isomorphism"],
        "algebra": probe_report["algebra"],
        "rank_argument": probe_report["rank_argument"],
    }
    return {
        "schema": 1,
        "kind": "d6_k7_rankone_star_2593240_exact_increment",
        "status": "COMPLETE_EXACT_REJECTION",
        "claim": (
            "Graph 2593240 is not the unit-edge graph of distinct points in R^6. "
            "Both required K7 seeds have exhaustive prior support boundaries; "
            "their isomorphic sole surviving families force a nonzero-star "
            "off-diagonal pattern, impossible for a rank-one PSD Schur complement."
        ),
        "semantics": dict(REPORT_SEMANTICS),
        "upstream_sha256": dict(sorted(UPSTREAM.items())),
        "frozen_exploratory_sha256": dict(sorted(FROZEN_EXPLORATORY.items())),
        "source_sha256": dict(sorted(source_hashes.items())),
        "input": {
            "artifact": "d6_k7_schur_3949382_increment_report.json",
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
            "positive_tetrad_star_control": {
                "point": ["2/3", "5", "1", "1", "1", "1"],
                "all_thirty_tetrads_zero": True,
                "rank_one_completion_exists": False,
                "failure": "g_01 and g_02 are nonzero while g_12 is zero",
            },
            "upstream_independent_replay": {
                "eligible_covers": 19_932,
                "current_passing_families": 88,
                "preceding_exact_increment_status": "PASS",
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
