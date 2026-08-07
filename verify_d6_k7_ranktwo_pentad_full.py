#!/usr/bin/env python3
"""Memory-bounded independent verifier for the full rank-two pentad scan."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import d6_k7_rank_reference as prior
import verify_d6_k7_ranktwo_pentad as checker
from verify_profile_d6 import lower_bound_18_graph


ROOT = Path(__file__).resolve().parent


def resolve(base: Path, value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else base / path


def atomic_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    temporary.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    os.replace(temporary, path)


def cover_subprocess(payload: tuple[int, dict, dict]) -> dict:
    ordinal, target, record = payload
    encoded = json.dumps(
        [ordinal, target, record], sort_keys=True, separators=(",", ":")
    )
    environment = os.environ.copy()
    for name in (
        "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
    ):
        environment[name] = "1"
    result = subprocess.run(
        [sys.executable, str(Path(__file__).resolve()), "--worker-payload", encoded],
        cwd=ROOT, env=environment,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
    )
    if result.returncode:
        raise RuntimeError(
            f"verification worker cover {ordinal}, graph {target['graph_index']} "
            f"exit {result.returncode}: {result.stderr.strip()}"
        )
    try:
        return json.loads(result.stdout)
    except json.JSONDecodeError as error:
        raise RuntimeError(
            f"verification worker {ordinal} invalid JSON; "
            f"stdout={result.stdout[:500]!r}, stderr={result.stderr[:500]!r}"
        ) from error


def target_gate_subprocess(report_path: Path) -> dict:
    result = subprocess.run(
        [
            sys.executable, str(Path(__file__).resolve()),
            "--target-gate", str(report_path),
        ],
        cwd=ROOT,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
    )
    if result.returncode:
        raise RuntimeError(
            f"independent target gate exit {result.returncode}: "
            f"{result.stderr.strip()}"
        )
    return json.loads(result.stdout)


def run_target_gate(report_path: Path) -> dict:
    report = json.loads(report_path.read_text(encoding="utf-8"))
    targets, profiles = checker.reconstruct_full_targets(report_path, report)
    manifest_record = report["provenance"]["target_manifest"]
    manifest_path = resolve(report_path.parent, manifest_record["path"])
    if checker.sha256(manifest_path) != manifest_record["sha256"]:
        raise ValueError("target manifest hash mismatch at independent gate")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest["targets"] != targets:
        raise ValueError("target manifest differs from independent exact replay")
    return {
        "status": "PASS",
        "targets": len(targets),
        "graphs": len({item["graph_index"] for item in targets}),
        "targets_sha256": checker.stable_hash(targets),
        "profiles_sha256": checker.stable_hash(profiles),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--report", type=Path,
        default=Path("d6_k7_ranktwo_pentad_full_report.json"),
    )
    parser.add_argument("--report-sha256")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--batch-size", type=int, default=16)
    parser.add_argument(
        "--checkpoint", type=Path,
        default=Path(".runs/d6_k7_ranktwo_pentad_verification_checkpoint.json"),
    )
    parser.add_argument("--resume", action="store_true")
    parser.add_argument(
        "--output", type=Path,
        default=Path("d6_k7_ranktwo_pentad_full_verification.json"),
    )
    parser.add_argument("--worker-payload", help=argparse.SUPPRESS)
    parser.add_argument("--target-gate", type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()

    if args.worker_payload is not None:
        ordinal, target, record = json.loads(args.worker_payload)
        print(json.dumps(
            checker.verify_full_cover_record(target, record, int(ordinal)),
            sort_keys=True, separators=(",", ":"),
        ))
        return
    if args.target_gate is not None:
        print(json.dumps(
            run_target_gate(args.target_gate),
            sort_keys=True, separators=(",", ":"),
        ))
        return
    if not args.report_sha256:
        parser.error("--report-sha256 is required")
    if args.workers <= 0 or args.batch_size <= 0:
        parser.error("workers and batch size must be positive")

    observed_hash = checker.sha256(args.report)
    if observed_hash != args.report_sha256:
        raise ValueError("full pentad report hash differs from explicit pin")
    report = json.loads(args.report.read_text(encoding="utf-8"))
    if (
        report.get("schema") != 1
        or report.get("kind")
        != "d6_k7_ranktwo_pentad_full_coefficientwise_scan"
        or report.get("status") != "COMPLETE"
    ):
        raise ValueError("unexpected full report schema/kind/status")
    for name, expected_hash in report["source_hashes"].items():
        if checker.sha256(ROOT / name) != expected_hash:
            raise ValueError(f"full report source hash mismatch: {name}")

    symbolic = {
        "elimination_derivation": checker.verify_elimination_derivation(),
        "generic_rank_two_gram_substitution": (
            checker.verify_generic_rank_two_gram_identity()
        ),
        "permutation_orbit": checker.verify_pentad_permutation_orbit(),
    }
    target_gate = target_gate_subprocess(args.report)
    if target_gate.get("status") != "PASS":
        raise ValueError("independent target gate did not pass")

    manifest_path = resolve(
        args.report.parent, report["provenance"]["target_manifest"]["path"]
    )
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    targets = manifest["targets"]
    records = report["covers"]
    if len(targets) != 53 or len(records) != 53:
        raise ValueError("full verification population is not 53 covers")

    configuration_hash = checker.stable_hash({
        "report_sha256": observed_hash,
        "checker_sha256": checker.sha256(ROOT / "verify_d6_k7_ranktwo_pentad.py"),
        "workers": args.workers,
        "batch_size": args.batch_size,
    })
    results: list[dict] = []
    if args.resume and args.checkpoint.exists():
        checkpoint = json.loads(args.checkpoint.read_text(encoding="utf-8"))
        if checkpoint.get("configuration_sha256") != configuration_hash:
            raise ValueError("verification checkpoint configuration mismatch")
        results = checkpoint["results"]
        if [item["ordinal"] for item in results] != list(range(len(results))):
            raise ValueError("verification checkpoint is not an ordinal prefix")
    jobs = list(zip(
        range(len(results), len(targets)),
        targets[len(results):], records[len(results):],
    ))
    for start in range(0, len(jobs), args.batch_size):
        batch = jobs[start:start + args.batch_size]
        with ThreadPoolExecutor(
            max_workers=min(args.workers, len(batch))
        ) as pool:
            results.extend(pool.map(cover_subprocess, batch))
        atomic_json(args.checkpoint, {
            "schema": 1,
            "kind": "d6_k7_ranktwo_pentad_verification_checkpoint",
            "configuration_sha256": configuration_hash,
            "completed": len(results),
            "results": results,
        })

    status_counts = Counter(item["status"] for item in results)
    classifications = Counter()
    rejected_by_graph: Counter[int] = Counter()
    cover_counts: Counter[int] = Counter()
    for item in results:
        classifications.update(item["classification_counts"])
        graph_index = int(item["graph_index"])
        cover_counts[graph_index] += 1
        if item["status"] == "REJECTED":
            rejected_by_graph[graph_index] += 1
    union_path = resolve(
        args.report.parent, report["provenance"]["union"]["path"]
    )
    union = json.loads(union_path.read_text(encoding="utf-8"))
    profiles = {
        int(item["index"]): item
        for item in union["cover_structure"]["residue_profiles"]
    }
    all_rejected = sorted(
        index for index, count in cover_counts.items()
        if rejected_by_graph[index] == count
    )
    marginal = sorted(
        index for index in all_rejected
        if not profiles[index]["has_saturating_cover"]
        and not profiles[index]["has_tetrad_resistant_near_clique_cover"]
    )
    exact_summary = {
        "cover_status_counts": dict(status_counts),
        "pentads_tested": sum(item["pentads_checked"] for item in results),
        "pentad_classification_counts": dict(classifications),
        "covers_with_multiple_one_sign_pentads": sum(
            item["one_sign_count"] > 1 for item in results
        ),
        "graphs_with_all_no_near_covers_rejected": len(all_rejected),
        "all_no_near_covers_rejected_indices": all_rejected,
        "marginal_graph_rejections": len(marginal),
        "marginal_graph_rejected_indices": marginal,
    }
    for name, value in exact_summary.items():
        if report["summary"][name] != value:
            raise ValueError(f"full report summary mismatch in {name}")

    positive = tuple(int(row) for row in lower_bound_18_graph())
    prior.validate_graph(positive)
    positive_k7 = sum(1 for _ in prior.clique_masks(positive, 7))
    if positive_k7 or report["positive_18_control"] != {
        "vertices": 18, "K7_seeds": 0,
        "status": "PASS_NOT_APPLICABLE_NO_K7",
    }:
        raise ValueError("known realizable 18-point positive control failed")

    output = {
        "schema": 1,
        "kind": "d6_k7_ranktwo_pentad_full_verification",
        "status": "PASS",
        "report": {"path": str(args.report), "sha256": observed_hash},
        "symbolic_checks": symbolic,
        "target_gate": target_gate,
        "checked": {
            "graphs": len(cover_counts),
            "covers": len(results),
            "pentads": sum(item["pentads_checked"] for item in results),
            "certificates": sum(item["certificate_checked"] for item in results),
            "marginal_graph_rejections": len(marginal),
            "positive_18_control": True,
        },
        "execution": {
            "workers": args.workers,
            "batch_size": args.batch_size,
            "fresh_subprocess_per_cover": True,
            "checkpoint": {
                "path": str(args.checkpoint),
                "sha256": checker.sha256(args.checkpoint),
            },
        },
        "source_hashes": {
            Path(__file__).name: checker.sha256(Path(__file__).resolve()),
            "verify_d6_k7_ranktwo_pentad.py": checker.sha256(
                ROOT / "verify_d6_k7_ranktwo_pentad.py"
            ),
            "verify_d6_k7_positive_polynomial_dual.py": checker.sha256(
                ROOT / "verify_d6_k7_positive_polynomial_dual.py"
            ),
        },
        "claim": (
            "The target corpus was independently replayed from frozen sources; "
            "all 1,113 canonical pentads and every stored positive identity "
            "were independently reconstructed with exact arithmetic."
        ),
    }
    atomic_json(args.output, output)
    print(json.dumps({
        "output": str(args.output),
        "sha256": checker.sha256(args.output),
        "status": output["status"],
        "checked": output["checked"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
