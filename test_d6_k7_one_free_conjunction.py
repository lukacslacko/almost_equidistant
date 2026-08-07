#!/usr/bin/env python3
"""Integrity checks for the produced K7 one-free conjunction artifacts."""

from __future__ import annotations

import hashlib
import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parent
REPORT = ROOT / "d6_k7_one_free_conjunction_report.json"
VERIFICATION = ROOT / "d6_k7_one_free_conjunction_verification.json"
PIPELINES = ("old", "singleton", "correlated", "combined")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_hash(value: object) -> str:
    return hashlib.sha256(
        json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
        ).encode("ascii")
    ).hexdigest()


class ArtifactIntegrityTests(unittest.TestCase):
    def test_report_and_verification_bind_every_boundary(self) -> None:
        if not REPORT.exists() and not VERIFICATION.exists():
            self.skipTest("official conjunction artifacts not generated yet")
        self.assertTrue(REPORT.exists())
        self.assertTrue(VERIFICATION.exists())
        report = json.loads(REPORT.read_text(encoding="utf-8"))
        verification = json.loads(VERIFICATION.read_text(encoding="utf-8"))
        self.assertEqual(report.get("schema"), 1)
        self.assertEqual(
            report.get("kind"), "d6_k7_one_free_edge_seed_conjunction"
        )
        self.assertEqual(report.get("status"), "COMPLETE")
        self.assertEqual(len(report["input_indices"]), 24)
        self.assertEqual(
            report["input_indices_sha256"], stable_hash(report["input_indices"])
        )
        self.assertEqual(
            [record["index"] for record in report["records"]],
            report["input_indices"],
        )
        summary = report["summary"]
        for pipeline in PIPELINES:
            rejected = summary["rejected_by_pipeline"][pipeline]
            survivors = summary["survivors_by_pipeline"][pipeline]
            self.assertEqual(
                summary["rejected_by_pipeline_sha256"][pipeline],
                stable_hash(rejected),
            )
            self.assertEqual(
                summary["survivors_by_pipeline_sha256"][pipeline],
                stable_hash(survivors),
            )
            self.assertEqual(
                rejected,
                [
                    record["index"]
                    for record in report["records"]
                    if record["decision_by_pipeline"][pipeline] == "REJECTED"
                ],
            )
            self.assertEqual(
                survivors,
                [
                    record["index"]
                    for record in report["records"]
                    if record["decision_by_pipeline"][pipeline] == "SURVIVOR"
                ],
            )
        self.assertFalse(summary["rejected_by_pipeline"]["old"])
        self.assertEqual(
            summary["survivors_by_pipeline"]["old"], report["input_indices"]
        )
        combined = set(summary["rejected_by_pipeline"]["combined"])
        self.assertLessEqual(
            set(summary["rejected_by_pipeline"]["singleton"]), combined
        )
        self.assertLessEqual(
            set(summary["rejected_by_pipeline"]["correlated"]), combined
        )
        for name, expected in report["source_sha256"].items():
            self.assertEqual(sha256(ROOT / name), expected)
        for name, expected in report["upstream_artifact_sha256"].items():
            self.assertEqual(sha256(ROOT / name), expected)

        self.assertEqual(verification.get("schema"), 1)
        self.assertEqual(
            verification.get("kind"),
            "d6_k7_one_free_edge_seed_conjunction_verification",
        )
        self.assertEqual(verification.get("status"), "PASS")
        self.assertEqual(verification["report"]["sha256"], sha256(REPORT))
        self.assertEqual(
            verification["upstream_artifact_sha256"],
            report["upstream_artifact_sha256"],
        )
        self.assertEqual(
            verification["exact_survivors_sha256"],
            summary["survivors_by_pipeline_sha256"]["combined"],
        )


if __name__ == "__main__":
    unittest.main()
