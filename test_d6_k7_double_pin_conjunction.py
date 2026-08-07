#!/usr/bin/env python3
"""Artifact controls for the exact K7 double-pin conjunction."""

import json
import hashlib
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parent


class DoublePinConjunctionArtifactTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not (ROOT / "d6_k7_double_pin_conjunction_report.json").exists():
            raise unittest.SkipTest("production report not generated yet")
        if not (ROOT / "d6_k7_double_pin_conjunction_verification.json").exists():
            raise unittest.SkipTest("production verification not generated yet")
        cls.report = json.loads((
            ROOT / "d6_k7_double_pin_conjunction_report.json"
        ).read_text(encoding="utf-8"))
        cls.verification = json.loads((
            ROOT / "d6_k7_double_pin_conjunction_verification.json"
        ).read_text(encoding="utf-8"))

    def test_exact_partition(self) -> None:
        summary = self.report["summary"]
        self.assertEqual(summary["rejected"], 130)
        self.assertEqual(summary["survivors"], 25)
        self.assertEqual(
            summary["survivor_indices_sha256"],
            "16872c94cce864ca17689714c05cb3a49e2cf6cb658d808c270d047c2794ac04",
        )
        self.assertEqual(len(summary["rejected_indices"]), 130)
        self.assertEqual(len(summary["survivor_indices"]), 25)
        self.assertFalse(
            set(summary["rejected_indices"]) & set(summary["survivor_indices"])
        )

    def test_full_quantifier_counts(self) -> None:
        totals = self.report["summary"]["totals"]
        self.assertEqual(totals["seeds"], 237)
        self.assertEqual(totals["eligible_covers"], 87_270)
        self.assertEqual(totals["current_covers"], 535)
        self.assertEqual(totals["labeled_z_families"], 19_716)
        self.assertEqual(totals["pre_double_pin_families"], 555)
        self.assertEqual(totals["double_pin_infeasible_families"], 336)
        self.assertEqual(totals["double_pin_passing_families"], 219)
        self.assertEqual(totals["double_pin_infeasible_covers"], 419)
        self.assertEqual(totals["double_pin_passing_covers"], 116)

    def test_interval_discovery_targets_and_control(self) -> None:
        decisions = {
            int(record["index"]): record["decision"]
            for record in self.report["records"]
        }
        self.assertEqual(decisions[379078], "REJECTED")
        self.assertEqual(decisions[2280137], "REJECTED")
        self.assertEqual(decisions[3936176], "SURVIVOR")

    def test_independent_verification_passed(self) -> None:
        self.assertEqual(self.verification["status"], "PASS")
        self.assertEqual(self.verification["checked"]["graphs"], 155)
        self.assertEqual(self.verification["checked"]["rejected_graphs"], 130)
        self.assertEqual(self.verification["checked"]["survivors"], 25)
        digest = hashlib.sha256((
            ROOT / "d6_k7_double_pin_conjunction_report.json"
        ).read_bytes()).hexdigest()
        self.assertEqual(self.verification["report"]["sha256"], digest)


if __name__ == "__main__":
    unittest.main()
