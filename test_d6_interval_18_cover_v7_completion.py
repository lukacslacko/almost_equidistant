#!/usr/bin/env python3
"""Focused tests for the completed d=6 n=18 campaign audit."""

from __future__ import annotations

import copy
import json
import unittest

import verify_d6_interval_18_cover_v7_completion as checker


class CompletionAuditTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.campaign = checker.load_json(checker.RAW_CAMPAIGN)
        cls.report = checker.load_json(checker.RAW_REPORT)
        cls.production, cls.selection = checker.rebuild_selection(
            cls.campaign["configuration"]
        )

    def test_independent_selection_reconstruction(self) -> None:
        self.assertEqual(len(self.production), 181)
        self.assertEqual(self.production[32]["index"], checker.KILLED_CLASS)
        self.assertEqual(self.production[-1]["index"], 11_957)
        self.assertEqual(
            checker.stable_hash([item["index"] for item in self.production]),
            checker.PRODUCTION_INDEX_HASH,
        )
        self.assertEqual(self.selection["backup_classes"], [126, 5_673])

    def test_live_completion_audit(self) -> None:
        result = checker.verify_completion(source_commit=None)
        self.assertEqual(result["status"], "PASS")
        audit = result["result_audit"]
        self.assertEqual(audit["checkpoint_files"], 181)
        self.assertEqual(audit["sole_certified_killed_class"], 2_100)
        self.assertEqual(audit["sole_unresolved_no_order_class"], 7_259)
        self.assertEqual(audit["nonrejecting_results"], 180)
        self.assertTrue(audit["decisions_exactly_reconstructed"])

    def test_status_semantics_fail_closed(self) -> None:
        altered = copy.deepcopy(self.report["results"])
        altered[0]["status"] = "KILLED"
        with self.assertRaises(AssertionError):
            checker.validate_result(altered[0], self.production[0])

    def test_decisions_reconstruction_detects_tamper(self) -> None:
        expected = checker.decision_text(self.report["results"])
        actual = checker.RAW_DECISIONS.read_text(encoding="utf-8")
        self.assertEqual(expected, actual)
        tampered = actual.replace("\tABORT\t", "\tKILLED\t", 1)
        self.assertNotEqual(expected, tampered)

    def test_checkpoint_identity_detects_tamper(self) -> None:
        path = (
            checker.RUN_DIR
            / "results"
            / "result_000032_0002100.json"
        )
        value = json.loads(path.read_text(encoding="utf-8"))
        value["result"]["index"] += 1
        with self.assertRaises(AssertionError):
            checker.validate_result(value["result"], self.production[32])


if __name__ == "__main__":
    unittest.main()
