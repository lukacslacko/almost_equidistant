#!/usr/bin/env python3
"""Focused controls for the independent v5 interval verifier."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import verify_d6_interval_v5 as verifier


class V5IntervalVerifierControls(unittest.TestCase):
    def test_reconstructs_exact_911_base(self) -> None:
        sample_report = verifier.ROOT / ".runs/d6_interval_v5_bench64_cap20000.json"
        if not sample_report.exists():
            self.skipTest("local benchmark is not present")
        configuration = json.loads(sample_report.read_text())["configuration"]
        graphs = verifier.reconstruct_v5(configuration)
        self.assertEqual(len(graphs), 64)
        self.assertEqual(
            verifier.object_sha256([graph["index"] for graph in graphs]),
            configuration["selection"]["indices_sha256"],
        )

    def test_completed_sample_audits_without_replay(self) -> None:
        sample_report = verifier.ROOT / ".runs/d6_interval_v5_bench64_cap20000.json"
        if not sample_report.exists():
            self.skipTest("local benchmark is not present")
        result = verifier.verify_report(sample_report, None, replay=False)
        self.assertEqual(result["status"], "PASS")
        self.assertEqual(result["selection"]["graphs"], 64)
        self.assertEqual(result["certified_killed"], 0)
        self.assertEqual(result["status_counts"]["ABORT"], 64)

    def test_report_hash_tamper_rejected_first(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "report.json"
            path.write_text("{}\n", encoding="utf-8")
            with self.assertRaisesRegex(AssertionError, "report hash"):
                verifier.verify_report(path, "0" * 64, replay=False)

    def test_does_not_import_production_wrapper(self) -> None:
        source = (verifier.ROOT / "verify_d6_interval_v5.py").read_text(
            encoding="utf-8"
        )
        self.assertNotIn("import run_d6_interval_v5", source)


if __name__ == "__main__":
    unittest.main()
