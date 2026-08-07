#!/usr/bin/env python3
"""Controls for the heuristic d=6 18-deletion numerical screen."""

from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

import run_d6_18_numerical_screen as screen
import verify_d6_18_numerical_screen as verifier


ROOT = Path(__file__).resolve().parent


class NumericalScreenTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.corpus = json.loads(
            (ROOT / "d6_residue_18_deletions.json").read_text(encoding="utf-8")
        )

    def test_standard_controls_are_complete_and_exactly_embedded(self) -> None:
        graph_text, controls = screen.graph_input(self.corpus)
        expected = [
            index
            for index, record in enumerate(self.corpus["unique_deletions"])
            if record["standard18_compatible"]
        ]
        self.assertEqual(sorted(controls), expected)
        self.assertEqual(len(controls), 14)
        self.assertEqual(len(graph_text.splitlines()), 12_712)
        self.assertEqual(
            {item["mod13_rigidity_rank_lower_bound"] for item in controls.values()},
            {86, 87},
        )

    def test_engine_seeded_positive_control(self) -> None:
        graph_text, controls = screen.graph_input(self.corpus)
        index = min(controls)
        with tempfile.TemporaryDirectory() as directory_name:
            directory = Path(directory_name)
            graph_path = directory / "graphs.txt"
            graph_path.write_text(graph_text, encoding="ascii")
            binary = directory / "screen"
            subprocess.run(
                [
                    "clang", "-O3", "-std=c11", "-Wall", "-Wextra", "-Werror",
                    str(ROOT / "d6_18_lm_screen.c"), "-lm", "-o", str(binary),
                ],
                cwd=ROOT,
                check=True,
            )
            completed = subprocess.run(
                [str(binary), str(graph_path), str(index), str(index + 1), "0", "1"],
                cwd=ROOT,
                capture_output=True,
                text=True,
                check=True,
            )
            row = screen.parse_result_line(completed.stdout.strip())
            self.assertEqual(row["index"], index)
            self.assertEqual(row["has_seed"], 1)
            self.assertLess(row["seed_residual"], 1e-20)
            self.assertGreater(row["seed_minimum_distance"], 1e-3)
            self.assertEqual(
                row["seed_numerical_rigidity_rank"],
                controls[index]["mod13_rigidity_rank_lower_bound"],
            )

    def test_completed_report_verifies_when_present(self) -> None:
        report = ROOT / "d6_18_numerical_screen_report.json"
        if not report.exists():
            self.skipTest("full numerical report has not been generated")
        result = verifier.verify(
            report, ROOT / "d6_residue_18_deletions.json"
        )
        self.assertEqual(result["status"], "PASS", result["failed_checks"])


if __name__ == "__main__":
    unittest.main()
