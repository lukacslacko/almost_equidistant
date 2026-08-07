#!/usr/bin/env python3
"""Controls for the exact standard-18 geometry package."""

from __future__ import annotations

import copy
import json
import tempfile
import unittest
from pathlib import Path

import d6_standard18_geometry as producer
import verify_d6_standard18_geometry as checker


class Standard18GeometryControls(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.report = producer.build_report()
        cls.temporary = tempfile.TemporaryDirectory()
        cls.path = Path(cls.temporary.name) / "report.json"
        cls.write(cls.report, cls.path)

    @classmethod
    def tearDownClass(cls) -> None:
        cls.temporary.cleanup()

    @staticmethod
    def write(value: object, path: Path) -> None:
        path.write_text(
            json.dumps(value, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )

    def test_exact_geometry_summary(self) -> None:
        self.assertEqual(self.report["unit_graph"]["edges"], 112)
        self.assertEqual(self.report["unit_graph"]["base_edges"], 80)
        self.assertEqual(self.report["unit_graph"]["K6_seeds"], 32)
        self.assertEqual(self.report["unit_graph"]["K7_seeds"], 0)
        self.assertEqual(self.report["rigidity"]["full_framework"]["exact_rank"], 87)
        self.assertEqual(self.report["affine_spanning"]["all_11_subsets_checked"], 4368)
        self.assertFalse(self.report["semantics"]["classification_claim"])

    def test_independent_checker(self) -> None:
        result = checker.verify(self.path, None)
        self.assertEqual(result["status"], "PASS")
        self.assertTrue(all(result["checks"].values()))
        self.assertTrue(result["exact_results"]["standard_coordinates_nonextendable"])

    def test_rank_or_scope_tamper_is_rejected(self) -> None:
        for mutate in ("rank", "scope"):
            report = copy.deepcopy(self.report)
            if mutate == "rank":
                report["rigidity"]["full_framework"]["exact_rank"] = 86
            else:
                report["semantics"]["classification_claim"] = True
            path = Path(self.temporary.name) / f"tampered-{mutate}.json"
            self.write(report, path)
            with self.assertRaises(ValueError):
                checker.verify(path, None)


if __name__ == "__main__":
    unittest.main(verbosity=2)
