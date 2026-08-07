#!/usr/bin/env python3
"""Controls for the exact rooted standard-18 attachment audit."""

from __future__ import annotations

import copy
import json
import tempfile
import unittest
from pathlib import Path

import d6_standard18_rooted_attachment as producer
import verify_d6_standard18_rooted_attachment as checker


class RootedStandard18AttachmentControls(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.report = producer.build_report()
        cls.temporary = tempfile.TemporaryDirectory()
        cls.report_path = Path(cls.temporary.name) / "report.json"
        cls.write(cls.report, cls.report_path)

    @classmethod
    def tearDownClass(cls) -> None:
        cls.temporary.cleanup()

    @staticmethod
    def write(value: object, path: Path) -> None:
        path.write_text(
            json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )

    def test_exact_summary_and_scope(self) -> None:
        summary = self.report["summary"]
        self.assertEqual(summary["compatible_classes"], 14)
        self.assertEqual(summary["compatible_occurrences"], 39)
        self.assertEqual(summary["rooted_witness_representatives"], 42)
        self.assertEqual(summary["rank_87_classes"], 9)
        self.assertEqual(summary["rank_86_classes"], 5)
        self.assertEqual(
            summary["parents_only_rank_86_deletions"], [3_950_119, 3_950_509]
        )
        self.assertTrue(self.report["semantics"]["conditional_on_standard_embedding"])
        self.assertEqual(self.report["semantics"]["parent_rejections"], 0)
        self.assertFalse(
            self.report["semantics"]["global_support_embedding_uniqueness_claim"]
        )

    def test_every_rooted_case_has_the_same_exact_obstruction(self) -> None:
        self.assertTrue(
            all(
                case["forced_common_scaled_squared_distance_Qsqrt3"]
                == [16, 3, 0, 1]
                and case["required_unit_scaled_squared_distance_Qsqrt3"]
                == [8, 1, 0, 1]
                and len(case["mapped_base_neighbours"]) in (11, 12)
                and case["mapped_pole_neighbour"] in (16, 17)
                for case in self.report["rooted_cases"]
            )
        )

    def test_independent_checker(self) -> None:
        result = checker.verify(self.report_path, None)
        self.assertEqual(result["status"], "PASS")
        self.assertTrue(all(result["checks"].values()))

    def test_overclaim_attachment_and_flex_tampering_are_rejected(self) -> None:
        variants = []

        overclaim = copy.deepcopy(self.report)
        overclaim["semantics"]["parent_rejections"] = 16
        variants.append(("overclaim", overclaim))

        attachment = copy.deepcopy(self.report)
        attachment["rooted_cases"][0][
            "forced_common_scaled_squared_distance_Qsqrt3"
        ] = [8, 1, 0, 1]
        attachment["summary"]["rooted_cases_sha256"] = producer.stable_hash(
            attachment["rooted_cases"]
        )
        variants.append(("attachment", attachment))

        flex = copy.deepcopy(self.report)
        rank86 = next(
            record
            for record in flex["rigidity"]["records"]
            if record["exact_rank"] == 86
        )
        rank86["rank_86_upper_bound_certificate"]["velocity_Qsqrt3"] = [
            [[0, 1, 0, 1] for _ in range(6)] for _ in range(18)
        ]
        flex["summary"]["rigidity_records_sha256"] = producer.stable_hash(
            flex["rigidity"]["records"]
        )
        variants.append(("flex", flex))

        for name, value in variants:
            with self.subTest(name=name):
                path = Path(self.temporary.name) / f"tampered-{name}.json"
                self.write(value, path)
                with self.assertRaises(ValueError):
                    checker.verify(path, None)


if __name__ == "__main__":
    unittest.main(verbosity=2)
