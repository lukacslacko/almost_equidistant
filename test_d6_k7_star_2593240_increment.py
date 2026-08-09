#!/usr/bin/env python3
"""Controls for the production graph-2593240 rank-one-star increment."""

from __future__ import annotations

import copy
import inspect
import json
import tempfile
import unittest
from pathlib import Path

import build_d6_k7_star_2593240_increment as builder
import verify_d6_k7_star_2593240_increment as checker


class K7Star2593240ProductionControls(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        source_hashes = {
            name: builder.sha256(builder.ROOT / name) for name in builder.SOURCE_FILES
        }
        cls.report = builder.build_report(
            {
                "commit": "0" * 40,
                "branch": "TEST",
                "source_sha256": source_hashes,
                "porcelain_lines": [],
                "porcelain_sha256": builder.stable_hash("TEST"),
                "proof_and_checker_sources_equal_committed_blobs": False,
            }
        )

    def test_union_ready_one_graph_increment(self) -> None:
        summary = self.report["summary"]
        self.assertEqual(summary["rejected_indices"], [2_593_240])
        self.assertEqual(summary["graphs_rejected"], 1)
        self.assertEqual(summary["graphs_surviving"], 14)
        self.assertNotIn(2_593_240, summary["ordered_survivor_indices"])
        self.assertEqual(
            builder.stable_hash(summary["ordered_survivor_indices"]),
            builder.EXPECTED_OUTPUT_SHA256,
        )
        self.assertEqual(
            self.report["input"]["artifact"],
            "d6_k7_schur_3949382_increment_report.json",
        )

    def test_independent_checker_imports_neither_builder_nor_probe(self) -> None:
        source = inspect.getsource(checker)
        self.assertNotIn("import build_d6_k7_star_2593240_increment", source)
        self.assertNotIn("import probe_d6_k7_2593240_star_boundary", source)
        self.assertNotIn("from build_d6_k7_star_2593240_increment", source)
        self.assertNotIn("from probe_d6_k7_2593240_star_boundary", source)

    def test_independent_algebra_uses_distinct_lex_route(self) -> None:
        algebra = checker.independent_algebra()
        self.assertEqual(algebra["status"], "EXACT_RANK_ONE_SUPPORT_CONTRADICTION")
        self.assertEqual(
            algebra["method"],
            "independent_one_stage_lex_elimination_and_support_closure",
        )
        self.assertEqual(algebra["tetrads_reconstructed"], 30)
        self.assertEqual(
            algebra["unique_positive_tetrad_point"],
            ["2/3", "5", "1", "1", "1", "1"],
        )
        self.assertEqual(algebra["star_nonzero_entries"]["0,2"], "8")

    def test_checker_replays_both_seed_quantifiers_and_optional_zeros(self) -> None:
        manifest, star, _star_check, _preceding, _preceding_check = (
            checker.load_upstream()
        )
        rebuilt = checker.independent_quantifier_and_patterns(manifest, star)
        self.assertEqual(
            [row["raw_eligible_covers"] for row in rebuilt["quantifier"]],
            [502, 502],
        )
        self.assertEqual(
            [row["current_covers"] for row in rebuilt["quantifier"]],
            [[0, 8], [0, 32]],
        )
        self.assertEqual(
            [row["pair_target_string"] for row in rebuilt["patterns"]],
            ["011111111110111", "111110111111101"],
        )
        self.assertTrue(rebuilt["system_isomorphism"]["checked_exactly"])
        self.assertEqual(rebuilt["entry_reasons"]["unresolved_optional_entries"], 0)

    def test_full_independent_report_replay_without_production_git_gate(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            report_path = Path(temporary) / "report.json"
            builder.atomic_json(report_path, self.report)
            result = checker.verify_report(
                report_path,
                checker.sha256(report_path),
                enforce_source_boundary=False,
            )
        self.assertEqual(result["status"], "PASS")
        self.assertEqual(result["conclusion"]["rejected_indices"], [2_593_240])
        self.assertFalse(result["conclusion"]["floating_point_used"])
        self.assertEqual(result["checked"]["required_k7_seeds"], 2)

    def test_checker_rejects_a_tampered_optional_zero(self) -> None:
        tampered = copy.deepcopy(self.report)
        tampered["certificate"]["exact_patterns"][0]["pair_targets"][0] = 1
        tampered["certificate_sha256"] = checker.stable_hash(tampered["certificate"])
        with tempfile.TemporaryDirectory() as temporary:
            report_path = Path(temporary) / "tampered.json"
            report_path.write_text(
                json.dumps(tampered, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "both exact patterns"):
                checker.verify_report(
                    report_path,
                    checker.sha256(report_path),
                    enforce_source_boundary=False,
                )

    def test_checker_rejects_a_missing_second_seed_branch(self) -> None:
        tampered = copy.deepcopy(self.report)
        tampered["certificate"]["quantifier"].pop()
        tampered["certificate_sha256"] = checker.stable_hash(tampered["certificate"])
        with tempfile.TemporaryDirectory() as temporary:
            report_path = Path(temporary) / "tampered.json"
            report_path.write_text(
                json.dumps(tampered, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "both-seed quantifier"):
                checker.verify_report(
                    report_path,
                    checker.sha256(report_path),
                    enforce_source_boundary=False,
                )


if __name__ == "__main__":
    unittest.main()
