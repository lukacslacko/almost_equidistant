#!/usr/bin/env python3
"""Controls for the production graph-3949382 exact Schur increment."""

from __future__ import annotations

import copy
import inspect
import json
import tempfile
import unittest
from pathlib import Path

import build_d6_k7_schur_3949382_increment as builder
import verify_d6_k7_schur_3949382_increment as checker


class K7Schur3949382ProductionControls(unittest.TestCase):
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
        self.assertEqual(summary["rejected_indices"], [3_949_382])
        self.assertEqual(summary["graphs_rejected"], 1)
        self.assertEqual(summary["graphs_surviving"], 15)
        self.assertNotIn(3_949_382, summary["ordered_survivor_indices"])
        self.assertEqual(
            builder.stable_hash(summary["ordered_survivor_indices"]),
            builder.EXPECTED_OUTPUT_SHA256,
        )

    def test_independent_checker_imports_neither_builder_nor_probe(self) -> None:
        source = inspect.getsource(checker)
        self.assertNotIn("import build_d6_k7_schur_3949382_increment", source)
        self.assertNotIn("import probe_d6_k7_schur_3949382", source)
        self.assertNotIn("from build_d6_k7_schur_3949382_increment", source)
        self.assertNotIn("from probe_d6_k7_schur_3949382", source)

    def test_independent_algebra_has_strict_contradiction_and_boundary_control(self) -> None:
        algebra = checker.independent_algebra()
        self.assertEqual(algebra["status"], "EXACT_STRICT_POSITIVITY_CONTRADICTION")
        self.assertEqual(algebra["method"], "independent_ideal_membership_A0_B0_elimination")
        self.assertTrue(algebra["boundary_control"]["all_ten_equations_zero"])
        self.assertFalse(algebra["boundary_control"]["strict_domain"])

    def test_checker_reconstructs_quantifier_and_optional_zero_justifications(self) -> None:
        manifest, star, _star_check = checker.load_upstream()
        rebuilt = checker.independent_quantifier_and_pattern(manifest, star)
        self.assertEqual(rebuilt["quantifier"]["raw_eligible_covers"], 128)
        self.assertEqual(rebuilt["quantifier"]["prior_layer_eliminated_raw_covers"], 125)
        self.assertEqual(rebuilt["quantifier"]["sole_new_branch_passing_families"], 1)
        self.assertEqual(rebuilt["pattern"]["pair_target_string"], "0111011011")
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
        self.assertEqual(result["conclusion"]["rejected_indices"], [3_949_382])
        self.assertFalse(result["conclusion"]["floating_point_used"])

    def test_checker_rejects_a_tampered_optional_zero(self) -> None:
        tampered = copy.deepcopy(self.report)
        tampered["certificate"]["exact_pattern"]["pair_targets"][0] = 1
        tampered["certificate_sha256"] = checker.stable_hash(tampered["certificate"])
        with tempfile.TemporaryDirectory() as temporary:
            report_path = Path(temporary) / "tampered.json"
            report_path.write_text(
                json.dumps(tampered, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "pattern field pair_targets"):
                checker.verify_report(
                    report_path,
                    checker.sha256(report_path),
                    enforce_source_boundary=False,
                )


if __name__ == "__main__":
    unittest.main()
