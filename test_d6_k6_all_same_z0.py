#!/usr/bin/env python3
"""Controls for the all-K6-systems same-Z0 production bundle."""

from __future__ import annotations

import ast
import json
import unittest

from d6_k6_all_same_z0 import ROOT, evaluate_graph, sha256
from test_d6_k6_lorentz import lower_bound_18_graph
from verify_d6_k6_all_same_z0 import IndependentZeroForcing


SOURCE = ROOT / "d6_k6_all_same_z0.py"
REPORT = ROOT / "d6_k6_all_same_z0_report.json"
ARCHIVE = ROOT / "d6_k6_all_same_z0_certificates.json"
VERIFIER = ROOT / "verify_d6_k6_all_same_z0.py"
VERIFICATION = ROOT / "d6_k6_all_same_z0_verification.json"


EXPECTED_HISTOGRAMS = {
    652_900: {
        "block+zero_forcing": 32,
        "nonbipartite": 22,
        "zero_forcing+nonbipartite": 10,
    },
    2_301_548: {
        "block+nonbipartite": 6,
        "nonbipartite": 25,
        "zero_forcing+nonbipartite": 1,
    },
    2_842_523: {
        "block+zero_forcing": 96,
        "nonbipartite": 15,
        "zero_forcing+nonbipartite": 16,
    },
    3_289_061: {
        "block+nonbipartite": 6,
        "nonbipartite": 25,
        "zero_forcing+nonbipartite": 1,
    },
}


class AllSameZ0ProductionControls(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.report = json.loads(REPORT.read_text(encoding="utf-8"))
        archive = json.loads(ARCHIVE.read_text(encoding="utf-8"))
        cls.certificates = {
            item["index"]: item["certificate"]
            for item in archive["certificates"]
        }

    def test_exact_full_result_and_artifact_boundary(self) -> None:
        expected = {
            "input_graphs": 990,
            "graphs_rejected": 4,
            "graphs_surviving": 986,
            "seeds_checked": 31_590,
            "impossible_seeds": 4,
            "z0_considered": 32_143,
            "z0_matchable": 32_143,
            "z0_zero_forcing_passed": 32_044,
            "z0_zero_forcing_failed": 99,
            "z0_block_passed_after_zero_forcing": 31_997,
            "z0_block_failed_after_zero_forcing": 47,
            "z0_nonbipartite_failed": 411,
            "zero_forcing_components_checked": 219_896,
            "block_components_checked": 219_663,
            "block_components_failed": 47,
            "generic_orientation_cases_checked": 439_326,
            "generic_orientation_support_new_failures": 2,
            "block_subsets_checked": 1_348_101,
            "pure_colorings_considered": 33_635,
            "pure_colorings_matchable": 31_586,
            "joint_support_searches": 31_586,
            "joint_support_dfs_nodes": 95_593,
            "zero_forcing_initial_sets_checked": 90_249,
        }
        self.assertEqual(self.report["schema"], "d6-k6-all-same-z0-v1")
        self.assertEqual(self.report["status"], "COMPLETE")
        for name, value in expected.items():
            with self.subTest(name=name):
                self.assertEqual(self.report[name], value)
        self.assertEqual(
            self.report["rejected_indices"],
            [652_900, 2_301_548, 2_842_523, 3_289_061],
        )
        self.assertEqual(
            self.report["prior_same_Z0_rejected_indices"],
            [652_900, 2_842_523],
        )
        self.assertEqual(
            self.report["new_rejected_indices"], [2_301_548, 3_289_061]
        )
        self.assertEqual(self.report["sources"][SOURCE.name], sha256(SOURCE))
        self.assertEqual(
            self.report["certificate_archive"]["sha256"], sha256(ARCHIVE)
        )

    def test_all_four_exhaustive_certificate_histograms(self) -> None:
        expected_seeds = {
            652_900: [1, 3, 10, 11, 14, 17],
            2_301_548: [2, 4, 10, 12, 14, 18],
            2_842_523: [0, 3, 4, 14, 15, 17],
            3_289_061: [4, 11, 12, 14, 15, 18],
        }
        expected_choices = {
            652_900: 64,
            2_301_548: 32,
            2_842_523: 127,
            3_289_061: 32,
        }
        self.assertEqual(set(self.certificates), set(EXPECTED_HISTOGRAMS))
        for index, histogram in EXPECTED_HISTOGRAMS.items():
            with self.subTest(index=index):
                certificate = self.certificates[index]
                self.assertEqual(certificate["seed"], expected_seeds[index])
                self.assertEqual(len(certificate["choices"]), expected_choices[index])
                self.assertEqual(certificate["histogram"], histogram)
                self.assertEqual(certificate["all_three_pass_count"], 0)
                self.assertNotIn(
                    "block+zero_forcing+nonbipartite", certificate["histogram"]
                )

    def test_two_new_certificates_are_pure_bipartite_quantifier_couplings(self) -> None:
        representatives = {
            2_301_548: {
                "block+nonbipartite": [],
                "nonbipartite": [1, 6],
                "zero_forcing+nonbipartite": [6, 11],
            },
            3_289_061: {
                "block+nonbipartite": [],
                "nonbipartite": [5, 7],
                "zero_forcing+nonbipartite": [7, 16],
            },
        }
        for index, expected_representatives in representatives.items():
            with self.subTest(index=index):
                certificate = self.certificates[index]
                self.assertEqual(
                    certificate["representative_Z0"], expected_representatives
                )
                choices = certificate["choices"]
                self.assertEqual(len(choices), 32)
                self.assertTrue(
                    all(
                        choice["nonbipartite_joint_support_passed"]
                        for choice in choices
                    )
                )
                self.assertEqual(
                    sum(choice["block_support_passed"] for choice in choices), 6
                )
                self.assertEqual(
                    sum(choice["zero_forcing_rank_passed"] for choice in choices),
                    1,
                )
                self.assertFalse(
                    any(
                        choice["block_support_passed"]
                        and choice["zero_forcing_rank_passed"]
                        for choice in choices
                    )
                )
                for choice in choices:
                    if not choice["zero_forcing_rank_passed"]:
                        self.assertTrue(choice["failed_zero_forcing_components"])
                        self.assertTrue(
                            all(
                                not item["passed"]
                                and item["required_dimension"] > 6
                                for item in choice[
                                    "failed_zero_forcing_components"
                                ]
                            )
                        )

    def test_known_realizable_eighteen_point_control(self) -> None:
        decision = evaluate_graph(lower_bound_18_graph())
        self.assertFalse(decision["rejected"])
        self.assertEqual(decision["seeds_checked"], 32)
        self.assertEqual(decision["impossible_seeds"], 0)
        self.assertTrue(self.report["positive_control"]["passed"])


class AllSameZ0IndependentControls(unittest.TestCase):
    def test_independent_zero_forcing_kernel(self) -> None:
        solver = IndependentZeroForcing()
        self.assertEqual(solver.number(()), 0)
        self.assertEqual(solver.number((0,)), 1)
        self.assertEqual(solver.number((0, 0, 0)), 3)
        # P4 has zero-forcing number one.
        self.assertEqual(solver.number((0b0010, 0b0101, 0b1010, 0b0100)), 1)
        # C4 has zero-forcing number two.
        self.assertEqual(solver.number((0b1010, 0b0101, 0b1010, 0b0101)), 2)

    def test_verifier_imports_no_production_engine(self) -> None:
        syntax = ast.parse(VERIFIER.read_text(encoding="utf-8"))
        imported = []
        for node in ast.walk(syntax):
            if isinstance(node, ast.Import):
                imported.extend(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module:
                imported.append(node.module)
        forbidden = [
            name
            for name in imported
            if name.startswith("d6_k6_all_same_z0")
            or name.startswith("d6_k6_bipartite_rank_reference")
        ]
        self.assertEqual(forbidden, [])
        self.assertIn("verify_d6_k6_same_z0", imported)

    def test_frozen_independent_verification(self) -> None:
        verification = json.loads(VERIFICATION.read_text(encoding="utf-8"))
        self.assertEqual(verification["status"], "PASS")
        self.assertEqual(verification["graphs_recomputed"], 990)
        self.assertEqual(verification["graphs_rejected"], 4)
        self.assertEqual(verification["graphs_surviving"], 986)
        self.assertEqual(
            verification["rejected_indices"],
            [652_900, 2_301_548, 2_842_523, 3_289_061],
        )
        self.assertEqual(
            verification["new_rejected_indices"], [2_301_548, 3_289_061]
        )
        self.assertTrue(verification["all_first_impossible_seeds_matched"])
        self.assertTrue(verification["all_archived_Z0_system_booleans_matched"])
        self.assertTrue(
            verification[
                "all_archived_zero_forcing_component_certificates_matched"
            ]
        )
        self.assertTrue(verification["positive_control"]["passed"])
        self.assertEqual(
            verification["inputs"][VERIFIER.name], sha256(VERIFIER)
        )
        observed = {
            item["index"]: item["histogram"]
            for item in verification["certificate_summaries"]
        }
        self.assertEqual(observed, EXPECTED_HISTOGRAMS)


if __name__ == "__main__":
    unittest.main(verbosity=2)
