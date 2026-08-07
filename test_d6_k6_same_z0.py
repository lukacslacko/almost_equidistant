#!/usr/bin/env python3
"""Controls for the exact K6 same-Z0 production and verification bundle."""

from __future__ import annotations

import ast
import json
import unittest

from d6_k6_same_z0 import ROOT, evaluate_graph, sha256
from test_d6_k6_lorentz import lower_bound_18_graph


SOURCE = ROOT / "d6_k6_same_z0.py"
REPORT = ROOT / "d6_k6_same_z0_report.json"
VERIFIER = ROOT / "verify_d6_k6_same_z0.py"
VERIFICATION = ROOT / "d6_k6_same_z0_verification.json"


class SameZ0ProductionControls(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.report = json.loads(REPORT.read_text(encoding="utf-8"))
        cls.rejected = {
            item["index"]: item["decision"]
            for item in cls.report["graph_results"]
            if item["decision"]["rejected"]
        }

    def test_exact_full_result_and_source_boundary(self) -> None:
        expected = {
            "input_graphs": 990,
            "graphs_rejected": 2,
            "graphs_surviving": 988,
            "seeds_checked": 31_602,
            "impossible_seeds": 2,
            "z0_considered": 32_093,
            "z0_matchable": 32_093,
            "z0_block_passed": 32_011,
            "z0_block_failed": 82,
            "z0_nonbipartite_failed": 411,
            "bipartite_components_checked": 219_831,
            "bipartite_components_failed": 82,
            "generic_orientation_cases_checked": 439_662,
            "generic_orientation_support_new_failures": 2,
            "block_subsets_checked": 1_358_802,
            "pure_colorings_considered": 33_649,
            "pure_colorings_matchable": 31_600,
            "joint_support_searches": 31_600,
            "joint_support_dfs_nodes": 95_607,
        }
        self.assertEqual(self.report["status"], "COMPLETE")
        self.assertEqual(self.report["schema"], "d6-k6-same-z0-v1")
        for name, value in expected.items():
            with self.subTest(name=name):
                self.assertEqual(self.report[name], value)
        self.assertEqual(self.report["rejected_indices"], [652_900, 2_842_523])
        self.assertEqual(self.report["sources"][SOURCE.name], sha256(SOURCE))

    def test_rejected_seed_certificates_partition_every_z0(self) -> None:
        expected = {
            652_900: {
                "seed": [1, 3, 10, 11, 14, 17],
                "choices": 64,
                "failure_histogram": {
                    "bipartite_block_support": 32,
                    "nonbipartite_joint_actual_support": 32,
                },
                "cross_histogram": {
                    "both": 0,
                    "block_only": 32,
                    "nonbipartite_only": 32,
                    "neither": 0,
                },
                "representatives": {
                    "block_only": [],
                    "nonbipartite_only": [4],
                },
            },
            2_842_523: {
                "seed": [0, 3, 4, 14, 15, 17],
                "choices": 127,
                "failure_histogram": {
                    "bipartite_block_support": 31,
                    "nonbipartite_joint_actual_support": 96,
                },
                "cross_histogram": {
                    "both": 0,
                    "block_only": 96,
                    "nonbipartite_only": 31,
                    "neither": 0,
                },
                "representatives": {
                    "block_only": [],
                    "nonbipartite_only": [6, 18],
                },
            },
        }
        self.assertEqual(set(self.rejected), set(expected))
        for index, wanted in expected.items():
            with self.subTest(index=index):
                decision = self.rejected[index]
                certificate = decision["certificate"]
                self.assertEqual(decision["first_impossible_seed"], wanted["seed"])
                self.assertEqual(certificate["seed"], wanted["seed"])
                failures = certificate["Z0_failures"]
                self.assertEqual(len(failures), wanted["choices"])
                histogram = {}
                for failure in failures:
                    reason = failure["reason"]
                    histogram[reason] = histogram.get(reason, 0) + 1
                    if reason == "nonbipartite_joint_actual_support":
                        self.assertEqual(
                            failure["allowed_mask_matchable_colorings"], 0
                        )
                        self.assertEqual(failure["coloring_failures"], [])
                    elif reason == "bipartite_block_support":
                        self.assertTrue(failure["failed_components"])
                        for component in failure["failed_components"]:
                            self.assertFalse(
                                any(case["passed"] for case in component["generic_cases"])
                            )
                            self.assertFalse(component["lightlike_case"]["passed"])
                self.assertEqual(histogram, wanted["failure_histogram"])
                cross = certificate["same_Z0_cross_classification"]
                self.assertEqual(cross["histogram"], wanted["cross_histogram"])
                self.assertEqual(
                    cross["representative_Z0"], wanted["representatives"]
                )
                self.assertTrue(cross["every_Z0_passes_exactly_one_system"])
                self.assertEqual(len(cross["choices"]), wanted["choices"])
                self.assertEqual(certificate["counts"]["joint_support_searches"], 0)
                self.assertEqual(
                    certificate["counts"]["generic_orientation_support_new_failures"],
                    0,
                )

    def test_known_realizable_eighteen_point_control(self) -> None:
        decision = evaluate_graph(lower_bound_18_graph())
        self.assertFalse(decision["rejected"])
        self.assertEqual(decision["seeds_checked"], 32)
        self.assertEqual(decision["impossible_seeds"], 0)
        self.assertTrue(self.report["positive_control"]["passed"])


class SameZ0IndependentVerifierControls(unittest.TestCase):
    def test_verifier_imports_no_production_or_k6_engine(self) -> None:
        syntax = ast.parse(VERIFIER.read_text(encoding="utf-8"))
        imported = []
        for node in ast.walk(syntax):
            if isinstance(node, ast.Import):
                imported.extend(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module:
                imported.append(node.module)
        self.assertFalse(
            [name for name in imported if name.startswith("d6_k6")],
            imported,
        )

    def test_frozen_independent_verification(self) -> None:
        verification = json.loads(VERIFICATION.read_text(encoding="utf-8"))
        self.assertEqual(verification["status"], "PASS")
        self.assertEqual(verification["graphs_recomputed"], 990)
        self.assertEqual(verification["graphs_rejected"], 2)
        self.assertEqual(verification["graphs_surviving"], 988)
        self.assertEqual(verification["rejected_indices"], [652_900, 2_842_523])
        self.assertTrue(verification["all_first_impossible_seeds_matched"])
        self.assertTrue(
            verification["all_rejected_seed_Z0_failure_partitions_matched"]
        )
        self.assertTrue(
            verification["all_rejected_seed_Z0_cross_classifications_matched"]
        )
        self.assertTrue(verification["positive_control"]["passed"])
        self.assertEqual(
            verification["inputs"][VERIFIER.name], sha256(VERIFIER)
        )
        expected_cross = {
            652_900: {"block_only": 32, "nonbipartite_only": 32},
            2_842_523: {"block_only": 96, "nonbipartite_only": 31},
        }
        observed = {
            item["index"]: item["cross_classification_histogram"]
            for item in verification["rejected_certificates"]
        }
        self.assertEqual(observed, expected_cross)


if __name__ == "__main__":
    unittest.main(verbosity=2)
