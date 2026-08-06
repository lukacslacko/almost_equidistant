#!/usr/bin/env python3
"""Controls for the exact K6 fused side-rank production bundle."""

from __future__ import annotations

import ast
import json
import unittest

from d6_k6_fused_side_rank import ROOT, sha256
from d6_k6_normal_block_support import OrthogonalBlock, check_block_subsets


SOURCE = ROOT / "d6_k6_fused_side_rank.py"
REPORT = ROOT / "d6_k6_fused_side_rank_report.json"
ARCHIVE = ROOT / "d6_k6_fused_side_rank_certificates.json"
VERIFIER = ROOT / "verify_d6_k6_fused_side_rank.py"
VERIFICATION = ROOT / "d6_k6_fused_side_rank_verification.json"

REJECTED = [
    122_871,
    129_414,
    652_900,
    665_962,
    719_091,
    1_401_532,
    1_948_946,
    2_301_548,
    2_842_523,
    3_274_548,
    3_289_061,
    3_959_774,
    3_962_868,
]
NEW_REJECTED = [
    122_871,
    129_414,
    665_962,
    719_091,
    1_401_532,
    1_948_946,
    3_274_548,
    3_959_774,
    3_962_868,
]
EXPECTED_HISTOGRAMS = {
    122_871: {"neither": 8, "nonbipartite_only": 24},
    129_414: {"neither": 4, "nonbipartite_only": 28},
    652_900: {"fused_only": 32, "nonbipartite_only": 32},
    665_962: {"neither": 4, "nonbipartite_only": 28},
    719_091: {"neither": 43, "nonbipartite_only": 84},
    1_401_532: {"nonbipartite_only": 64},
    1_948_946: {"nonbipartite_only": 32},
    2_301_548: {"nonbipartite_only": 32},
    2_842_523: {"fused_only": 96, "nonbipartite_only": 31},
    3_274_548: {"nonbipartite_only": 16},
    3_289_061: {"nonbipartite_only": 32},
    3_959_774: {"nonbipartite_only": 32},
    3_962_868: {"neither": 4, "nonbipartite_only": 28},
}


class FusedKernelControls(unittest.TestCase):
    def test_sidewise_max_is_stronger_than_two_separate_totals(self) -> None:
        z0_dimension = 2
        zero_forcing = (3, 1)
        inertia = (1, 3)
        self.assertEqual(z0_dimension + sum(zero_forcing), 6)
        self.assertEqual(z0_dimension + sum(inertia), 6)

        old_inertia_blocks = (
            OrthogonalBlock("A", 1, 0b000111),
            OrthogonalBlock("B", 3, 0b111000),
            OrthogonalBlock("Z0:x", 1, 0b000001),
            OrthogonalBlock("Z0:y", 1, 0b000010),
        )
        self.assertTrue(check_block_subsets(old_inertia_blocks).passed)
        fused_blocks = (
            OrthogonalBlock("A", max(zero_forcing[0], inertia[0]), 0b000111),
            OrthogonalBlock("B", max(zero_forcing[1], inertia[1]), 0b111000),
            OrthogonalBlock("Z0:x", 1, 0b000001),
            OrthogonalBlock("Z0:y", 1, 0b000010),
        )
        result = check_block_subsets(fused_blocks)
        self.assertFalse(result.passed)


class FusedProductionControls(unittest.TestCase):
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
            "graphs_rejected": 13,
            "graphs_surviving": 977,
            "seeds_checked": 31_399,
            "impossible_seeds": 13,
            "z0_considered": 32_337,
            "z0_matchable": 32_337,
            "z0_fused_passed": 31_792,
            "z0_fused_failed": 545,
            "z0_nonbipartite_failed": 406,
            "bipartite_components_checked": 219_765,
            "bipartite_components_failed": 545,
            "components_old_separate_passed_fused_failed": 15,
            "generic_orientation_cases_checked": 439_530,
            "orientation_fusion_new_failures": 485,
            "old_block_subsets_checked": 1_403_884,
            "fused_block_subsets_checked": 1_398_992,
            "pure_colorings_considered": 33_426,
            "pure_colorings_matchable": 31_386,
            "joint_support_searches": 31_386,
            "joint_support_dfs_nodes": 95_166,
            "zero_forcing_initial_sets_checked": 89_386,
        }
        self.assertEqual(self.report["schema"], "d6-k6-fused-side-rank-v1")
        self.assertEqual(self.report["status"], "COMPLETE")
        for name, value in expected.items():
            with self.subTest(name=name):
                self.assertEqual(self.report[name], value)
        self.assertEqual(self.report["rejected_indices"], REJECTED)
        self.assertEqual(self.report["new_rejected_indices"], NEW_REJECTED)
        self.assertEqual(self.report["sources"][SOURCE.name], sha256(SOURCE))
        self.assertEqual(
            self.report["certificate_archive"]["sha256"], sha256(ARCHIVE)
        )

    def test_fixed_corpus_witness_passes_old_rules_but_fails_fusion(self) -> None:
        certificate = self.certificates[122_871]
        empty = next(choice for choice in certificate["choices"] if not choice["Z0"])
        component = next(
            item
            for item in empty["failed_fused_components"]
            if item["old_separate_systems_passed"]
        )
        self.assertEqual(
            component["component"], [0, 1, 2, 5, 8, 10, 14, 15, 16, 17]
        )
        self.assertEqual(component["zero_forcing"]["rank_lower_A"], 2)
        self.assertEqual(component["zero_forcing"]["rank_lower_B"], 3)
        self.assertEqual(component["zero_forcing"]["required_dimension"], 5)
        self.assertTrue(component["zero_forcing"]["passed"])
        positive = component["orientations"][0]
        self.assertEqual(positive["case"], "A_positive")
        self.assertEqual(
            (positive["inertia_rank_lower_A"], positive["inertia_rank_lower_B"]),
            (4, 2),
        )
        self.assertTrue(positive["old_inertia_support_passed"])
        self.assertEqual(
            (positive["fused_rank_lower_A"], positive["fused_rank_lower_B"]),
            (4, 3),
        )
        self.assertFalse(positive["fused_support_passed"])
        self.assertTrue(positive["old_separate_tests_passed_but_fused_failed"])
        self.assertEqual(positive["fused_first_failure"]["rank_lower"], 7)
        self.assertEqual(positive["fused_first_failure"]["coordinate_capacity"], 6)
        self.assertFalse(component["lightlike"]["passed"])

    def test_all_exhaustive_certificate_histograms(self) -> None:
        self.assertEqual(set(self.certificates), set(REJECTED))
        self.assertEqual(
            sum(len(certificate["choices"]) for certificate in self.certificates.values()),
            654,
        )
        for index, histogram in EXPECTED_HISTOGRAMS.items():
            with self.subTest(index=index):
                certificate = self.certificates[index]
                self.assertEqual(certificate["histogram"], histogram)
                self.assertEqual(certificate["common_pass_count"], 0)
        for index in NEW_REJECTED:
            self.assertTrue(
                all(
                    not choice["fused_side_rank_passed"]
                    for choice in self.certificates[index]["choices"]
                )
            )

    def test_known_realizable_eighteen_point_control(self) -> None:
        self.assertTrue(self.report["positive_control"]["passed"])
        self.assertEqual(self.report["positive_control"]["K6_seeds"], 32)
        self.assertEqual(self.report["positive_control"]["impossible_seeds"], 0)


class FusedIndependentVerifierControls(unittest.TestCase):
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
            if name.startswith("d6_k6_fused_side_rank")
            or name.startswith("d6_k6_normal_inertia")
            or name.startswith("d6_k6_bipartite_rank_reference")
        ]
        self.assertEqual(forbidden, [])
        self.assertIn("verify_d6_k6_all_same_z0", imported)

    def test_frozen_independent_verification(self) -> None:
        verification = json.loads(VERIFICATION.read_text(encoding="utf-8"))
        self.assertEqual(verification["status"], "PASS")
        self.assertEqual(verification["graphs_recomputed"], 990)
        self.assertEqual(verification["graphs_rejected"], 13)
        self.assertEqual(verification["graphs_surviving"], 977)
        self.assertEqual(verification["rejected_indices"], REJECTED)
        self.assertEqual(verification["new_rejected_indices"], NEW_REJECTED)
        self.assertTrue(verification["all_first_impossible_seeds_matched"])
        self.assertTrue(
            verification["all_archived_fused_component_certificates_matched"]
        )
        self.assertTrue(
            verification["all_archived_nonbipartite_booleans_matched"]
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
