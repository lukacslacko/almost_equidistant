#!/usr/bin/env python3
"""Controls for the exact K6 PSD Z-matrix refinement."""

from __future__ import annotations

import ast
import hashlib
import json
import unittest
from pathlib import Path

import d6_k6_psd_zmatrix as production
import verify_d6_k6_psd_zmatrix as independent
from d6_k6_bipartite_rank_reference import ZeroForcingSolver
from d6_k6_normal_inertia import InertiaCache
from test_d6_k6_lorentz import lower_bound_18_graph


ROOT = Path(__file__).resolve().parent
EXPECTED_HASHES = {
    "d6_k6_psd_zmatrix.py": (
        "24c2b3abfce5074560183d17f77778f96298bef6e38ece5bc937ccb2a4cba736"
    ),
    "d6_k6_psd_zmatrix_report.json": (
        "abb920fbd9eb848502cca37e32c2aa86bf220ab1773e6006ea1fbee4bc2f4e30"
    ),
    "d6_k6_psd_zmatrix_certificates.json": (
        "62fc7170d97056941bcaf553c45cfa1cabdf20117ffc69fbd8429b3764ec999c"
    ),
    "verify_d6_k6_psd_zmatrix.py": (
        "ca4ca62f75d02846f12e59e16eb9eab4473fe5438fc496e72a25a39db77f9042"
    ),
    "d6_k6_psd_zmatrix_verification.json": (
        "8261b58e788f120e3d1eee6dc481a4989e3eb75f9f03babbd552474b4a5302eb"
    ),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


class PsdZMatrixTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.report = json.loads(
            (ROOT / "d6_k6_psd_zmatrix_report.json").read_text()
        )
        cls.archive = json.loads(
            (ROOT / "d6_k6_psd_zmatrix_certificates.json").read_text()
        )
        cls.verification = json.loads(
            (ROOT / "d6_k6_psd_zmatrix_verification.json").read_text()
        )

    def test_small_connected_graph_controls(self) -> None:
        rows = {
            "K1": (0,),
            "K2": (2, 1),
            "P3": (2, 5, 2),
            "C5": (18, 5, 10, 20, 9),
        }
        expected = {
            "K1": {"positive": 1, "negative": 0},
            "K2": {"positive": 2, "negative": 1},
            "P3": {"positive": 2, "negative": 2},
            "C5": {"positive": 3, "negative": 4},
        }
        observed = {}
        for name, graph in rows.items():
            observed[name] = {}
            for sign in ("positive", "negative"):
                result = production.side_rank_lower(
                    graph,
                    tuple(range(len(graph))),
                    sign,
                    InertiaCache(),
                    ZeroForcingSolver(),
                )
                observed[name][sign] = result.rank_lower
        self.assertEqual(observed, expected)
        self.assertEqual(independent.kernel_controls(), expected)
        self.assertIsNone(independent.pf_nullity_upper(rows["C5"], "positive"))
        self.assertEqual(independent.pf_nullity_upper(rows["C5"], "negative"), 1)
        self.assertEqual(independent.pf_nullity_upper(rows["P3"], "positive"), 1)

    def test_frozen_production_totals(self) -> None:
        self.assertEqual(self.report["status"], "COMPLETE")
        self.assertEqual(self.report["input_graphs"], 977)
        self.assertEqual(self.report["graphs_rejected"], 116)
        self.assertEqual(self.report["graphs_surviving"], 861)
        self.assertEqual(
            production.stable_hash(self.report["rejected_indices"]),
            "adb32fa4c00f6209ce7e246fd49337a6c030b20ddd49193765cf3451afb4c36c",
        )
        expected = {
            "seeds_checked": 28_597,
            "impossible_seeds": 116,
            "z0_considered": 35_844,
            "z0_matchable": 35_844,
            "z0_psd_zmatrix_passed": 28_733,
            "z0_psd_zmatrix_failed": 7_111,
            "z0_nonbipartite_failed": 252,
            "bipartite_components_checked": 221_715,
            "bipartite_components_failed": 7_111,
            "prior_fused_components_newly_failed": 695,
            "generic_orientation_cases_checked": 443_430,
            "prior_fused_orientations_newly_failed": 1_959,
            "side_graph_connected_blocks_checked": 541_380,
            "psd_zmatrix_blocks_applicable": 517_485,
            "psd_zmatrix_strict_component_rank_improvements": 14_579,
            "prior_fused_block_subsets_checked": 1_918_388,
            "psd_zmatrix_block_subsets_checked": 1_870_804,
        }
        for name, value in expected.items():
            self.assertEqual(self.report[name], value, name)

    def test_every_rejection_has_exhaustive_psd_failure(self) -> None:
        self.assertEqual(self.archive["status"], "COMPLETE")
        self.assertEqual(
            self.archive["rejected_indices"], self.report["rejected_indices"]
        )
        self.assertEqual(len(self.archive["certificates"]), 116)
        choices = [
            row
            for item in self.archive["certificates"]
            for row in item["certificate"]["choices"]
        ]
        self.assertEqual(len(choices), 4_987)
        self.assertTrue(all(not row["psd_zmatrix_side_rank_passed"] for row in choices))
        self.assertTrue(
            all(item["certificate"]["common_pass_count"] == 0 for item in self.archive["certificates"])
        )

    def test_fixed_strict_witness(self) -> None:
        item = next(
            item for item in self.archive["certificates"] if item["index"] == 3278
        )
        certificate = item["certificate"]
        self.assertEqual(certificate["seed"], [0, 1, 3, 8, 12, 16])
        empty = next(row for row in certificate["choices"] if row["Z0"] == [])
        failure = empty["first_failed_psd_zmatrix_component"]
        self.assertEqual(
            failure["component"], [2, 5, 6, 9, 10, 11, 13, 14, 17, 18]
        )
        orientations = failure["orientations"]
        self.assertEqual(
            [
                (
                    case["side_A"]["prior_fused_rank_lower"],
                    case["side_A"]["rank_lower"],
                    case["side_B"]["prior_fused_rank_lower"],
                    case["side_B"]["rank_lower"],
                )
                for case in orientations
            ],
            [(3, 3, 3, 5), (2, 3, 4, 4)],
        )
        self.assertTrue(all(case["prior_fused_support_passed"] for case in orientations))
        self.assertTrue(all(not case["psd_zmatrix_support_passed"] for case in orientations))
        self.assertFalse(failure["lightlike"]["passed"])

    def test_known_realizable_positive_control(self) -> None:
        result = production.evaluate_graph(lower_bound_18_graph())
        self.assertFalse(result["rejected"])
        self.assertEqual(result["seeds_checked"], 32)
        self.assertEqual(result["impossible_seeds"], 0)

    def test_independent_verifier_and_import_boundary(self) -> None:
        self.assertEqual(self.verification["status"], "PASS")
        self.assertEqual(self.verification["graphs_recomputed"], 977)
        self.assertEqual(self.verification["graphs_rejected"], 116)
        self.assertTrue(self.verification["all_exhaustive_Z0_certificates_matched"])
        tree = ast.parse(
            (ROOT / "verify_d6_k6_psd_zmatrix.py").read_text(encoding="utf-8")
        )
        imported = {
            node.module
            for node in ast.walk(tree)
            if isinstance(node, ast.ImportFrom)
        }
        imported.update(
            alias.name
            for node in ast.walk(tree)
            if isinstance(node, ast.Import)
            for alias in node.names
        )
        self.assertNotIn("d6_k6_psd_zmatrix", imported)

    def test_frozen_hashes(self) -> None:
        self.assertEqual(
            {name: sha256(ROOT / name) for name in EXPECTED_HASHES},
            EXPECTED_HASHES,
        )


if __name__ == "__main__":
    unittest.main()
