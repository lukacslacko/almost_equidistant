#!/usr/bin/env python3
"""Focused tests for the fixed-maximum-clique basis-extension pilot."""

from __future__ import annotations

import json
import unittest
from pathlib import Path

import d6_k7_arbitrary_basis_psd_dual as dual
import run_d6_k7_arbitrary_basis_extension_pilot as pilot
import verify_d6_k7_arbitrary_basis_extension_pilot as verifier


ROOT = Path(__file__).resolve().parent


class FixedCliqueBasisExtensionPilotTest(unittest.TestCase):
    def test_zero_outside_column_certificate(self) -> None:
        system = dual.basis_affine_system((0, 0), (0,))
        certificate = pilot.easy_basis_certificate(system)
        self.assertEqual(certificate["kind"], "zero_outside_column")
        pilot.verify_easy_basis_certificate(system, certificate)

    def test_identical_outside_columns_certificate(self) -> None:
        # Vertices 1 and 2 have the same column (1) to core {0}.
        system = dual.basis_affine_system((6, 1, 1), (0,))
        certificate = pilot.easy_basis_certificate(system)
        self.assertEqual(certificate["kind"], "identical_outside_columns")
        pilot.verify_easy_basis_certificate(system, certificate)
        tampered = json.loads(json.dumps(certificate))
        tampered["column"][0] = 0
        with self.assertRaises(ValueError):
            pilot.verify_easy_basis_certificate(system, tampered)

    def test_extension_universe_counts(self) -> None:
        target = {
            "fixed_maximum_clique": [0, 2],
            "graph_n": [0, 0, 0, 0, 0],
            "rank_upper": 4,
        }
        cores = pilot.extension_cores(target)
        self.assertEqual(len(cores), 1 + 3 + 3)
        self.assertEqual(cores[0], (0, 2))
        self.assertEqual(cores[-1], (0, 2, 3, 4))
        self.assertEqual(len(cores), len(set(cores)))
        self.assertTrue(all({0, 2}.issubset(core) for core in cores))

    def test_committed_report_verifies(self) -> None:
        report_path = ROOT / "d6_k7_arbitrary_basis_extension_pilot_report.json"
        if not report_path.exists():
            self.skipTest("pilot report has not been generated yet")
        result = verifier.verify_report(
            report_path, verifier.sha256(report_path)
        )
        self.assertEqual(result["status"], "PASS")


if __name__ == "__main__":
    unittest.main()
