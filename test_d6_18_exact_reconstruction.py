#!/usr/bin/env python3
"""Controls for the exact nonstandard d=6 18-set certificates."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import verify_d6_18_exact as verifier


ROOT = Path(__file__).resolve().parent
CERTIFICATE = ROOT / "d6_18_exact_reconstruction.json"
CORPUS = ROOT / "d6_residue_18_deletions.json"


class ExactReconstructionTests(unittest.TestCase):
    def test_independent_exact_verifier(self) -> None:
        result = verifier.verify(CERTIFICATE, CORPUS)
        self.assertEqual(result["status"], "PASS")
        checked = result["configurations_checked"]
        self.assertEqual([item["unit_edges"] for item in checked], [110, 111])
        self.assertTrue(all(item["globally_nonextendable"] for item in checked))
        self.assertTrue(all(item["switching_formula_recomputed"] for item in checked))

    def test_certificate_tamper_is_detected(self) -> None:
        certificate = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
        matrix = certificate["configurations"][0]["exact_squared_distance_matrix"]
        matrix[0][1] = matrix[1][0] = "1/3"
        with tempfile.TemporaryDirectory() as directory_name:
            path = Path(directory_name) / "tampered.json"
            path.write_text(json.dumps(certificate), encoding="utf-8")
            with self.assertRaises(AssertionError):
                verifier.verify(path, CORPUS)

    def test_three_nonisometric_switching_types_are_recorded(self) -> None:
        certificate = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
        theorem = certificate["switching_family_theorem"]
        self.assertEqual(
            theorem["clebsch_distance2_graph"]["strongly_regular_parameters"],
            [16, 5, 0, 2],
        )
        self.assertEqual(
            theorem["classification_within_this_switching_operation"][
                "number_of_types"
            ],
            3,
        )
        self.assertEqual(
            certificate["pairwise_nonisometry_witness"],
            {
                "invariant": "number of unit-distance pairs",
                "nonstandard_counts": [110, 111],
                "standard18_count": 112,
            },
        )


if __name__ == "__main__":
    unittest.main()
