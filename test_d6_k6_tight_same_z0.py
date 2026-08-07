#!/usr/bin/env python3
"""Focused controls for the strongest exact K6 same-Z0 conjunction."""

from __future__ import annotations

import unittest
from collections import Counter

import d6_k6_tight_same_z0 as production
import verify_d6_k6_tight_same_z0 as verifier
from d6_k6_support_reference import support_domains as production_support_domains


class TightSameZ0Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        records, indices, _ = production.load_input()
        cls.production_records = {int(row["index"]): row for row in records}
        cls.indices = indices
        independent_records, independent_indices = verifier.load_input()
        cls.independent_records = {
            int(row["index"]): row for row in independent_records
        }
        cls.independent_indices = independent_indices

    def test_frozen_625_parent_boundary(self) -> None:
        self.assertEqual(len(self.production_records), 625)
        self.assertEqual(self.indices, self.independent_indices)
        self.assertEqual(
            production.stable_hash(self.indices),
            production.EXPECTED_INPUT_SHA256,
        )
        self.assertNotIn(552_851, self.indices)  # repeated-arm rejection
        self.assertIn(1_442_098, self.indices)  # first new rejection
        self.assertIn(1_495_055, self.indices)  # second new rejection

    def test_fixed_rejections_and_exact_Z0_partitions(self) -> None:
        expected = {
            1_442_098: {
                "seeds_checked": 10,
                "seed_mask": 110_597,
                "seed": [0, 2, 12, 13, 15, 16],
                "failures": {
                    "tight_bipartite_system": 125,
                    "nonbipartite_two_light_ray_actual_support": 2,
                },
                "cross": {
                    "bipartite_only": 2,
                    "neither": 79,
                    "nonbipartite_only": 46,
                },
            },
            1_495_055: {
                "seeds_checked": 30,
                "seed_mask": 346_136,
                "seed": [3, 4, 11, 14, 16, 18],
                "failures": {
                    "tight_bipartite_system": 126,
                    "nonbipartite_two_light_ray_actual_support": 1,
                },
                "cross": {
                    "bipartite_only": 1,
                    "neither": 58,
                    "nonbipartite_only": 68,
                },
            },
        }
        for index, wanted in expected.items():
            with self.subTest(index=index):
                result = production.evaluate_record(
                    self.production_records[index]
                )
                self.assertTrue(result["rejected"])
                self.assertEqual(result["seeds_checked"], wanted["seeds_checked"])
                self.assertEqual(
                    result["first_impossible_seed_mask"], wanted["seed_mask"]
                )
                self.assertEqual(result["first_impossible_seed"], wanted["seed"])
                certificate = result["certificate"]
                self.assertEqual(
                    Counter(row["reason"] for row in certificate["Z0_failures"]),
                    wanted["failures"],
                )
                cross = certificate["same_Z0_cross_classification"]
                self.assertEqual(cross["histogram"], wanted["cross"])
                self.assertNotIn("both", cross["histogram"])
                self.assertGreater(cross["histogram"]["bipartite_only"], 0)
                self.assertGreater(cross["histogram"]["nonbipartite_only"], 0)
                self.assertEqual(len(cross["choices"]), 127)

    def test_independent_replay_of_every_rejected_Z0(self) -> None:
        for index in (1_442_098, 1_495_055):
            with self.subTest(index=index):
                produced = production.evaluate_record(
                    self.production_records[index]
                )
                checked = verifier.evaluate_record(
                    self.independent_records[index]
                )
                self.assertTrue(checked["rejected"])
                self.assertEqual(
                    checked["seed_mask"], produced["first_impossible_seed_mask"]
                )
                self.assertEqual(checked["seed"], produced["first_impossible_seed"])
                self.assertEqual(
                    checked["failures"], produced["certificate"]["Z0_failures"]
                )
                self.assertEqual(
                    checked["cross_classification"],
                    produced["certificate"]["same_Z0_cross_classification"],
                )

    def test_actual_support_domains_preserve_optional_zeros(self) -> None:
        nonzero_expected = {1, 2, 3, 4, 5, 6, 7}
        self.assertEqual(
            set(production_support_domains(0b111, False)),
            nonzero_expected,
        )
        self.assertEqual(
            set(verifier.base.actual_support_domains(0b111, False)),
            nonzero_expected,
        )
        # A zero-factor point needs at least three actual coordinates, but an
        # allowed fourth coordinate may still vanish.
        zero_expected = {0b0111, 0b1011, 0b1101, 0b1110, 0b1111}
        self.assertEqual(
            set(production_support_domains(0b1111, True)),
            zero_expected,
        )
        self.assertEqual(
            set(verifier.base.actual_support_domains(0b1111, True)),
            zero_expected,
        )

    def test_positive_control_and_import_independence(self) -> None:
        self.assertEqual(
            production.positive_control(), {"passed": True, "K6_seeds": 32}
        )
        self.assertEqual(
            verifier.positive_control(), {"passed": True, "K6_seeds": 32}
        )
        verifier.assert_import_independence()


if __name__ == "__main__":
    unittest.main(verbosity=2)
