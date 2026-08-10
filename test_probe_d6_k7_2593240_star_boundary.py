#!/usr/bin/env python3
"""Tests for the exact graph-2593240 rank-one star probe."""

from __future__ import annotations

import unittest

import probe_d6_k7_2593240_star_boundary as probe


class Graph2593240StarProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.report = probe.build_report()

    def test_complete_exact_contradiction(self) -> None:
        self.assertEqual(self.report["status"], "EXACT_CONTRADICTION")
        self.assertEqual(self.report["target"]["index"], 2_593_240)
        self.assertEqual(len(self.report["target"]["required_k7_seeds"]), 2)
        self.assertEqual(
            [row["raw_eligible_covers"] for row in self.report["quantifier"]],
            [502, 502],
        )
        self.assertEqual(
            [row["current_covers"] for row in self.report["quantifier"]],
            [[0, 8], [0, 32]],
        )
        self.assertTrue(
            all(
                row["sole_new_branch_passing_families"] == 1
                for row in self.report["quantifier"]
            )
        )

    def test_exact_patterns_and_regression_strings(self) -> None:
        strings = [
            row["pair_target_string"] for row in self.report["exact_patterns"]
        ]
        self.assertEqual(strings, ["011111111110111", "111110111111101"])
        # These two strings appeared in an early numerical handoff and were
        # independently found to have one transposed bit each.
        self.assertNotEqual(strings[0], "011111111111011")
        self.assertNotEqual(strings[1], "110111111111101")
        self.assertTrue(self.report["system_isomorphism"]["checked_exactly"])

    def test_elimination_and_star_support(self) -> None:
        algebra = self.report["algebra"]
        self.assertEqual(algebra["tetrads_checked"], 30)
        self.assertEqual(
            algebra["unique_positive_tetrad_point"],
            ["2/3", "5", "1", "1", "1", "1"],
        )
        self.assertEqual(
            algebra["star_nonzero_entries"],
            {
                "0,1": "-8/3",
                "0,2": "8",
                "0,3": "-8/3",
                "0,4": "-8/3",
                "0,5": "-28/3",
            },
        )
        self.assertIn("g_12", algebra["rank_one_contradiction"])

    def test_optional_nonedge_is_not_silently_zeroed(self) -> None:
        graph = (0, 0)
        with self.assertRaisesRegex(ValueError, "optional K entry remains"):
            probe.exact_k_entry(graph, (1, 1), 0, 1)
        self.assertEqual(
            probe.exact_k_entry(graph, (1, 2), 0, 1),
            (0, "disjoint_propagated_support_supersets"),
        )

    def test_bad_pair_target_breaks_certificate(self) -> None:
        pattern = {
            "basis_masks": probe.EXPECTED_BASIS_MASKS[0],
            "pair_targets": (1,) + probe.EXPECTED_PAIR_TARGETS[0][1:],
        }
        with self.assertRaises(ValueError):
            probe.exact_algebra_certificate(pattern)


if __name__ == "__main__":
    unittest.main()
