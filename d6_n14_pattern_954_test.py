#!/usr/bin/env python3
"""Positive and negative controls for the pattern-954 obstruction."""

from __future__ import annotations

import copy
import json
import unittest
from pathlib import Path

import d6_n14_pattern_954_verify as verifier


ROOT = Path(__file__).resolve().parent


def remove_edge(adjacency: list[int], first: int, second: int) -> None:
    adjacency[first] &= ~(1 << second)
    adjacency[second] &= ~(1 << first)


def add_edge(adjacency: list[int], first: int, second: int) -> None:
    adjacency[first] |= 1 << second
    adjacency[second] |= 1 << first


class AlgebraicCoreControls(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = json.loads(
            (ROOT / "d6_n14_pattern_954_input.json").read_text(encoding="utf-8")
        )
        cls.core = list(cls.payload["obstruction_core"]["adjacency"])

    def test_full_independent_verifier(self) -> None:
        report = verifier.verify(
            ROOT / "d6_n14_pattern_954_input.json",
            ROOT / "aeq_d6_n14.txt",
        )
        self.assertEqual(report["status"], "PASS")
        self.assertTrue(report["required_edge_only"])
        self.assertEqual(report["sign_cases_checked"], 4)

    def test_exact_core_hypotheses(self) -> None:
        defects = verifier.verify_obstruction_hypotheses(
            self.core, verifier.EXPECTED_SEED, verifier.EXPECTED_ROLES
        )
        self.assertEqual(
            defects,
            {
                role: sorted(mask)
                for role, mask in verifier.EXPECTED_DEFECT_BOUNDS.items()
            },
        )

    def test_source_pattern_is_an_edge_supergraph_of_core_certificate(self) -> None:
        source = verifier.source_adjacency(ROOT / "aeq_d6_n14.txt")
        original_roles = {"A": 0, "B": 3, "C": 5, "D": 8, "E": 2, "F": 6}
        defects = verifier.verify_obstruction_hypotheses(
            source, (4, 7, 9, 10, 11, 12, 13), original_roles
        )
        self.assertEqual(defects["E"], [2, 4])
        self.assertEqual(defects["F"], [2, 3])

    def test_adding_required_edges_preserves_certificate(self) -> None:
        complete = copy.copy(self.core)
        for first in range(len(complete)):
            for second in range(first):
                add_edge(complete, first, second)
        defects = verifier.verify_obstruction_hypotheses(
            complete, verifier.EXPECTED_SEED, verifier.EXPECTED_ROLES
        )
        self.assertTrue(all(not mask for mask in defects.values()))

    def test_missing_bridge_blocks_certificate(self) -> None:
        mutated = copy.copy(self.core)
        first = verifier.EXPECTED_ROLES["E"]
        second = verifier.EXPECTED_ROLES["F"]
        remove_edge(mutated, first, second)
        with self.assertRaisesRegex(AssertionError, "E-F bridge"):
            verifier.verify_obstruction_hypotheses(
                mutated, verifier.EXPECTED_SEED, verifier.EXPECTED_ROLES
            )

    def test_missing_triangle_edge_blocks_certificate(self) -> None:
        mutated = copy.copy(self.core)
        first = verifier.EXPECTED_ROLES["A"]
        second = verifier.EXPECTED_ROLES["C"]
        remove_edge(mutated, first, second)
        with self.assertRaisesRegex(AssertionError, "triangle"):
            verifier.verify_obstruction_hypotheses(
                mutated, verifier.EXPECTED_SEED, verifier.EXPECTED_ROLES
            )

    def test_expanded_defect_mask_blocks_certificate(self) -> None:
        mutated = copy.copy(self.core)
        role = verifier.EXPECTED_ROLES["E"]
        seed_coordinate_zero = verifier.EXPECTED_SEED[0]
        self.assertTrue(mutated[role] & (1 << seed_coordinate_zero))
        remove_edge(mutated, role, seed_coordinate_zero)
        with self.assertRaisesRegex(AssertionError, "unallowed defect"):
            verifier.verify_obstruction_hypotheses(
                mutated, verifier.EXPECTED_SEED, verifier.EXPECTED_ROLES
            )

    def test_all_four_sign_cases_are_exactly_impossible(self) -> None:
        records = verifier.verify_algebra()
        self.assertEqual(
            {(row["sigma"], row["tau"]) for row in records},
            {(-1, -1), (-1, 1), (1, -1), (1, 1)},
        )
        self.assertTrue(
            all(row["bridge_product_minus_one"] != [0, 1, 0, 1]
                for row in records)
        )


if __name__ == "__main__":
    unittest.main()
