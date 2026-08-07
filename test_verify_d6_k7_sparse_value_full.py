#!/usr/bin/env python3
"""Controls for the independent full sparse-value verifier."""

from __future__ import annotations

import unittest

import verify_d6_k7_sparse_value_full as verifier


def edge(first: int, second: int) -> int:
    return (1 << first) | (1 << second)


class ExactAlgebraControls(unittest.TestCase):
    def test_mobius_order_discriminants_and_one_defect_orbit(self) -> None:
        result = verifier.validate_mobius_lemmas()
        self.assertEqual(
            result["fixed_point_discriminants"],
            [[-3, 0], [-27, 0], [-108, 0], [-243, 0], [-243, 0]],
        )
        self.assertEqual(result["one_defect_orbit_sqrt7_coefficients"][-1], "0")


class IndependentValueRuleControls(unittest.TestCase):
    def test_three_cycle_fails_but_six_cycle_survives(self) -> None:
        triangle = tuple(sorted((edge(0, 1), edge(1, 2), edge(0, 2))))
        self.assertEqual(
            verifier.independent_value_failure(triangle),
            "forbidden_two_defect_cycle",
        )
        six_cycle = tuple(sorted(edge(i, (i + 1) % 6) for i in range(6)))
        self.assertIsNone(verifier.independent_value_failure(six_cycle))

    def test_one_defect_count_and_uniqueness(self) -> None:
        self.assertEqual(
            verifier.independent_value_failure((1, 2, 4)), "three_one_defects"
        )
        self.assertEqual(
            verifier.independent_value_failure((1, 1)), "duplicate_one_defect"
        )

    def test_parallel_type_isolated_rule(self) -> None:
        supports = tuple(sorted((edge(0, 1), edge(0, 1), edge(1, 2))))
        self.assertEqual(
            verifier.independent_value_failure(supports),
            "parallel_two_defect_type_touches_another",
        )

    def test_one_defect_path_and_branch_rules(self) -> None:
        joined = tuple(sorted((1 << 0, 1 << 2, edge(0, 1), edge(1, 2))))
        self.assertEqual(
            verifier.independent_value_failure(joined),
            "one_defects_joined_by_two_defects",
        )
        branched = tuple(sorted((1 << 0, edge(0, 1), edge(1, 2),
                                 edge(1, 3), edge(1, 4))))
        self.assertEqual(
            verifier.independent_value_failure(branched),
            "one_defect_component_branches",
        )

    def test_incompatible_branch_distance(self) -> None:
        supports = tuple(sorted((
            edge(0, 1), edge(0, 2), edge(0, 3),
            edge(1, 4), edge(1, 5),
        )))
        self.assertEqual(
            verifier.independent_value_failure(supports),
            "incompatible_two_defect_branch_distance",
        )


class SmallSupportSearchControls(unittest.TestCase):
    def test_disjoint_required_edge_is_rejected(self) -> None:
        result = verifier.independent_small_support_check(
            (0b10, 0b01), (0b0001, 0b0010)
        )
        self.assertFalse(result.feasible)
        self.assertEqual(dict(result.prunes), {"disjoint_required_small_support": 1})

    def test_two_bit_domains_include_singletons_and_pair(self) -> None:
        result = verifier.independent_small_support_check((0,), (0b0011,))
        self.assertTrue(result.feasible)
        self.assertEqual(result.witness, ((0, 1),))

    def test_three_forced_one_defects_are_exhaustively_rejected(self) -> None:
        result = verifier.independent_small_support_check(
            (0, 0, 0), (1 << 0, 1 << 1, 1 << 2)
        )
        self.assertFalse(result.feasible)
        self.assertGreater(result.search_nodes, 0)
        self.assertEqual(dict(result.prunes), {"three_one_defects": 1})


if __name__ == "__main__":
    unittest.main()
