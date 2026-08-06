#!/usr/bin/env python3
"""Unit controls for the independent full K7 support verifier."""

from __future__ import annotations

import unittest

import verify_d6_k7_support_full as verifier


class IndependentSupportKernelTests(unittest.TestCase):
    def test_domains_are_all_subsets_of_size_at_least_three(self) -> None:
        self.assertEqual(verifier.independent_support_domains(0b111), (0b111,))
        domain = verifier.independent_support_domains(0b1111)
        self.assertEqual(len(domain), 5)
        self.assertEqual(domain[-1], 0b1111)

    def test_injective_representatives(self) -> None:
        self.assertTrue(verifier.has_injective_representatives((0b111, 0b111)))
        self.assertFalse(verifier.has_injective_representatives((1, 1)))

    def test_labeled_enumeration_keeps_equal_domain_labels(self) -> None:
        families = tuple(
            verifier.independent_labeled_support_families((0b1111, 0b1111))
        )
        self.assertEqual(len(families), 25)
        self.assertIn((0b0111, 0b1011), families)
        self.assertIn((0b1011, 0b0111), families)
        self.assertGreater(len(families), len(set(tuple(sorted(x)) for x in families)))
        self.assertTrue(all(verifier.independently_valid_support_family(x)
                            for x in families))

    def test_pair_intersection_one_is_invalid(self) -> None:
        self.assertFalse(
            verifier.independently_valid_support_family((0b00111, 0b11001))
        )
        self.assertTrue(
            verifier.independently_valid_support_family((0b01101, 0b11001))
        )

    def test_singleton_closure_cascades_to_empty(self) -> None:
        supports = (0b0011001, 0b0011010)
        closed, deletions = verifier.independent_singleton_closure(
            0b0000011, supports
        )
        self.assertEqual(closed, 0)
        self.assertEqual(deletions, 2)

    def test_empty_and_disjoint_failures_are_distinguished(self) -> None:
        empty_supports = (0b0011001, 0b0011010)
        failure, _ = verifier.independently_propagate_family(
            (0,), (0b11,), empty_supports
        )
        self.assertEqual(failure, "empty_propagated_mask")

        # A required edge remains, but propagation makes its endpoint masks
        # disjoint without emptying either one.
        graph_n = (0b10, 0b01)
        failure, _ = verifier.independently_propagate_family(
            graph_n,
            (0b0011, 0b0110),
            (0b0001110,),
        )
        self.assertEqual(failure, "disjoint_required_edge")

    def test_surviving_family_is_not_a_rejection(self) -> None:
        failure, deletions = verifier.independently_propagate_family(
            (0b10, 0b01), (0b0111, 0b0111), (0b1110000,)
        )
        self.assertIsNone(failure)
        self.assertEqual(deletions, 0)


if __name__ == "__main__":
    unittest.main()
