#!/usr/bin/env python3
"""Tests for correlated pinned products in the K7 one-free-edge lemma."""

from __future__ import annotations

import unittest
from dataclasses import replace
from fractions import Fraction
from itertools import product

import d6_k7_correlated_one_free_edge as correlated
import d6_k7_one_free_edge as base


def pinned_case(
    first_mask: int,
    second_mask: int,
    *,
    bridge_shared_coordinates: tuple[int, ...] = (),
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """Make independent pinning triangles, optionally joined by bridges."""

    free_bit = 1
    if not first_mask & free_bit or not second_mask & free_bit:
        raise ValueError("synthetic masks must share free coordinate zero")
    masks = [first_mask, second_mask]
    edges = [(0, 1)]
    auxiliaries = {}
    for endpoint, mask in enumerate((first_mask, second_mask)):
        for coordinate in base.bits(mask ^ free_bit):
            first_auxiliary = len(masks)
            masks.extend((1 << coordinate, 1 << coordinate))
            second_auxiliary = first_auxiliary + 1
            auxiliaries[(endpoint, coordinate)] = (
                first_auxiliary,
                second_auxiliary,
            )
            edges.extend(
                (
                    (endpoint, first_auxiliary),
                    (endpoint, second_auxiliary),
                    (first_auxiliary, second_auxiliary),
                )
            )
    for coordinate in bridge_shared_coordinates:
        edges.append(
            (
                auxiliaries[(0, coordinate)][0],
                auxiliaries[(1, coordinate)][0],
            )
        )
    graph = [0] * len(masks)
    for first, second in edges:
        graph[first] |= 1 << second
        graph[second] |= 1 << first
    return tuple(graph), tuple(masks)


class ExactTableTests(unittest.TestCase):
    def test_complete_integer_product_table(self) -> None:
        hits = set()
        for first_size in range(1, 8):
            for first_sum in base.sign_sums(first_size):
                for second_size in range(first_size, 8):
                    for second_sum in base.sign_sums(second_size):
                        rational, irrational = base.quadratic_product(
                            base.one_free_value(first_size, first_sum),
                            base.one_free_value(second_size, second_sum),
                        )
                        if not irrational and rational.denominator == 1:
                            hits.add(
                                (
                                    first_size,
                                    first_sum,
                                    second_size,
                                    second_sum,
                                    int(rational),
                                )
                            )
        self.assertEqual(
            hits,
            {
                (1, 0, 4, -3, 4),
                (1, 0, 4, 3, 4),
                (3, -2, 3, 2, -1),
                (3, 2, 3, -2, -1),
                (3, 0, 4, -3, 3),
                (3, 0, 4, 3, 3),
                (4, -3, 4, -3, 7),
                (4, -3, 4, 3, 7),
                (4, 3, 4, -3, 7),
                (4, 3, 4, 3, 7),
                (4, -3, 5, 0, 2),
                (4, 3, 5, 0, 2),
                (4, -3, 7, 0, 1),
                (4, 3, 7, 0, 1),
                (5, -2, 5, -2, 4),
                (5, -2, 5, 2, -4),
                (5, 2, 5, -2, -4),
                (5, 2, 5, 2, 4),
            },
        )

    def test_coarse_independent_sign_compatibility_types(self) -> None:
        compatible = set()
        for first_size in range(1, 8):
            for second_size in range(first_size, 8):
                for shared_pinned_count in range(
                    1, min(first_size, second_size)
                ):
                    # The masks share the free coordinate and m pinned ones.
                    if first_size + second_size - (shared_pinned_count + 1) > 7:
                        continue
                    for first_signs in product(
                        (-1, 1), repeat=first_size - 1
                    ):
                        first_sum = sum(first_signs)
                        for second_signs in product(
                            (-1, 1), repeat=second_size - 1
                        ):
                            second_sum = sum(second_signs)
                            shared_sum = sum(
                                first_signs[position] * second_signs[position]
                                for position in range(shared_pinned_count)
                            )
                            value = base.quadratic_product(
                                base.one_free_value(first_size, first_sum),
                                base.one_free_value(second_size, second_sum),
                            )
                            if (
                                value[0] + shared_sum,
                                value[1],
                            ) == (Fraction(1), Fraction(0)):
                                compatible.add(
                                    (
                                        first_size,
                                        second_size,
                                        shared_pinned_count,
                                        first_sum,
                                        second_sum,
                                        shared_sum,
                                    )
                                )
        self.assertEqual(
            compatible,
            {
                (4, 5, 1, -3, 0, -1),
                (4, 5, 1, 3, 0, -1),
                (4, 5, 3, -3, 0, -1),
                (4, 5, 3, 3, 0, -1),
            },
        )


class CorrelatedCertificateTests(unittest.TestCase):
    def test_incompatible_independent_shared_signs_reject(self) -> None:
        graph, masks = pinned_case(0b11, 0b11)
        certificate = correlated.find_correlated_one_free_edge(graph, masks)
        self.assertIsNotNone(certificate)
        assert certificate is not None
        correlated.verify_certificate(graph, masks, certificate)
        restored = correlated.certificate_from_json(
            correlated.certificate_json(certificate)
        )
        self.assertEqual(restored, certificate)
        correlated.verify_certificate(graph, masks, restored)

    def test_four_five_exception_has_a_correlated_witness(self) -> None:
        # Sizes 4 and 5, with one pinned shared coordinate, admit P=+/-3,
        # Q=0 and shared product -1: F(4,P)F(5,0)-1 = 1.
        graph, masks = pinned_case(0b0001111, 0b1110011)
        self.assertIsNone(correlated.find_correlated_one_free_edge(graph, masks))

    def test_component_bridge_removes_four_five_witness(self) -> None:
        # Joining the two coordinate-1 pinning triangles forces that shared
        # sign product to +1, eliminating the only compatible assignment.
        graph, masks = pinned_case(
            0b0001111,
            0b1110011,
            bridge_shared_coordinates=(1,),
        )
        certificate = correlated.find_correlated_one_free_edge(graph, masks)
        self.assertIsNotNone(certificate)
        assert certificate is not None
        correlated.verify_certificate(graph, masks, certificate)

    def test_tampering_fails(self) -> None:
        graph, masks = pinned_case(0b11, 0b11)
        certificate = correlated.find_correlated_one_free_edge(graph, masks)
        assert certificate is not None
        with self.assertRaises(ValueError):
            correlated.verify_certificate(
                graph,
                masks,
                replace(certificate, assignments_checked=1),
            )
        with self.assertRaises(ValueError):
            correlated.verify_certificate(
                graph,
                masks,
                replace(certificate, first_pin_components=()),
            )


if __name__ == "__main__":
    unittest.main()
