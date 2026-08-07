#!/usr/bin/env python3
"""Tests for the exact K7 singleton one-free-edge obstruction."""

from __future__ import annotations

import unittest
from dataclasses import replace
from fractions import Fraction

import d6_k7_one_free_edge as obstruction


def graph_from_edges(count: int, edges: list[tuple[int, int]]) -> tuple[int, ...]:
    graph = [0] * count
    for first, second in edges:
        graph[first] |= 1 << second
        graph[second] |= 1 << first
    return tuple(graph)


def synthetic_case() -> tuple[tuple[int, ...], tuple[int, ...]]:
    # Vertices 0 and 1 have masks {0,1} and {0,2}.  Their required edge has
    # singleton intersection 0.  Triangles (0,2,3) and (1,4,5) pin their
    # respective remaining coordinates 1 and 2.
    masks = (0b0000011, 0b0000101, 0b0000010, 0b0000010, 0b0000100, 0b0000100)
    edges = [
        (0, 1),
        (0, 2),
        (0, 3),
        (2, 3),
        (1, 4),
        (1, 5),
        (4, 5),
    ]
    return graph_from_edges(len(masks), edges), masks


class ExactValueTests(unittest.TestCase):
    def test_reciprocal_classification(self) -> None:
        hits = set()
        for first_size in range(1, 8):
            for second_size in range(1, 8):
                for first_sum in obstruction.sign_sums(first_size):
                    for second_sum in obstruction.sign_sums(second_size):
                        product = obstruction.quadratic_product(
                            obstruction.one_free_value(first_size, first_sum),
                            obstruction.one_free_value(second_size, second_sum),
                        )
                        if product == (Fraction(1), Fraction(0)):
                            hits.add(
                                (first_size, first_sum, second_size, second_sum)
                            )
        self.assertEqual(
            hits,
            {
                (4, -3, 7, 0),
                (4, 3, 7, 0),
                (7, 0, 4, -3),
                (7, 0, 4, 3),
            },
        )

    def test_exception_cannot_have_singleton_intersection(self) -> None:
        # Two subsets of a seven-set with sizes k,l can meet in one point only
        # when k+l<=8.  Every such size pair has no reciprocal values.
        for first_size in range(1, 8):
            for second_size in range(1, 8):
                if first_size + second_size <= 8:
                    self.assertFalse(
                        obstruction.reciprocal_possible(first_size, second_size)
                    )

    def test_one_free_values_are_nonzero(self) -> None:
        for mask_size in range(1, 8):
            for sign_sum in obstruction.sign_sums(mask_size):
                self.assertNotEqual(
                    obstruction.one_free_value(mask_size, sign_sum),
                    (Fraction(0), Fraction(0)),
                )


class CertificateTests(unittest.TestCase):
    def test_find_verify_and_round_trip(self) -> None:
        graph, masks = synthetic_case()
        certificate = obstruction.find_one_free_edge(graph, masks)
        self.assertIsNotNone(certificate)
        assert certificate is not None
        self.assertEqual((certificate.first, certificate.second), (0, 1))
        self.assertEqual(certificate.shared_coordinate, 0)
        obstruction.verify_certificate(graph, masks, certificate)
        restored = obstruction.certificate_from_json(
            obstruction.certificate_json(certificate)
        )
        self.assertEqual(restored, certificate)
        obstruction.verify_certificate(graph, masks, restored)

    def test_missing_pin_does_not_reject(self) -> None:
        graph, masks = synthetic_case()
        graph = graph_from_edges(
            len(graph),
            [
                (0, 1),
                (0, 2),
                (0, 3),
                # Delete (2,3), making the coordinate-1 component bipartite.
                (1, 4),
                (1, 5),
                (4, 5),
            ],
        )
        self.assertIsNone(obstruction.find_one_free_edge(graph, masks))

    def test_tampered_certificate_fails(self) -> None:
        graph, masks = synthetic_case()
        certificate = obstruction.find_one_free_edge(graph, masks)
        assert certificate is not None
        with self.assertRaises(ValueError):
            obstruction.verify_certificate(
                graph,
                masks,
                replace(certificate, shared_coordinate=1),
            )
        with self.assertRaises(ValueError):
            obstruction.verify_certificate(
                graph,
                masks,
                replace(certificate, first_pin_components=()),
            )


if __name__ == "__main__":
    unittest.main()
