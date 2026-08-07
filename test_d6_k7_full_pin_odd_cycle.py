#!/usr/bin/env python3
"""Controls for the generalized odd-cycle/full-pin K7 obstruction."""

import unittest

from d6_k7_full_pin_odd_cycle import (
    certificate_from_json,
    certificate_json,
    find_full_pin,
    verify_certificate,
)


def graph(count: int, edges):
    adjacency = [0] * count
    for first, second in edges:
        adjacency[first] |= 1 << second
        adjacency[second] |= 1 << first
    return tuple(adjacency)


class FullPinOddCycleTests(unittest.TestCase):
    def test_three_coordinate_full_pin(self) -> None:
        # Central 0 lies in a singleton-intersection triangle at each bit.
        edges = []
        for first, second in ((1, 2), (3, 4), (5, 6)):
            edges += [(0, first), (0, second), (first, second)]
        masks = (
            0b111,
            0b1001, 0b10001,
            0b1010, 0b100010,
            0b100100, 0b1000100,
        )
        certificate = find_full_pin(graph(7, edges), masks)
        self.assertIsNotNone(certificate)
        self.assertEqual(certificate.central, 0)
        verify_certificate(graph(7, edges), masks, certificate)
        verify_certificate(
            graph(7, edges), masks,
            certificate_from_json(certificate_json(certificate)),
        )

    def test_odd_five_cycle_pins(self) -> None:
        edges = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)]
        masks = (1, 1, 1, 1, 1)
        self.assertIsNotNone(find_full_pin(graph(5, edges), masks))

    def test_even_cycle_does_not_pin(self) -> None:
        edges = [(0, 1), (1, 2), (2, 3), (3, 0)]
        masks = (1, 1, 1, 1)
        self.assertIsNone(find_full_pin(graph(4, edges), masks))

    def test_only_some_coordinates_pinned_survives(self) -> None:
        edges = [(0, 1), (0, 2), (1, 2)]
        masks = (0b11, 0b101, 0b1001)
        self.assertIsNone(find_full_pin(graph(3, edges), masks))


if __name__ == "__main__":
    unittest.main()
