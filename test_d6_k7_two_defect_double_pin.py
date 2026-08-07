#!/usr/bin/env python3
"""Controls for the exact two-defect double-pin obstruction."""

import unittest

from d6_k7_two_defect_double_pin import (
    certificate_from_json,
    certificate_json,
    find_double_pin,
    verify_certificate,
)


def graph(count: int, triangles: tuple[tuple[int, int, int], ...]):
    adjacency = [0] * count
    for triangle in triangles:
        for offset, first in enumerate(triangle):
            for second in triangle[offset + 1:]:
                adjacency[first] |= 1 << second
                adjacency[second] |= 1 << first
    return tuple(adjacency)


class DoublePinTests(unittest.TestCase):
    def test_positive_pattern(self) -> None:
        adjacency = graph(5, ((0, 1, 2), (0, 3, 4)))
        masks = (0b11, 0b101, 0b1001, 0b110, 0b10010)
        certificate = find_double_pin(adjacency, masks)
        self.assertIsNotNone(certificate)
        self.assertEqual(certificate.central, 0)
        verify_certificate(adjacency, masks, certificate)
        verify_certificate(
            adjacency, masks, certificate_from_json(certificate_json(certificate))
        )

    def test_missing_triangle_edge_survives(self) -> None:
        adjacency = list(graph(5, ((0, 1, 2), (0, 3, 4))))
        adjacency[3] &= ~(1 << 4)
        adjacency[4] &= ~(1 << 3)
        masks = (0b11, 0b101, 0b1001, 0b110, 0b10010)
        self.assertIsNone(find_double_pin(adjacency, masks))

    def test_non_singleton_pair_intersection_survives(self) -> None:
        adjacency = graph(5, ((0, 1, 2), (0, 3, 4)))
        masks = (0b11, 0b101, 0b1001, 0b111, 0b10011)
        self.assertIsNone(find_double_pin(adjacency, masks))

    def test_central_three_mask_survives(self) -> None:
        adjacency = graph(5, ((0, 1, 2), (0, 3, 4)))
        masks = (0b111, 0b1001, 0b10001, 0b110, 0b10010)
        self.assertIsNone(find_double_pin(adjacency, masks))


if __name__ == "__main__":
    unittest.main()
