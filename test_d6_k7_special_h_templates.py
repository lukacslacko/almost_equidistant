#!/usr/bin/env python3

import unittest
from itertools import combinations

from d6_k7_special_h_templates import classify_template


def graph_from_masks(masks: list[int], nonedges: set[tuple[int, int]]) -> tuple[int, ...]:
    """Build a six-clique plus five outside masks for focused controls."""

    n = 11
    adj = [0] * n
    for first, second in combinations(range(6), 2):
        adj[first] |= 1 << second
        adj[second] |= 1 << first
    for outside, mask in enumerate(masks, start=6):
        for basis in range(6):
            if mask & (1 << basis):
                adj[outside] |= 1 << basis
                adj[basis] |= 1 << outside
    normalized_nonedges = {tuple(sorted(edge)) for edge in nonedges}
    for first, second in combinations(range(5), 2):
        if (first, second) in normalized_nonedges:
            continue
        u, v = 6 + first, 6 + second
        adj[u] |= 1 << v
        adj[v] |= 1 << u
    return tuple(adj)


class SpecialHTemplateTests(unittest.TestCase):
    clique = (1 << 6) - 1

    def test_type_one(self) -> None:
        masks = [
            0b110100,  # A={2,4,5}
            0b111110,  # B={1,2,3,4,5}
            0b111000,  # C={3,4,5}
            0b010011,  # D={0,1,4}
            0b011100,  # E={2,3,4}
        ]
        witness = classify_template(
            graph_from_masks(masks, {(1, 3)}), self.clique
        )
        self.assertIsNotNone(witness)
        self.assertEqual(witness["template"], "type_I_one_nonedge")
        self.assertEqual((witness["B"], witness["D"]), (7, 9))

    def test_type_two(self) -> None:
        masks = [
            0b110101,  # A={0,2,4,5}
            0b111110,  # B={1,2,3,4,5}
            0b111001,  # C={0,3,4,5}
            0b010010,  # D={1,4}
            0b011101,  # E={0,2,3,4}
        ]
        witness = classify_template(
            graph_from_masks(masks, {(0, 1), (1, 2), (1, 4)}),
            self.clique,
        )
        self.assertIsNotNone(witness)
        self.assertEqual(
            witness["template"], "type_II_three_nonedge_star"
        )
        self.assertEqual((witness["B"], witness["D"]), (7, 9))

    def test_type_one_requires_the_nonedge(self) -> None:
        masks = [0b110100, 0b111110, 0b111000, 0b010011, 0b011100]
        self.assertIsNone(classify_template(
            graph_from_masks(masks, set()), self.clique
        ))

    def test_type_two_rejects_near_miss_mask(self) -> None:
        masks = [0b110101, 0b111110, 0b111001, 0b010010, 0b011100]
        self.assertIsNone(classify_template(
            graph_from_masks(masks, {(0, 1), (1, 2), (1, 4)}),
            self.clique,
        ))

    def test_non_six_basis_is_not_applicable(self) -> None:
        masks = [0b110100, 0b111110, 0b111000, 0b010011, 0b011100]
        graph = graph_from_masks(masks, {(1, 3)})
        self.assertIsNone(classify_template(graph, (1 << 5) - 1))


if __name__ == "__main__":
    unittest.main()
