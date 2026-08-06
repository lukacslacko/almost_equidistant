#!/usr/bin/env python3
"""Controls for exact propagation through K7--K6 overlaps."""

from __future__ import annotations

import copy
import unittest
from fractions import Fraction
from itertools import combinations

import d6_k7_reflection_overlap as reflection
import d6_k7_reflection_verify as independent_verifier


def graph_from_cliques(
    vertex_count: int,
    cliques: list[tuple[int, ...]],
    extra_edges: tuple[tuple[int, int], ...] = (),
) -> list[int]:
    adjacency = [0] * vertex_count
    for clique in cliques:
        for first, second in combinations(clique, 2):
            adjacency[first] |= 1 << second
            adjacency[second] |= 1 << first
    for first, second in extra_edges:
        adjacency[first] |= 1 << second
        adjacency[second] |= 1 << first
    return adjacency


class ReflectionCoordinateTests(unittest.TestCase):
    def test_exact_root_metric_and_first_reflections(self) -> None:
        coordinates = reflection.root_coordinates(tuple(range(7)))
        first = reflection.reflect_coordinate(coordinates, 6, tuple(range(6)))
        second = reflection.reflect_coordinate(
            coordinates, 5, (0, 1, 2, 3, 4, 6)
        )
        self.assertEqual(
            first,
            tuple([Fraction(1, 3)] * 6 + [Fraction(-1)]),
        )
        for facet_vertex in range(6):
            self.assertEqual(
                reflection.squared_distance(first, coordinates[facet_vertex]), 1
            )
        self.assertEqual(
            reflection.squared_distance(first, coordinates[6]), Fraction(7, 3)
        )
        self.assertEqual(
            reflection.squared_distance(first, second), Fraction(16, 9)
        )

    def test_single_K7_survives(self) -> None:
        adjacency = graph_from_cliques(7, [tuple(range(7))])
        result = reflection.analyze_graph(adjacency)
        self.assertFalse(result.rejected)
        self.assertEqual(result.k7_cliques, 1)
        self.assertEqual(result.overlap_edges, 0)

    def test_honest_facet_reflection_survives(self) -> None:
        adjacency = graph_from_cliques(
            8,
            [tuple(range(7)), (0, 1, 2, 3, 4, 5, 7)],
        )
        result = reflection.analyze_graph(adjacency)
        self.assertFalse(result.rejected)
        self.assertEqual(result.k7_cliques, 2)
        self.assertEqual(result.overlap_edges, 1)
        self.assertEqual(result.nontrivial_components, 1)
        self.assertEqual(result.maximum_component_vertices, 8)

    def test_realized_two_step_chain_survives(self) -> None:
        # Vertex 7 reflects 6 through [0,...,5].  Vertex 8 then reflects 5
        # through [0,...,4,7].  Taking only the two required K7 edge sets is
        # an exact realization, so this is a nontrivial positive control.
        adjacency = graph_from_cliques(
            9,
            [
                tuple(range(7)),
                (0, 1, 2, 3, 4, 5, 7),
                (0, 1, 2, 3, 4, 7, 8),
            ],
        )
        result = reflection.analyze_graph(adjacency)
        self.assertFalse(result.rejected)
        self.assertGreaterEqual(result.overlap_edges, 2)
        self.assertGreaterEqual(result.maximum_component_vertices, 9)


class ReflectionContradictionTests(unittest.TestCase):
    def assert_replays(self, adjacency: list[int], expected: str) -> dict:
        result = reflection.analyze_graph(adjacency)
        self.assertTrue(result.rejected)
        self.assertEqual(result.reason, expected)
        self.assertIsNotNone(result.certificate)
        reflection.verify_certificate(adjacency, result.certificate)
        independent_verifier.verify(adjacency, result.certificate)
        return result.certificate

    def test_required_edge_between_opposite_apices_is_wrong(self) -> None:
        adjacency = graph_from_cliques(
            8,
            [tuple(range(7)), (0, 1, 2, 3, 4, 5, 7)],
            ((6, 7),),
        )
        certificate = self.assert_replays(adjacency, "wrong_required_distance")
        self.assertEqual(
            certificate["terminal"]["squared_distance"],
            {"numerator": 7, "denominator": 3},
        )

    def test_three_extensions_of_one_facet_force_collision(self) -> None:
        adjacency = graph_from_cliques(
            9,
            [
                (0, 1, 2, 3, 4, 5, 6),
                (0, 1, 2, 3, 4, 5, 7),
                (0, 1, 2, 3, 4, 5, 8),
            ],
        )
        certificate = self.assert_replays(adjacency, "collision")
        self.assertNotEqual(
            certificate["terminal"]["new_vertex"],
            certificate["terminal"]["existing_vertex"],
        )

    def test_inconsistent_overlap_cycle_rejects(self) -> None:
        adjacency = graph_from_cliques(
            9,
            [
                (0, 1, 2, 3, 4, 5, 6),
                (0, 1, 2, 3, 4, 5, 7),
                (0, 1, 2, 3, 4, 7, 8),
                (0, 1, 2, 3, 4, 6, 8),
            ],
        )
        result = reflection.analyze_graph(adjacency)
        self.assertTrue(result.rejected)
        self.assertIn(
            result.reason,
            {"coordinate_conflict", "collision", "wrong_required_distance"},
        )
        reflection.verify_certificate(adjacency, result.certificate)
        independent_verifier.verify(adjacency, result.certificate)

    def test_tampered_certificate_fails(self) -> None:
        adjacency = graph_from_cliques(
            8,
            [tuple(range(7)), (0, 1, 2, 3, 4, 5, 7)],
            ((6, 7),),
        )
        certificate = reflection.analyze_graph(adjacency).certificate
        self.assertIsNotNone(certificate)
        tampered = copy.deepcopy(certificate)
        tampered["terminal"]["squared_distance"]["numerator"] += 1
        with self.assertRaises(ValueError):
            reflection.verify_certificate(adjacency, tampered)
        with self.assertRaises(ValueError):
            independent_verifier.verify(adjacency, tampered)

    def test_nonedges_are_not_required_to_be_nonunit(self) -> None:
        # Removing a required edge cannot create a rejection.  In particular,
        # the two apices in this valid reflection pair are a graph nonedge;
        # the checker only uses their actual nonunit distance to accept them.
        adjacency = graph_from_cliques(
            8,
            [tuple(range(7)), (0, 1, 2, 3, 4, 5, 7)],
        )
        result = reflection.analyze_graph(adjacency)
        self.assertFalse(result.rejected)
        self.assertFalse(adjacency[6] & (1 << 7))


if __name__ == "__main__":
    unittest.main()
