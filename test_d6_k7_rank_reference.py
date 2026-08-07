#!/usr/bin/env python3
"""Controls for the independent K7 existential-cover rank reference."""

from __future__ import annotations

import unittest
from fractions import Fraction
from itertools import combinations, product

from d6_k7_rank_reference import (
    SUPPORT_DOMAINS,
    CliqueStructureSolver,
    SupportSolver,
    ZeroForcingSolver,
    basis_kernel_compatibility,
    component_inertia_nullity_cap,
    component_inertia_nullity_caps,
    complement_graph,
    components,
    eligible_covers,
    direct_cover_cap_failure,
    matching_size,
    rational_nullspace,
    rank_upper_bounds,
    saturating_clique_mask_compatibility,
    support_family_valid,
)


def graph_from_edges(n: int, edges: list[tuple[int, int]]) -> tuple[int, ...]:
    adj = [0] * n
    for u, v in edges:
        adj[u] |= 1 << v
        adj[v] |= 1 << u
    return tuple(adj)


def graph_with_clique_masks(
    clique_size: int,
    masks: list[int],
    remainder_edges: list[tuple[int, int]] | None = None,
) -> tuple[int, ...]:
    edges = list(combinations(range(clique_size), 2))
    for local, mask in enumerate(masks):
        vertex = clique_size + local
        edges.extend(
            (vertex, coordinate)
            for coordinate in range(clique_size)
            if mask & (1 << coordinate)
        )
    edges.extend(
        (clique_size + first, clique_size + second)
        for first, second in (remainder_edges or [])
    )
    return graph_from_edges(clique_size + len(masks), edges)


class CoverControls(unittest.TestCase):
    def test_all_covers_include_strict_supercovers(self) -> None:
        # One edge 0--1 and an eligible isolated vertex 2.
        ladj = graph_from_edges(3, [(0, 1)])
        observed = eligible_covers(ladj, 0b111, cap=2)
        self.assertEqual(
            observed,
            [0b001, 0b010, 0b011, 0b101, 0b110],
        )
        self.assertIn(0b101, observed)  # nonminimal cover with isolated 2


class SupportControls(unittest.TestCase):
    def test_dense_householder_support_passes(self) -> None:
        self.assertTrue(support_family_valid([0x7F] * 7))
        self.assertIsNotNone(SupportSolver().solve([0x7F] * 7))

    def test_single_coordinate_intersection_fails(self) -> None:
        domains = [0b0000111, 0b0011001]
        self.assertIsNone(SupportSolver().solve(domains))

    def test_matching_failure(self) -> None:
        self.assertEqual(matching_size([0b111] * 4), 3)
        self.assertIsNone(SupportSolver().solve([0b111] * 4))

    def test_row_support_size_two_fails(self) -> None:
        supports = [
            0b0011011,
            0b0011010,
            0b0111110,
            0b1011101,
            0b1011010,
            0b0011101,
            0b0111100,
        ]
        self.assertFalse(support_family_valid(supports))

    def test_row_single_intersection_fails(self) -> None:
        supports = [
            0b1010011,
            0b0111110,
            0b1101011,
            0b1101110,
            0b0101111,
            0b1100011,
            0b1011001,
        ]
        # These pass all column-side checks, isolating the row condition.
        self.assertTrue(
            all((a & b).bit_count() != 1 for a, b in combinations(supports, 2))
        )
        self.assertEqual(matching_size(supports), 7)
        self.assertFalse(support_family_valid(supports))

    def test_solver_agrees_with_direct_product_search(self) -> None:
        cases = [
            (0b001111, 0b011110, 0b111100),
            (0b000111, 0b001110, 0b011100, 0b111000),
            (0b001111, 0b001111, 0b110011, 0b110011),
            (0b000111, 0b000111, 0b000111, 0b000111),
        ]
        solver = SupportSolver()
        for allowed in cases:
            brute = any(
                support_family_valid(chosen)
                for chosen in product(*(SUPPORT_DOMAINS[d] for d in allowed))
            )
            self.assertEqual(solver.solve(allowed) is not None, brute)


class ZeroForcingControls(unittest.TestCase):
    def setUp(self) -> None:
        self.solver = ZeroForcingSolver()

    def test_known_families(self) -> None:
        for n in range(1, 8):
            empty = (0,) * n
            complete = complement_graph(empty)
            path = graph_from_edges(n, [(i, i + 1) for i in range(n - 1)])
            self.assertEqual(self.solver.solve(empty).number, n)
            self.assertEqual(
                self.solver.solve(complete).number,
                1 if n == 1 else n - 1,
            )
            self.assertEqual(self.solver.solve(path).number, 1)
        for n in range(3, 8):
            cycle = graph_from_edges(
                n, [(i, (i + 1) % n) for i in range(n)]
            )
            self.assertEqual(self.solver.solve(cycle).number, 2)

    def test_disjoint_union_additivity(self) -> None:
        graph = graph_from_edges(4, [(0, 1), (1, 2)])  # P3 plus isolate
        self.assertEqual(self.solver.solve(graph).number, 2)
        self.assertEqual(sorted(map(len, components(graph))), [1, 3])

    def test_component_inertia_bounds(self) -> None:
        isolates = (0,) * 5
        path4 = graph_from_edges(4, [(0, 1), (1, 2), (2, 3)])
        two_cycles = graph_from_edges(
            8,
            [
                (0, 1), (1, 2), (2, 3), (3, 0),
                (4, 5), (5, 6), (6, 7), (7, 4),
            ],
        )
        self.assertEqual(
            component_inertia_nullity_cap(isolates, self.solver), 0
        )
        self.assertEqual(
            component_inertia_nullity_cap(path4, self.solver), 1
        )
        self.assertEqual(self.solver.solve(two_cycles).number, 4)
        self.assertEqual(
            component_inertia_nullity_cap(two_cycles, self.solver), 3
        )
        caps, checks = component_inertia_nullity_caps(two_cycles, self.solver)
        self.assertEqual(
            caps,
            {"trivial": 4, "tier1": 4, "tier2": 3, "tier3": 3, "exact": 3},
        )
        self.assertEqual(checks["trivial"], 0)
        self.assertGreater(checks["tier2"], checks["tier1"])

    def test_cover_consumes_support_dimensions(self) -> None:
        self.assertEqual(rank_upper_bounds(5, 7, 7, 7), (0, 1))
        self.assertEqual(rank_upper_bounds(9, 5, 7, 3), (4, 5))
        self.assertEqual(rank_upper_bounds(9, 5, 6, 3), (3, 4))
        self.assertEqual(rank_upper_bounds(9, 5, 2, 3), (0, 1))
        self.assertEqual(rank_upper_bounds(12, 7, 7, 0), (7, 8))
        self.assertEqual(
            [direct_cover_cap_failure(size) for size in range(8)],
            [None, None, None, None, "cap4_to3", "cap5_to4", "cap6_to5", "cap7_to6"],
        )


class CliqueStructureControls(unittest.TestCase):
    def test_positive_definite_clique_and_perpendicular_degree(self) -> None:
        solver = CliqueStructureSolver()
        complete4 = graph_from_edges(4, list(combinations(range(4), 2)))
        clique_failure = solver.solve(complete4, 3)
        self.assertEqual(clique_failure.clique_number, 4)
        self.assertGreater(clique_failure.clique_number, 3)

        # G is K3 plus an isolate, so F is the four-vertex star.
        complement_star = graph_from_edges(
            4, list(combinations((1, 2, 3), 2))
        )
        degree_failure = solver.solve(complement_star, 3)
        self.assertEqual(degree_failure.clique_number, 3)
        self.assertEqual(degree_failure.f_maximum_degree, 3)
        self.assertGreater(degree_failure.f_maximum_degree, 3 - 1)

    def test_exact_rational_nullspace(self) -> None:
        basis = rational_nullspace([[1, 1, 0], [0, 1, 1]], 3)
        self.assertEqual(basis, ((Fraction(1), Fraction(-1), Fraction(1)),))

    def test_basis_kernel_positive_and_negative_controls(self) -> None:
        # P has columns 100, 011, 010, 101 and kernel (1,1,-1,-1).
        # A_R=K2,2 sends that vector to -2 times itself, so D=2 works.
        positive = graph_with_clique_masks(
            3,
            [0b001, 0b110, 0b010, 0b101],
            [(0, 2), (0, 3), (1, 2), (1, 3)],
        )
        compatible, witness, nullity = basis_kernel_compatibility(
            positive, 0b111
        )
        self.assertTrue(compatible)
        self.assertIsNone(witness)
        self.assertEqual(nullity, 1)

        # A bow tie has duplicate P columns and forces D=1, not D>1.
        bow_tie = (30, 17, 9, 5, 3)
        compatible, witness, nullity = basis_kernel_compatibility(
            bow_tie, (1 << 0) | (1 << 1) | (1 << 4)
        )
        self.assertFalse(compatible)
        self.assertEqual(nullity, 1)
        self.assertIsNotNone(witness)
        self.assertEqual(
            witness["failure_kind"], "diagonal_not_greater_than_one"
        )
        self.assertEqual(witness["forced_diagonal"], [1, 1])

    def test_each_saturating_clique_mask_failure(self) -> None:
        cases = [
            ([0b000], [], "empty_basis_neighbour_mask"),
            ([0b111], [], "full_basis_neighbour_mask"),
            (
                [0b011, 0b011],
                [],
                "duplicate_basis_neighbour_masks",
            ),
            ([0b001, 0b010], [], "disjoint_basis_neighbour_masks"),
            (
                [0b001, 0b011],
                [],
                "nonedge_comparable_basis_masks",
            ),
            (
                [0b011, 0b101],
                [(0, 1)],
                "edge_basis_masks_cover_clique",
            ),
        ]
        for masks, edges, failure_kind in cases:
            graph = graph_with_clique_masks(3, masks, edges)
            compatible, witness = saturating_clique_mask_compatibility(
                graph, 0b111
            )
            self.assertFalse(compatible)
            self.assertIsNotNone(witness)
            self.assertEqual(witness["failure_kind"], failure_kind)

    def test_saturating_clique_mask_positive_control(self) -> None:
        graph = graph_with_clique_masks(
            4, [0b0011, 0b0101], [(0, 1)]
        )
        compatible, witness = saturating_clique_mask_compatibility(
            graph, 0b1111
        )
        self.assertTrue(compatible)
        self.assertIsNone(witness)


if __name__ == "__main__":
    unittest.main()
