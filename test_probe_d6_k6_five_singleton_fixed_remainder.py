#!/usr/bin/env python3
"""Exact controls for the five-singleton fixed-remainder package."""

from __future__ import annotations

import itertools
import random
import unittest
from types import SimpleNamespace

import sympy as sp

import probe_d6_k6_five_singleton_fixed_remainder as production
import verify_probe_d6_k6_five_singleton_fixed_remainder as verifier


def graph(size, edges):
    rows = [0] * size
    for left, right in edges:
        rows[left] |= 1 << right
        rows[right] |= 1 << left
    return tuple(rows)


def instance(neighbourhoods, fixed, remaining_edges=(), forbidden=None):
    if len(neighbourhoods) != 8 or len(fixed) != 8:
        raise ValueError("five basis vertices leave eight outside points")
    edges = list(itertools.combinations(range(5), 2))
    defects = [1 << coordinate for coordinate in range(5)]
    for offset, neighbourhood in enumerate(neighbourhoods):
        local = 5 + offset
        for coordinate in range(5):
            if neighbourhood & (1 << coordinate):
                edges.append((coordinate, local))
        defects.append(neighbourhood | (0 if fixed[offset] else 1 << 5))
    edges.extend((5 + left, 5 + right) for left, right in remaining_edges)
    if forbidden is not None:
        offset, coordinate = forbidden
        defects[5 + offset] &= ~(1 << coordinate)
    return (
        graph(13, edges),
        SimpleNamespace(outside=tuple(range(13)), defects=tuple(defects)),
        {coordinate: 1 << coordinate for coordinate in range(5)},
    )


class FixedRemainderTests(unittest.TestCase):
    def test_diagonal_gives_two_fixed_roots(self):
        t, m = sp.symbols("t m")
        c = 1 + m * t
        z = c - 3 * t
        diagonal = sp.expand(z * z - (c * c + 6 - 6 * m * t * t))
        self.assertEqual(sp.expand(diagonal - 3 * (3 * t * t - 2 * t - 2)), 0)
        for epsilon in (-1, 1):
            root = (1 + epsilon * sp.sqrt(7)) / 3
            self.assertEqual(sp.simplify(3 * root * root - 2 * root - 2), 0)

    def test_pair_rule_is_exact_for_all_five_bit_shapes(self):
        for left in range(32):
            for right in range(32):
                common = (left & right).bit_count()
                coefficient = left.bit_count() + right.bit_count() - 3 - 2 * common
                for left_bit in (False, True):
                    for right_bit in (False, True):
                        t = verifier.fixed_root(left_bit)
                        s = verifier.fixed_root(right_bit)
                        exact_unit = verifier.canonical(t + s + coefficient * t * s) == 0
                        discrete_unit = (
                            (left ^ right).bit_count() == 4
                            and left_bit != right_bit
                        )
                        self.assertEqual(
                            exact_unit, discrete_unit,
                            (left, right, left_bit, right_bit),
                        )

    def test_optional_nonedge_has_no_distance_equation(self):
        points = (
            production.FixedPoint(0, 0, 0b00001, 0b00001, None, None),
            production.FixedPoint(1, 1, 0b00010, 0b00010, None, None),
        )
        no_edge = graph(2, [])
        # All four root choices pass although most do not make a unit pair.
        for signs in range(4):
            reason, _ = production.branch_reason(points, signs, no_edge)
            self.assertEqual(reason, "relaxed_free_points_ignored")

        # Distinctness is mandatory even across a candidate nonedge.
        repeated = (points[0], production.FixedPoint(
            1, 1, points[0].neighbourhood, points[0].allowed, None, None
        ))
        self.assertEqual(production.branch_reason(repeated, 0, no_edge)[0], "collision")
        self.assertEqual(
            production.branch_reason(repeated, 1, no_edge)[0],
            "relaxed_free_points_ignored",
        )

    def test_producer_and_independent_extensions_agree(self):
        rng = random.Random(0xF1ED)
        for _ in range(40):
            neighbourhoods = tuple(rng.randrange(32) for _ in range(8))
            fixed = tuple(bool(rng.randrange(2)) for _ in range(8))
            edges = [
                pair for pair in itertools.combinations(range(8), 2)
                if rng.randrange(4) == 0
            ]
            adjacency, synthetic, assignments = instance(
                neighbourhoods, fixed, edges
            )
            left, left_detail = production.fixed_remainder_extension(
                adjacency, synthetic, tuple(range(5)), assignments
            )
            right, right_detail = verifier.independent_extension(
                adjacency, synthetic, tuple(range(5)), assignments
            )
            self.assertEqual(left, right)
            self.assertEqual(left_detail.get("counts"), right_detail.get("counts"))
            self.assertEqual(
                left_detail.get("branch_outcomes"),
                right_detail.get("branch_outcomes"),
            )

    def test_forbidden_required_coordinate_kills_every_sign(self):
        adjacency, synthetic, assignments = instance(
            (1, 2, 4, 8, 16, 3, 5, 9),
            (True,) * 8,
            forbidden=(0, 0),
        )
        for extension in (
            production.fixed_remainder_extension,
            verifier.independent_extension,
        ):
            feasible, detail = extension(
                adjacency, synthetic, tuple(range(5)), assignments
            )
            self.assertFalse(feasible)
            self.assertEqual(
                detail["counts"]["fixed_sign_branch_fixed_support_not_allowed"],
                1 << 8,
            )

    def test_free_points_are_explicitly_ignored(self):
        adjacency, synthetic, assignments = instance(
            (1, 2, 4, 8, 16, 3, 5, 9),
            (False,) * 8,
        )
        for extension in (
            production.fixed_remainder_extension,
            verifier.independent_extension,
        ):
            feasible, detail = extension(
                adjacency, synthetic, tuple(range(5)), assignments
            )
            self.assertTrue(feasible)
            self.assertEqual(len(detail["fixed_points"]), 0)
            self.assertEqual(len(detail["free_vertices_ignored"]), 8)
            self.assertEqual(
                detail["counts"]["fixed_sign_branch_relaxed_free_points_ignored"],
                1,
            )

    def test_frozen_boundary_and_one_rejection(self):
        left_records, left_indices = production.load_records()
        right_records, right_indices = verifier.load_records()
        self.assertEqual(left_indices, right_indices)
        self.assertEqual(len(left_indices), 249)
        self.assertEqual(
            production.base.stable_hash(left_indices), production.EXPECTED_INPUT_SHA256
        )
        target = 428414
        production.worker_initializer()
        verifier.worker_initializer()
        left = next(row for row in left_records if row["index"] == target)
        right = next(row for row in right_records if row["index"] == target)
        self.assertTrue(production.evaluate_record(left)["rejected"])
        self.assertTrue(verifier.evaluate_record(right)["rejected"])

    def test_source_independence(self):
        verifier.assert_import_independence()


if __name__ == "__main__":
    unittest.main()
