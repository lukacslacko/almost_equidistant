#!/usr/bin/env python3
"""Focused exact controls for the K7 sparse-value refinement."""

from __future__ import annotations

import unittest
from fractions import Fraction

import d6_k7_small_support_value as value


Q7 = tuple[Fraction, Fraction]


def q7_add(x: Q7, y: Q7) -> Q7:
    return x[0] + y[0], x[1] + y[1]


def q7_neg(x: Q7) -> Q7:
    return -x[0], -x[1]


def q7_mul(x: Q7, y: Q7) -> Q7:
    return x[0] * y[0] + 7 * x[1] * y[1], x[0] * y[1] + x[1] * y[0]


def q7_div(x: Q7, y: Q7) -> Q7:
    norm = y[0] * y[0] - 7 * y[1] * y[1]
    if not norm:
        raise ZeroDivisionError
    numerator = q7_mul(x, (y[0], -y[1]))
    return numerator[0] / norm, numerator[1] / norm


def matrix_mul(first, second):
    return tuple(
        tuple(
            q7_add(
                q7_mul(first[i][0], second[0][j]),
                q7_mul(first[i][1], second[1][j]),
            )
            for j in range(2)
        )
        for i in range(2)
    )


def matrix_power(matrix, exponent):
    answer = (((Fraction(1), Fraction(0)), (Fraction(0), Fraction(0))),
              ((Fraction(0), Fraction(0)), (Fraction(1), Fraction(0))))
    for _ in range(exponent):
        answer = matrix_mul(answer, matrix)
    return answer


def mobius(matrix, x: Q7) -> Q7:
    numerator = q7_add(q7_mul(matrix[0][0], x), matrix[0][1])
    denominator = q7_add(q7_mul(matrix[1][0], x), matrix[1][1])
    return q7_div(numerator, denominator)


def intersection_graph(masks):
    adjacency = [0] * len(masks)
    for first in range(len(masks)):
        for second in range(first):
            if masks[first] & masks[second]:
                adjacency[first] |= 1 << second
                adjacency[second] |= 1 << first
    return tuple(adjacency)


class ExactMobiusControls(unittest.TestCase):
    def test_projective_order_and_no_real_short_fixed_point(self):
        # T(x)=(x-rho)/(rho*x-4), rho^2=7.
        zero = (Fraction(0), Fraction(0))
        one = (Fraction(1), Fraction(0))
        rho = (Fraction(0), Fraction(1))
        matrix = ((one, q7_neg(rho)), (rho, (Fraction(-4), Fraction(0))))
        sixth = matrix_power(matrix, 6)
        self.assertEqual(
            sixth,
            (((Fraction(-27), Fraction(0)), zero),
             (zero, (Fraction(-27), Fraction(0)))),
        )
        expected_discriminants = (-3, -27, -108, -243, -243)
        for exponent, expected in enumerate(expected_discriminants, 1):
            power = matrix_power(matrix, exponent)
            a, b = power[0]
            c, d = power[1]
            difference = q7_add(d, q7_neg(a))
            discriminant = q7_add(
                q7_mul(difference, difference),
                (4 * q7_mul(b, c)[0], 4 * q7_mul(b, c)[1]),
            )
            self.assertEqual(discriminant, (Fraction(expected), Fraction(0)))

    def test_one_defect_path_orbit(self):
        zero = (Fraction(0), Fraction(0))
        one = (Fraction(1), Fraction(0))
        rho = (Fraction(0), Fraction(1))
        matrix = ((one, q7_neg(rho)), (rho, (Fraction(-4), Fraction(0))))
        expected = (
            (Fraction(0), Fraction(1, 4)),
            (Fraction(0), Fraction(1, 3)),
            (Fraction(0), Fraction(2, 5)),
            (Fraction(0), Fraction(1, 2)),
            (Fraction(0), Fraction(1)),
            zero,
            (Fraction(0), Fraction(1, 4)),
        )
        current = expected[0]
        observed = [current]
        for _ in range(6):
            current = mobius(matrix, current)
            observed.append(current)
        self.assertEqual(tuple(observed), expected)

    def test_full_branch_sign_orbits(self):
        one = (Fraction(1), Fraction(0))
        rho = (Fraction(0), Fraction(1))
        matrix = ((one, q7_neg(rho)), (rho, (Fraction(-4), Fraction(0))))
        expected = {
            one: (
                one,
                (Fraction(1, 3), Fraction(1, 3)),
                (Fraction(1, 2), Fraction(1, 2)),
                (Fraction(-1), Fraction(0)),
                (Fraction(-1, 3), Fraction(1, 3)),
                (Fraction(-1, 2), Fraction(1, 2)),
                one,
            ),
            (Fraction(-1), Fraction(0)): (
                (Fraction(-1), Fraction(0)),
                (Fraction(-1, 3), Fraction(1, 3)),
                (Fraction(-1, 2), Fraction(1, 2)),
                one,
                (Fraction(1, 3), Fraction(1, 3)),
                (Fraction(1, 2), Fraction(1, 2)),
                (Fraction(-1), Fraction(0)),
            ),
        }
        signs = {one, (Fraction(-1), Fraction(0))}
        for start, wanted in expected.items():
            current = start
            observed = [current]
            for _ in range(6):
                current = mobius(matrix, current)
                observed.append(current)
            self.assertEqual(tuple(observed), wanted)
            self.assertEqual(
                [step for step, point in enumerate(observed) if point in signs],
                [0, 3, 6],
            )


class ValueRuleControls(unittest.TestCase):
    def test_cycle_lengths(self):
        triangle = {(0, 1), (1, 2), (0, 2)}
        six_cycle = {(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (0, 5)}
        self.assertEqual(value.simple_cycle_lengths(triangle), {3})
        self.assertEqual(value.simple_cycle_lengths(six_cycle), {6})
        self.assertEqual(value.value_failure((3, 6, 5)),
                         "forbidden_two_defect_cycle")
        self.assertIsNone(value.value_failure((3, 6, 12, 24, 48, 33)))

    def test_known_one_and_two_defect_caps(self):
        self.assertEqual(value.value_failure((1, 2, 4)),
                         "three_one_defects")
        self.assertEqual(value.value_failure((1, 1)),
                         "duplicate_one_defect")
        self.assertEqual(value.value_failure((3, 3, 3)),
                         "three_same_two_defects")
        self.assertEqual(value.value_failure((1, 3, 5)),
                         "two_defects_repeat_at_one_defect")

    def test_one_defect_path_and_branch(self):
        self.assertEqual(value.value_failure((1, 3, 6, 4)),
                         "one_defects_joined_by_two_defects")
        self.assertEqual(value.value_failure((1, 3, 6, 10)),
                         "one_defect_component_branches")

    def test_parallel_type_cannot_touch_another_type(self):
        self.assertIsNone(value.value_failure((3, 3)))
        self.assertEqual(
            value.value_failure((3, 3, 5)),
            "parallel_two_defect_type_touches_another",
        )

    def test_two_nearby_branch_coordinates_are_impossible(self):
        # Branches 0 and 1 are joined by a two-edge path through coordinate 2.
        # Each has two additional leaves, using all seven seed coordinates.
        supports = (5, 9, 17, 6, 34, 66)
        self.assertEqual(
            value.value_failure(supports),
            "incompatible_two_defect_branch_distance",
        )

    def test_mask_csp_separates_support_and_value(self):
        masks = (3, 6, 5)
        graph = intersection_graph(masks)
        support_only = value.check_small_support_masks(
            graph, masks, apply_value_constraints=False
        )
        exact_value = value.check_small_support_masks(graph, masks)
        self.assertTrue(support_only.feasible)
        self.assertFalse(exact_value.feasible)
        self.assertIn(
            "forbidden_two_defect_cycle",
            dict(exact_value.failures),
        )

    def test_six_cycle_remains_a_survivor(self):
        masks = (3, 6, 12, 24, 48, 33)
        result = value.check_small_support_masks(
            intersection_graph(masks), masks
        )
        self.assertTrue(result.feasible)


if __name__ == "__main__":
    unittest.main()
