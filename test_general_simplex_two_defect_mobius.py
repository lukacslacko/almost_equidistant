#!/usr/bin/env python3
"""Exact symbolic controls for general simplex two-defect dynamics."""

from __future__ import annotations

import unittest
from fractions import Fraction
from math import isqrt


QRoot = tuple[Fraction, Fraction]
Matrix = tuple[tuple[QRoot, QRoot], tuple[QRoot, QRoot]]


def normalize(x: QRoot, radicand: int) -> QRoot:
    """Use the rational embedding when the nominal square root is rational."""
    root = isqrt(radicand)
    if root * root == radicand:
        return x[0] + root * x[1], Fraction(0)
    return x


def add(x: QRoot, y: QRoot, radicand: int) -> QRoot:
    return normalize((x[0] + y[0], x[1] + y[1]), radicand)


def negate(x: QRoot, radicand: int) -> QRoot:
    return normalize((-x[0], -x[1]), radicand)


def mul(x: QRoot, y: QRoot, radicand: int) -> QRoot:
    return normalize((
        x[0] * y[0] + radicand * x[1] * y[1],
        x[0] * y[1] + x[1] * y[0],
    ), radicand)


def matrix(dimension: int) -> Matrix:
    radicand = dimension + 1
    one = (Fraction(1), Fraction(0))
    root = normalize((Fraction(0), Fraction(1)), radicand)
    return (
        (one, negate(root, radicand)),
        (root, (Fraction(-(dimension + 2), 2), Fraction(0))),
    )


def identity() -> Matrix:
    zero = (Fraction(0), Fraction(0))
    one = (Fraction(1), Fraction(0))
    return ((one, zero), (zero, one))


def matrix_mul(first: Matrix, second: Matrix, radicand: int) -> Matrix:
    return tuple(
        tuple(
            add(
                mul(first[i][0], second[0][j], radicand),
                mul(first[i][1], second[1][j], radicand),
                radicand,
            )
            for j in range(2)
        )
        for i in range(2)
    )  # type: ignore[return-value]


def matrix_power(value: Matrix, exponent: int, radicand: int) -> Matrix:
    answer = identity()
    for _ in range(exponent):
        answer = matrix_mul(answer, value, radicand)
    return answer


def scalar_value(value: Matrix, radicand: int) -> QRoot | None:
    zero = (Fraction(0), Fraction(0))
    entries = tuple(
        tuple(normalize(value[i][j], radicand) for j in range(2))
        for i in range(2)
    )
    if entries[0][1] != zero or entries[1][0] != zero:
        return None
    return entries[0][0] if entries[0][0] == entries[1][1] else None


def rational_value(value: QRoot, radicand: int) -> Fraction:
    value = normalize(value, radicand)
    if value[1] != 0:
        raise AssertionError("entry is irrational")
    return value[0]


class GeneralMobiusControls(unittest.TestCase):
    def test_trace_determinant_and_discriminant(self):
        for dimension in range(1, 13):
            radicand = dimension + 1
            value = matrix(dimension)
            trace_exact = add(value[0][0], value[1][1], radicand)
            determinant_exact = add(
                mul(value[0][0], value[1][1], radicand),
                negate(mul(value[0][1], value[1][0], radicand),
                       radicand),
                radicand,
            )
            self.assertEqual(trace_exact, (Fraction(-dimension, 2), 0))
            self.assertEqual(determinant_exact,
                             (Fraction(dimension, 2), 0))
            trace = trace_exact[0]
            determinant = determinant_exact[0]
            discriminant = trace * trace - 4 * determinant
            self.assertEqual(
                discriminant, Fraction(dimension * (dimension - 8), 4)
            )
            # Constant in (A-r)(B-r)=d/2 agrees with the expanded
            # normalized two-defect equation.
            self.assertEqual(
                Fraction(dimension + 1) - Fraction(dimension, 2),
                Fraction(dimension + 2, 2),
            )

    def test_finite_projective_orders(self):
        cases = ((2, 3, Fraction(1)), (4, 4, Fraction(-4)),
                 (6, 6, Fraction(-27)))
        for dimension, order, expected_scalar in cases:
            value = matrix(dimension)
            radicand = dimension + 1
            for exponent in range(1, order):
                self.assertIsNone(
                    scalar_value(matrix_power(value, exponent, radicand),
                                 radicand)
                )
            self.assertEqual(
                scalar_value(matrix_power(value, order, radicand), radicand),
                (expected_scalar, Fraction(0)),
            )

    def test_odd_elliptic_cases_have_no_short_identity_power(self):
        # The algebraic-integer argument in the proof excludes every power;
        # this exact finite sweep is a regression control on the matrices.
        for dimension in (1, 3, 5, 7):
            value = matrix(dimension)
            self.assertTrue(all(
                scalar_value(
                    matrix_power(value, exponent, dimension + 1),
                    dimension + 1,
                ) is None
                for exponent in range(1, 33)
            ))

    def test_dimension_eight_is_nontrivial_parabolic(self):
        value = matrix(8)
        radicand = 9
        plus_two_identity: Matrix = (
            ((Fraction(2), Fraction(0)), (Fraction(0), Fraction(0))),
            ((Fraction(0), Fraction(0)), (Fraction(2), Fraction(0))),
        )
        nilpotent = tuple(
            tuple(add(value[i][j], plus_two_identity[i][j], radicand)
                  for j in range(2))
            for i in range(2)
        )
        self.assertEqual(matrix_mul(nilpotent, nilpotent, radicand),
                         (((Fraction(0), Fraction(0)),) * 2,) * 2)
        self.assertIsNone(scalar_value(value, radicand))

        # For T(x)=(a*x+b)/(c*x+e), the fixed polynomial is
        # c*x^2+(e-a)*x-b.  Extract it from M_8, rather than restating it.
        a = rational_value(value[0][0], radicand)
        b = rational_value(value[0][1], radicand)
        c = rational_value(value[1][0], radicand)
        e = rational_value(value[1][1], radicand)
        fixed_polynomial = (c, e - a, -b)
        self.assertEqual(fixed_polynomial,
                         (Fraction(3), Fraction(-6), Fraction(3)))
        self.assertEqual(sum(coefficient for coefficient in fixed_polynomial),
                         0)
        self.assertEqual(fixed_polynomial[1] ** 2
                         - 4 * fixed_polynomial[0] * fixed_polynomial[2], 0)
        # phi(1)=(3-5)/(1-3)=1, and conversion gives t=3,u_i=u_j=1.
        self.assertEqual(Fraction(3 - 5, 1 - 3), 1)
        self.assertEqual(Fraction(3) * 1 / 3, 1)


if __name__ == "__main__":
    unittest.main()
