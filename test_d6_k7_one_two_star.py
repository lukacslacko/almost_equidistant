#!/usr/bin/env python3
"""Exact controls for the K7 one-free-neighbour/two-free-center star."""

from __future__ import annotations

import copy
import json
import unittest
from fractions import Fraction
from pathlib import Path

import d6_k7_one_two_star as star
import d6_k7_rank_reference as reference


def fixture(index: int, seed: list[int], masks: list[int]) -> tuple[tuple[int, ...], list[int]]:
    manifest = json.loads(
        (Path(__file__).resolve().parent / "d6_current_residue_manifest_v6.json").read_text(
            encoding="utf-8"
        )
    )
    graph = next(
        item for item in manifest["classes"]["K7"]["graphs"]
        if int(item["index"]) == index
    )
    adjacency = tuple(map(int, graph["adjacency"]))
    seed_mask = sum(1 << vertex for vertex in seed)
    rebuilt_seed, outside, _defects, _ladj, _eligible = reference.seed_instance(
        adjacency, seed_mask
    )
    if rebuilt_seed != seed:
        raise AssertionError("fixture seed differs")
    return reference.induced_graph(adjacency, outside), masks


class QuadraticFieldControls(unittest.TestCase):
    def test_positive_embedding_sign_is_exact(self) -> None:
        values = (
            ((Fraction(-3), Fraction(1)), -1),
            ((Fraction(-2), Fraction(1)), 1),
            ((Fraction(3), Fraction(-1)), 1),
            ((Fraction(2), Fraction(-1)), -1),
            ((Fraction(0), Fraction(1)), 1),
            ((Fraction(0), Fraction(-1)), -1),
        )
        for value, expected in values:
            self.assertEqual(star.qsign(value), expected)

    def test_two_free_factorization_has_no_hidden_zero_denominator(self) -> None:
        for size in range(2, 8):
            for pinned_sum in range(-(size - 2), size - 1, 2):
                center, radius = star.diagonal_two_free(size, pinned_sum)
                self.assertNotEqual(radius, star.ZERO)
                first = star.q(2)
                if first == center:
                    first = star.q(3)
                second = star.qadd(
                    center, star.qdiv(radius, star.qsub(first, center))
                )
                self.assertEqual(
                    star.qmul(
                        star.qsub(first, center),
                        star.qsub(second, center),
                    ),
                    radius,
                )


class LineHyperbolaControls(unittest.TestCase):
    def setUp(self) -> None:
        self.center = star.ZERO
        self.radius = star.ONE

    def test_secant_tangent_and_negative_discriminant(self) -> None:
        secant = star.lines_hyperbola_decision(
            ((star.ONE, star.ONE, star.q(3)),), self.center, self.radius
        )
        tangent = star.lines_hyperbola_decision(
            ((star.ONE, star.ONE, star.q(2)),), self.center, self.radius
        )
        negative = star.lines_hyperbola_decision(
            ((star.ONE, star.ONE, star.ZERO),), self.center, self.radius
        )
        self.assertTrue(secant.feasible)
        self.assertEqual(secant.reason, "rank1_secant")
        self.assertTrue(tangent.feasible)
        self.assertEqual(tangent.reason, "rank1_tangent")
        self.assertEqual(star.qsign(tangent.discriminant), 0)
        self.assertFalse(negative.feasible)
        self.assertEqual(negative.reason, "rank1_negative_discriminant")

    def test_coincident_and_parallel_lines(self) -> None:
        coincident = star.lines_hyperbola_decision(
            (
                (star.ONE, star.ONE, star.q(3)),
                (star.q(2), star.q(2), star.q(6)),
            ),
            self.center,
            self.radius,
        )
        parallel = star.lines_hyperbola_decision(
            (
                (star.ONE, star.ZERO, star.ZERO),
                (star.ONE, star.ZERO, star.ONE),
            ),
            self.center,
            self.radius,
        )
        self.assertTrue(coincident.feasible)
        self.assertEqual(coincident.rank, 1)
        self.assertFalse(parallel.feasible)
        self.assertEqual(parallel.reason, "parallel_inconsistent")

    def test_rank_two_on_and_off_diagonal(self) -> None:
        feasible = star.lines_hyperbola_decision(
            (
                (star.ONE, star.ZERO, star.ONE),
                (star.ZERO, star.ONE, star.ONE),
            ),
            self.center,
            self.radius,
        )
        mismatch = star.lines_hyperbola_decision(
            (
                (star.ONE, star.ZERO, star.ONE),
                (star.ZERO, star.ONE, star.q(2)),
            ),
            self.center,
            self.radius,
        )
        self.assertTrue(feasible.feasible)
        self.assertEqual(feasible.reason, "rank2_on_diagonal")
        self.assertFalse(mismatch.feasible)
        self.assertEqual(mismatch.reason, "rank2_diagonal_mismatch")

    def test_zero_line_and_rank_two_inconsistency(self) -> None:
        consistent_zero = star.lines_hyperbola_decision(
            ((star.ZERO, star.ZERO, star.ZERO),), self.center, self.radius
        )
        inconsistent_zero = star.lines_hyperbola_decision(
            ((star.ZERO, star.ZERO, star.ONE),), self.center, self.radius
        )
        inconsistent_rank_two = star.lines_hyperbola_decision(
            (
                (star.ONE, star.ZERO, star.ZERO),
                (star.ZERO, star.ONE, star.ZERO),
                (star.ONE, star.ZERO, star.ONE),
            ),
            self.center,
            self.radius,
        )
        self.assertTrue(consistent_zero.feasible)
        self.assertEqual(consistent_zero.rank, 0)
        self.assertFalse(inconsistent_zero.feasible)
        self.assertEqual(inconsistent_zero.reason, "zero_line_inconsistent")
        self.assertFalse(inconsistent_rank_two.feasible)
        self.assertEqual(inconsistent_rank_two.reason, "rank2_inconsistent")

    def test_optional_zero_coordinate_is_accepted(self) -> None:
        # (x-1)(y-1)=-1 and x=0 have the real solution (x,y)=(0,2).
        decision = star.lines_hyperbola_decision(
            ((star.ONE, star.ZERO, star.ZERO),),
            star.ONE,
            star.q(-1),
        )
        self.assertTrue(decision.feasible)
        self.assertEqual(decision.reason, "rank1_fixed_coordinate")


class CorpusCertificateControls(unittest.TestCase):
    def setUp(self) -> None:
        self.graph, self.masks = fixture(
            2592657,
            [0, 3, 4, 11, 14, 17, 18],
            [36, 66, 94, 33, 91, 13, 39, 72, 53, 80, 50, 127],
        )

    def test_exact_fixture_has_exhaustive_four_assignment_certificate(self) -> None:
        certificate = star.find_one_two_star(self.graph, self.masks)
        self.assertIsNotNone(certificate)
        assert certificate is not None
        self.assertEqual(certificate.center, 10)
        self.assertEqual(certificate.one_free_neighbours, (0, 1, 3, 9))
        self.assertEqual(certificate.sign_variables, 2)
        self.assertEqual(certificate.assignments_checked, 4)
        self.assertEqual(
            {failure.reason for failure in certificate.failures},
            {"rank2_diagonal_mismatch"},
        )
        star.verify_certificate(self.graph, self.masks, certificate)

    def test_json_roundtrip_and_tamper_rejection(self) -> None:
        certificate = star.find_one_two_star(self.graph, self.masks)
        assert certificate is not None
        payload = star.certificate_json(certificate)
        rebuilt = star.certificate_from_json(payload)
        self.assertEqual(rebuilt, certificate)
        star.verify_certificate(self.graph, self.masks, rebuilt)

        tampered = copy.deepcopy(payload)
        tampered["failures"][0]["reason"] = "rank1_negative_discriminant"
        with self.assertRaisesRegex(ValueError, "sign replay differs"):
            star.verify_certificate(
                self.graph,
                self.masks,
                star.certificate_from_json(tampered),
            )

    def test_known_passing_family_is_not_rejected(self) -> None:
        graph, masks = fixture(
            316173,
            [2, 5, 7, 9, 13, 15, 18],
            [
                36,
                105,
                127,
                18,
                42,
                68,
                35,
                81,
                88,
                22,
                11,
                116,
            ],
        )
        self.assertIsNone(star.find_one_two_star(graph, masks))


if __name__ == "__main__":
    unittest.main()
