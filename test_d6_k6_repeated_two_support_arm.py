#!/usr/bin/env python3
"""Focused exact controls for the repeated-two-support-arm K6 layer."""

from __future__ import annotations

import unittest
from fractions import Fraction

import d6_k6_repeated_two_support_arm as production
import verify_d6_k6_repeated_two_support_arm as verifier


class RepeatedTwoSupportArmTests(unittest.TestCase):
    def test_synthetic_forbidden_repeated_arm(self) -> None:
        # Only the two anchor--arm pairs are required.  The arm pair is a
        # candidate nonedge here and is deliberately unused by the theorem.
        graph = (0b110, 0b001, 0b001)
        absolute = (10, 11, 12)
        defects = {10: 0b0001, 11: 0b0011, 12: 0b0011}
        certificate = production.repeated_two_support_arm_obstruction(
            graph, absolute, defects, removed=0
        )
        self.assertEqual(
            certificate["kind"],
            "singleton_two_identical_two_support_arms",
        )
        self.assertEqual(certificate["supports"], [[0], [0, 1], [0, 1]])
        self.assertEqual(
            {tuple(pair) for pair in certificate["required_anchor_pairs"]},
            {(10, 11), (10, 12)},
        )
        self.assertEqual(certificate["arm_pair_distance"], "unused")

    def test_exact_normalized_reconstruction(self) -> None:
        # The normalized arm is w=(alpha,a), with c=1/(1-alpha-a).
        # Every equality below is rational; r itself need not be adjoined.
        values = (
            Fraction(0), Fraction(1, 2), Fraction(2), Fraction(5),
            Fraction(8), Fraction(21, 2), Fraction(29),
        )
        for R in values:
            self.assertNotIn(R, (1, 7))  # generic, non-lightlike branch
            kappa = (1 - R) / 6
            anchor_normalized = (7 - R) / 12
            alpha = kappa / anchor_normalized
            self.assertEqual(alpha, 2 * (1 - R) / (7 - R))
            beta = 1 - alpha
            delta = kappa - alpha * alpha
            arm_second = (delta + beta * beta) / (2 * beta)

            # Optional-zero branch resolved by the equation itself.
            self.assertEqual(arm_second, (R + 5) / 12)
            self.assertNotEqual(arm_second, 0)
            self.assertEqual(
                arm_second * arm_second,
                delta + (beta - arm_second) ** 2,
            )

            denominator = beta - arm_second
            self.assertEqual(denominator, (R + 5) ** 2 / (12 * (7 - R)))
            self.assertNotEqual(denominator, 0)
            scale = 1 / denominator
            self.assertEqual(scale, 12 * (7 - R) / (R + 5) ** 2)

            u_i = scale * alpha
            u_j = scale * arm_second
            self.assertEqual(scale, 1 + u_i + u_j)
            # This is K6-diag after replacing z^2 by R*c^2.
            self.assertEqual(
                R * scale * scale,
                scale * scale + 6 - 6 * (u_i * u_i + u_j * u_j),
            )

    def test_R7_and_lightlike_are_outside_the_rule_derivation(self) -> None:
        # R=7 makes the nonempty singleton anchor's defect zero.  R=1 is the
        # lightlike direction, which the inherited solver checks separately.
        R = Fraction(7)
        self.assertEqual((7 - R) / 12, 0)
        self.assertEqual((1 - R) / 6, -1)
        R = Fraction(1)
        self.assertEqual((1 - R) / 6, 0)

    def test_optional_zero_and_required_edge_guards(self) -> None:
        anchor_v = (0b110, 0b001, 0b001)
        absolute = (10, 11, 12)

        # A third allowed coordinate leaves a genuinely two-variable arm, so
        # the cancellation does not determine the point.
        defects = {10: 0b0001, 11: 0b0011, 12: 0b0111}
        self.assertIsNone(
            production.repeated_two_support_arm_obstruction(
                anchor_v, absolute, defects, removed=0
            )
        )

        # Different second coordinates are not the repeated-arm motif (they
        # belong to the separately frozen singleton-fan theorem).
        defects = {10: 0b0001, 11: 0b0011, 12: 0b0101}
        self.assertIsNone(
            production.repeated_two_support_arm_obstruction(
                anchor_v, absolute, defects, removed=0
            )
        )

        # One missing anchor edge removes the reconstruction equation.  The
        # mutual arm pair remains irrelevant either way.
        missing_anchor_edge = (0b010, 0b001, 0b000)
        defects = {10: 0b0001, 11: 0b0011, 12: 0b0011}
        self.assertIsNone(
            production.repeated_two_support_arm_obstruction(
                missing_anchor_edge, absolute, defects, removed=0
            )
        )

    def test_fixed_marginal_rejection(self) -> None:
        records, _, _ = production.load_input()
        record = next(row for row in records if int(row["index"]) == 552851)
        result = production.evaluate_record(record)
        self.assertTrue(result["rejected"])
        self.assertEqual(result["first_impossible_seed_mask"], 139339)
        self.assertEqual(result["first_impossible_seed"], [0, 1, 3, 6, 13, 17])
        collisions = []
        for row in result["certificate"]["rows"]:
            failure = row.get("failure") or {}
            for branch in failure.get("branches", []):
                certificate = branch.get("certificate") or {}
                obstruction = certificate.get("obstruction") or {}
                if obstruction.get("kind") == (
                    "singleton_two_identical_two_support_arms"
                ):
                    collisions.append(obstruction)
        self.assertTrue(collisions)
        self.assertEqual(collisions[0]["vertices"], [2, 12, 18])
        self.assertEqual(collisions[0]["supports"], [[1], [1, 5], [1, 5]])

    def test_parent_residue_and_positive_control(self) -> None:
        records, indices, _ = production.load_input()
        self.assertEqual(len(records), 634)
        self.assertEqual(len(indices), 634)
        self.assertNotIn(58458, indices)  # frozen singleton-fan rejection
        self.assertIn(552851, indices)  # first new marginal rejection
        control = production.positive_control()
        self.assertTrue(control["passed"])
        self.assertEqual(control["K6_seeds"], 32)

    def test_independent_transcription_and_import_boundary(self) -> None:
        graph = (0b110, 0b001, 0b001)
        absolute = (10, 11, 12)
        defects = {10: 0b0001, 11: 0b0011, 12: 0b0011}
        self.assertTrue(
            verifier.repeated_two_support_arm_obstructed(
                graph, absolute, defects, removed=0
            )
        )
        verifier.assert_import_independence()


if __name__ == "__main__":
    unittest.main()
