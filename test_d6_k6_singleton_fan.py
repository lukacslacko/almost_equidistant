#!/usr/bin/env python3
"""Focused exact controls for the incremental K6 singleton-fan layer."""

from __future__ import annotations

import unittest
from fractions import Fraction

import d6_k6_singleton_fan as production
import verify_d6_k6_singleton_fan as verifier


class SingletonFanTests(unittest.TestCase):
    def test_synthetic_forbidden_fan(self) -> None:
        graph = (0b110, 0b101, 0b011)  # required K3
        absolute = (10, 11, 12)
        defects = {10: 0b0001, 11: 0b0011, 12: 0b0101}
        certificate = production.singleton_fan_obstruction(
            graph, absolute, defects, removed=0
        )
        self.assertEqual(
            certificate["kind"],
            "singleton_two_distinct_two_support_arms",
        )
        self.assertEqual(certificate["supports"], [[0], [0, 1], [0, 2]])

    def test_optional_zero_semantics_need_forced_two_supports(self) -> None:
        graph = (0b110, 0b101, 0b011)
        absolute = (10, 11, 12)
        # The last arm still has a third allowed coordinate.  It need not be
        # a genuine two-support vector, so the fan lemma is not asserted.
        defects = {10: 0b0001, 11: 0b0011, 12: 0b1101}
        self.assertIsNone(
            production.singleton_fan_obstruction(
                graph, absolute, defects, removed=0
            )
        )
        # Two arms with the same second coordinate also do not give the
        # distinct-arm algebraic equation.
        defects = {10: 0b0001, 11: 0b0011, 12: 0b0011}
        self.assertIsNone(
            production.singleton_fan_obstruction(
                graph, absolute, defects, removed=0
            )
        )

    def test_exact_R_equation(self) -> None:
        # alpha=2(1-R)/(7-R).  The arm--arm edge requires
        # alpha^2=(1-R)/6.  Clearing denominators gives
        # (R-1)(R+5)^2=0.  Generic non-lightlike orientation excludes R=1;
        # R=-5 is impossible because R=r^2>=0.
        for R in range(0, 31):
            if R == 7:
                continue
            alpha = Fraction(2 * (1 - R), 7 - R)
            difference = alpha * alpha - Fraction(1 - R, 6)
            cleared = difference * 6 * (R - 7) ** 2
            self.assertEqual(cleared, (R - 1) * (R + 5) ** 2)
            if R != 1:
                self.assertNotEqual(difference, 0)

    def test_fixed_marginal_rejection(self) -> None:
        records, _, _ = production.load_input()
        record = next(row for row in records if int(row["index"]) == 58458)
        result = production.evaluate_record(record)
        self.assertTrue(result["rejected"])
        self.assertEqual(result["first_impossible_seed"], [3, 8, 9, 11, 14, 16])
        fan_certificates = []
        for row in result["certificate"]["rows"]:
            failure = row.get("failure") or {}
            for branch in failure.get("branches", []):
                certificate = branch.get("certificate") or {}
                obstruction = certificate.get("obstruction") or {}
                if obstruction.get("kind") == "singleton_two_distinct_two_support_arms":
                    fan_certificates.append(obstruction)
        self.assertTrue(fan_certificates)

    def test_parent_residue_and_positive_control(self) -> None:
        records, indices, _ = production.load_input()
        self.assertEqual(len(records), 649)
        self.assertEqual(len(indices), 649)
        self.assertNotIn(32033, indices)  # already frozen out by the 107 layer
        control = production.positive_control()
        self.assertTrue(control["passed"])
        self.assertEqual(control["K6_seeds"], 32)

    def test_independent_transcription_and_import_boundary(self) -> None:
        graph = (0b110, 0b101, 0b011)
        absolute = (10, 11, 12)
        defects = {10: 0b0001, 11: 0b0011, 12: 0b0101}
        self.assertTrue(
            verifier.singleton_fan_obstructed(
                graph, absolute, defects, removed=0
            )
        )
        verifier.assert_import_independence()


if __name__ == "__main__":
    unittest.main()
