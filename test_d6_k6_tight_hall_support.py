#!/usr/bin/env python3
"""Focused controls for exact K6 tight-Hall support propagation."""

from __future__ import annotations

import ast
import unittest
from pathlib import Path

import d6_k6_psd_z_hereditary as hall
import d6_k6_tight_hall_support as production
import verify_d6_k6_psd_z_hereditary as independent_prior
import verify_d6_k6_tight_hall_support as verifier


class TightHallSupportTests(unittest.TestCase):
    def test_tight_union_is_detected_with_allowed_supersets(self) -> None:
        groups = (
            hall.Group(
                "first",
                (
                    hall.Choice(0, 0, (), "empty"),
                    hall.Choice(1, 0b0011, (10,), "selected"),
                ),
                True,
            ),
            hall.Group(
                "second",
                (
                    hall.Choice(0, 0, (), "empty"),
                    hall.Choice(1, 0b0011, (11,), "selected"),
                ),
                True,
            ),
        )
        # Two orthogonal rank-one spans contained in the same two-coordinate
        # plane fill it, even though either vector may use a strict subset of
        # its allowed support.
        tight = production.tight_states(groups)
        self.assertIn(0b0011, tight)
        self.assertEqual(sum(choice.rank for _, choice in tight[0b0011]), 2)

    def test_required_gram_entry_after_support_deletion(self) -> None:
        graph = (0b10, 0b01)
        absolute = (7, 8)
        defects = {7: 0b0011, 8: 0b0110}
        certificate = production.support_obstruction(
            graph, absolute, defects, removed=0b0010
        )
        self.assertEqual(
            certificate["kind"],
            "required_nonzero_gram_entry_lost_common_support",
        )

    def test_same_line_singleton_collision_and_exact_roots(self) -> None:
        graph = (0b10, 0b01)
        absolute = (4, 9)
        defects = {4: 0b0011, 9: 0b0101}
        certificate = production.support_obstruction(
            graph, absolute, defects, removed=0b0110
        )
        self.assertEqual(certificate["kind"], "same_line_singleton_collision")
        # With u=a e_i and common Lorentz ratio r=z/c, the diagonal equation
        # factors exactly as
        #   (a+1)((r^2+5)a+(r^2-7))=0.
        # The first root is the seed point and has zero Lorentz factor; hence
        # at most one non-colliding point remains on this coordinate line.
        for r_squared in (0, 1, 2, 7, 19):
            for a in (-1, (7 - r_squared) / (r_squared + 5)):
                left = (r_squared - 1) * (a + 1) ** 2 + 6 * a**2 - 6
                self.assertAlmostEqual(left, 0.0)

    def test_fixed_corpus_rejection(self) -> None:
        records, _, _ = production.load_input()
        record = next(row for row in records if int(row["index"]) == 32033)
        result = production.evaluate_record(record)
        self.assertTrue(result["rejected"])
        self.assertEqual(result["first_impossible_seed"], [0, 1, 5, 8, 13, 18])

    def test_positive_18_control(self) -> None:
        result = production.positive_control()
        self.assertTrue(result["passed"])
        self.assertEqual(result["K6_seeds"], 32)

    def test_independent_container_enumeration_agrees_synthetically(self) -> None:
        spans = (
            independent_prior.Span(
                "first",
                (
                    independent_prior.Option(0, 0, (), "empty"),
                    independent_prior.Option(1, 0b0011, (10,), "selected"),
                ),
                True,
            ),
            independent_prior.Span(
                "second",
                (
                    independent_prior.Option(0, 0, (), "empty"),
                    independent_prior.Option(1, 0b0011, (11,), "selected"),
                ),
                True,
            ),
        )
        tight = dict(verifier.tight_containers(spans))
        self.assertIn(0b0011, tight)

    def test_verifier_import_independence(self) -> None:
        verifier.assert_import_independence()
        tree = ast.parse(Path(verifier.__file__).read_text(encoding="utf-8"))
        imports = {
            node.module
            for node in ast.walk(tree)
            if isinstance(node, ast.ImportFrom) and node.module
        }
        self.assertNotIn("d6_k6_tight_hall_support", imports)


if __name__ == "__main__":
    unittest.main()
