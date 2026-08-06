#!/usr/bin/env python3
"""Tests for the exact K7 degree-two/preordering pilot."""

from __future__ import annotations

import json
import unittest
from fractions import Fraction as Q
from pathlib import Path

import d6_k7_dual_degree2_pilot as pilot
import d6_k7_dual_degree2_verify as independent
import d6_k7_positive_polynomial_dual as degree_one
import d6_k7_rank_reference as prior


ROOT = Path(__file__).resolve().parent


class DegreeTwoPreorderingTest(unittest.TestCase):
    def test_outside_diagonal_formula(self) -> None:
        # For A={0}, T*x0-x0^2-T = x0*x1-x1-1.
        self.assertEqual(
            pilot.outside_diagonal_polynomial((1, 0)),
            {(1, 1): Q(1), (0, 1): Q(-1), (0, 0): Q(-1)},
        )

    def test_synthetic_strict_inequality_certificate(self) -> None:
        # f=x-1=0 and g=x-2>0 are incompatible because f=g+1.
        equations = [{(1,): Q(1), (0,): Q(-1)}]
        generators = [
            pilot.Generator((7,), {(1,): Q(1), (0,): Q(-2)})
        ]
        certificate = pilot.find_raw_certificate(
            equations, generators, variables=1, output_degree=2
        )
        self.assertIsNotNone(certificate)
        self.assertTrue(certificate["generator_multipliers"])
        pilot.verify_raw_certificate(
            equations, generators, variables=1, certificate=certificate
        )

    def test_feasible_strict_system_has_no_certificate(self) -> None:
        # f=x-2=0 and g=x-1>0 have the feasible point x=2.
        equations = [{(1,): Q(1), (0,): Q(-2)}]
        generators = [
            pilot.Generator((7,), {(1,): Q(1), (0,): Q(-1)})
        ]
        self.assertIsNone(pilot.find_raw_certificate(
            equations, generators, variables=1, output_degree=2
        ))

    def test_real_clique_certificate_is_checked_independently(self) -> None:
        sample = json.loads(
            (ROOT / "d6_k7_rank_sample.json").read_text(encoding="utf-8")
        )
        graph = next(
            item for item in sample["graphs"] if item["index"] == 2592168
        )
        target = json.loads(
            (ROOT / "d6_k7_h_mps_triage_report.json").read_text(
                encoding="utf-8"
            )
        )["triage"]["targets"][0]
        graph_n = tuple(prior.induced_graph(
            tuple(graph["adjacency"]), target["nvertices"]
        ))
        clique_mask = sum(1 << vertex for vertex in target["clique"])
        # This target already has a degree-one certificate, so the richer
        # cone must also be able to locate a valid exact certificate.
        baseline = degree_one.find_clique_certificate(
            graph_n, clique_mask, multiplier_degree=1
        )
        self.assertIsNotNone(baseline)
        certificate = pilot.find_clique_certificate(
            graph_n, clique_mask, output_degree=4,
            maximum_generator_order=2,
        )
        self.assertIsNotNone(certificate)
        pilot.verify_clique_certificate(graph_n, clique_mask, certificate)
        independent.verify_clique_certificate(
            graph_n, clique_mask, certificate
        )
        tampered = json.loads(json.dumps(certificate))
        target_items = (
            tampered["positive_polynomial"]
            or tampered["generator_multipliers"]
        )
        target_items[0]["coefficient"][0] += 1
        with self.assertRaises(ValueError):
            independent.verify_clique_certificate(
                graph_n, clique_mask, tampered
            )


if __name__ == "__main__":
    unittest.main()
