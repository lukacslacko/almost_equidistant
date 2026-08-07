#!/usr/bin/env python3
"""Tests for the exact K7 positive-polynomial dual certificate."""

from __future__ import annotations

import json
import unittest
from fractions import Fraction as Q
from pathlib import Path

import d6_k7_rank_reference as prior
import d6_k7_positive_polynomial_dual as dual
import verify_d6_k7_positive_polynomial_dual as independent


ROOT = Path(__file__).resolve().parent


class PositivePolynomialDualTest(unittest.TestCase):
    def test_synthetic_infeasible_identity_is_certified(self) -> None:
        # The equality -1=0 has the certificate (-1)*(-1)=1.
        equations = [{(0,): Q(-1)}]
        certificate = dual.find_raw_certificate(
            equations, variables=1, multiplier_degree=0
        )
        self.assertIsNotNone(certificate)
        positive = dual.verify_raw_certificate(equations, 1, certificate)
        self.assertEqual(positive, {(0,): Q(1)})

    def test_positive_feasible_equation_has_no_certificate(self) -> None:
        # x-1=0 has the strictly positive solution x=1, so no sound
        # positive-polynomial certificate can exist at any degree.
        equations = [{(1,): Q(1), (0,): Q(-1)}]
        self.assertIsNone(dual.find_raw_certificate(
            equations, variables=1, multiplier_degree=1
        ))

    def test_exact_certificate_on_frozen_sample_target(self) -> None:
        sample = json.loads(
            (ROOT / "d6_k7_rank_sample.json").read_text(encoding="utf-8")
        )
        graph = next(
            graph for graph in sample["graphs"] if graph["index"] == 2592168
        )
        target = json.loads(
            (ROOT / "d6_k7_h_mps_triage_report.json").read_text(
                encoding="utf-8"
            )
        )["triage"]["targets"][0]
        graph_n = prior.induced_graph(
            tuple(graph["adjacency"]), target["nvertices"]
        )
        clique_mask = sum(1 << vertex for vertex in target["clique"])
        certificate = dual.find_clique_certificate(
            graph_n, clique_mask, multiplier_degree=1
        )
        self.assertIsNotNone(certificate)
        positive = dual.verify_clique_certificate(
            graph_n, clique_mask, certificate
        )
        self.assertTrue(positive)
        self.assertTrue(all(value > 0 for value in positive.values()))
        independent.verify_clique_certificate(
            tuple(graph_n), clique_mask, certificate
        )

        # The checker must reject a changed rational identity.
        tampered = json.loads(json.dumps(certificate))
        tampered["multipliers"][0]["coefficient"][0] += 1
        with self.assertRaises(ValueError):
            dual.verify_clique_certificate(graph_n, clique_mask, tampered)
        with self.assertRaises(ValueError):
            independent.verify_clique_certificate(
                tuple(graph_n), clique_mask, tampered
            )


if __name__ == "__main__":
    unittest.main()
