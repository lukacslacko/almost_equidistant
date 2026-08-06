#!/usr/bin/env python3
"""Tests for arbitrary-principal-basis PSD-atom certificates."""

from __future__ import annotations

import json
import unittest
from pathlib import Path

import d6_k7_arbitrary_basis_psd_dual as dual
import d6_k7_rank_reference as prior


ROOT = Path(__file__).resolve().parent


class ArbitraryBasisPsdDualTest(unittest.TestCase):
    def test_zero_outside_column_is_impossible(self) -> None:
        # For core {0}, the proposed outside column is zero.  The necessary
        # strict Schur diagonal 0^T H 0 > 1 is immediately impossible.
        system = dual.basis_affine_system((0, 0), (0,))
        certificate = dual.find_certificate(system)
        self.assertIsNotNone(certificate)
        dual.verify_certificate(system, certificate)
        self.assertEqual(certificate["psd_atoms"], [])
        self.assertTrue(certificate["diagonal_multipliers"])

    def test_feasible_one_dimensional_basis_is_not_rejected(self) -> None:
        # With K_01=1, any scalar H>1 satisfies the only outside condition.
        system = dual.basis_affine_system((2, 1), (0,))
        self.assertIsNone(dual.find_certificate(system))

    def test_exact_frozen_sample_core_certificate(self) -> None:
        sample = json.loads(
            (ROOT / "d6_k7_rank_sample.json").read_text(encoding="utf-8")
        )
        graph = next(
            graph for graph in sample["graphs"] if graph["index"] == 58186
        )
        # Frozen K7 seed [0,2,5,8,9,12,17], cover zmask=128.
        nvertices = [1, 3, 4, 6, 7, 10, 11, 14, 15, 16, 18]
        graph_n = prior.induced_graph(tuple(graph["adjacency"]), nvertices)
        system = dual.basis_affine_system(graph_n, (0, 1, 2, 3, 4))
        certificate = dual.find_certificate(system, maximum_atom_support=3)
        self.assertIsNotNone(certificate)
        dual.verify_certificate(system, certificate)

        tampered = json.loads(json.dumps(certificate))
        if tampered["psd_atoms"]:
            tampered["psd_atoms"][0]["weight"][0] += 1
        else:
            tampered["diagonal_multipliers"][0]["value"][0] += 1
        with self.assertRaises(ValueError):
            dual.verify_certificate(system, tampered)


if __name__ == "__main__":
    unittest.main()
