#!/usr/bin/env python3
"""Focused controls for the independent support-capacity pilot checker."""

from __future__ import annotations

import json
import unittest
from pathlib import Path

import verify_d6_k7_support_capacity_pilot as verifier


ROOT = Path(__file__).resolve().parent


def graph_from_edges(order: int, edges):
    graph = [0] * order
    for first, second in edges:
        graph[first] |= 1 << second
        graph[second] |= 1 << first
    return tuple(graph)


class IndependentCapacityControls(unittest.TestCase):
    def test_direct_witness_retains_optional_nonedge(self):
        graph = graph_from_edges(3, ((0, 1),))
        masks = (0b11, 0b111, 0b110)
        # Vertices 1 and 2 are a candidate nonedge; their actual supports may
        # intersect or not.  Only the required edge 0--1 is checked.
        verifier.direct_capacity_witness_check(
            graph, masks, (), (0b11, 0b11, 0b110)
        )

    def test_direct_witness_rejects_capacity_overflow(self):
        with self.assertRaisesRegex(ValueError, "multiplicity"):
            verifier.direct_capacity_witness_check(
                graph_from_edges(2, ()), (0b1, 0b1), (), (0b1, 0b1)
            )

    def test_direct_witness_rejects_disjoint_required_edge(self):
        with self.assertRaisesRegex(ValueError, "disjoint"):
            verifier.direct_capacity_witness_check(
                graph_from_edges(2, ((0, 1),)),
                (0b1, 0b10), (), (0b1, 0b10),
            )

    def test_stratified_selection_reconstructs_committed_sample(self):
        union = json.loads(
            (ROOT / "d6_k7_rankone_pattern_union.json").read_text()
        )
        report = json.loads(
            (ROOT / "d6_k7_support_capacity_pilot_report.json").read_text()
        )
        profiles = union["cover_structure"]["residue_profiles"]
        expected = report["selection"]["selected_indices"]
        self.assertEqual(
            verifier.independent_sample(profiles, len(expected)), expected
        )


if __name__ == "__main__":
    unittest.main()
