#!/usr/bin/env python3
"""Focused exact controls for the opposite-light-ray rank package."""

from __future__ import annotations

import hashlib
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import d6_k6_opposite_ray_rank as producer
import verify_d6_k6_opposite_ray_rank as verifier
from d6_k6_lorentz_reference import build_instance
from test_d6_k6_lorentz import lower_bound_18_graph, opposite_color_triangles


def edge(rows: list[int], left: int, right: int) -> None:
    rows[left] |= 1 << right
    rows[right] |= 1 << left


class OppositeRayRankControls(unittest.TestCase):
    def test_psd_z_component_rank_lower(self) -> None:
        self.assertEqual(producer.psd_z_component_rank_lower((0,)), (0, []))
        self.assertEqual(producer.psd_z_component_rank_lower((2, 1)), (1, [[0, 1]]))
        path = (2, 5, 2)
        self.assertEqual(
            producer.psd_z_component_rank_lower(path), (2, [[0, 1, 2]])
        )
        two_edges = (2, 1, 8, 4, 0)
        self.assertEqual(
            producer.psd_z_component_rank_lower(two_edges),
            (2, [[0, 1], [2, 3]]),
        )

    def test_zero_cross_degree_is_rank_one_not_zero(self) -> None:
        instance = SimpleNamespace(outside=tuple(range(7)))
        adjacency = [0] * 7
        # Two source vertices, five target vertices, and no cross edges.
        decision = producer.directed_ray_check(
            adjacency, instance, 0, 0b11, sum(1 << i for i in range(2, 7))
        )
        self.assertEqual(decision["component_rank_lower"], 0)
        self.assertEqual(decision["zero_cross_degree_source"], [0, 1])
        self.assertEqual(decision["rank_lower"], 2)
        self.assertEqual(decision["rank_capacity"], 1)
        self.assertFalse(decision["passed"])

    def test_b5_common_neighbour_rank_boundary(self) -> None:
        instance = SimpleNamespace(outside=tuple(range(8)))
        adjacency = [0] * 8
        # Three source vertices all meet target vertex 3.  Their common-
        # neighbour graph is K3, hence rank lower 2 > 6-5.
        for source in range(3):
            edge(adjacency, source, 3)
        target = sum(1 << i for i in range(3, 8))
        decision = producer.directed_ray_check(
            adjacency, instance, 0, 0b111, target
        )
        self.assertEqual(decision["common_neighbour_components"], [[0, 1, 2]])
        self.assertEqual(decision["rank_lower"], 2)
        self.assertEqual(decision["rank_capacity"], 1)
        self.assertFalse(decision["passed"])

    def test_all_component_colorings_are_needed(self) -> None:
        adjacency = opposite_color_triangles()
        instance = build_instance(adjacency, tuple(range(6)))
        decision, counts = producer.ray_assignment_decision(adjacency, instance, 0)
        self.assertTrue(decision["passed"])
        self.assertEqual(decision["expected_colorings_mod_ray_swap"], 2)
        self.assertEqual(decision["colorings_considered"], 2)
        self.assertEqual(counts["ray_coloring_bin_hall_fail"], 1)
        self.assertEqual(decision["witness"]["colors"], [0, 1])

    def test_fixed_corpus_rejection_and_independent_agreement(self) -> None:
        records, _, _ = producer.load_input()
        record = next(row for row in records if row["index"] == 207582)
        first = producer.evaluate_record(record)
        second = verifier.evaluate_record(record)
        self.assertEqual(first, second)
        self.assertTrue(first["rejected"])
        self.assertEqual(first["first_impossible_seed"], [0, 2, 9, 14, 16, 17])
        self.assertEqual(len(first["certificate"]["Z0_failures"]), 64)

    def test_standard_18_passes_all_32_k6_seeds(self) -> None:
        adjacency = lower_bound_18_graph()
        result = producer.evaluate_record({"index": -18, "adjacency": adjacency})
        independent = verifier.evaluate_record(
            {"index": -18, "adjacency": verifier.standard_18_graph()}
        )
        self.assertFalse(result["rejected"])
        self.assertEqual(result["seeds_checked"], 32)
        self.assertEqual(result, independent)

    def test_verifier_import_independence(self) -> None:
        verifier.assert_import_independence()

    def test_porcelain_binding_allows_only_untracked_dirt(self) -> None:
        lines = ["?? discovery.json", "?? scratch.py"]
        status = "\n".join(lines)
        cleanly_bound = {
            "porcelain_lines": lines,
            "porcelain_sha256": hashlib.sha256(status.encode("utf-8")).hexdigest(),
            "dirty": True,
        }
        observed = verifier.validate_porcelain_binding(cleanly_bound)
        self.assertEqual(observed["untracked_entries"], 2)
        tracked = dict(cleanly_bound)
        tracked["porcelain_lines"] = [" M d6_k6_opposite_ray_rank.py"]
        tracked_status = "\n".join(tracked["porcelain_lines"])
        tracked["porcelain_sha256"] = hashlib.sha256(
            tracked_status.encode("utf-8")
        ).hexdigest()
        with self.assertRaisesRegex(ValueError, "tracked source was dirty"):
            verifier.validate_porcelain_binding(tracked)
        tampered = dict(cleanly_bound, porcelain_sha256="0" * 64)
        with self.assertRaisesRegex(ValueError, "porcelain binding differs"):
            verifier.validate_porcelain_binding(tampered)

    def test_deterministic_gzip_certificates(self) -> None:
        payload = {"schema": "control", "rows": [1, 2, 3]}
        with tempfile.TemporaryDirectory() as directory:
            first = Path(directory) / "first.json.gz"
            second = Path(directory) / "second.json.gz"
            raw_first = producer.atomic_gzip_json(first, payload)
            raw_second = producer.atomic_gzip_json(second, payload)
            self.assertEqual(raw_first, raw_second)
            self.assertEqual(first.read_bytes(), second.read_bytes())


if __name__ == "__main__":
    unittest.main()
