#!/usr/bin/env python3
"""Focused controls for the independent pattern-954 report verifier."""

from __future__ import annotations

import copy
import json
import unittest
from pathlib import Path

import d6_n14_pattern_954_containment_report_verify as verifier


ROOT = Path(__file__).resolve().parent


class IndependentReportControls(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.verification = verifier.verify(ROOT / verifier.REPORT_NAME)
        cls.pattern = tuple(cls.verification["pattern"]["adjacency"])
        cls.selected, cls.targets, _ = verifier.rebuild_targets()
        report = json.loads(
            (ROOT / verifier.REPORT_NAME).read_text(encoding="utf-8")
        )
        cls.first_hit = copy.deepcopy(report["certified_hits"][0])

    def test_full_independent_verification(self) -> None:
        self.assertEqual(self.verification["status"], "PASS")
        self.assertTrue(self.verification["required_edge_only"])
        self.assertEqual(
            self.verification["mapping_verification"],
            {
                "mappings_checked": 3_403,
                "required_edges_checked_per_mapping": 54,
                "required_edge_incidence_checks": 183_762,
                "all_injective": True,
            },
        )

    def test_core_is_independently_reconstructed(self) -> None:
        core = verifier.reconstruct_core(ROOT / "aeq_d6_n14.txt")
        self.assertEqual(core, verifier.EXPECTED_CORE)
        self.assertEqual(sum(row.bit_count() for row in core) // 2, 54)

    def test_residue_is_exact_ordered_hit_complement(self) -> None:
        hits = set(self.verification["certified_hit_indices"])
        residue = [index for index in self.selected if index not in hits]
        self.assertEqual(residue, self.verification["residue_indices"])
        self.assertEqual(len(residue), 9_436)
        self.assertEqual(verifier.stable_hash(residue), verifier.RESIDUE_INDICES_SHA256)

    def test_all_extra_target_edges_are_permitted(self) -> None:
        full = (1 << 19) - 1
        complete = tuple(full ^ (1 << vertex) for vertex in range(19))
        self.assertEqual(
            verifier.verify_mapping(self.pattern, complete, list(range(13))), 54
        )

    def test_mapping_collision_is_rejected(self) -> None:
        mapping = copy.copy(self.first_hit["mapping"])
        mapping[1] = mapping[0]
        target = self.targets[self.first_hit["index"]]
        with self.assertRaisesRegex(AssertionError, "not injective"):
            verifier.verify_mapping(self.pattern, target, mapping)

    def test_lost_required_edge_is_rejected(self) -> None:
        mapping = self.first_hit["mapping"]
        target = list(self.targets[self.first_hit["index"]])
        first, second = next(
            (first, second)
            for first in range(13) for second in range(first)
            if self.pattern[first] & (1 << second)
        )
        left, right = mapping[first], mapping[second]
        target[left] &= ~(1 << right)
        target[right] &= ~(1 << left)
        with self.assertRaisesRegex(AssertionError, "loses required edge"):
            verifier.verify_mapping(self.pattern, target, mapping)

    def test_malformed_timeout_is_not_accepted(self) -> None:
        row = {
            "index": self.selected[0],
            "nodes": 50_000,
            "status": "TIMEOUT",
            "wall_seconds": 0.0,
        }
        with self.assertRaisesRegex(AssertionError, "did not exhaust"):
            verifier.validate_result_record(
                row, self.selected[0], self.pattern,
                self.targets[self.selected[0]],
            )

    def test_immutable_chunk_manifest_is_pinned(self) -> None:
        report = json.loads(
            (ROOT / verifier.REPORT_NAME).read_text(encoding="utf-8")
        )
        manifest = verifier.chunk_manifest(Path(report["checkpoint_directory"]))
        self.assertEqual(len(manifest), verifier.CHUNK_COUNT)
        self.assertEqual(sum(row["bytes"] for row in manifest), verifier.CHUNK_BYTES)
        self.assertEqual(
            verifier.stable_hash(manifest), verifier.CHUNK_MANIFEST_SHA256
        )

    def test_config_mutation_changes_pinned_hash(self) -> None:
        mutated = copy.deepcopy(verifier.EXPECTED_CONFIG)
        mutated["node_limit"] += 1
        self.assertNotEqual(verifier.stable_hash(mutated), verifier.CONFIG_SHA256)


if __name__ == "__main__":
    unittest.main()
