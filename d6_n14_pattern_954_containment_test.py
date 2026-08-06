#!/usr/bin/env python3
"""Controls for the exact pattern-954 non-induced containment kernel."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import d6_n14_pattern_954_containment as containment


ROOT = Path(__file__).resolve().parent


class ContainmentControls(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.temporary = tempfile.TemporaryDirectory(dir=ROOT / ".runs")
        cls.library = Path(cls.temporary.name) / "pattern954.dylib"
        cls.payload, cls.pattern, cls.rows = containment.load_inputs(
            ROOT / "d6_n14_pattern_954_input.json"
        )
        cls.compilation = containment.compile_kernel(
            ROOT / "d6_n14_pattern_954_containment.c", cls.library, "cc"
        )
        cls.matcher = containment.Matcher(cls.library, cls.pattern)

    @classmethod
    def tearDownClass(cls) -> None:
        cls.temporary.cleanup()

    def test_strict_compile_is_clean_and_hash_pinned(self) -> None:
        self.assertEqual(self.compilation["stderr"], "")
        self.assertEqual(
            containment.sha256(ROOT / "d6_n14_pattern_954_containment.c"),
            containment.EXPECTED_C_SOURCE_SHA256,
        )

    def test_identity_embedding_returns_checked_mapping(self) -> None:
        target = containment.embedded_pattern_target(self.pattern)
        result = self.matcher.call(target, 50_000)
        self.assertEqual(result["status"], "HIT")
        containment.verify_hit_mapping(
            self.pattern, target, result["mapping"]
        )

    def test_edge_supergraph_returns_checked_mapping(self) -> None:
        target = containment.complete_target()
        result = self.matcher.call(target, 50_000)
        self.assertEqual(result["status"], "HIT")
        containment.verify_hit_mapping(
            self.pattern, target, result["mapping"]
        )

    def test_empty_graph_is_no_hit(self) -> None:
        result = self.matcher.call((0,) * containment.TARGET_ORDER, 50_000)
        self.assertEqual(result, {"status": "NO_HIT", "nodes": 0})

    def test_cap_exhaustion_is_timeout_not_no_hit(self) -> None:
        target = containment.embedded_pattern_target(self.pattern)
        result = self.matcher.call(target, 1)
        self.assertEqual(result["status"], "TIMEOUT")
        self.assertGreater(result["nodes"], 1)

    def test_invalid_graph_is_not_silently_discarded(self) -> None:
        target = [0] * containment.TARGET_ORDER
        target[0] = 1
        result = self.matcher.call(target, 50_000)
        self.assertEqual(result["status"], "INVALID")

    def test_known_current_residue_hit(self) -> None:
        by_index = {row["index"]: row for row in self.rows}
        target = by_index[containment.KNOWN_HIT_INDEX]["adjacency"]
        result = self.matcher.call(target, 50_000)
        self.assertEqual(result["status"], "HIT")
        containment.verify_hit_mapping(
            self.pattern, target, result["mapping"]
        )

    def test_hit_mapping_checker_rejects_collisions_and_lost_edges(self) -> None:
        target = containment.embedded_pattern_target(self.pattern)
        collision = list(range(containment.PATTERN_ORDER))
        collision[1] = collision[0]
        with self.assertRaisesRegex(ValueError, "not injective"):
            containment.verify_hit_mapping(self.pattern, target, collision)
        lost_edge_target = list(target)
        first = 0
        second = (self.pattern[first] & -self.pattern[first]).bit_length() - 1
        lost_edge_target[first] &= ~(1 << second)
        lost_edge_target[second] &= ~(1 << first)
        with self.assertRaisesRegex(ValueError, "loses pattern edge"):
            containment.verify_hit_mapping(
                self.pattern, lost_edge_target, list(range(13))
            )


if __name__ == "__main__":
    unittest.main()
