#!/usr/bin/env python3
"""Controls for the exact K6 empty-essential / virtual-K7 package."""

from __future__ import annotations

import ast
import json
import tempfile
import unittest
from pathlib import Path

import d6_k6_empty_essential_profile as profiler
import d6_k6_empty_essential_virtual_k7 as conjunction
import verify_d6_k6_empty_essential_profile as profile_verifier
import verify_d6_k6_empty_essential_virtual_k7 as conjunction_verifier


ROOT = Path(__file__).resolve().parent


class EmptyStateControls(unittest.TestCase):
    def test_opposite_pair_excludes_required_edge_only_in_pair_state(self) -> None:
        adjacency = [0, 1 << 2, 1 << 1]
        states = {()}
        first = profiler.extend_states(states, (1,), adjacency, True)
        self.assertEqual(first, {(1,)})
        second = profiler.extend_states(first, (2,), adjacency, True)
        self.assertEqual(second, set())
        relaxed = profiler.extend_states(first, (2,), adjacency, False)
        self.assertEqual(relaxed, {(1, 2)})

    def test_all_nonempty_choice_is_retained(self) -> None:
        adjacency = (0, 0, 0)
        states = profiler.extend_states({()}, (None, 1), adjacency, True)
        self.assertEqual(states, {(), (1,)})


class VirtualK7Controls(unittest.TestCase):
    def test_augmentation_adds_only_apex_star(self) -> None:
        adjacency = list(conjunction.complete_graph(6)) + [0]
        augmented, added, seed_mask = conjunction.augment_with_apex(
            adjacency, list(range(6)), 6
        )
        self.assertEqual(seed_mask, (1 << 7) - 1)
        self.assertEqual(len(added), 6)
        self.assertEqual(augmented, conjunction.complete_graph(7))

    def test_positive_and_negative_controls(self) -> None:
        observed = conjunction.controls()
        self.assertTrue(observed["positive_regular_K7"]["passed"])
        self.assertEqual(
            observed["positive_regular_K7"]["passing_small_covers"], [0]
        )
        self.assertTrue(observed["negative_K8"]["passed"])

    def test_independent_cover_enumerator_includes_supercovers(self) -> None:
        # One L-edge 0--1 with eligible isolate 2.
        ladj = (0b010, 0b001, 0)
        observed = conjunction_verifier.independent_covers(ladj, 0b111)
        self.assertIn(0b101, observed)
        self.assertIn(0b110, observed)
        self.assertNotIn(0, observed)


class ArtifactControls(unittest.TestCase):
    def test_profile_and_conjunction_totals(self) -> None:
        profile = json.loads(
            (ROOT / "d6_k6_empty_essential_profile.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(
            profile["aggregate"]["seed_classification_counts"],
            {"ALL_NONEMPTY_WITNESS": 25243, "EMPTY_ESSENTIAL": 111},
        )
        self.assertEqual(
            profile["aggregate"]["totals"]["graphs_with_empty_essential_seed"],
            49,
        )
        result = json.loads(
            (ROOT / "d6_k6_empty_essential_virtual_k7.json").read_text(
                encoding="utf-8"
            )
        )
        totals = result["aggregate"]["totals"]
        self.assertEqual(
            totals,
            {
                "direct_cap_covers": 2237,
                "eligible_covers": 2386,
                "impossible_seeds": 111,
                "promotions": 291,
                "rejected_graphs": 49,
                "small_covers": 149,
                "target_graphs": 49,
                "target_seeds": 111,
            },
        )
        self.assertEqual(
            result["aggregate"]["promotion_status_counts"],
            {"ELIMINATED": 291},
        )
        self.assertEqual(len(result["aggregate"]["ordered_residue_indices"]), 756)

    def test_independent_verifications_pass_and_bind_reports(self) -> None:
        profile_check = json.loads(
            (ROOT / "d6_k6_empty_essential_profile_verification.json").read_text(
                encoding="utf-8"
            )
        )
        conjunction_check = json.loads(
            (
                ROOT
                / "d6_k6_empty_essential_virtual_k7_verification.json"
            ).read_text(encoding="utf-8")
        )
        self.assertEqual(profile_check["status"], "PASS")
        self.assertEqual(conjunction_check["status"], "PASS")
        self.assertEqual(
            profile_check["report_sha256"],
            profiler.sha256(ROOT / "d6_k6_empty_essential_profile.json"),
        )
        self.assertEqual(
            conjunction_check["report_sha256"],
            conjunction.sha256(ROOT / "d6_k6_empty_essential_virtual_k7.json"),
        )
        self.assertEqual(conjunction_check["rejected_graphs"], 49)
        self.assertEqual(conjunction_check["ordered_residue_graphs"], 756)

    def test_verifiers_do_not_import_production_kernels(self) -> None:
        profile_verifier.assert_import_independence()
        conjunction_verifier.assert_import_independence()
        tree = ast.parse(
            (ROOT / "verify_d6_k6_empty_essential_virtual_k7.py").read_text(
                encoding="utf-8"
            )
        )
        imported = {
            node.module
            for node in ast.walk(tree)
            if isinstance(node, ast.ImportFrom) and node.module
        }
        self.assertNotIn("d6_k6_empty_essential_virtual_k7", imported)

    def test_explicit_hash_boundary_rejects_tamper_before_replay(self) -> None:
        source = ROOT / "d6_k6_empty_essential_virtual_k7.json"
        with tempfile.TemporaryDirectory() as directory:
            tampered = Path(directory) / "tampered.json"
            data = bytearray(source.read_bytes())
            data[-2:-1] = b" "
            tampered.write_bytes(data)
            with self.assertRaisesRegex(ValueError, "explicit report SHA-256"):
                conjunction_verifier.verify(
                    tampered,
                    conjunction.sha256(source),
                    Path(directory) / "verification.json",
                    1,
                )


if __name__ == "__main__":
    unittest.main()
