#!/usr/bin/env python3
"""Focused controls for the independent graph-3936177 replay verifier."""

from __future__ import annotations

import copy
import tempfile
import unittest
from pathlib import Path

import run_d6_interval_k7_3936177 as wrapper
import verify_d6_interval_k7_3936177 as verifier


TARGET_ADJACENCY_SHA256 = (
    "264bc830ee02673f329bfdf52659d75672b084fe0e3128e597103855a868fb35"
)


def production_selection_configuration() -> dict:
    all_graphs, provenance = wrapper.v6.load_v6_graphs(
        wrapper.ROOT / "d6_current_residue_manifest_v6.json",
        wrapper.ROOT / "d6_current_residue_manifest_v6_verification.json",
    )
    graphs = wrapper.select_target(all_graphs)
    return {
        "selection": {
            "graphs": 1,
            "target_index": wrapper.TARGET_INDEX,
            "population_counts": {"K7": 1},
            "indices_sha256": wrapper.engine.stable_hash(
                [graph["index"] for graph in graphs]
            ),
            "selection_layer": provenance,
        },
        "selector": copy.deepcopy(verifier.EXPECTED_SELECTOR),
    }


class FocusedReplayVerifierControls(unittest.TestCase):
    def test_independent_reconstruction_returns_exact_target(self) -> None:
        graphs = verifier.reconstruct_target(production_selection_configuration())
        self.assertEqual(len(graphs), 1)
        self.assertEqual(graphs[0]["index"], verifier.TARGET_INDEX)
        self.assertEqual(graphs[0]["population"], "K7")
        self.assertEqual(
            verifier.object_sha256(list(graphs[0]["adjacency"])),
            TARGET_ADJACENCY_SHA256,
        )

    def test_selector_tamper_is_rejected(self) -> None:
        configuration = production_selection_configuration()
        configuration["selector"]["index"] += 1
        with self.assertRaisesRegex(AssertionError, "exact selector"):
            verifier.reconstruct_target(configuration)

    def test_selection_trust_tier_tamper_is_rejected(self) -> None:
        configuration = production_selection_configuration()
        configuration["selection"]["selection_layer"]["certificate_accounting"][
            "interval_non_KILLED_used"
        ] = True
        with self.assertRaisesRegex(AssertionError, "selection trust tiers"):
            verifier.reconstruct_target(configuration)

    def test_report_hash_tamper_is_rejected_first(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "report.json"
            path.write_text("{}\n", encoding="utf-8")
            with self.assertRaisesRegex(AssertionError, "report hash"):
                verifier.verify_report(path, "0" * 64, replay=False)

    def test_production_parameters_are_frozen(self) -> None:
        self.assertEqual(
            verifier.EXPECTED_SEARCH,
            {
                "orders": 4,
                "cap": 2_000_000,
                "slices": 24,
                "zero_circle_only": False,
                "include_bulk_order": False,
            },
        )
        self.assertEqual(verifier.EXPECTED_SELECTOR["index"], 3_936_177)

    def test_focused_wrapper_hash_is_frozen_in_checker(self) -> None:
        filename, expected = verifier.EXPECTED_SOURCES["focused_wrapper"]
        self.assertEqual(filename, "run_d6_interval_k7_3936177.py")
        self.assertEqual(expected, verifier.file_sha256(verifier.ROOT / filename))

    def test_checker_imports_neither_launch_wrapper(self) -> None:
        source = Path(verifier.__file__).read_text(encoding="utf-8")
        self.assertNotIn("import run_d6_interval_k7_3936177", source)
        self.assertNotIn("import run_d6_interval_v6", source)

    def test_all_checker_dependency_hashes_match(self) -> None:
        for _key, (filename, expected) in verifier.EXPECTED_SOURCES.items():
            self.assertEqual(verifier.file_sha256(verifier.ROOT / filename), expected)
        for _key, (filename, expected) in verifier.EXPECTED_INPUTS.items():
            self.assertEqual(verifier.file_sha256(verifier.ROOT / filename), expected)


if __name__ == "__main__":
    unittest.main()
