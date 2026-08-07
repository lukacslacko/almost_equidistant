#!/usr/bin/env python3
"""Focused controls for the independent v6 K7 interval verifier."""

from __future__ import annotations

import copy
import inspect
import tempfile
import unittest
from collections import Counter
from pathlib import Path

import run_d6_interval_v6 as wrapper
import verify_d6_interval_v6 as verifier


EXPECTED_K7_INDICES = [
    316173,
    2581209,
    2592657,
    2593240,
    3595554,
    3648882,
    3729907,
    3785980,
    3888410,
    3935560,
    3936177,
    3936310,
    3936435,
    3945490,
    3945555,
    3945557,
    3945564,
    3947605,
    3949382,
]


def production_selection_configuration() -> dict:
    all_graphs, provenance = wrapper.load_v6_graphs(
        wrapper.ROOT / "d6_current_residue_manifest_v6.json",
        wrapper.ROOT / "d6_current_residue_manifest_v6_verification.json",
    )
    graphs = wrapper.select_k7(
        all_graphs,
        shards=1,
        shard=0,
        sample=None,
        sample_seed=600_019_006,
    )
    return {
        "selection": {
            "graphs": len(graphs),
            "population_counts": dict(
                Counter(graph["population"] for graph in graphs)
            ),
            "indices_sha256": wrapper.object_hash(
                [graph["index"] for graph in graphs]
            ),
            "selection_layer": provenance,
        },
        "sampling": {
            "sample": None,
            "sample_seed": 600_019_006,
            "population": "K7",
            "shards": 1,
            "shard": 0,
        },
    }


class V6IntervalVerifierControls(unittest.TestCase):
    def test_independent_reconstruction_returns_exact_frozen_19(self) -> None:
        graphs = verifier.reconstruct_v6(production_selection_configuration())
        self.assertEqual([graph["index"] for graph in graphs], EXPECTED_K7_INDICES)
        self.assertEqual(
            verifier.object_sha256([graph["index"] for graph in graphs]),
            verifier.EXPECTED_CLASS_HASHES["K7"],
        )
        self.assertTrue(all(graph["population"] == "K7" for graph in graphs))

    def test_selection_trust_tier_tamper_is_rejected(self) -> None:
        configuration = copy.deepcopy(production_selection_configuration())
        configuration["selection"]["selection_layer"]["certificate_accounting"][
            "interval_non_KILLED_used"
        ] = True
        with self.assertRaisesRegex(AssertionError, "selection trust tiers"):
            verifier.reconstruct_v6(configuration)

    def test_non_k7_population_selector_is_rejected(self) -> None:
        configuration = copy.deepcopy(production_selection_configuration())
        configuration["sampling"]["population"] = "K6"
        with self.assertRaisesRegex(AssertionError, "K7-only selector"):
            verifier.reconstruct_v6(configuration)

    def test_report_hash_tamper_is_rejected_first(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "report.json"
            path.write_text("{}\n", encoding="utf-8")
            with self.assertRaisesRegex(AssertionError, "report hash"):
                verifier.verify_report(path, "0" * 64, replay=False)

    def test_production_cap_and_worker_defaults_are_frozen(self) -> None:
        signature = inspect.signature(verifier.verify_report)
        self.assertEqual(signature.parameters["expected_cap"].default, 500_000)
        self.assertEqual(signature.parameters["expected_workers"].default, 9)

    def test_wrapper_hash_is_frozen_in_checker(self) -> None:
        filename, expected = verifier.EXPECTED_SOURCES["wrapper"]
        self.assertEqual(filename, "run_d6_interval_v6.py")
        self.assertEqual(expected, verifier.file_sha256(verifier.ROOT / filename))

    def test_checker_does_not_import_production_wrapper(self) -> None:
        source = (verifier.ROOT / "verify_d6_interval_v6.py").read_text(
            encoding="utf-8"
        )
        self.assertNotIn("import run_d6_interval_v6", source)


if __name__ == "__main__":
    unittest.main()
