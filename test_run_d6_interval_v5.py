#!/usr/bin/env python3
"""Focused controls for the v5-bound interval selection wrapper."""

from __future__ import annotations

import copy
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import run_d6_interval_v5 as wrapper


class V5SelectionControls(unittest.TestCase):
    def test_exact_manifest_reconstructs_all_current_graphs(self) -> None:
        graphs, provenance = wrapper.load_v5_graphs(
            wrapper.ROOT / "d6_current_residue_manifest_v5.json",
            wrapper.ROOT / "d6_current_residue_manifest_v5_verification.json",
        )
        self.assertEqual(len(graphs), 911)
        self.assertEqual(
            {name: sum(graph["population"] == name for graph in graphs) for name in ("K7", "K6")},
            {"K7": 155, "K6": 756},
        )
        self.assertEqual(len({graph["index"] for graph in graphs}), 911)
        self.assertEqual(provenance["verification"]["status"], "PASS")

    def test_unpinned_manifest_is_rejected_before_parsing(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            bad = Path(directory) / "manifest.json"
            bad.write_text("{}\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "manifest hash"):
                wrapper.load_v5_graphs(
                    bad,
                    wrapper.ROOT / "d6_current_residue_manifest_v5_verification.json",
                )

    def test_embedded_graph_tamper_fails_structural_hash(self) -> None:
        manifest_path = wrapper.ROOT / "d6_current_residue_manifest_v5.json"
        payload = json.loads(manifest_path.read_text(encoding="utf-8"))
        payload = copy.deepcopy(payload)
        payload["classes"]["K7"]["graphs"][0]["adjacency"][0] ^= 2
        with tempfile.TemporaryDirectory() as directory:
            tampered = Path(directory) / "manifest.json"
            tampered.write_text(json.dumps(payload), encoding="utf-8")
            verification_payload = json.loads(
                (wrapper.ROOT / "d6_current_residue_manifest_v5_verification.json").read_text(
                    encoding="utf-8"
                )
            )
            verification_payload["manifest"]["sha256"] = wrapper.engine.sha256(tampered)
            verification = Path(directory) / "verification.json"
            verification.write_text(json.dumps(verification_payload), encoding="utf-8")
            with (
                patch.object(wrapper, "EXPECTED_MANIFEST_SHA256", wrapper.engine.sha256(tampered)),
                patch.object(
                    wrapper,
                    "EXPECTED_VERIFICATION_SHA256",
                    wrapper.engine.sha256(verification),
                ),
                self.assertRaisesRegex(ValueError, "class graphs hash"),
            ):
                wrapper.load_v5_graphs(tampered, verification)


if __name__ == "__main__":
    unittest.main()
