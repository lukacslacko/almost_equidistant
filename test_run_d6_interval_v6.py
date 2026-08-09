#!/usr/bin/env python3
"""Focused controls for the provenance-pinned v6 K7 interval wrapper."""

from __future__ import annotations

import copy
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import run_d6_interval_v6 as wrapper


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


class V6SelectionControls(unittest.TestCase):
    def load(self) -> tuple[list[dict], dict]:
        return wrapper.load_v6_graphs(
            wrapper.ROOT / "d6_current_residue_manifest_v6.json",
            wrapper.ROOT / "d6_current_residue_manifest_v6_verification.json",
        )

    def test_exact_manifest_reconstructs_all_644_current_graphs(self) -> None:
        graphs, provenance = self.load()
        self.assertEqual(len(graphs), 644)
        self.assertEqual(
            {
                name: sum(graph["population"] == name for graph in graphs)
                for name in ("K7", "K6")
            },
            {"K7": 19, "K6": 625},
        )
        self.assertEqual(len({graph["index"] for graph in graphs}), 644)
        self.assertEqual(provenance["verification"]["status"], "PASS")
        self.assertEqual(
            provenance["certificate_accounting"]["interval_trust_assumptions"],
            wrapper.EXPECTED_TRUST,
        )

    def test_full_k7_selection_is_exactly_the_frozen_19(self) -> None:
        graphs, _provenance = self.load()
        selected = wrapper.select_k7(
            graphs,
            shards=1,
            shard=0,
            sample=None,
            sample_seed=600_019_006,
        )
        indices = [graph["index"] for graph in selected]
        self.assertEqual(indices, EXPECTED_K7_INDICES)
        self.assertEqual(wrapper.object_hash(indices), wrapper.EXPECTED_CLASS_HASHES["K7"])
        self.assertEqual(
            [graph["ordinal"] for graph in selected], list(range(len(selected)))
        )

    def test_sharding_partitions_and_sampling_is_deterministic(self) -> None:
        graphs, _provenance = self.load()
        shards = [
            wrapper.select_k7(
                graphs,
                shards=3,
                shard=shard,
                sample=None,
                sample_seed=600_019_006,
            )
            for shard in range(3)
        ]
        shard_sets = [{graph["index"] for graph in selected} for selected in shards]
        self.assertEqual(set.union(*shard_sets), set(EXPECTED_K7_INDICES))
        self.assertTrue(all(left.isdisjoint(right) for left, right in zip(shard_sets, shard_sets[1:])))
        first = wrapper.select_k7(
            graphs,
            shards=1,
            shard=0,
            sample=7,
            sample_seed=12345,
        )
        second = wrapper.select_k7(
            graphs,
            shards=1,
            shard=0,
            sample=7,
            sample_seed=12345,
        )
        self.assertEqual(first, second)
        self.assertEqual(len(first), 7)

    def test_explicit_index_selection_is_exact_and_validated(self) -> None:
        graphs, _provenance = self.load()
        requested = frozenset({3945557, 316173, 2593240})
        selected = wrapper.select_k7(
            graphs,
            shards=1,
            shard=0,
            sample=None,
            sample_seed=600_019_006,
            indices=requested,
        )
        self.assertEqual(
            [graph["index"] for graph in selected],
            [316173, 2593240, 3945557],
        )
        self.assertEqual([graph["ordinal"] for graph in selected], [0, 1, 2])
        with self.assertRaisesRegex(ValueError, "absent from v6 K7 class"):
            wrapper.select_k7(
                graphs,
                shards=1,
                shard=0,
                sample=None,
                sample_seed=600_019_006,
                indices=frozenset({123456789}),
            )

    def test_unpinned_manifest_is_rejected_before_parsing(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            bad = Path(directory) / "manifest.json"
            bad.write_text("{}\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "manifest hash"):
                wrapper.load_v6_graphs(
                    bad,
                    wrapper.ROOT
                    / "d6_current_residue_manifest_v6_verification.json",
                )

    def test_embedded_graph_tamper_fails_structural_hash(self) -> None:
        manifest_path = wrapper.ROOT / "d6_current_residue_manifest_v6.json"
        payload = copy.deepcopy(json.loads(manifest_path.read_text(encoding="utf-8")))
        payload["classes"]["K7"]["graphs"][0]["adjacency"][0] ^= 2
        with tempfile.TemporaryDirectory() as directory:
            tampered = Path(directory) / "manifest.json"
            tampered.write_text(json.dumps(payload), encoding="utf-8")
            verification_payload = json.loads(
                (
                    wrapper.ROOT
                    / "d6_current_residue_manifest_v6_verification.json"
                ).read_text(encoding="utf-8")
            )
            verification_payload["manifest"]["sha256"] = wrapper.engine.sha256(
                tampered
            )
            verification = Path(directory) / "verification.json"
            verification.write_text(json.dumps(verification_payload), encoding="utf-8")
            with (
                patch.object(
                    wrapper,
                    "EXPECTED_MANIFEST_SHA256",
                    wrapper.engine.sha256(tampered),
                ),
                patch.object(
                    wrapper,
                    "EXPECTED_VERIFICATION_SHA256",
                    wrapper.engine.sha256(verification),
                ),
                self.assertRaisesRegex(ValueError, "class graphs hash"),
            ):
                wrapper.load_v6_graphs(tampered, verification)

    def test_changed_interval_trust_statement_is_rejected(self) -> None:
        manifest_path = wrapper.ROOT / "d6_current_residue_manifest_v6.json"
        payload = copy.deepcopy(json.loads(manifest_path.read_text(encoding="utf-8")))
        payload["certificate_accounting"]["interval"]["trust_assumptions"][
            "transcendentals"
        ] = "weaker unrecorded libm assumption"
        with tempfile.TemporaryDirectory() as directory:
            tampered = Path(directory) / "manifest.json"
            tampered.write_text(json.dumps(payload), encoding="utf-8")
            verification_payload = json.loads(
                (
                    wrapper.ROOT
                    / "d6_current_residue_manifest_v6_verification.json"
                ).read_text(encoding="utf-8")
            )
            verification_payload["manifest"]["sha256"] = wrapper.engine.sha256(
                tampered
            )
            verification = Path(directory) / "verification.json"
            verification.write_text(json.dumps(verification_payload), encoding="utf-8")
            with (
                patch.object(
                    wrapper,
                    "EXPECTED_MANIFEST_SHA256",
                    wrapper.engine.sha256(tampered),
                ),
                patch.object(
                    wrapper,
                    "EXPECTED_VERIFICATION_SHA256",
                    wrapper.engine.sha256(verification),
                ),
                self.assertRaisesRegex(ValueError, "trust assumptions"),
            ):
                wrapper.load_v6_graphs(tampered, verification)


if __name__ == "__main__":
    unittest.main()
