#!/usr/bin/env python3
"""Focused, accelerator-free controls for the 181-class all-K6-gauge screen."""

from __future__ import annotations

import gzip
import tempfile
import unittest
from pathlib import Path

import numpy as np

import run_d6_18_gpu_mps_cover_v7 as runner
import run_d6_18_gpu_mps_screen as screen


class CoverV7AllK6GaugeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.dependencies = runner.dependency_boundary()
        cls.campaign, cls.provenance, cls.deletions = runner.load_campaign()
        cls.metadata = runner.gauge_metadata(cls.campaign)

    def test_exact_frozen_campaign_and_all_k6_manifest(self) -> None:
        self.assertEqual(len(self.campaign), runner.EXPECTED_CLASS_COUNT)
        audit = runner.audit_full_gauge_manifest(self.metadata)
        self.assertEqual(audit["gauge_count"], runner.EXPECTED_GAUGE_COUNT)
        self.assertEqual(
            audit["gauge_manifest_sha256"],
            runner.EXPECTED_GAUGE_MANIFEST_SHA256,
        )
        self.assertEqual(
            audit["variant_ids_sha256"], runner.EXPECTED_VARIANT_IDS_SHA256
        )
        self.assertEqual(
            audit["population_counts"], runner.EXPECTED_GAUGE_POPULATIONS
        )
        self.assertEqual(
            audit["known_positive_gauge_count"],
            runner.EXPECTED_POSITIVE_GAUGES,
        )

    def test_every_variant_is_a_distinct_required_edge_k6(self) -> None:
        by_class = {item["class_index"]: item for item in self.campaign}
        seen = set()
        for info in self.metadata:
            self.assertNotIn(info["variant_id"], seen)
            seen.add(info["variant_id"])
            self.assertEqual(
                info["variant_id"],
                info["class_index"] * runner.VARIANT_STRIDE
                + info["gauge_ordinal"],
            )
            seed = info["seed"]
            self.assertEqual(len(seed), 6)
            self.assertEqual(seed, sorted(seed))
            rows = by_class[info["class_index"]]["raw"]["adjacency"]
            for position, first in enumerate(seed):
                for second in seed[:position]:
                    self.assertTrue(rows[first] & (1 << second))

    def test_known_standard_positive_aligns_in_every_gauge(self) -> None:
        positive = [
            item
            for item in self.campaign
            if item["class_index"] == runner.EXPECTED_POSITIVE_CLASS_INDEX
        ]
        self.assertEqual(len(positive), 1)
        variants, metadata = runner.make_variants(positive)
        self.assertEqual(len(variants), runner.EXPECTED_POSITIVE_GAUGES)
        self.assertEqual(len(metadata), runner.EXPECTED_POSITIVE_GAUGES)
        for record in variants:
            self.assertTrue(record["standard_compatible"])
            aligned = record["standard_aligned"]
            self.assertIsNotNone(aligned)
            rms, maximum, minimum = screen.numpy_metrics(
                aligned, record["adjacency"]
            )
            self.assertLess(rms, 1e-14)
            self.assertLess(maximum, 1e-13)
            self.assertGreater(minimum, 0.1)

    def test_optional_nonedges_do_not_enter_unit_edge_metrics(self) -> None:
        # Only {0,1} is a required unit edge.  Moving vertex 2 from a collision
        # to an arbitrary non-unit location changes distinctness but cannot
        # change either required-edge residual metric.
        adjacency = np.zeros((18, 18), dtype=np.float32)
        adjacency[0, 1] = adjacency[1, 0] = 1.0
        points = np.zeros((18, 6), dtype=np.float64)
        points[1, 0] = 1.0
        for vertex in range(3, 18):
            points[vertex, 0] = 10.0 * vertex
        collided = screen.numpy_metrics(points, adjacency)
        points[2, 0] = 7.25
        separated = screen.numpy_metrics(points, adjacency)
        self.assertEqual(collided[:2], separated[:2])
        self.assertEqual(collided[:2], (0.0, 0.0))
        self.assertLess(collided[2], separated[2])
        self.assertEqual(
            self.provenance["candidate_nonedges"],
            "unconstrained and may also have distance one",
        )

    def test_variant_random_stream_is_chunk_independent(self) -> None:
        selected = self.campaign[:2]
        variants, _ = runner.make_variants(selected)
        variants = [variants[0], variants[-1]]
        together = screen.initial_outside(
            variants, 3, 2, 6181814099, 0.7, True, "seed_sphere"
        )
        separate = np.concatenate(
            [
                screen.initial_outside(
                    [record], 3, 2, 6181814099, 0.7, True, "seed_sphere"
                )
                for record in variants
            ],
            axis=0,
        )
        self.assertTrue(np.array_equal(together, separate))

    def test_defaults_use_mps_and_only_one_cpu_refiner(self) -> None:
        args = runner.parser().parse_args(["--selection-only"])
        self.assertEqual(args.device, "mps")
        self.assertEqual(args.cpu_workers, 1)
        self.assertEqual(args.total_restarts, 8)
        self.assertEqual(args.lm_top_k, 2)
        config = runner.configuration(args, self.campaign, "source")
        self.assertEqual(config["mathematical_rejections"], 0)
        self.assertEqual(config["mathematical_realizability_conclusions"], 0)
        self.assertIn("unconstrained", config["candidate_nonedges"])

    def test_checkpoint_round_trip_checks_identity_and_hash(self) -> None:
        metadata = [self.metadata[0]]
        payload = {
            "schema": "d6-18-gpu-mps-cover-v7-checkpoint-v1",
            "proof_status": "HEURISTIC_ONLY_NO_REJECTIONS",
            "config_sha256": "config",
            "class_range": [0, 1],
            "class_indices": [metadata[0]["class_index"]],
            "gauge_count": 1,
            "gauge_metadata": metadata,
            "rows": [{"variant_id": metadata[0]["variant_id"]}],
            "witnesses": [],
            "retained_endpoints": [],
            "lm_attempts": [],
            "wave_profiles": [],
            "lm_profile": {},
            "mathematical_rejections": 0,
            "mathematical_realizability_conclusions": 0,
        }
        payload["checkpoint_sha256"] = screen.stable_hash(payload)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "checkpoint.json.gz"
            screen.atomic_gzip_json(path, payload)
            loaded = runner.load_checkpoint(
                path,
                "config",
                0,
                1,
                [metadata[0]["class_index"]],
            )
            self.assertEqual(loaded["checkpoint_sha256"], payload["checkpoint_sha256"])
            with gzip.open(path, "rt", encoding="ascii") as stream:
                self.assertIn("checkpoint_sha256", stream.read())


if __name__ == "__main__":
    unittest.main()
