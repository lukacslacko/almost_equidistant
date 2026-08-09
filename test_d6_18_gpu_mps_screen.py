#!/usr/bin/env python3
"""Fast deterministic controls for the heuristic d=6 MPS screen."""

from __future__ import annotations

import gzip
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

import run_d6_18_gpu_mps_screen as screen
import verify_d6_18_gpu_mps_screen as verifier


ROOT = Path(__file__).resolve().parent


class GPUMPSScreenTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.corpus = json.loads(
            (ROOT / "d6_residue_18_deletions.json").read_text(encoding="utf-8")
        )
        cls.raw_records = cls.corpus["unique_deletions"]

    def test_all_graphs_have_deterministic_k6_gauge(self) -> None:
        compatible = 0
        k7 = 0
        for index, raw in enumerate(self.raw_records):
            record = screen.reorder_record(index, raw)
            seed = record["seed_old_vertices"]
            self.assertIn(len(seed), (6, 7))
            self.assertEqual(len(seed), record["gauge_size"])
            self.assertEqual(seed, sorted(seed))
            for first in seed:
                for second in seed:
                    if first > second:
                        self.assertTrue(raw["adjacency"][first] & (1 << second))
            compatible += record["standard_compatible"]
            k7 += record["gauge_size"] == 7
        self.assertEqual(compatible, 14)
        self.assertEqual(k7, 1616)

    def test_standard_controls_align_and_are_exact_embeddings(self) -> None:
        simplex = screen.simplex_coordinates()
        controls = 0
        for index, raw in enumerate(self.raw_records):
            if not raw["standard18_compatible"]:
                continue
            controls += 1
            record = screen.reorder_record(index, raw)
            aligned = record["standard_aligned"]
            self.assertIsNotNone(aligned)
            np.testing.assert_allclose(aligned[:6], simplex, atol=1e-10, rtol=0)
            rms, maximum, minimum = screen.numpy_metrics(aligned, record["adjacency"])
            self.assertLess(rms, 1e-14)
            self.assertLess(maximum, 1e-13)
            self.assertGreater(minimum, 0.1)
            self.assertTrue(verifier.exact_standard_control(raw))
        self.assertEqual(controls, 14)

    def test_random_streams_are_batch_independent(self) -> None:
        selected = [
            screen.reorder_record(index, self.raw_records[index])
            for index in (17, 352, 1000)
        ]
        together = screen.initial_outside(
            selected, 7, 3, 600180777, 0.7, True, "seed_sphere"
        )
        separate = np.concatenate(
            [
                screen.initial_outside(
                    [record], 7, 3, 600180777, 0.7, True,
                    "seed_sphere"
                )
                for record in selected
            ],
            axis=0,
        )
        self.assertTrue(np.array_equal(together, separate))

    def test_atomic_gzip_round_trip_and_payload_hash(self) -> None:
        payload = {
            "schema": "d6-18-gpu-mps-checkpoint-v1",
            "proof_status": "HEURISTIC_ONLY_NO_REJECTIONS",
            "config_sha256": "test",
            "range": [0, 1],
            "rows": [{"index": 0}],
        }
        payload["checkpoint_sha256"] = screen.stable_hash(payload)
        with tempfile.TemporaryDirectory() as directory_name:
            path = Path(directory_name) / "checkpoint.json.gz"
            screen.atomic_gzip_json(path, payload)
            with gzip.open(path, "rt", encoding="ascii") as stream:
                loaded = json.load(stream)
        claimed = loaded.pop("checkpoint_sha256")
        self.assertEqual(screen.stable_hash(loaded), claimed)

    def test_mps_metric_replay_allows_only_binary32_reduction_noise(self) -> None:
        # Worst observed minimum-distance discrepancy in the complete sweep:
        # identical serialized float32 coordinates, evaluated by MPS float32
        # versus the independent Python float64 checker.
        self.assertTrue(
            verifier.close_mps_metric(0.11541603273264017, 0.11541681736707687)
        )
        self.assertTrue(
            verifier.close_mps_metric(1.643707678855577e-7, 1.1920928955078125e-7)
        )
        self.assertFalse(verifier.close_mps_metric(0.1154160, 0.1155160))
        self.assertFalse(
            verifier.close_metric(1.643707678855577e-7, 1.1920928955078125e-7)
        )

    def test_cpu_seeded_control_smoke(self) -> None:
        index = next(
            index for index, raw in enumerate(self.raw_records)
            if raw["standard18_compatible"]
        )
        record = screen.reorder_record(index, self.raw_records[index])
        result = screen.optimize_wave(
            [record], restarts=2, wave=0, steps=1, learning_rate=0.03,
            seedbase=600180777, initialization_scale=0.7,
            collision_distance=0.12, collision_weight=0.05,
            device_name="cpu", seed_standard_controls=True,
            initialization="seed_sphere",
        )
        control = result["initial_controls"][index]
        self.assertLess(control["edge_rms"], 2e-6)
        self.assertGreater(control["minimum_distance"], 0.01)

    def test_analytic_lm_keeps_exact_distinct_control(self) -> None:
        index = next(
            index for index, raw in enumerate(self.raw_records)
            if raw["standard18_compatible"]
        )
        record = screen.reorder_record(index, self.raw_records[index])
        result = screen.cpu_lm_refine(
            record["standard_aligned"], record["adjacency"],
            record["gauge_size"], 50, 0.05, 100.0,
        )
        self.assertLess(result["edge_rms"], 1e-12)
        self.assertGreater(result["minimum_distance"], 0.01)

    def test_seed_sphere_initialization_makes_seed_edges_exact(self) -> None:
        selected = []
        for wanted in (6, 7):
            index, raw = next(
                (index, raw) for index, raw in enumerate(self.raw_records)
                if screen.reorder_record(index, raw)["gauge_size"] == wanted
            )
            selected.append(screen.reorder_record(index, raw))
        values = screen.initial_outside(
            selected, 5, 4, 1234567, 0.7, False, "seed_sphere"
        )
        fixed = screen.simplex_coordinates()
        for local, record in enumerate(selected):
            for restart in range(5):
                points = np.concatenate((fixed, values[local, restart]), axis=0)
                if record["gauge_size"] == 7:
                    np.testing.assert_allclose(
                        points[6], screen.full_simplex_coordinates()[6], atol=1e-7
                    )
                for vertex in range(record["gauge_size"], 18):
                    for seed in range(record["gauge_size"]):
                        if record["adjacency"][vertex, seed]:
                            squared = np.sum((points[vertex] - points[seed]) ** 2)
                            self.assertAlmostEqual(float(squared), 1.0, places=5)

    def test_completed_report_verifies_when_present(self) -> None:
        report = ROOT / "d6_18_gpu_mps_screen_report.json"
        if not report.exists():
            self.skipTest("full GPU report has not been assembled")
        result = verifier.verify(report, ROOT / "d6_residue_18_deletions.json")
        self.assertEqual(result["status"], "PASS", result["failed_checks"])


if __name__ == "__main__":
    unittest.main()
