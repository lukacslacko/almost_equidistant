#!/usr/bin/env python3
"""Accelerator-free boundary and geometry tests for the n=19 MPS screen."""

from __future__ import annotations

import unittest
from unittest import mock
import tempfile
from pathlib import Path

import numpy as np

import run_d6_19_gpu_mps_residue_v1 as runner


class D619MPSResidueV1Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.dependencies = runner.dependency_boundary()
        cls.campaign, cls.provenance = runner.load_campaign()
        cls.metadata = runner.gauge_metadata(cls.campaign)

    def test_published_physical_and_active_boundaries(self) -> None:
        self.assertEqual(len(self.campaign), 261)
        active = [x for x in self.campaign if x["campaign_role"] == "active_unresolved"]
        negative = [
            x for x in self.campaign if x["campaign_role"] == "exact_negative_calibration"
        ]
        self.assertEqual(len(active), 260)
        self.assertEqual([x["candidate_index"] for x in negative], [3936435])
        self.assertEqual(
            runner.kernel.stable_hash([x["candidate_index"] for x in active]),
            runner.EXPECTED_CAMPAIGN_INDICES_SHA256,
        )
        self.assertEqual(self.provenance["active_unresolved_graph_count"], 260)

    def test_frozen_all_gauge_and_intrinsic_dof_census(self) -> None:
        audit = runner.audit_full_gauge_manifest(self.metadata)
        self.assertEqual(audit["gauge_count"], 7721)
        self.assertEqual(
            audit["population_counts"],
            {"K7_active": 34, "K7_exact_negative_calibration": 3, "K6_only_active": 7684},
        )
        self.assertEqual(audit["active_intrinsic_dof_ranges"]["K7"], [25, 31])
        self.assertEqual(audit["active_intrinsic_dof_ranges"]["K6_only"], [27, 36])
        self.assertEqual(audit["gauges_with_forced_S0_reflections"], 6)

    def test_every_frozen_variant_builds_without_collision(self) -> None:
        built = 0
        for item in self.campaign:
            variants, metadata = runner.make_variants([item])
            self.assertEqual(len(variants), len(metadata))
            built += len(variants)
        self.assertEqual(built, 7721)

    def test_manifold_enforces_every_seed_link_exactly(self) -> None:
        examples = [
            next(item for item in self.campaign if item["source_class"] == "K7"),
            next(item for item in self.campaign if item["source_class"] == "K6_only"),
        ]
        for item in examples:
            variant = runner.make_variants([item])[0][0]
            raw = runner.initial_manifold_parameters(
                [variant], 3, 2, 1234567, 0.7, False
            )[0]
            points = runner.manifold_points_numpy(variant, raw)
            for realization in points:
                for outside in range(variant["gauge_size"], runner.N):
                    for fixed in range(variant["gauge_size"]):
                        if variant["adjacency"][outside, fixed]:
                            squared = float(
                                np.sum((realization[outside] - realization[fixed]) ** 2)
                            )
                            self.assertAlmostEqual(squared, 1.0, places=11)

    def test_cpu_adam_projection_preserves_seed_links(self) -> None:
        item = next(x for x in self.campaign if x["source_class"] == "K6_only")
        variant = runner.make_variants([item])[0][0]
        result = runner.optimize_wave_manifold_n19(
            [variant],
            restarts=2,
            wave=0,
            steps=3,
            learning_rate=0.03,
            seedbase=99117,
            initialization_scale=0.7,
            collision_distance=0.12,
            collision_weight=0.05,
            device_name="cpu",
            seed_standard_controls=False,
            initialization="seed_sphere",
        )
        for points in result["coordinates"][0]:
            for outside in range(variant["gauge_size"], runner.N):
                for fixed in range(variant["gauge_size"]):
                    if variant["adjacency"][outside, fixed]:
                        squared = float(np.sum((points[outside] - points[fixed]) ** 2))
                        self.assertLess(abs(squared - 1.0), 2e-6)

    def test_k7_s0_uses_distinct_reflected_root(self) -> None:
        info = next(x for x in self.metadata if x["s0_reflection_count"])
        item = next(x for x in self.campaign if x["candidate_index"] == info["candidate_index"])
        variant = runner.reorder_record(item, info)
        reflected = [
            p
            for p in variant["manifold_parameters"]
            if p.get("reason") == "distinct_reflection_of_missing_K7_seed"
        ]
        self.assertEqual(len(reflected), info["s0_reflection_count"])
        seed = runner.kernel.full_simplex_coordinates()
        for parameters in reflected:
            missing = seed[parameters["missing_seed_slot"]]
            self.assertGreater(np.linalg.norm(parameters["point"] - missing), 0.5)

    @staticmethod
    def synthetic_item(rows: list[int], label: str) -> dict:
        return {
            "candidate_index": 9_000_001,
            "source_class": "K7",
            "campaign_role": "test_control",
            "class_id": label,
            "adjacency": rows,
            "edges": runner.edge_count(rows),
            "control_coordinates": None,
        }

    @staticmethod
    def add_edge(rows: list[int], first: int, second: int) -> None:
        rows[first] |= 1 << second
        rows[second] |= 1 << first

    def test_forced_twin_collision_fails_closed(self) -> None:
        rows = [0] * runner.N
        for first in range(7):
            for second in range(first):
                self.add_edge(rows, first, second)
        for twin in (7, 8):
            for fixed in range(6):
                self.add_edge(rows, twin, fixed)
        info = {
            "variant_id": 1,
            "candidate_index": 9_000_001,
            "source_class": "K7",
            "campaign_role": "test_control",
            "gauge_ordinal": 0,
            "gauge_size": 7,
            "seed": list(range(7)),
        }
        with self.assertRaisesRegex(ValueError, "forced-twin collision"):
            runner.reorder_record(self.synthetic_item(rows, "forced_twins"), info)

    def test_k8_fails_closed(self) -> None:
        rows = [0] * runner.N
        for first in range(8):
            for second in range(first):
                self.add_edge(rows, first, second)
        info = {
            "variant_id": 2,
            "candidate_index": 9_000_001,
            "source_class": "K7",
            "campaign_role": "test_control",
            "gauge_ordinal": 0,
            "gauge_size": 7,
            "seed": list(range(7)),
        }
        with self.assertRaisesRegex(ValueError, "K8"):
            runner.reorder_record(self.synthetic_item(rows, "K8"), info)

    def test_optional_nonedge_never_enters_unit_edge_metric(self) -> None:
        adjacency = np.zeros((runner.N, runner.N), dtype=np.float32)
        adjacency[0, 1] = adjacency[1, 0] = 1.0
        points = np.zeros((runner.N, runner.D), dtype=np.float64)
        points[1, 0] = 1.0
        for vertex in range(3, runner.N):
            points[vertex, 0] = 10.0 * vertex
        with runner.n19_kernel_context():
            collided = runner.kernel.numpy_metrics(points, adjacency)
            points[2, 0] = 7.25
            separated = runner.kernel.numpy_metrics(points, adjacency)
        self.assertEqual(collided[:2], (0.0, 0.0))
        self.assertEqual(collided[:2], separated[:2])
        self.assertLess(collided[2], separated[2])

    def test_exact_k6_and_k7_positive_control_inputs(self) -> None:
        variants, metadata = runner.positive_control_variants()
        self.assertEqual(len(variants), len(metadata), 3)
        self.assertEqual([x["gauge_size"] for x in metadata], [6, 6, 7])
        raw = runner.initial_manifold_parameters(variants, 1, 0, 778899, 0.7, True)
        for local, variant in enumerate(variants):
            points = runner.manifold_points_numpy(variant, raw[local])[0]
            with runner.n19_kernel_context():
                rms, maximum, minimum = runner.kernel.numpy_metrics(
                    points, variant["adjacency"]
                )
            # The exact coordinates are encoded through the final float32
            # manifold parameter tensor before this evaluation.
            self.assertLess(rms, 2e-7)
            self.assertLess(maximum, 5e-7)
            self.assertGreater(minimum, 0.1)

    def test_defaults_fill_gpu_and_limit_lm_process_threads(self) -> None:
        args = runner.parser().parse_args(["--selection-only"])
        self.assertEqual(args.device, "mps")
        self.assertEqual(args.total_restarts, 16)
        self.assertEqual(args.steps, 1000)
        self.assertEqual(args.lm_top_k, 2)
        self.assertEqual(args.cpu_refine_max_nfev, 1500)
        self.assertEqual(args.cpu_workers, 10)
        self.assertEqual(args.class_chunk_size, 8)

    def test_positive_control_gate_requires_final_kernel_recovery(self) -> None:
        rows = [
            {
                "initial_standard_control_edge_rms": 1e-8,
                "initial_standard_control_minimum_distance": 0.8,
                "lm_candidate": recovered,
            }
            for recovered in (True, True, False)
        ]
        config = {
            "distinctness_threshold": 0.01,
            "lm_top_k": 2,
            "cpu_refine_max_nfev": 1500,
        }
        audit = runner.audit_positive_controls({"rows": rows}, config)
        self.assertEqual(audit["status"], "FAIL_FINAL_KERNEL_RECOVERY")
        self.assertEqual(audit["lm_candidate_recoveries"], 2)
        rows[-1]["lm_candidate"] = True
        self.assertEqual(
            runner.audit_positive_controls({"rows": rows}, config)["status"],
            "PASS",
        )
        config["lm_top_k"] = 0
        self.assertEqual(
            runner.audit_positive_controls({"rows": rows}, config)["status"],
            "PASS_SEEDED_INPUT_ONLY",
        )

    def test_lm_thread_environment_is_set_before_spawn(self) -> None:
        with mock.patch.dict(runner.os.environ, {}, clear=True):
            runner.set_lm_thread_environment()
            self.assertEqual(
                {name: runner.os.environ.get(name) for name in runner.LM_THREAD_ENVIRONMENT},
                runner.LM_THREAD_ENVIRONMENT,
            )

    def test_completed_output_cannot_cross_configuration_boundary(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "report.json"
            runner.kernel.atomic_json(path, {"config_sha256": "first"})
            runner.guard_output_target(path, "first", compressed=False)
            with self.assertRaisesRegex(ValueError, "another configuration"):
                runner.guard_output_target(path, "second", compressed=False)


if __name__ == "__main__":
    unittest.main()
