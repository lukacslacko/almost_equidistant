#!/usr/bin/env python3
"""Focused non-kernel controls for the v7 dense 18-deletion wrapper."""

from __future__ import annotations

import json
import subprocess
import unittest
from collections import Counter, defaultdict
from unittest.mock import patch

import run_d6_interval_18_cover_v7 as wrapper


EXPECTED_K7 = [
    316173,
    2581209,
    3648882,
    3729907,
    3935560,
    3936310,
    3936435,
    3945490,
    3945555,
    3945557,
    3945564,
    3947605,
]


class DenseDeletionCoverControls(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        k6_report = json.loads(
            (
                wrapper.ROOT
                / "d6_k6_opposite_ray_rank_conjunction_report.json"
            ).read_text(encoding="utf-8")
        )
        cls.parents = [
            {"class": "K7", "index": index} for index in EXPECTED_K7
        ] + [
            {"class": "K6_only", "index": index}
            for index in k6_report["ordered_residue_indices"]
        ]
        cls.deletions, cls.deletion_provenance = wrapper.load_deletion_boundary(
            wrapper.ROOT / "d6_residue_18_deletions_v2.json",
            wrapper.ROOT / "d6_residue_18_deletions_v2_verification.json",
        )
        cls.base_graphs, cls.cover = wrapper.build_dense_deletion_cover(
            cls.parents, cls.deletions
        )
        cls.graphs, cls.augmentation = wrapper.build_actionable_campaign(
            cls.base_graphs, cls.cover, cls.deletions
        )

    def test_exact_frozen_179_class_cover(self) -> None:
        self.assertEqual(len(self.parents), 263)
        self.assertEqual(len(self.base_graphs), wrapper.EXPECTED_SELECTED_COUNT)
        self.assertEqual(
            self.cover["greedy_gain_histogram"],
            wrapper.EXPECTED_GREEDY_GAIN_HISTOGRAM,
        )
        self.assertEqual(
            self.cover["selected_population_counts"],
            wrapper.EXPECTED_SELECTED_POPULATIONS,
        )
        self.assertEqual(
            self.cover["selected_class_indices_sha256"],
            wrapper.EXPECTED_SELECTED_CLASS_INDICES_SHA256,
        )
        self.assertEqual(
            self.cover["selected_class_ids_sha256"],
            wrapper.EXPECTED_SELECTED_CLASS_IDS_SHA256,
        )
        self.assertEqual(
            self.cover["greedy_steps_sha256"],
            wrapper.EXPECTED_GREEDY_STEPS_SHA256,
        )
        self.assertEqual(
            sum(
                int(gain) * count
                for gain, count in self.cover["greedy_gain_histogram"].items()
            ),
            263,
        )

    def test_actionable_181_class_augmentation(self) -> None:
        self.assertEqual(len(self.graphs), wrapper.EXPECTED_PRODUCTION_COUNT)
        self.assertEqual(
            self.augmentation["production_population_counts"],
            wrapper.EXPECTED_PRODUCTION_POPULATIONS,
        )
        self.assertEqual(
            self.augmentation["production_class_indices_sha256"],
            wrapper.EXPECTED_PRODUCTION_CLASS_INDICES_SHA256,
        )
        self.assertEqual(
            self.augmentation["blocked_parents"],
            ["K6_only:364827", "K6_only:3335955"],
        )
        frozen = [
            {
                "blocked_parent": record["blocked_parent"],
                "class_index": record["class_index"],
                "class_id": record["class_id"],
            }
            for record in self.augmentation["backup_records"]
        ]
        self.assertEqual(frozen, wrapper.EXPECTED_ACTIONABLE_BACKUPS)
        selected_indices = {
            record["class_index"] for record in self.cover["selected_records"]
        }
        production_indices = {graph["index"] for graph in self.graphs}
        self.assertIn(wrapper.EXPECTED_KNOWN_POSITIVE_CLASS_INDEX, production_indices)
        self.assertEqual(
            production_indices - selected_indices,
            {record["class_index"] for record in frozen},
        )
        for record in self.augmentation["backup_records"]:
            self.assertFalse(record["known_realizable_positive_control"])
            self.assertFalse(record["standard18_compatible"])
            adjacency = self.deletions["unique_deletions"][record["class_index"]][
                "adjacency"
            ]
            orders = wrapper.engine.cdriver6.gen_orders(adjacency, 18, kmax=4)
            self.assertTrue(orders)
            self.assertEqual(
                {wrapper.engine.ncircle(adjacency, seed, order) for seed, order in orders},
                {2},
            )

    def test_every_parent_is_covered_by_a_maximum_edge_deletion(self) -> None:
        dense_by_parent = {
            record["parent"]: record for record in self.cover["parents"]["dense_records"]
        }
        covered = set()
        selected = set(self.cover["selected_class_indices_greedy_order"])
        for record in self.cover["selected_records"]:
            class_index = record["class_index"]
            self.assertIn(class_index, selected)
            for parent in record["active_parents"]:
                dense = dense_by_parent[parent]
                eligible = {
                    item["class_index"] for item in dense["dense_classes"]
                }
                self.assertIn(class_index, eligible)
                self.assertEqual(record["edges"], dense["maximum_deletion_edges"])
                covered.add(parent)
        self.assertEqual(covered, set(dense_by_parent))
        self.assertEqual(len(covered), 263)

    def test_tie_break_is_explicitly_largest_class_id(self) -> None:
        # Recompute only the final greedy loop from the committed dense records.
        coverage: dict[int, set[str]] = defaultdict(set)
        class_ids = {
            position: record["class_id"]
            for position, record in enumerate(self.deletions["unique_deletions"])
        }
        for record in self.cover["parents"]["dense_records"]:
            for option in record["dense_classes"]:
                coverage[option["class_index"]].add(record["parent"])

        def greedy(largest: bool) -> list[int]:
            uncovered = set(self.cover["parents"]["ordered_tagged"])
            answer = []
            while uncovered:
                candidates = [index for index, cover in coverage.items() if cover & uncovered]
                primary = max(len(coverage[index] & uncovered) for index in candidates)
                tied = [
                    index
                    for index in candidates
                    if len(coverage[index] & uncovered) == primary
                ]
                chosen = (max if largest else min)(tied, key=lambda index: class_ids[index])
                answer.append(chosen)
                uncovered.difference_update(coverage[chosen])
            return answer

        largest = greedy(True)
        smallest = greedy(False)
        self.assertEqual(
            largest, self.cover["selected_class_indices_greedy_order"]
        )
        self.assertEqual(len(largest), 179)
        self.assertEqual(len(smallest), 180)

    def test_known_standard_positive_is_selected_and_protected(self) -> None:
        positives = self.cover["known_positive_selected"]
        self.assertEqual(len(positives), 1)
        self.assertEqual(
            positives[0]["class_index"],
            wrapper.EXPECTED_KNOWN_POSITIVE_CLASS_INDEX,
        )
        self.assertEqual(
            positives[0]["class_id"], wrapper.EXPECTED_KNOWN_POSITIVE_CLASS_ID
        )
        self.assertEqual(
            positives[0]["active_parents"],
            ["K6_only:364827", "K6_only:3335955"],
        )

    def test_selection_and_sharding_are_deterministic(self) -> None:
        shards = [
            wrapper.select_run_graphs(
                self.graphs,
                shards=3,
                shard=shard,
                sample=None,
                sample_seed=17,
                indices=None,
            )
            for shard in range(3)
        ]
        shard_sets = [{graph["index"] for graph in shard} for shard in shards]
        self.assertEqual(set.union(*shard_sets), {graph["index"] for graph in self.graphs})
        self.assertTrue(
            all(
                left.isdisjoint(right)
                for position, left in enumerate(shard_sets)
                for right in shard_sets[position + 1 :]
            )
        )
        first = wrapper.select_run_graphs(
            self.graphs,
            shards=1,
            shard=0,
            sample=11,
            sample_seed=123,
            indices=None,
        )
        second = wrapper.select_run_graphs(
            self.graphs,
            shards=1,
            shard=0,
            sample=11,
            sample_seed=123,
            indices=None,
        )
        self.assertEqual(first, second)

    def test_placement_profile_exposes_no_order_exception_without_kernel(self) -> None:
        by_index = {graph["index"]: graph for graph in self.graphs}
        sample = [
            by_index[wrapper.EXPECTED_NO_ORDER_CLASS_INDEX],
            by_index[wrapper.EXPECTED_KNOWN_POSITIVE_CLASS_INDEX],
        ]
        with patch.object(
            wrapper.engine.cdriver6,
            "decide6",
            side_effect=AssertionError("placement profiling called the kernel"),
        ):
            profile = wrapper.profile_placement_orders(sample, max_orders=4)
        self.assertEqual(
            profile["no_order_class_indices"],
            [wrapper.EXPECTED_NO_ORDER_CLASS_INDEX],
        )
        self.assertTrue(profile["all_generated_orders_use_one_or_two_circles"])
        self.assertTrue(
            set(map(int, profile["generated_circle_count_histogram"])) <= {1, 2}
        )

    def test_engine_order_scope_restores_global_and_controls_are_mockable(self) -> None:
        previous = wrapper.engine.N
        with wrapper.engine_order_18():
            self.assertEqual(wrapper.engine.N, 18)
        self.assertEqual(wrapper.engine.N, previous)

        fake_result = wrapper.engine.Result(
            0,
            -1,
            "POSITIVE_CONTROL",
            "ABORT",
            -1,
            -1,
            1,
            1,
            1,
            0,
            0,
            1,
            1000,
            1,
            0.01,
        )
        fake_base = {
            "status": "PASS",
            "positive": {"description": "old"},
            "negative": {},
        }
        with (
            patch.object(wrapper.engine, "run_kernel_controls", return_value=fake_base),
            patch.object(wrapper.engine, "analyze_graph", return_value=fake_result) as analyze,
        ):
            controls = wrapper.run_n18_kernel_controls(self.deletions, slices=3)
        self.assertEqual(controls["status"], "PASS")
        self.assertEqual(analyze.call_count, 4)
        self.assertTrue(
            all(
                "elapsed_seconds" not in control["result"]
                for control in controls["known_realizable_exact18"]
            )
        )
        for call in analyze.call_args_list:
            self.assertFalse(call.kwargs["zero_circle_only"])
            self.assertTrue(call.kwargs["include_bulk_order"])

    def test_pinned_independently_verified_v7_parent_boundary(self) -> None:
        parents, provenance = wrapper.load_v7_parents(
            wrapper.ROOT / "d6_current_residue_manifest_v7.json",
            wrapper.ROOT / "d6_current_residue_manifest_v7_verification.json",
        )
        self.assertEqual(
            [(parent["class"], parent["index"]) for parent in parents],
            [(parent["class"], parent["index"]) for parent in self.parents],
        )
        self.assertEqual(provenance["combined_count"], 263)
        self.assertEqual(provenance["verification"]["status"], "PASS")

    def test_deletion_semantics_keep_nonedges_optional(self) -> None:
        self.assertEqual(
            self.deletions["semantics"]["nonedges"],
            "unconstrained and may also have distance one",
        )
        self.assertIn("unconstrained", wrapper.TRUST_ASSUMPTIONS["candidate_nonedges"])
        self.assertEqual(self.deletion_provenance["verification"]["status"], "PASS")

    def test_source_boundary_ignores_unrelated_head_and_porcelain(self) -> None:
        source_names = (
            "run_d6_interval_18_cover_v7.py",
            "run_d6_interval_residue.py",
            "cdriver6.py",
            "ckernel6.c",
            "ival.py",
        )
        blobs = {name: (wrapper.ROOT / name).read_bytes() for name in source_names}
        source_commit = "a" * 40
        calls = []

        def fake_run(arguments, **kwargs):
            command = tuple(arguments[1:])
            calls.append(command)
            if command == ("branch", "--show-current"):
                output = "codex/dimension6\n"
            elif command[:4] == ("rev-list", "-1", "HEAD", "--"):
                self.assertEqual(command[4:], source_names)
                output = source_commit + "\n"
            elif command[:1] == ("show",):
                revision, name = command[1].split(":", 1)
                self.assertIn(revision, {"HEAD", source_commit})
                output = blobs[name]
            else:
                raise AssertionError(f"unexpected git call {command}")
            return subprocess.CompletedProcess(arguments, 0, stdout=output, stderr=b"")

        with patch.object(wrapper.subprocess, "run", side_effect=fake_run):
            first = wrapper.committed_source_boundary()
            second = wrapper.committed_source_boundary()
        self.assertEqual(first, second)
        self.assertEqual(first["source_commit"], source_commit)
        self.assertNotIn("commit", first)
        self.assertNotIn("porcelain_sha256", first)
        self.assertFalse(any(call[:1] in {("status",), ("rev-parse",)} for call in calls))

    def test_source_boundary_fails_closed_on_executable_dirt(self) -> None:
        source_names = (
            "run_d6_interval_18_cover_v7.py",
            "run_d6_interval_residue.py",
            "cdriver6.py",
            "ckernel6.c",
            "ival.py",
        )
        blobs = {name: (wrapper.ROOT / name).read_bytes() for name in source_names}
        source_commit = "b" * 40

        def fake_run(arguments, **kwargs):
            command = tuple(arguments[1:])
            if command == ("branch", "--show-current"):
                output = "codex/dimension6\n"
            elif command[:4] == ("rev-list", "-1", "HEAD", "--"):
                output = source_commit + "\n"
            elif command[:1] == ("show",):
                _revision, name = command[1].split(":", 1)
                output = b"dirty mismatch" if name == source_names[0] else blobs[name]
            else:
                raise AssertionError(f"unexpected git call {command}")
            return subprocess.CompletedProcess(arguments, 0, stdout=output, stderr=b"")

        with patch.object(wrapper.subprocess, "run", side_effect=fake_run):
            with self.assertRaisesRegex(ValueError, "production source"):
                wrapper.committed_source_boundary()


if __name__ == "__main__":
    unittest.main(verbosity=2)
