#!/usr/bin/env python3
"""Focused controls for the source-bound graph-3936177 replay wrapper."""

from __future__ import annotations

import copy
import unittest

import run_d6_interval_k7_3936177 as wrapper


TARGET_ADJACENCY_SHA256 = (
    "264bc830ee02673f329bfdf52659d75672b084fe0e3128e597103855a868fb35"
)


class FocusedReplayWrapperControls(unittest.TestCase):
    def load(self) -> tuple[list[dict], dict]:
        return wrapper.v6.load_v6_graphs(
            wrapper.ROOT / "d6_current_residue_manifest_v6.json",
            wrapper.ROOT / "d6_current_residue_manifest_v6_verification.json",
        )

    def test_selection_is_exactly_the_frozen_target(self) -> None:
        graphs, _provenance = self.load()
        selected = wrapper.select_target(graphs)
        self.assertEqual(len(selected), 1)
        self.assertEqual(selected[0]["index"], wrapper.TARGET_INDEX)
        self.assertEqual(selected[0]["population"], "K7")
        self.assertEqual(selected[0]["ordinal"], 0)
        self.assertEqual(
            wrapper.engine.stable_hash([selected[0]["index"]]),
            wrapper.TARGET_INDICES_SHA256,
        )
        self.assertEqual(
            wrapper.engine.stable_hash(list(selected[0]["adjacency"])),
            TARGET_ADJACENCY_SHA256,
        )

    def test_missing_duplicate_and_wrong_class_targets_are_rejected(self) -> None:
        graphs, _provenance = self.load()
        target = next(graph for graph in graphs if graph["index"] == wrapper.TARGET_INDEX)
        without = [graph for graph in graphs if graph["index"] != wrapper.TARGET_INDEX]
        with self.assertRaisesRegex(ValueError, "exactly once"):
            wrapper.select_target(without)
        with self.assertRaisesRegex(ValueError, "exactly once"):
            wrapper.select_target([*graphs, copy.deepcopy(target)])
        wrong = [copy.deepcopy(target)]
        wrong[0]["population"] = "K6"
        with self.assertRaisesRegex(ValueError, "not in the K7"):
            wrapper.select_target(wrong)

    def test_search_parameters_and_worker_count_are_not_cli_variables(self) -> None:
        self.assertEqual(
            (wrapper.ORDERS, wrapper.CAP, wrapper.SLICES, wrapper.WORKERS),
            (4, 2_000_000, 24, 1),
        )

    def test_generated_order_list_contains_the_discovery_order_without_using_its_result(self) -> None:
        graphs, _provenance = self.load()
        target = wrapper.select_target(graphs)[0]
        orders = wrapper.engine.cdriver6.gen_orders(
            target["adjacency"], wrapper.engine.N, kmax=wrapper.ORDERS
        )
        self.assertEqual(len(orders), 2)
        seed, placement = orders[1]
        self.assertEqual(list(seed), [1, 6, 8, 11, 13, 16, 18])
        self.assertEqual(
            list(placement), [2, 15, 12, 7, 5, 10, 9, 4, 3, 0, 17, 14]
        )

    def test_all_frozen_dependency_hashes_match(self) -> None:
        for _key, (filename, expected) in wrapper.EXPECTED_TRACKED_FILES.items():
            self.assertEqual(wrapper.engine.sha256(wrapper.ROOT / filename), expected)
        for _key, (filename, expected) in wrapper.EXPECTED_INPUT_FILES.items():
            self.assertEqual(wrapper.engine.sha256(wrapper.ROOT / filename), expected)

    def test_modified_tracked_launch_is_rejected_before_blob_lookup(self) -> None:
        fake_git = {
            "available": True,
            "branch": "codex/dimension6",
            "commit": "0" * 40,
            "dirty": True,
            "porcelain_lines": [" M ckernel6.c"],
            "porcelain_sha256": "irrelevant",
        }
        with self.assertRaisesRegex(ValueError, "modified tracked file"):
            wrapper.validate_committed_launch(fake_git, {})

    def test_wrong_branch_is_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "wrong launch branch"):
            wrapper.validate_committed_launch(
                {"available": True, "branch": "main"}, {}
            )


if __name__ == "__main__":
    unittest.main()
