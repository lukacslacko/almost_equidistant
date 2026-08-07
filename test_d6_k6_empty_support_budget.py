#!/usr/bin/env python3
"""Exact controls for the K6 empty-defect support budget."""

from __future__ import annotations

import ast
import json
import tempfile
import unittest
from pathlib import Path

import d6_k6_empty_support_budget as target
import d6_k6_psd_z_hereditary as hall
import d6_k6_psd_zmatrix as ranks
from d6_k6_bipartite_rank_reference import ZeroForcingSolver
from d6_k6_normal_inertia import InertiaCache
from test_d6_k6_lorentz import lower_bound_18_graph


class EmptySupportBudgetControls(unittest.TestCase):
    def test_positive_singleton_is_nonempty_negative_singleton_may_be_empty(self) -> None:
        graph = (0,)
        absolute = (7,)
        inertia = InertiaCache()
        forcing = ZeroForcingSolver()
        positive = ranks.side_rank_lower(
            graph, absolute, "positive", inertia, forcing
        )
        negative = ranks.side_rank_lower(
            graph, absolute, "negative", inertia, forcing
        )
        self.assertEqual(positive.rank_lower, 1)
        self.assertEqual(negative.rank_lower, 0)

    @staticmethod
    def singleton_system(count: int) -> tuple[
        tuple[hall.Group, ...], tuple[target.EmptySingleton, ...]
    ]:
        groups = tuple(
            hall.Group(
                f"negative:F{i}",
                (
                    hall.Choice(0, 0, (), "omit"),
                    hall.Choice(0, 1, (i,), "old_zero_rank"),
                ),
                True,
            )
            for i in range(count)
        )
        zero = tuple(target.EmptySingleton(i, i, 1) for i in range(count))
        return groups, zero

    def test_one_empty_on_one_negative_line_is_allowed(self) -> None:
        groups, zero = self.singleton_system(2)
        strengthened = target.force_singleton_status(groups, zero, 0)
        self.assertTrue(hall.check_groups(strengthened).passed)
        # With neither point empty, two rank-one vectors confined to the same
        # coordinate fail Hall exactly.
        self.assertFalse(
            hall.check_groups(
                target.force_singleton_status(groups, zero, None)
            ).passed
        )

    def test_two_empty_points_on_one_negative_line_are_forbidden(self) -> None:
        groups, zero = self.singleton_system(3)
        # Permitting at most one empty leaves two nonzero orthogonal spans in
        # one coordinate, hence every allowed choice fails.
        for chosen in (None, 0, 1, 2):
            self.assertFalse(
                hall.check_groups(
                    target.force_singleton_status(groups, zero, chosen)
                ).passed
            )

    def test_global_two_point_budget(self) -> None:
        self.assertTrue(target.global_component_costs_pass([1, 1]))
        self.assertTrue(target.global_component_costs_pass([0, 1, 0, 1]))
        self.assertFalse(target.global_component_costs_pass([1, 1, 1]))
        self.assertFalse(target.global_component_costs_pass([0, None]))

    def test_known_realizable_eighteen_point_control(self) -> None:
        decision = target.evaluate_graph(lower_bound_18_graph())
        self.assertFalse(decision["rejected"])
        self.assertEqual(decision["seeds_checked"], 32)

    def test_fixed_new_pilot_witness(self) -> None:
        _, records, _, _ = target.verify_parent_input()
        record = next(row for row in records if row["index"] == 202_556)
        decision = target.evaluate_graph(record["adjacency"])
        self.assertTrue(decision["rejected"])
        self.assertEqual(
            decision["first_impossible_seed"], [3, 8, 11, 12, 13, 18]
        )
        encoded = json.dumps(decision["certificate"], sort_keys=True)
        self.assertNotIn('"failures"', encoded)

    def test_checkpoint_is_source_input_dependency_and_config_bound(self) -> None:
        indices = [17, 23]
        dependencies = {"parent": "abc"}
        completed = [{"index": 17, "decision": {"rejected": False}}]
        payload = target.checkpoint_payload(
            "source-a", indices, dependencies, completed
        )
        self.assertEqual(
            target.validate_checkpoint_payload(
                payload, "source-a", indices, dependencies
            ),
            completed,
        )
        with self.assertRaises(ValueError):
            target.validate_checkpoint_payload(
                payload, "source-b", indices, dependencies
            )
        with self.assertRaises(ValueError):
            target.validate_checkpoint_payload(
                payload, "source-a", [17, 29], dependencies
            )
        with self.assertRaises(ValueError):
            target.validate_checkpoint_payload(
                payload, "source-a", indices, {"parent": "changed"}
            )
        reordered = dict(payload)
        reordered["completed"] = [
            {"index": 23, "decision": {"rejected": False}}
        ]
        with self.assertRaises(ValueError):
            target.validate_checkpoint_payload(
                reordered, "source-a", indices, dependencies
            )
        with tempfile.TemporaryDirectory() as directory:
            checkpoint = Path(directory) / "checkpoint.json"
            target.atomic_json(checkpoint, payload)
            self.assertEqual(
                target.load_checkpoint(
                    checkpoint, "source-a", indices, dependencies
                ),
                completed,
            )

    def test_independent_checker_does_not_import_production_kernel(self) -> None:
        checker = Path(__file__).with_name(
            "verify_d6_k6_empty_support_budget.py"
        )
        tree = ast.parse(checker.read_text(encoding="utf-8"))
        imports = {
            node.module
            for node in ast.walk(tree)
            if isinstance(node, ast.ImportFrom) and node.module
        }
        imports.update(
            alias.name
            for node in ast.walk(tree)
            if isinstance(node, ast.Import)
            for alias in node.names
        )
        self.assertNotIn("d6_k6_empty_support_budget", imports)

    def test_infrastructure_abort_is_never_a_rejection(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            checkpoint = Path(directory) / "checkpoint.json"
            target.atomic_json(
                checkpoint,
                target.checkpoint_payload("source", [1], {}, []),
            )
            try:
                raise RuntimeError("synthetic worker failure")
            except RuntimeError as error:
                payload = target.infrastructure_abort_payload(
                    "source", {}, [1], [], checkpoint, error
                )
        self.assertEqual(payload["status"], "INFRA_ABORT")
        self.assertNotIn("rejected", payload)


if __name__ == "__main__":
    unittest.main()
