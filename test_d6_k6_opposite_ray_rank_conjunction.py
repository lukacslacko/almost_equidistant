#!/usr/bin/env python3
"""Focused controls for opposite-light-ray rank and support conjunction."""

from __future__ import annotations

import hashlib
import unittest
from types import SimpleNamespace

import d6_k6_opposite_ray_rank_conjunction as production
import verify_d6_k6_opposite_ray_rank_conjunction as verifier


def graph_with_edges(size: int, edges: list[tuple[int, int]]) -> tuple[int, ...]:
    rows = [0] * size
    for left, right in edges:
        rows[left] |= 1 << right
        rows[right] |= 1 << left
    return tuple(rows)


class OppositeRayRankConjunctionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        records, indices, _ = production.load_input()
        cls.records = {int(row["index"]): row for row in records}
        cls.indices = indices
        independent, independent_indices = verifier.load_input()
        cls.independent = {int(row["index"]): row for row in independent}
        cls.independent_indices = independent_indices

    def test_certified_623_parent_boundary(self) -> None:
        self.assertEqual(len(self.records), 623)
        self.assertEqual(self.indices, self.independent_indices)
        self.assertEqual(
            production.stable_hash(self.indices),
            production.EXPECTED_INPUT_SHA256,
        )
        self.assertNotIn(1_442_098, self.indices)
        self.assertIn(87_757, self.indices)
        self.assertEqual(production.EXPECTED_REJECTIONS, 372)
        self.assertEqual(verifier.EXPECTED_REJECTIONS, 372)

    def test_PSD_Z_component_rank_lower(self) -> None:
        triangle = (0b110, 0b101, 0b011)
        self.assertEqual(
            production.psd_z_component_rank_lower(triangle),
            (2, [[0, 1, 2]]),
        )
        two_edges = (0b0010, 0b0001, 0b1000, 0b0100)
        self.assertEqual(
            production.psd_z_component_rank_lower(two_edges),
            (2, [[0, 1], [2, 3]]),
        )
        self.assertEqual(
            production.psd_z_component_rank_lower((0, 0, 0)),
            (0, []),
        )

    def test_saturated_and_zero_cross_projection_failures(self) -> None:
        instance = SimpleNamespace(outside=tuple(range(8)))
        source = 0b00000011
        basis6 = 0b11111100

        # Two source rows sharing one target give a nontrivial residual
        # support component.  A six-vector target basis leaves rank capacity
        # zero, so this is impossible.
        adjacency = graph_with_edges(8, [(0, 2), (1, 2)])
        check = production.directed_projection_check(
            adjacency, instance, 0, source, basis6
        )
        self.assertFalse(check["passed"])
        self.assertEqual(check["rank_lower"], 1)
        self.assertEqual(check["rank_capacity"], 0)
        self.assertEqual(check["common_neighbour_components"], [[0, 1]])

        # A row with no target-light required neighbour has zero cross row,
        # hence residual diagonal exactly one and another forced rank unit.
        adjacency = graph_with_edges(8, [(1, 2)])
        check = production.directed_projection_check(
            adjacency, instance, 0, source, basis6
        )
        self.assertFalse(check["passed"])
        self.assertEqual(check["zero_cross_degree_source"], [0])
        self.assertEqual(check["zero_cross_degree_rank"], 1)

    def test_both_saturated_integer_type_equation(self) -> None:
        self.assertTrue(
            production.both_saturated_with_light(
                0, (0b00111111, 0b11111100)
            )
        )
        self.assertFalse(
            production.both_saturated_with_light(
                0b00111111, (0b00111111, 0b00111111)
            )
        )
        # Irrationality forces p=0 in
        # 6=k+p(4-2sqrt(3))+4n, p+n=6-k.  The sole integer solution is the
        # no-light case k=6, p=n=0.
        solutions = []
        for k in range(7):
            p = 0
            n = 6 - k
            if 6 == k + 4 * n:
                solutions.append((k, p, n))
        self.assertEqual(solutions, [(6, 0, 0)])

    def test_tight_actual_support_rule_and_optional_zeros(self) -> None:
        assignments = {0: 0b0011, 1: 0b0011, 2: 0b0111}
        # Rows 0,1 have a tight two-coordinate union.  Row 2 intersects that
        # union in two coordinates, which the older pair rule allowed, but it
        # cannot be orthogonal to their full span.
        self.assertFalse(
            production.tight_orthogonal_supports_pass(0b111, assignments)
        )
        assignments[2] = 0b0100
        self.assertTrue(
            production.tight_orthogonal_supports_pass(0b111, assignments)
        )
        self.assertEqual(
            set(production.support_domains(0b111, False)),
            {1, 2, 3, 4, 5, 6, 7},
        )

    def test_opposite_required_edge_local_support_rules(self) -> None:
        instance = SimpleNamespace(outside=(0, 1))
        adjacency = graph_with_edges(2, [(0, 1)])
        bins = (1, 2)
        checks = (
            production.opposite_required_edge_supports_pass,
            verifier.independent_cross_edge_local_pass,
        )
        for check in checks:
            # A required opposite-ray edge has nonzero dot product, so the
            # actual supports cannot be disjoint.
            self.assertFalse(
                check(adjacency, instance, 0, bins, {0: 1, 1: 2})
            )
            # Equal singleton supports force both unit rows to be +e_i and
            # c=2, giving dot 1 instead of the required 4/3.
            self.assertFalse(
                check(adjacency, instance, 0, bins, {0: 1, 1: 1})
            )
            # Overlapping nonsingleton supports remain a necessary
            # possibility; no coefficient existence is claimed.
            self.assertTrue(
                check(adjacency, instance, 0, bins, {0: 0b11, 1: 0b101})
            )

    def test_tight_selected_isolated_cross_edge_line_energy(self) -> None:
        instance = SimpleNamespace(outside=tuple(range(7)))
        z0 = sum(1 << local for local in range(2, 7))
        bins = (z0 | 1, z0 | 2)
        adjacency = graph_with_edges(7, [(0, 1)])
        assignments = {
            0: 0b100000,
            1: 0b100000,
            **{local: 0b011111 for local in range(2, 7)},
        }
        # The one edge line by itself is coordinate-tight.  Parseval and the
        # exact opposite-ray scalar values exclude such a tight family.
        self.assertFalse(
            production.cross_component_intersection_supports_pass(
                adjacency, instance, z0, bins, assignments
            )
        )
        self.assertFalse(
            verifier.independent_cross_intersection_pass(
                adjacency, instance, z0, bins, assignments
            )
        )

        # A two-coordinate equal support is not tight by itself, but the five
        # displayed Z0 masks make a larger tight line family, which the same
        # all-ones/Parseval identity excludes.
        assignments[0] = assignments[1] = 0b110000
        self.assertFalse(
            production.cross_component_intersection_supports_pass(
                adjacency, instance, z0, bins, assignments
            )
        )
        self.assertFalse(
            verifier.independent_cross_intersection_pass(
                adjacency, instance, z0, bins, assignments
            )
        )

        # With no Z0 lines, a two-coordinate edge-line mask is not a tight
        # subfamily, so this theorem alone leaves the block feasible.
        no_z0_instance = SimpleNamespace(outside=(0, 1))
        no_z0_bins = (1, 2)
        no_z0_assignments = {0: 0b110000, 1: 0b110000}
        self.assertTrue(
            production.cross_component_intersection_supports_pass(
                adjacency[:2], no_z0_instance, 0,
                no_z0_bins, no_z0_assignments,
            )
        )
        self.assertTrue(
            verifier.independent_cross_intersection_pass(
                adjacency[:2], no_z0_instance, 0,
                no_z0_bins, no_z0_assignments,
            )
        )

    def test_fixed_rejection_and_independent_full_partition(self) -> None:
        produced = production.evaluate_record(self.records[87_757])
        checked = verifier.evaluate_record(self.independent[87_757])
        self.assertTrue(produced["rejected"])
        self.assertTrue(checked["rejected"])
        self.assertEqual(produced["seeds_checked"], 17)
        self.assertEqual(produced["first_impossible_seed_mask"], 278_740)
        self.assertEqual(
            produced["first_impossible_seed"], [2, 4, 6, 7, 14, 18]
        )
        self.assertEqual(len(produced["certificate"]["Z0_failures"]), 127)
        self.assertEqual(checked["seed_mask"], 278_740)
        self.assertEqual(
            checked["failures"], produced["certificate"]["Z0_failures"]
        )

        # This graph survives the 285-rejection component-intersection
        # baseline and is rejected only after the isolated-edge line-energy
        # theorem is conjoined on the same actual-support assignment.
        produced = production.evaluate_record(self.records[204_844])
        checked = verifier.evaluate_record(self.independent[204_844])
        self.assertTrue(produced["rejected"])
        self.assertEqual(produced["seeds_checked"], 9)
        self.assertEqual(produced["first_impossible_seed_mask"], 164_641)
        self.assertEqual(
            produced["first_impossible_seed"], [0, 5, 8, 9, 15, 17]
        )
        self.assertEqual(len(produced["certificate"]["Z0_failures"]), 64)
        self.assertEqual(checked["seed_mask"], 164_641)
        self.assertEqual(
            checked["failures"], produced["certificate"]["Z0_failures"]
        )

    def test_positive_controls_and_import_independence(self) -> None:
        self.assertEqual(
            production.positive_control(), {"passed": True, "K6_seeds": 32}
        )
        self.assertEqual(
            verifier.positive_control(), {"passed": True, "K6_seeds": 32}
        )
        verifier.assert_import_independence()

    def test_porcelain_binding_allows_only_untracked_dirt(self) -> None:
        lines = ["?? discovery.json", "?? scratch.py"]
        status = "\n".join(lines)
        binding = {
            "porcelain_lines": lines,
            "porcelain_sha256": hashlib.sha256(
                status.encode("utf-8")
            ).hexdigest(),
            "dirty": True,
        }
        self.assertEqual(
            verifier.validate_porcelain_binding(binding)["untracked_entries"],
            2,
        )
        tracked = dict(binding)
        tracked["porcelain_lines"] = [
            " M d6_k6_opposite_ray_rank_conjunction.py"
        ]
        tracked["porcelain_sha256"] = hashlib.sha256(
            tracked["porcelain_lines"][0].encode("utf-8")
        ).hexdigest()
        with self.assertRaisesRegex(ValueError, "tracked source was dirty"):
            verifier.validate_porcelain_binding(tracked)
        tampered = dict(binding, porcelain_sha256="0" * 64)
        with self.assertRaisesRegex(ValueError, "porcelain binding differs"):
            verifier.validate_porcelain_binding(tampered)


if __name__ == "__main__":
    unittest.main(verbosity=2)
