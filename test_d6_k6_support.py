#!/usr/bin/env python3
"""Controls for the K6 joint actual-support intersection prototype."""

from __future__ import annotations

import json
import unittest
from pathlib import Path

from d6_k6_lorentz_reference import K6LorentzInstance, support_matching
from d6_k6_support_reference import (
    actual_supports_for_bins,
    evaluate_support_graph,
    intersection_compatible,
)
from test_d6_k6_lorentz import lower_bound_18_graph


ROOT = Path(__file__).resolve().parent
RESIDUE = ROOT / "d6_k6_support_residue.json"


def artificial_instance(defects: list[int], eligible: int = 0) -> K6LorentzInstance:
    outside = tuple(range(6, 6 + len(defects)))
    return K6LorentzInstance(
        tuple(range(6)),
        outside,
        tuple(defects),
        tuple(0 for _ in outside),
        eligible,
    )


class K6SupportControls(unittest.TestCase):
    def test_single_shared_coordinate_cannot_cancel(self) -> None:
        self.assertFalse(intersection_compatible(0b000111, 0b000100))
        self.assertTrue(intersection_compatible(0b000011, 0b001100))
        self.assertTrue(intersection_compatible(0b000011, 0b000011))

    def test_z0_support_is_shared_across_both_bins(self) -> None:
        # z has D={0,1,2,3}; the left and right light vectors have singleton
        # supports {0} and {1}.  Each bin separately can choose a three-slot z
        # support avoiding its singleton.  No *shared* size-three support can
        # avoid both, so the correct joint search rejects.
        instance = artificial_instance([0b001111, 0b000001, 0b000010], eligible=1)
        z0 = 0b001
        left_only = actual_supports_for_bins(instance, z0, (0b011, 0b001))
        right_only = actual_supports_for_bins(instance, z0, (0b001, 0b101))
        joint = actual_supports_for_bins(instance, z0, (0b011, 0b101))
        self.assertTrue(left_only.feasible)
        self.assertTrue(right_only.feasible)
        self.assertFalse(joint.feasible)

    def test_allowed_mask_matching_can_be_strictly_too_weak(self) -> None:
        # The Z0 support is forced to {0,1,2}; the light support is forced to
        # {0}.  Allowed masks match to distinct coordinates, but actual support
        # intersection is exactly one for every choice.
        instance = artificial_instance([0b000111, 0b000001], eligible=1)
        z0 = 0b01
        bins = (0b11, 0b01)
        self.assertIsNotNone(support_matching(bins[0], instance.defects))
        self.assertIsNotNone(support_matching(bins[1], instance.defects))
        self.assertFalse(actual_supports_for_bins(instance, z0, bins).feasible)

    def test_intersection_two_is_conservatively_accepted(self) -> None:
        # This is a support screen, not a coefficient solver.  Two shared
        # coordinates may permit cancellation and therefore must pass here.
        instance = artificial_instance([0b000111, 0b000011], eligible=1)
        result = actual_supports_for_bins(instance, 0b01, (0b11, 0b01))
        self.assertTrue(result.feasible)
        assert result.supports is not None
        self.assertEqual(result.supports[0], 0b000111)
        self.assertEqual(result.supports[1], 0b000011)

    def test_empty_light_support_domain_fails(self) -> None:
        instance = artificial_instance([0])
        self.assertFalse(actual_supports_for_bins(instance, 0, (1, 0)).feasible)

    def test_realizable_18_points_pass_over_every_k6_seed(self) -> None:
        decision = evaluate_support_graph(lower_bound_18_graph(), scan_all=True)
        self.assertTrue(decision.applicable)
        self.assertFalse(decision.rejected)
        self.assertEqual(decision.seeds_checked, 32)
        self.assertEqual(decision.impossible_seeds, 0)

    def test_fixed_real_residue_witnesses(self) -> None:
        if not RESIDUE.exists():
            self.skipTest("run extract_d6_k6_support_residue.py first")
        expected = {
            132_876: 1,
            550_637: 2,
            1_323_686: 1,
            2_640_615: 1,
            2_642_489: 1,
            3_107_353: 1,
            3_819_932: 2,
            3_952_058: 1,
        }
        with RESIDUE.open(encoding="utf-8") as stream:
            residue = json.load(stream)
        selected = {
            graph["index"]: graph["adjacency"]
            for graph in residue["graphs"]
            if graph["index"] in expected
        }
        self.assertEqual(set(selected), set(expected))
        for index, impossible_seeds in expected.items():
            with self.subTest(index=index):
                decision = evaluate_support_graph(selected[index], scan_all=True)
                self.assertTrue(decision.applicable)
                self.assertTrue(decision.rejected)
                self.assertEqual(decision.impossible_seeds, impossible_seeds)


if __name__ == "__main__":
    unittest.main(verbosity=2)
