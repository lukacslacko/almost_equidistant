#!/usr/bin/env python3
"""Synthetic controls for the exact near-clique active/sign layer."""

from __future__ import annotations

import copy
import unittest
from itertools import product

import d6_k7_rankone_active_sign as locator
import verify_d6_k7_rankone_active_sign as verifier


def targets(
    vertices: list[int], entries: dict[tuple[int, int], int]
) -> dict[tuple[int, int], int]:
    expected = {
        (left, right)
        for position, left in enumerate(vertices)
        for right in vertices[position + 1:]
    }
    if set(entries) != expected:
        raise AssertionError("test target table is incomplete")
    return entries


class PairClassificationControls(unittest.TestCase):
    def kind(self, left: int, right: int, target: int, width: int = 3) -> str:
        return locator.classify_pair(left, right, target, width).classification

    def test_nonedge_disjoint_nonempty_is_positive(self) -> None:
        self.assertEqual(self.kind(0b001, 0b010, 0), locator.FORCED_POSITIVE)

    def test_nonedge_with_empty_mask_is_identically_zero(self) -> None:
        self.assertEqual(self.kind(0, 0b001, 0), locator.IDENTICALLY_ZERO)

    def test_nonedge_comparable_overlap_is_negative(self) -> None:
        self.assertEqual(self.kind(0b001, 0b011, 0), locator.FORCED_NEGATIVE)

    def test_nonedge_three_nonempty_regions_is_flexible(self) -> None:
        self.assertEqual(self.kind(0b011, 0b110, 0), locator.FLEXIBLE)

    def test_edge_disjoint_is_positive(self) -> None:
        self.assertEqual(self.kind(0b001, 0b010, 1), locator.FORCED_POSITIVE)

    def test_edge_union_full_is_positive(self) -> None:
        self.assertEqual(self.kind(0b011, 0b110, 1), locator.FORCED_POSITIVE)

    def test_edge_with_overlap_and_outside_is_flexible(self) -> None:
        self.assertEqual(
            self.kind(0b0011, 0b0110, 1, width=4),
            locator.FLEXIBLE,
        )


class ActiveSignControls(unittest.TestCase):
    def test_empty_and_singleton_active_sets_are_safe(self) -> None:
        vertices = [0, 1]
        pair_targets = targets(vertices, {(0, 1): 0})
        masks = [0b0011, 0b0110]
        pairs = locator.classify_system(4, vertices, masks, pair_targets)
        result = locator.assess_system(4, vertices, masks, pair_targets)
        self.assertFalse(result.rejected)
        self.assertEqual(result.relaxed_active, ())
        self.assertTrue(locator.branch_is_pairwise_compatible(pairs, []))
        self.assertTrue(locator.branch_is_pairwise_compatible(pairs, [0]))
        self.assertTrue(locator.branch_is_pairwise_compatible(pairs, [1]))
        self.assertTrue(locator.branch_is_pairwise_compatible(
            pairs, [0, 1], {0: 1, 1: -1}
        ))

    def test_empty_and_full_masks_are_diagonally_forced_active(self) -> None:
        for mask in (0, 0b111):
            with self.subTest(mask=mask):
                result = locator.assess_system(3, [7], [mask], {})
                self.assertFalse(result.rejected)
                self.assertEqual(result.diagonal_forced_active, (7,))
                self.assertEqual(result.relaxed_active, (7,))

    def test_extreme_mask_diagonal_zero_contradiction(self) -> None:
        vertices = [0, 1]
        pair_targets = targets(vertices, {(0, 1): 0})
        masks = [0, 0b111]
        result = locator.assess_system(3, vertices, masks, pair_targets)
        self.assertTrue(result.rejected)
        self.assertEqual(result.reason, "ZERO_INSIDE_FORCED_ACTIVE")
        self.assertEqual(result.diagonal_forced_active, (0, 1))
        certificate = locator.make_certificate(
            3, vertices, masks, pair_targets, result
        )
        self.assertEqual(
            verifier.verify_certificate(certificate)["status"], "VERIFIED"
        )

    def test_flexible_pair_never_forces_activity(self) -> None:
        vertices = [4, 9]
        pair_targets = targets(vertices, {(4, 9): 0})
        pairs = locator.classify_system(
            4, vertices, [0b0011, 0b0110], pair_targets
        )
        result = locator.assess_system(
            4, vertices, [0b0011, 0b0110], pair_targets
        )
        self.assertEqual(pairs[0].classification, locator.FLEXIBLE)
        self.assertEqual(result.forced_active, ())
        self.assertFalse(result.rejected)

    def test_no_saturating_alpha_two_pair_condition_is_no_go(self) -> None:
        # Exhaust the local hypotheses of the proof-note proposition for
        # three clique coordinates and three remainder vertices.  No mask is
        # full, and a zero target is allowed only when the two masks cover the
        # clique (the consequence of alpha(G)<=2).
        vertices = [0, 1, 2]
        full = 0b111
        for masks in product(range(full), repeat=len(vertices)):
            pairs = [(0, 1), (0, 2), (1, 2)]
            domains = [
                (0, 1) if masks[left] | masks[right] == full else (1,)
                for left, right in pairs
            ]
            for values in product(*domains):
                pair_targets = dict(zip(pairs, values, strict=True))
                result = locator.assess_system(
                    3, vertices, masks, pair_targets
                )
                self.assertFalse(
                    result.rejected,
                    (masks, pair_targets, result),
                )

    def test_identically_zero_forced_active_pair_rejects(self) -> None:
        vertices = [0, 1, 2]
        pair_targets = targets(vertices, {
            (0, 1): 0,
            (0, 2): 1,
            (1, 2): 1,
        })
        result = locator.assess_system(2, vertices, [0, 1, 2], pair_targets)
        self.assertTrue(result.rejected)
        self.assertEqual(result.reason, "ZERO_INSIDE_FORCED_ACTIVE")
        certificate = locator.make_certificate(
            2, vertices, [0, 1, 2], pair_targets, result
        )
        checked = verifier.verify_certificate(certificate)
        self.assertEqual(checked["status"], "VERIFIED")

    def test_unbalanced_forced_sign_triangle_rejects(self) -> None:
        vertices = [0, 1, 2]
        pair_targets = targets(vertices, {
            (0, 1): 0,
            (0, 2): 1,
            (1, 2): 1,
        })
        masks = [0b01, 0b11, 0b10]
        result = locator.assess_system(2, vertices, masks, pair_targets)
        self.assertTrue(result.rejected)
        self.assertEqual(result.reason, "FORCED_SIGN_PARITY")
        certificate = locator.make_certificate(
            2, vertices, masks, pair_targets, result
        )
        checked = verifier.verify_certificate(certificate)
        self.assertEqual(checked["reason"], "FORCED_SIGN_PARITY")

    def test_balanced_forced_sign_triangle_survives_relaxation(self) -> None:
        vertices = [0, 1, 2]
        pair_targets = targets(vertices, {
            (0, 1): 0,
            (0, 2): 1,
            (1, 2): 0,
        })
        masks = [0b01, 0b11, 0b10]
        result = locator.assess_system(2, vertices, masks, pair_targets)
        self.assertFalse(result.rejected)
        self.assertEqual(result.relaxed_active, (0, 1, 2))
        pairs = locator.classify_system(2, vertices, masks, pair_targets)
        self.assertTrue(locator.branch_is_pairwise_compatible(
            pairs,
            result.relaxed_active,
            dict(result.relaxed_signs or ()),
        ))

    def test_tampered_certificate_is_rejected_independently(self) -> None:
        vertices = [0, 1, 2]
        pair_targets = targets(vertices, {
            (0, 1): 0,
            (0, 2): 1,
            (1, 2): 1,
        })
        result = locator.assess_system(2, vertices, [0, 1, 2], pair_targets)
        certificate = locator.make_certificate(
            2, vertices, [0, 1, 2], pair_targets, result
        )
        tampered = copy.deepcopy(certificate)
        tampered["assessment"]["witness"]["zero_pair"] = [1, 2]
        with self.assertRaises(ValueError):
            verifier.verify_certificate(tampered)
        tampered_diagonal = copy.deepcopy(certificate)
        tampered_diagonal["assessment"]["diagonal_forced_active"] = []
        with self.assertRaises(ValueError):
            verifier.verify_certificate(tampered_diagonal)


if __name__ == "__main__":
    unittest.main()
