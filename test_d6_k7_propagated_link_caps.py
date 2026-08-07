#!/usr/bin/env python3
"""Exact controls for propagated K7 clique-link caps."""

from __future__ import annotations

import unittest

import d6_k7_propagated_link_caps as caps


class PropagatedLinkCapControls(unittest.TestCase):
    def test_all_coordinate_subsets_are_present(self):
        self.assertEqual(
            sum(map(len, caps.SUBSETS_BY_SIZE.values())), 119
        )
        self.assertEqual(caps.OUTSIDE_CAPS, {1: 1, 2: 2, 3: 7, 4: 8, 5: 11})

    def test_each_cap_boundary_and_first_failure(self):
        for size, outside_cap in caps.OUTSIDE_CAPS.items():
            subset = (1 << size) - 1
            boundary = [subset] * outside_cap
            self.assertIsNone(caps.first_link_cap_failure((), boundary))
            failure = caps.first_link_cap_failure((), (*boundary, subset))
            self.assertIsNotNone(failure)
            self.assertEqual(failure.subset_mask, subset)
            self.assertEqual(failure.subset_size, size)
            self.assertEqual(failure.outside_cap, outside_cap)
            self.assertEqual(failure.forced_count, outside_cap + 1)

    def test_fixed_zero_factor_supports_are_counted(self):
        # Seven points of type contained in {0,1,2} meet the K4-link outside
        # cap; adding one propagated N mask exceeds it.
        self.assertIsNone(caps.first_link_cap_failure((7, 7, 7), (7, 7, 7, 7)))
        failure = caps.first_link_cap_failure((7, 7, 7), (7, 7, 7, 7, 7))
        self.assertIsNotNone(failure)
        self.assertEqual(failure.subset_mask, 7)
        self.assertEqual(failure.zero_factor_count, 3)
        self.assertEqual(failure.nonzero_factor_count, 5)

    def test_uncontained_masks_are_conservatively_ignored(self):
        self.assertIsNone(caps.first_link_cap_failure((), (127,) * 12))

    def test_invalid_masks_fail_closed(self):
        for bad in (0, 128, -1):
            with self.assertRaises(ValueError):
                caps.first_link_cap_failure((bad,), ())
            with self.assertRaises(ValueError):
                caps.first_link_cap_failure((), (bad,))


if __name__ == "__main__":
    unittest.main()

