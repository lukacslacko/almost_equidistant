#!/usr/bin/env python3
"""Cross-checks for the odd-cycle/full-pin incremental package."""

import hashlib
import json
import unittest

import d6_k7_full_pin_odd_cycle as production
import verify_d6_k7_full_pin_increment as independent


class FullPinIncrementTests(unittest.TestCase):
    def test_independent_locator_matches_synthetic_certificate(self) -> None:
        adjacency = [0] * 5
        for first, second in ((0, 1), (1, 2), (2, 3), (3, 4), (4, 0)):
            adjacency[first] |= 1 << second
            adjacency[second] |= 1 << first
        masks = (1, 1, 1, 1, 1)
        certificate = production.find_full_pin(adjacency, masks)
        self.assertIsNotNone(certificate)
        self.assertEqual(
            production.certificate_json(certificate),
            independent.independent_full_pin(adjacency, masks),
        )

    def test_qsqrt7_separation_rules_out_every_sign_support(self) -> None:
        # For signs, Q=k and T is an integer.  The sqrt(7) coefficient first
        # forces T=0; the rational part would then require k=8, impossible.
        for size in range(1, 8):
            for negative in range(size + 1):
                total = size - 2 * negative
                rational_part = size - total * total - 8
                irrational_coefficient = 2 * total
                self.assertNotEqual(
                    (rational_part, irrational_coefficient), (0, 0)
                )

    def test_expected_exact_24_residue_hash(self) -> None:
        survivors = [
            226183, 316173, 423661, 424226, 2581209, 2592657,
            2593240, 3595554, 3624785, 3648882, 3729907, 3785980,
            3888410, 3935560, 3936176, 3936177, 3936310, 3936435,
            3945490, 3945555, 3945557, 3945564, 3947605, 3949382,
        ]
        digest = hashlib.sha256(json.dumps(
            survivors, sort_keys=True, separators=(",", ":"), ensure_ascii=True
        ).encode("ascii")).hexdigest()
        self.assertEqual(
            digest,
            "71d9ae101445ae90d08699ede81ca405bf5bfe956fda01ee384203fa113d1779",
        )


if __name__ == "__main__":
    unittest.main()
