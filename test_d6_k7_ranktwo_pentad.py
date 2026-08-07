#!/usr/bin/env python3
"""Focused exact tests for the rank-two Schur pentad layer."""

from __future__ import annotations

import json
import unittest
from fractions import Fraction
from pathlib import Path

import d6_k7_ranktwo_pentad as production
import d6_k7_rank_reference as prior
import run_d6_k7_ranktwo_pentad_full as full_runner
import verify_d6_k7_ranktwo_pentad as verifier
from verify_profile_d6 import lower_bound_18_graph


ROOT = Path(__file__).resolve().parent


def first_target() -> dict:
    payload = json.loads(
        (ROOT / "d6_k7_ranktwo_pentad_targets.json").read_text(
            encoding="utf-8"
        )
    )
    return payload["targets"][0]


class RankTwoPentadTest(unittest.TestCase):
    def test_exact_target_manifest_shape(self) -> None:
        payload = json.loads(
            (ROOT / "d6_k7_ranktwo_pentad_targets.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(payload["status"], "COMPLETE")
        self.assertEqual(len(payload["targets"]), 53)
        self.assertEqual(
            len({item["graph_index"] for item in payload["targets"]}), 50
        )
        self.assertTrue(all(
            len(item["graph_n"]) == 12
            and item["rank_upper"] == 7
            and item["maximum_clique_size"] == 5
            and item["zmask"] == 0
            for item in payload["targets"]
        ))

    def test_elimination_derivation(self) -> None:
        result = verifier.verify_elimination_derivation()
        self.assertEqual(result, {
            "formal_variables": 10, "terms": 12, "degree": 5,
        })

    def test_generic_two_dimensional_gram_identity(self) -> None:
        result = verifier.verify_generic_rank_two_gram_identity()
        self.assertEqual(result["residual_terms"], 0)

    def test_all_labelings_are_the_same_up_to_sign(self) -> None:
        result = verifier.verify_pentad_permutation_orbit()
        self.assertEqual(result["same_sign"], 60)
        self.assertEqual(result["opposite_sign"], 60)
        self.assertEqual(result["genuinely_distinct_up_to_sign"], 1)

    def test_numeric_rank_two_gram_identity(self) -> None:
        vectors = ((1, 2), (3, -1), (2, 5), (-2, 4), (7, 3))
        pair_polynomials = {}
        for first in range(5):
            for second in range(first + 1, 5):
                value = sum(
                    vectors[first][coordinate] * vectors[second][coordinate]
                    for coordinate in range(2)
                )
                pair_polynomials[(first, second)] = {(0,): Fraction(value)}
        result = production.pentad_polynomial(
            range(5), pair_polynomials, variables=1
        )
        self.assertEqual(result, {})

    def test_first_core_certificate(self) -> None:
        target = first_target()
        graph_n = tuple(target["graph_n"])
        clique_mask = sum(
            1 << vertex for vertex in target["fixed_maximum_clique"]
        )
        certificate = production.find_coefficientwise_certificate(
            graph_n, clique_mask
        )
        self.assertIsNotNone(certificate)
        positive = production.verify_coefficientwise_certificate(
            graph_n, clique_mask, certificate
        )
        self.assertEqual(len(positive), 201)

    def test_first_full_record_independent_reconstruction(self) -> None:
        target = first_target()
        record = full_runner.evaluate_cover((0, target))
        checked = verifier.verify_full_cover_record(target, record, 0)
        self.assertEqual(checked["status"], "REJECTED")
        self.assertEqual(checked["pentads_checked"], 21)
        self.assertTrue(checked["certificate_checked"])

    def test_realizable_18_point_control_has_no_k7(self) -> None:
        adjacency = tuple(lower_bound_18_graph())
        prior.validate_graph(adjacency)
        self.assertIsNone(next(prior.clique_masks(adjacency, 7), None))

    def test_tampered_positive_coefficient_is_rejected(self) -> None:
        target = first_target()
        graph_n = tuple(target["graph_n"])
        clique_mask = sum(
            1 << vertex for vertex in target["fixed_maximum_clique"]
        )
        certificate = production.find_coefficientwise_certificate(
            graph_n, clique_mask
        )
        certificate = json.loads(json.dumps(certificate))
        certificate["positive_polynomial"][0]["coefficient"][0] += 1
        with self.assertRaises(ValueError):
            production.verify_coefficientwise_certificate(
                graph_n, clique_mask, certificate
            )


if __name__ == "__main__":
    unittest.main()
