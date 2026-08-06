#!/usr/bin/env python3
"""Controls for the plain-Python dimension-six reference filters."""

from __future__ import annotations

import json
import unittest
from itertools import combinations
from pathlib import Path

from d6_reference_filters import (
    add_clique,
    add_edge,
    evaluate_sample,
    is_clique,
    k6_clique_hall,
    k6_clique_hall_seed,
    k7_bounded_cover_seed,
    k7_clique_hall,
    k7_clique_hall_seed,
    k7_disjoint_edge_bounded_cover,
    k7_tight_cover_matching,
    k7_tight_cover_matching_seed,
    validate_graph,
)


ROOT = Path(__file__).resolve().parent
SAMPLE = ROOT / "d6_reference_sample.json"


def reflected_facet_graph(count: int) -> list[int]:
    """A regular K7 plus zero, one, or two distinct facet reflections."""

    if count not in (0, 1, 2):
        raise ValueError("this control supplies at most two reflections")
    adj = [0] * (7 + count)
    add_clique(adj, range(7))
    for offset in range(count):
        x = 7 + offset
        missed_seed_vertex = offset
        for q in range(7):
            if q != missed_seed_vertex:
                add_edge(adj, x, q)
    # Distinct facet reflections are at squared distance 16/9, so there is
    # deliberately no required edge between vertices 7 and 8.
    validate_graph(adj)
    return adj


def lower_bound_18_graph() -> list[int]:
    """The realizable 16-point half-cube graph plus two unit apices."""

    base = [x for x in range(32) if x.bit_count() % 2 == 1]
    adj = [0] * 18
    for i, x in enumerate(base):
        for j, y in enumerate(base[:i]):
            if (x ^ y).bit_count() == 2:
                add_edge(adj, i, j)
    for apex in (16, 17):
        for v in range(16):
            add_edge(adj, apex, v)
    validate_graph(adj)
    return adj


def k7_hall_failure() -> list[int]:
    """A K7 seed and a 3-clique using only two allowed coordinates."""

    adj = [0] * 10
    add_clique(adj, range(7))
    add_clique(adj, range(7, 10))
    for x in range(7, 10):
        for q in range(2, 7):
            add_edge(adj, x, q)
    validate_graph(adj)
    return adj


def k6_hall_failure() -> list[int]:
    """A K6 seed and a 4-clique needing two defects plus the normal slot."""

    adj = [0] * 10
    add_clique(adj, range(6))
    add_clique(adj, range(6, 10))
    for x in range(6, 10):
        for q in range(2, 6):
            add_edge(adj, x, q)
    validate_graph(adj)
    return adj


def ineligible_disjoint_edge_failure() -> list[int]:
    """A disjoint-defect edge whose two endpoints have at most two defects."""

    adj = [0] * 9
    add_clique(adj, range(7))
    for q in range(2, 7):
        add_edge(adj, 7, q)  # D_7 = {0,1}
    for q in (0, 1, 4, 5, 6):
        add_edge(adj, 8, q)  # D_8 = {2,3}
    add_edge(adj, 7, 8)
    validate_graph(adj)
    return adj


def cardinality_cover_failure() -> list[int]:
    """An alpha-at-most-two control whose relative-seed L is K_{1,8}.

    Vertex 7 is the ineligible center with D={0}.  Vertices 8,...,15 are
    eligible leaves, two of each of four three-defect types.  The center is
    adjacent to every leaf.  Leaves are adjacent exactly when their defect
    masks intersect, so no leaf-leaf edge enters L.  Consequently an eligible
    cover must contain all eight leaves: cap 8 passes and cap 7 fails.
    """

    adj = [0] * 16
    add_clique(adj, range(7))
    center = 7
    for q in range(1, 7):
        add_edge(adj, center, q)
    defect_types = (
        {1, 2, 3},
        {1, 2, 3},
        {4, 5, 6},
        {4, 5, 6},
        {1, 2, 4},
        {1, 2, 4},
        {3, 5, 6},
        {3, 5, 6},
    )
    for offset, defects in enumerate(defect_types):
        leaf = 8 + offset
        for q in range(7):
            if q not in defects:
                add_edge(adj, leaf, q)
        add_edge(adj, center, leaf)
    for i, left_defects in enumerate(defect_types):
        for j, right_defects in enumerate(defect_types[:i]):
            if left_defects.intersection(right_defects):
                add_edge(adj, 8 + i, 8 + j)
    validate_graph(adj)
    return adj


def tight_cover_matching_failure() -> list[int]:
    """An alpha-at-most-two tight K1,7 cover with no coordinate matching."""

    adj = [0] * 15
    add_clique(adj, range(7))
    center = 7
    for q in range(1, 7):
        add_edge(adj, center, q)  # D_center = {0}
    defect_types = (
        {1, 2, 3},
        {1, 2, 3},
        {4, 5, 6},
        {4, 5, 6},
        {1, 2, 4},
        {1, 2, 4},
        {3, 5, 6},
    )
    for offset, defects in enumerate(defect_types):
        leaf = 8 + offset
        for q in range(7):
            if q not in defects:
                add_edge(adj, leaf, q)
        add_edge(adj, center, leaf)
    for i, left_defects in enumerate(defect_types):
        for j, right_defects in enumerate(defect_types[:i]):
            if left_defects.intersection(right_defects):
                add_edge(adj, 8 + i, 8 + j)
    validate_graph(adj)
    return adj


def has_independent_triple(adj: list[int]) -> bool:
    """Return whether the graph violates the almost-equidistant graph rule."""

    return any(
        not (adj[u] & (1 << v))
        and not (adj[u] & (1 << w))
        and not (adj[v] & (1 << w))
        for u, v, w in combinations(range(len(adj)), 3)
    )


class ReferenceFilterControls(unittest.TestCase):
    def test_reflected_facet_positive_controls(self) -> None:
        for count in (1, 2):
            with self.subTest(reflections=count):
                adj = reflected_facet_graph(count)
                hall = k7_clique_hall(adj)
                cover = k7_disjoint_edge_bounded_cover(adj)
                self.assertTrue(hall.applicable)
                self.assertTrue(cover.applicable)
                self.assertFalse(hall.rejected)
                self.assertFalse(cover.rejected)
                self.assertFalse(k7_tight_cover_matching(adj).rejected)

    def test_lower_bound_18_passes_k6_rule(self) -> None:
        decision = k6_clique_hall(lower_bound_18_graph())
        self.assertTrue(decision.applicable)
        self.assertGreater(decision.seeds_checked, 0)
        self.assertFalse(decision.rejected)

    def test_k7_hall_failure_and_explicit_k8_witness(self) -> None:
        adj = k7_hall_failure()
        witness = k7_clique_hall_seed(adj, tuple(range(7)))
        self.assertIsNotNone(witness)
        assert witness is not None
        self.assertEqual(witness["hall_limit"], 2)
        self.assertEqual(len(witness["outside_clique"]), 3)
        self.assertEqual(len(witness["forced_K8"]), 8)
        self.assertTrue(is_clique(adj, witness["forced_K8"]))
        self.assertTrue(k7_clique_hall(adj).rejected)

    def test_k6_hall_failure_and_explicit_k8_witness(self) -> None:
        adj = k6_hall_failure()
        witness = k6_clique_hall_seed(adj, tuple(range(6)))
        self.assertIsNotNone(witness)
        assert witness is not None
        self.assertEqual(witness["hall_limit"], 3)
        self.assertEqual(len(witness["outside_clique"]), 4)
        self.assertEqual(len(witness["forced_K8"]), 8)
        self.assertTrue(is_clique(adj, witness["forced_K8"]))
        # Every K6 Hall failure contains the displayed K8 and therefore is
        # outside the K6-only population.  The graph-level filter correctly
        # reports not-applicable; the seed-level routine is the negative
        # control for the algebra itself.
        self.assertFalse(k6_clique_hall(adj).applicable)

    def test_ineligible_disjoint_edge_failure(self) -> None:
        adj = ineligible_disjoint_edge_failure()
        witness = k7_bounded_cover_seed(adj, tuple(range(7)))
        self.assertIsNotNone(witness)
        assert witness is not None
        self.assertEqual(witness["failure_kind"], "ineligible_disjoint_edge")
        self.assertEqual(witness["ineligible_edge"], [7, 8])
        self.assertTrue(k7_disjoint_edge_bounded_cover(adj).rejected)

    def test_tight_cover_matching_isolating_failure(self) -> None:
        adj = tight_cover_matching_failure()
        seed = tuple(range(7))
        self.assertFalse(has_independent_triple(adj))
        # This is an isolating control relative to the distinguished seed:
        # Hall and the old cap-seven cover rule both pass, while cap six fails
        # with exact minimum cover size seven.  (The prescribed graph has
        # other K7 seeds, so helper-level isolation is the meaningful check.)
        self.assertIsNone(k7_clique_hall_seed(adj, seed))
        self.assertIsNone(k7_bounded_cover_seed(adj, seed, cover_bound=7))
        old_graph_decision = k7_disjoint_edge_bounded_cover(adj)
        self.assertTrue(old_graph_decision.rejected)
        self.assertNotEqual(old_graph_decision.witness["seed"], list(seed))
        cap_six = k7_bounded_cover_seed(adj, seed, cover_bound=6)
        self.assertIsNotNone(cap_six)
        assert cap_six is not None
        self.assertEqual(cap_six["minimum_eligible_cover_size"], 7)

        witness = k7_tight_cover_matching_seed(adj, seed)
        self.assertIsNotNone(witness)
        assert witness is not None
        self.assertEqual(
            witness["failure_kind"], "tight_cover_no_perfect_matching"
        )
        self.assertEqual(witness["size7_covers_checked"], 1)
        self.assertEqual(
            witness["L_edges"], [[7, leaf] for leaf in range(8, 15)]
        )
        self.assertEqual(witness["eligible_vertices"], list(range(8, 15)))
        failed = witness["failed_cover_examples"]
        self.assertEqual(len(failed), 1)
        self.assertEqual(failed[0]["cover"], list(range(8, 15)))
        self.assertEqual(failed[0]["coordinate_union"], list(range(1, 7)))
        self.assertTrue(k7_tight_cover_matching(adj).rejected)

    def test_cover_cardinality_failure_is_not_a_forbidden_edge(self) -> None:
        adj = cardinality_cover_failure()
        self.assertFalse(has_independent_triple(adj))
        witness = k7_bounded_cover_seed(adj, tuple(range(7)))
        self.assertIsNotNone(witness)
        assert witness is not None
        self.assertEqual(witness["failure_kind"], "cover_cardinality")
        self.assertIsNone(witness["ineligible_edge"])
        self.assertEqual(witness["minimum_eligible_cover_size"], 8)
        self.assertEqual(
            sorted(witness["L_edges"]), [[7, leaf] for leaf in range(8, 16)]
        )
        self.assertEqual(witness["eligible_vertices"], list(range(8, 16)))
        eligible = set(witness["eligible_vertices"])
        self.assertTrue(
            all(eligible.intersection(edge) for edge in witness["L_edges"])
        )
        self.assertIsNone(
            k7_bounded_cover_seed(adj, tuple(range(7)), cover_bound=8)
        )
        self.assertTrue(k7_disjoint_edge_bounded_cover(adj).rejected)

    def test_fixed_real_sample(self) -> None:
        if not SAMPLE.exists():
            self.skipTest("run extract_d6_reference_sample.py first")
        with SAMPLE.open(encoding="utf-8") as stream:
            sample = json.load(stream)
        self.assertEqual(sample["schema"], 2)
        counts = sample["stratum_counts"]
        self.assertGreater(counts["K7_certified"], 0)
        self.assertGreater(counts["K7_deferred"], 0)
        self.assertGreater(counts["K6_only_deferred"], 0)
        self.assertEqual(counts["K6_only_certified"], 0)
        self.assertEqual(counts["K7_deferred_tight_targeted"], 4)
        result = evaluate_sample(sample, include_graph_results=True)
        self.assertEqual(result["sample_graphs"], len(sample["graphs"]))
        # Both Hall rules are K8 checks in disguise; the source corpus was
        # generated K8-free, so neither may reject a sampled graph.
        self.assertTrue(
            all(
                stratum["rejected"]["K7_clique_Hall"] == 0
                and stratum["rejected"]["K6_clique_Hall"] == 0
                for stratum in result["strata"].values()
            )
        )
        for name in ("K7_certified", "K7_deferred", "K6_only_deferred"):
            self.assertEqual(
                result["strata"][name]["rejected"]
                ["K7_tight_cover_matching"],
                0,
            )
        self.assertEqual(
            result["strata"]["K7_certified"]["rejected"]
            ["K7_disjoint_edge_bounded_cover"],
            16,
        )
        self.assertEqual(
            result["strata"]["K7_deferred"]["rejected"]
            ["K7_disjoint_edge_bounded_cover"],
            14,
        )
        self.assertEqual(
            result["strata"]["K6_only_deferred"]["rejected"]
            ["K7_disjoint_edge_bounded_cover"],
            0,
        )
        targeted = result["strata"]["K7_deferred_tight_targeted"]
        self.assertEqual(
            targeted["rejected"]["K7_disjoint_edge_bounded_cover"], 2
        )
        self.assertEqual(
            targeted["rejected"]["K7_tight_cover_matching"], 4
        )
        expected_old_cover = {2_751_055, 2_887_126}
        expected_tight = expected_old_cover | {3_368_019, 3_936_501}
        observed_old_cover = set()
        observed_tight = set()
        for graph in result["graph_results"]:
            if graph["stratum"] != "K7_deferred_tight_targeted":
                continue
            decisions = graph["decisions"]
            if decisions["K7_disjoint_edge_bounded_cover"]["rejected"]:
                observed_old_cover.add(graph["index"])
            if decisions["K7_tight_cover_matching"]["rejected"]:
                observed_tight.add(graph["index"])
        self.assertEqual(observed_old_cover, expected_old_cover)
        self.assertEqual(observed_tight, expected_tight)


if __name__ == "__main__":
    unittest.main(verbosity=2)
