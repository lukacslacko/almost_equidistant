#!/usr/bin/env python3
"""Controls for the exact labeled-support propagation refinement."""

from __future__ import annotations

import csv
import unittest
from itertools import product
from pathlib import Path
from tempfile import TemporaryDirectory

import d6_k7_rank_reference as reference
from d6_k7_support_propagation import (
    DECISION_FIELDS,
    analyze_support_assignment,
    atomic_decisions,
    labeled_support_families,
    load_completed,
    propagate_n_masks,
    propagate_one_mask,
    refine_baseline_passing_cover,
)


def graph_from_edges(n: int, edges: list[tuple[int, int]]) -> tuple[int, ...]:
    adjacency = [0] * n
    for first, second in edges:
        adjacency[first] |= 1 << second
        adjacency[second] |= 1 << first
    return tuple(adjacency)


class LabeledEnumerationControls(unittest.TestCase):
    def test_matches_direct_labeled_product(self) -> None:
        cases = (
            (0b000111,),
            (0b001111, 0b011110),
            (0b001111, 0b001111, 0b110011),
            (0b011111, 0b101111, 0b110111),
        )
        for allowed in cases:
            observed = set(labeled_support_families(allowed))
            expected = {
                chosen
                for chosen in product(
                    *(reference.SUPPORT_DOMAINS[mask] for mask in allowed)
                )
                if reference.support_family_valid(chosen)
            }
            self.assertEqual(observed, expected)

    def test_equal_domains_remain_labeled(self) -> None:
        families = set(labeled_support_families((0b1111, 0b1111)))
        first = (0b0111, 0b1011)
        second = tuple(reversed(first))
        self.assertIn(first, families)
        self.assertIn(second, families)


class PropagationControls(unittest.TestCase):
    def test_fixed_point_cascade(self) -> None:
        # The first support removes coordinate 0.  Only then does the second
        # support have singleton intersection and remove coordinate 1.
        supports = (0b0011001, 0b0001011)
        closed, deletions = propagate_one_mask(0b0000111, supports)
        self.assertEqual(closed, 0b0000100)
        self.assertEqual(deletions, 2)

    def test_order_independent_simultaneous_closure(self) -> None:
        supports = (0b0011001, 0b0001011)
        forward = propagate_one_mask(0b0000111, supports)
        backward = propagate_one_mask(0b0000111, tuple(reversed(supports)))
        self.assertEqual(forward, backward)

    def test_multiple_masks_and_deletion_total(self) -> None:
        masks, deletions = propagate_n_masks(
            (0b0000111, 0b1110000),
            (0b0011001, 0b0001011),
        )
        self.assertEqual(masks, (0b0000100, 0b1100000))
        self.assertEqual(deletions, 3)

    def test_empty_mask_is_an_exact_failure(self) -> None:
        graph = (0,)
        analysis = analyze_support_assignment(
            graph,
            (0b0000111,),
            (0b0000001,),
            reference.ZeroForcingSolver(),
            reference.CliqueStructureSolver(),
        )
        self.assertEqual(analysis.failure, "empty_propagated_mask")

    def test_required_edge_disjoint_is_an_exact_failure(self) -> None:
        graph = graph_from_edges(2, [(0, 1)])
        analysis = analyze_support_assignment(
            graph,
            (),
            (0b0000001, 0b0000010),
            reference.ZeroForcingSolver(),
            reference.CliqueStructureSolver(),
        )
        self.assertEqual(analysis.failure, "disjoint_required_edge")

    def test_no_fixed_support_is_identity(self) -> None:
        graph = graph_from_edges(2, [(0, 1)])
        masks = (0b0000011, 0b0000110)
        analysis = analyze_support_assignment(
            graph,
            (),
            masks,
            reference.ZeroForcingSolver(),
            reference.CliqueStructureSolver(),
        )
        self.assertEqual(analysis.propagated_masks, masks)
        self.assertEqual(analysis.deletions, 0)

    def test_propagated_term_rank_triggers_existing_clique_rule(self) -> None:
        # Both allowed masks initially meet through coordinate 0 and support a
        # required edge.  Orthogonality to S_z deletes 0 from each, collapsing
        # both masks onto coordinate 3.  Thus nu_N=1 and the required K2 is a
        # positive-definite clique of size two, an exact contradiction.
        graph = graph_from_edges(2, [(0, 1)])
        analysis = analyze_support_assignment(
            graph,
            (0b0000111,),
            (0b0001001, 0b0001001),
            reference.ZeroForcingSolver(),
            reference.CliqueStructureSolver(),
        )
        self.assertEqual(analysis.propagated_masks, (0b0001000, 0b0001000))
        self.assertEqual(analysis.n_term_rank, 1)
        self.assertEqual(analysis.k_rank_upper, 1)
        self.assertEqual(analysis.failure, "clique")

    def test_cover_refinement_positive_control(self) -> None:
        # One Z support and one N mask overlap in two coordinates, so no
        # deletion occurs.  A one-vertex K has rank upper bound one and passes.
        result = refine_baseline_passing_cover(
            (0,),
            (0b0000111,),
            (0b0000011,),
            reference.ZeroForcingSolver(),
            reference.CliqueStructureSolver(),
        )
        self.assertTrue(result.passes)
        self.assertEqual(result.passing_assignment, (0b0000111,))


class CheckpointControls(unittest.TestCase):
    @staticmethod
    def decision_row(index: int) -> list[object]:
        row: list[object] = [0] * len(DECISION_FIELDS)
        row[0] = index
        row[1] = 1
        row[2] = "SURVIVOR"
        row[3] = "SURVIVOR"
        return row

    def test_uncommitted_torn_tail_is_discarded(self) -> None:
        with TemporaryDirectory() as directory:
            path = Path(directory) / "decisions.partial"
            atomic_decisions(path, (self.decision_row(101),))
            with path.open("a", encoding="ascii") as stream:
                stream.write("102\t1\tSURV")
            completed = load_completed(path, (101, 102), committed_rows=1)
            self.assertEqual(len(completed), 1)
            with path.open(newline="", encoding="ascii") as stream:
                rows = list(csv.reader(stream, delimiter="\t"))
            self.assertEqual(len(rows), 2)
            self.assertEqual(len(rows[1]), len(DECISION_FIELDS))

    def test_truncated_committed_row_is_rejected(self) -> None:
        with TemporaryDirectory() as directory:
            path = Path(directory) / "decisions.partial"
            with path.open("w", encoding="ascii") as stream:
                stream.write("\t".join(DECISION_FIELDS) + "\n")
                stream.write("101\t1\tSURV")
            with self.assertRaisesRegex(ValueError, "has 3 fields"):
                load_completed(path, (101,), committed_rows=1)


if __name__ == "__main__":
    unittest.main()
