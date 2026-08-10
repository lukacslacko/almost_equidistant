#!/usr/bin/env python3
"""Focused controls for the exact saturated-singleton K6 basis layer."""

from __future__ import annotations

import hashlib
import unittest
from types import SimpleNamespace
from unittest.mock import patch

import d6_k6_saturated_singleton_basis as production
import verify_d6_k6_saturated_singleton_basis as verifier


def graph_with_edges(
    size: int, edges: list[tuple[int, int]],
) -> tuple[int, ...]:
    rows = [0] * size
    for left, right in edges:
        rows[left] |= 1 << right
        rows[right] |= 1 << left
    return tuple(rows)


def encoded_graph(size: int, encoded: int) -> tuple[int, ...]:
    edges = []
    position = 0
    for right in range(size):
        for left in range(right):
            if encoded & (1 << position):
                edges.append((left, right))
            position += 1
    return graph_with_edges(size, edges)


def singleton_basis_instance(
    neighbourhoods: tuple[int, ...],
    defects: tuple[int, ...] | None = None,
) -> tuple[tuple[int, ...], SimpleNamespace, int, dict[int, int]]:
    """Make a six-basis/seven-remaining synthetic extension instance."""

    if len(neighbourhoods) != 7:
        raise ValueError("a K6 seed has exactly thirteen outside vertices")
    edges = []
    for offset, support in enumerate(neighbourhoods):
        remaining = 6 + offset
        for coordinate in range(6):
            if support & (1 << coordinate):
                edges.append((coordinate, remaining))
    if defects is None:
        defect_rows = (0b111111,) * 13
    else:
        if len(defects) != 13:
            raise ValueError("defect rows must cover every outside vertex")
        defect_rows = defects
    instance = SimpleNamespace(outside=tuple(range(13)), defects=defect_rows)
    assignments = {coordinate: 1 << coordinate for coordinate in range(6)}
    return graph_with_edges(13, edges), instance, 0b111111, assignments


class SaturatedSingletonBasisTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        records, indices, dependencies = production.load_input()
        cls.records = {int(row["index"]): row for row in records}
        cls.indices = indices
        cls.dependencies = dependencies
        independent, independent_indices = verifier.load_input()
        cls.independent = {int(row["index"]): row for row in independent}
        cls.independent_indices = independent_indices

    def test_frozen_251_parent_boundary(self) -> None:
        self.assertEqual(len(self.records), production.EXPECTED_INPUT)
        self.assertEqual(production.EXPECTED_INPUT, 251)
        self.assertEqual(self.indices, self.independent_indices)
        self.assertEqual(
            production.stable_hash(self.indices),
            production.EXPECTED_INPUT_SHA256,
        )
        self.assertEqual(self.dependencies, production.EXPECTED_DEPENDENCIES)
        self.assertEqual(
            verifier.EXPECTED_PRODUCTION_DEPENDENCIES,
            production.EXPECTED_DEPENDENCIES,
        )
        self.assertEqual(
            production.EXPECTED_REJECTED_INDICES,
            tuple(verifier.EXPECTED_REJECTED_INDICES),
        )
        self.assertEqual(production.EXPECTED_REJECTIONS, 2)
        self.assertIn(405_497, self.indices)
        self.assertIn(3_138_618, self.indices)
        self.assertIn(3_673_988, self.indices)
        for name, expected in production.EXPECTED_DEPENDENCIES.items():
            with self.subTest(parent_artifact=name):
                self.assertEqual(production.sha256(production.ROOT / name), expected)

    def test_enumeration_and_independent_parity_solver_agree(self) -> None:
        # Exhaust all pairs of the 64 possible reconstructed basis supports,
        # both with and without a required edge.
        for left in range(64):
            for right in range(64):
                for required in range(2):
                    adjacency = encoded_graph(2, required)
                    enumerated, _ = production.sign_csp(
                        (left, right), adjacency
                    )
                    propagated, _ = verifier.reconstructed_sign_system(
                        (left, right), adjacency
                    )
                    self.assertEqual(
                        enumerated,
                        propagated,
                        (left, right, required),
                    )

        # Exercise parity cycles and all support degrees/intersection shapes
        # on triples, while keeping the unit test compact.
        representatives = (0, 1, 2, 3, 5, 7, 10, 15, 21, 31, 42, 63)
        for first in representatives:
            for second in representatives:
                for third in representatives:
                    neighbourhoods = (first, second, third)
                    for encoded in range(8):
                        adjacency = encoded_graph(3, encoded)
                        enumerated, _ = production.sign_csp(
                            neighbourhoods, adjacency
                        )
                        propagated, _ = verifier.reconstructed_sign_system(
                            neighbourhoods, adjacency
                        )
                        self.assertEqual(
                            enumerated,
                            propagated,
                            (neighbourhoods, encoded),
                        )

    def test_exact_empty_nonempty_edge_rules(self) -> None:
        edge = graph_with_edges(2, [(0, 1)])

        # An empty apex and a degree-four reconstructed point are adjacent
        # exactly when their two algebraic signs agree.
        enumerated, enumeration = production.sign_csp((0, 0b001111), edge)
        propagated, parity = verifier.reconstructed_sign_system(
            (0, 0b001111), edge
        )
        self.assertTrue(enumerated)
        self.assertTrue(propagated)
        self.assertEqual(enumeration["signs"][0], enumeration["signs"][1])
        self.assertIn([0, 1, 0], parity["constraints"])

        # Degree three cannot meet the exact empty/nonempty edge equation.
        enumerated, _ = production.sign_csp((0, 0b000111), edge)
        propagated, detail = verifier.reconstructed_sign_system(
            (0, 0b000111), edge
        )
        self.assertFalse(enumerated)
        self.assertFalse(propagated)
        self.assertEqual(
            detail["reason"], "empty_apex_edge_has_wrong_basis_degree"
        )

        # Two empty apices may coexist only with opposite normal signs and
        # can never be joined by a required unit edge.
        no_edge = graph_with_edges(2, [])
        enumerated, enumeration = production.sign_csp((0, 0), no_edge)
        propagated, parity = verifier.reconstructed_sign_system((0, 0), no_edge)
        self.assertTrue(enumerated)
        self.assertTrue(propagated)
        self.assertNotEqual(enumeration["signs"][0], enumeration["signs"][1])
        self.assertIn([0, 1, 1], parity["constraints"])
        self.assertFalse(production.sign_csp((0, 0), edge)[0])
        self.assertEqual(
            verifier.reconstructed_sign_system((0, 0), edge)[1]["reason"],
            "required_edge_between_empty_apices",
        )

    def test_exact_nonempty_edge_and_collision_rules(self) -> None:
        edge = graph_with_edges(2, [(0, 1)])
        disjoint_pairs = (0b000011, 0b001100)
        enumerated, enumeration = production.sign_csp(disjoint_pairs, edge)
        propagated, parity = verifier.reconstructed_sign_system(
            disjoint_pairs, edge
        )
        self.assertTrue(enumerated)
        self.assertTrue(propagated)
        self.assertNotEqual(enumeration["signs"][0], enumeration["signs"][1])
        self.assertIn([0, 1, 1], parity["constraints"])

        # Here m+n-2k=2 rather than 4, so no choice of the two roots works.
        overlapping_pairs = (0b000011, 0b000110)
        self.assertFalse(production.sign_csp(overlapping_pairs, edge)[0])
        propagated, detail = verifier.reconstructed_sign_system(
            overlapping_pairs, edge
        )
        self.assertFalse(propagated)
        self.assertEqual(
            detail["reason"], "nonempty_edge_neighbourhood_equation"
        )

        # A repeated reconstructed neighbourhood needs opposite signs.  A
        # third copy makes pairwise distinctness inconsistent.
        repeated = (0b10101, 0b10101)
        enumerated, enumeration = production.sign_csp(
            repeated, graph_with_edges(2, [])
        )
        propagated, parity = verifier.reconstructed_sign_system(
            repeated, graph_with_edges(2, [])
        )
        self.assertTrue(enumerated)
        self.assertTrue(propagated)
        self.assertNotEqual(enumeration["signs"][0], enumeration["signs"][1])
        self.assertIn([0, 1, 1], parity["constraints"])
        triple = (0b10101,) * 3
        self.assertFalse(
            production.sign_csp(triple, graph_with_edges(3, []))[0]
        )
        self.assertFalse(
            verifier.reconstructed_sign_system(
                triple, graph_with_edges(3, [])
            )[0]
        )

    def test_reconstructed_support_allowed_and_assigned_rules(self) -> None:
        neighbourhoods = (0b11, 0b1, 0b10, 0b100, 0b101, 0b110, 0b111)
        adjacency, instance, basis, assignments = singleton_basis_instance(
            neighbourhoods
        )

        # Allowed seed defects are supersets: unused allowed coordinates may
        # still be zero.  Both implementations accept this strict-subset
        # reconstruction, including one empty actual support.
        for extension in (
            production.saturated_singleton_extension,
            verifier.saturated_singleton_extension,
        ):
            with self.subTest(extension=extension.__module__, case="optional"):
                feasible, detail = extension(
                    adjacency, instance, basis, assignments
                )
                self.assertTrue(feasible, detail)
                self.assertTrue(detail["applicable"])

        # A basis-required coordinate outside the seed's allowed defect mask
        # is an exact contradiction.
        defect_rows = [0b111111] * 13
        defect_rows[6] = 0b1
        forbidden_instance = SimpleNamespace(
            outside=instance.outside, defects=tuple(defect_rows)
        )
        for extension in (
            production.saturated_singleton_extension,
            verifier.saturated_singleton_extension,
        ):
            with self.subTest(extension=extension.__module__, case="allowed"):
                feasible, detail = extension(
                    adjacency, forbidden_instance, basis, assignments
                )
                self.assertFalse(feasible)
                self.assertEqual(
                    detail["reason"],
                    "reconstructed_support_not_allowed_by_seed_defects",
                )
                self.assertEqual(detail["vertex"], 6)

        # If the other light bin already assigned this point's actual
        # support, reconstruction must equal that complete leaf assignment.
        mismatched = dict(assignments)
        mismatched[6] = 0b1
        for extension in (
            production.saturated_singleton_extension,
            verifier.saturated_singleton_extension,
        ):
            with self.subTest(extension=extension.__module__, case="assigned"):
                feasible, detail = extension(
                    adjacency, instance, basis, mismatched
                )
                self.assertFalse(feasible)
                self.assertEqual(
                    detail["reason"],
                    "assigned_light_support_differs_from_reconstruction",
                )
                self.assertEqual(detail["vertex"], 6)

    def test_failed_singleton_leaves_do_not_discard_alternatives(self) -> None:
        # Four completed singleton-basis leaves fail for this graph.  Every
        # affected actual-support DFS nevertheless finds another leaf, all
        # 35 K6 seeds pass, and both engines retain the graph.
        produced = production.evaluate_record(self.records[405_497])
        checked = verifier.evaluate_record(self.independent[405_497])
        self.assertFalse(produced["rejected"])
        self.assertFalse(checked["rejected"])
        self.assertEqual(produced["seeds_checked"], 35)
        self.assertEqual(produced["counts"]["actual_support_searches"], 35)
        self.assertEqual(produced["counts"]["joint_coloring_passed"], 35)
        self.assertEqual(
            produced["singleton_basis_counts"],
            {
                "saturated_singleton_branches_checked": 4,
                "saturated_singleton_branches_failed": 4,
            },
        )
        self.assertEqual(len(produced["singleton_basis_witnesses"]), 4)
        self.assertTrue(all(
            witness["reason"]
            == "reconstructed_support_not_allowed_by_seed_defects"
            for witness in produced["singleton_basis_witnesses"]
        ))

    def test_both_fixed_rejections_and_independent_certificates(self) -> None:
        expected = {
            3_138_618: {
                "seeds_checked": 22,
                "seed_mask": 204_888,
                "seed": [3, 4, 6, 13, 16, 17],
                "branches": 1,
            },
            3_673_988: {
                "seeds_checked": 16,
                "seed_mask": 279_716,
                "seed": [2, 5, 7, 10, 14, 18],
                "branches": 2,
            },
        }
        for index, wanted in expected.items():
            with self.subTest(index=index):
                produced = production.evaluate_record(self.records[index])
                checked = verifier.evaluate_record(self.independent[index])
                self.assertTrue(produced["rejected"])
                self.assertTrue(checked["rejected"])
                self.assertEqual(produced["seeds_checked"], wanted["seeds_checked"])
                self.assertEqual(
                    produced["first_impossible_seed_mask"], wanted["seed_mask"]
                )
                self.assertEqual(produced["first_impossible_seed"], wanted["seed"])
                self.assertEqual(checked["seed_mask"], wanted["seed_mask"])
                self.assertEqual(checked["seed"], wanted["seed"])
                self.assertEqual(
                    produced["singleton_basis_counts"],
                    {
                        "saturated_singleton_branches_checked": wanted["branches"],
                        "saturated_singleton_branches_failed": wanted["branches"],
                    },
                )
                self.assertEqual(
                    checked["failures"],
                    produced["certificate"]["Z0_failures"],
                )

    def test_all_positive_seeds_and_import_independence(self) -> None:
        self.assertEqual(
            production.positive_control(), {"passed": True, "K6_seeds": 32}
        )
        self.assertEqual(
            verifier.positive_control(), {"passed": True, "K6_seeds": 32}
        )
        verifier.assert_import_independence()

    def test_package_source_and_launch_porcelain_binding(self) -> None:
        self.assertEqual(production.PACKAGE_SOURCES, verifier.PACKAGE_SOURCES)
        self.assertEqual(
            set(production.PACKAGE_SOURCES),
            {
                "d6_k6_saturated_singleton_basis.py",
                "verify_d6_k6_saturated_singleton_basis.py",
                "test_d6_k6_saturated_singleton_basis.py",
                "d6_k6_saturated_singleton_basis.md",
            },
        )
        source_hashes = {
            name: verifier.sha256(verifier.ROOT / name)
            for name in verifier.PACKAGE_SOURCES
        }
        lines = ["?? scratch-proof-search.json"]
        porcelain = "\n".join(lines)
        report = {
            "source_sha256": source_hashes,
            "execution": {
                "git": {
                    "available": True,
                    "commit": "1" * 40,
                    "branch": "codex/dimension6",
                    "tracked_clean": True,
                    "dirty": True,
                    "porcelain_lines": lines,
                    "porcelain_sha256": hashlib.sha256(
                        porcelain.encode("utf-8")
                    ).hexdigest(),
                    "committed_source_sha256": source_hashes,
                }
            },
        }

        def committed_blob(command, **_kwargs):
            name = command[-1].split(":", 1)[1]
            return SimpleNamespace(stdout=(verifier.ROOT / name).read_bytes())

        with patch.object(verifier.subprocess, "run", side_effect=committed_blob):
            binding = verifier.validate_git_provenance(report)
        self.assertTrue(binding["tracked_sources_clean_at_launch"])
        self.assertEqual(binding["untracked_entries"], 1)

        tracked_lines = [" M d6_k6_saturated_singleton_basis.py"]
        tracked = dict(report["execution"]["git"])
        tracked.update({
            "porcelain_lines": tracked_lines,
            "porcelain_sha256": hashlib.sha256(
                tracked_lines[0].encode("utf-8")
            ).hexdigest(),
        })
        with self.assertRaisesRegex(ValueError, "tracked source was dirty"):
            verifier.validate_porcelain_binding(tracked)

        tampered = dict(report)
        tampered["source_sha256"] = dict(source_hashes)
        tampered["source_sha256"][production.PACKAGE_SOURCES[0]] = "0" * 64
        with self.assertRaisesRegex(ValueError, "source boundary differs"):
            verifier.validate_git_provenance(tampered)


if __name__ == "__main__":
    unittest.main(verbosity=2)
