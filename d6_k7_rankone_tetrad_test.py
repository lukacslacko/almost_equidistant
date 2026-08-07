#!/usr/bin/env python3
"""Tests for the exact K7 rank-one Schur tetrad pilot."""

from __future__ import annotations

import json
import tempfile
import unittest
from dataclasses import asdict
from fractions import Fraction as Q
from pathlib import Path

import d6_k7_rankone_tetrad_pilot as pilot
import d6_k7_rankone_tetrad_verify as independent
import d6_k7_rank_reference as prior
import run_d6_k7_rankone_tetrad_full as full
import verify_d6_k7_rankone_tetrad_full as full_verify


ROOT = Path(__file__).resolve().parent


class RankOneTetradTest(unittest.TestCase):
    def test_constant_tetrad_contradiction_is_certified(self) -> None:
        equations = [{(0,): Q(1)}]
        certificate = pilot.find_raw_certificate(
            equations, (), variables=1, output_degree=4
        )
        self.assertIsNotNone(certificate)
        self.assertEqual(
            pilot.verify_raw_certificate(equations, (), 1, certificate),
            {(0,): Q(1)},
        )

    def test_equation_with_positive_root_has_no_certificate(self) -> None:
        equations = [{(1,): Q(1), (0,): Q(-1)}]
        self.assertIsNone(pilot.find_raw_certificate(
            equations, (), variables=1, output_degree=4
        ))

    def test_triangle_sign_identity_is_checked(self) -> None:
        # E=1+x, t=x gives the exact sound form E=1+t for x>0.
        equations = [{(0,): Q(1), (1,): Q(1)}]
        signs = [{(1,): Q(1)}]
        certificate = {
            "equation_multipliers": [{
                "equation": 0, "monomial": [0], "coefficient": [1, 1],
            }],
            "sign_multipliers": [{
                "sign": 0, "monomial": [0], "coefficient": [1, 1],
            }],
            "positive_polynomial": [{
                "monomial": [0], "coefficient": [1, 1],
            }],
        }
        pilot.verify_raw_certificate(equations, signs, 1, certificate)
        tampered = json.loads(json.dumps(certificate))
        tampered["positive_polynomial"][0]["coefficient"][0] = 2
        with self.assertRaises(ValueError):
            pilot.verify_raw_certificate(equations, signs, 1, tampered)

    def test_real_near_saturating_system_matches_independent_rebuild(self) -> None:
        sample = json.loads(
            (ROOT / "d6_k7_rank_sample.json").read_text(encoding="utf-8")
        )
        graph = next(
            item for item in sample["graphs"] if item["index"] == 369959
        )
        adj = tuple(graph["adjacency"])
        support_solver = prior.SupportSolver()
        zero_forcing = prior.ZeroForcingSolver()
        clique_solver = prior.CliqueStructureSolver()
        selected = None
        for seed_mask in prior.clique_masks(adj, 7):
            _seed, outside, defects, ladj, eligible = prior.seed_instance(
                adj, seed_mask
            )
            total_term_rank = prior.matching_size(defects)
            for zmask in prior.eligible_covers(ladj, eligible):
                analysis = prior.analyze_cover(
                    adj, outside, defects, zmask,
                    support_solver, zero_forcing, total_term_rank,
                    clique_solver,
                )
                if analysis.enhanced_joint_failed:
                    continue
                nvertices = [
                    outside[index] for index in range(len(outside))
                    if not (zmask & (1 << index))
                ]
                graph_n = tuple(prior.induced_graph(adj, nvertices))
                if tuple(prior.clique_masks(
                    graph_n, analysis.k_rank_upper
                )):
                    continue
                near = tuple(prior.clique_masks(
                    graph_n, analysis.k_rank_upper - 1
                ))
                if near:
                    selected = (graph_n, near[0])
                    break
            if selected is not None:
                break
        self.assertIsNotNone(selected)
        graph_n, clique_mask = selected
        located = pilot.rank_one_system(graph_n, clique_mask)
        rebuilt = independent.reconstruct_rank_one_system(
            graph_n, clique_mask
        )
        self.assertEqual(located.clique, rebuilt[0])
        self.assertEqual(located.remainder, rebuilt[1])
        self.assertEqual(located.masks, rebuilt[2])
        self.assertEqual(located.pair_labels, rebuilt[3])
        self.assertEqual(located.tetrad_labels, rebuilt[4])
        self.assertEqual(located.tetrad_equations, rebuilt[5])
        self.assertEqual(located.triangle_labels, rebuilt[6])
        self.assertEqual(located.triangle_signs, rebuilt[7])

    def test_current_and_frozen_pilot_reports_verify(self) -> None:
        self.assertEqual(
            independent.verify_report(
                ROOT / "d6_k7_rankone_tetrad_sample.json",
                ROOT / "d6_k7_rankone_tetrad_pilot_report.json",
            ),
            {"graphs": 12, "certificates": 12},
        )
        self.assertEqual(
            independent.verify_report(
                ROOT / "d6_k7_rank_sample.json",
                ROOT / "d6_k7_rankone_tetrad_frozen_control_report.json",
            ),
            {"graphs": 9, "certificates": 20},
        )

    def test_full_selection_and_one_production_graph(self) -> None:
        selection_path = ROOT / "d6_k7_rankone_tetrad_full_selection.json"
        selection = json.loads(selection_path.read_text(encoding="utf-8"))
        self.assertEqual(len(selection["selected_indices"]), 12_839)
        selected, provenance, rows, certificates = full.load_inputs(
            ROOT / ".runs/d6_k7_rank_survivors.json",
            selection_path,
            full.sha256(selection_path),
        )
        self.assertEqual(len(selected), 12_839)
        self.assertEqual(provenance["selected"], 12_839)
        graph = next(item for item in selected if item["index"] == 1466007)
        result = full.analyze_graph(
            graph, rows[1466007], certificates.get(1466007, [])
        )
        self.assertTrue(result.marginal_tetrad_rejected)
        self.assertEqual(result.prior_dual_passing_covers, 1)
        self.assertEqual(result.tetrad_failed_covers, 1)
        serialized = asdict(result)
        decision = {
            name: str(serialized[name])
            for name in full.DECISION_FIELDS
            if name in serialized
        }
        decision.update({
            "ordinal": str(graph["ordinal"]),
            "index": str(graph["index"]),
            "status": "REJECTED",
            "error_type": "",
        })
        verified = full_verify.verify_graph((
            graph,
            rows[1466007],
            certificates.get(1466007, []),
            decision,
            result.tetrad_failure_witnesses,
        ))
        self.assertEqual(verified["status"], "PASS")
        self.assertEqual(verified["production_status"], "REJECTED")

    def test_checkpoint_roundtrip_is_config_pinned(self) -> None:
        selected = [{"ordinal": 0, "index": 7}]
        record = {
            "schema": 1, "ordinal": 0, "index": 7,
            "status": "SURVIVOR", "result": {},
        }
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            path = full.checkpoint_path(directory, 0, 7)
            full.atomic_json(path, {
                "schema": 1, "config_sha256": "abc", "record": record,
            })
            self.assertEqual(
                full.load_checkpoints(directory, "abc", selected),
                {0: record},
            )
            with self.assertRaises(ValueError):
                full.load_checkpoints(directory, "changed", selected)


if __name__ == "__main__":
    unittest.main()
