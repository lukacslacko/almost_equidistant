#!/usr/bin/env python3
"""Independent controls for the exact graph-3949382 Schur probe."""

from __future__ import annotations

import inspect
import unittest

import sympy as sp

import probe_d6_k7_schur_3949382 as probe


class K7Schur3949382Controls(unittest.TestCase):
    def setUp(self) -> None:
        self.boundary = probe.load_boundary()
        self.audit = probe.quantifier_audit(self.boundary)
        self.pattern = probe.reconstruct_exact_pattern(self.boundary, self.audit)

    def test_full_seed_cover_quantifier_is_frozen(self) -> None:
        self.assertEqual(len(self.audit["raw_covers"]), 128)
        self.assertEqual(
            probe.stable_hash(list(self.audit["raw_covers"])),
            probe.EXPECTED_RAW_COVERS_SHA256,
        )
        self.assertEqual(self.audit["current_covers"], (0, 2048, 3072))
        self.assertEqual(self.audit["prior_layer_eliminated_raw_covers"], 125)
        self.assertEqual(self.audit["current_infeasible_covers"], (2048, 3072))
        self.assertEqual(self.audit["sole_new_branch"], 0)

    def test_every_zero_entry_uses_support_disjointness_not_a_bare_nonedge(self) -> None:
        graph_n = self.pattern["graph_n"]
        propagated = self.audit["propagated_masks"]
        for row, vertex in enumerate(probe.REMAINDER):
            for coordinate, basis_vertex in enumerate(probe.CLIQUE):
                target = (self.pattern["basis_masks"][row] >> coordinate) & 1
                if target == 0:
                    self.assertFalse(propagated[vertex] & propagated[basis_vertex])
                    self.assertFalse(graph_n[vertex] & (1 << basis_vertex))
        for pair, target in zip(
            probe.PAIR_ORDER, self.pattern["pair_targets"], strict=True
        ):
            first = probe.REMAINDER[pair[0]]
            second = probe.REMAINDER[pair[1]]
            if target == 0:
                self.assertFalse(propagated[first] & propagated[second])
                self.assertFalse(graph_n[first] & (1 << second))
            else:
                self.assertTrue(graph_n[first] & (1 << second))

    def test_exact_pattern_is_reconstructed_not_archived(self) -> None:
        self.assertEqual(self.pattern["basis_masks"], (83, 47, 98, 86, 90))
        self.assertEqual(
            "".join(map(str, self.pattern["pair_targets"])), "0111011011"
        )
        self.assertEqual(
            [self.audit["outside"][index] for index in probe.CLIQUE],
            [0, 3, 9, 10, 13, 15, 17],
        )

    def test_alternative_exact_elimination_certificate(self) -> None:
        """Cross-check by a second elimination, distinct from the probe's."""

        t, b, e, f, g = sp.symbols("t b e f g")
        total = 1 + 3 * t + b + e + f + g
        p = t + b + e + g
        q = 3 * t + b + f
        r = b + f + g
        equations = (
            sp.expand((t + b) * total - p * q),
            sp.expand((p - t - 1) * total - p**2),
            sp.expand((b + g - 1) * total - p * r),
            sp.expand((b + f - 1) * total - q * r),
        )
        e0, e1, e2, e3 = equations
        a0 = b * (f + t) - f - t**2 - 2 * t
        b0 = b * f + b * t + b - f * t - 2 * f - g * t - 4 * t - 2
        self.assertEqual(
            sp.expand((f + 2 * t + 1) * a0 - (f + t) * e0 - (f + 2 * t) * e1),
            0,
        )
        self.assertEqual(
            sp.expand(
                (b + t - 1) * b0
                - (b - 2) * (e0 + e1)
                - t * (e2 + e3)
            ),
            0,
        )

        b_solution = (f + t**2 + 2 * t) / (f + t)
        g_solution = (
            -f**2 * t - f**2 - 3 * f * t - f + t**3 - t**2
        ) / (t * (f + t))
        e_solution = (f * t + f - t**2 + 2 * t + 1) / t
        substitution = {b: b_solution, e: e_solution, g: g_solution}
        self.assertEqual(sp.factor(e0.subs(substitution)), 0)
        self.assertEqual(sp.factor(e1.subs(substitution)), 0)
        final = t**2 - t * f - t - f
        expected = (f + 2 * t + 1) ** 2 * final / (t * (f + t))
        self.assertEqual(sp.factor(e2.subs(substitution) - expected), 0)
        self.assertEqual(sp.factor(e3.subs(substitution) + expected), 0)

        f_solution = t * (t - 1) / (t + 1)
        final_g = sp.factor(g_solution.subs(f, f_solution))
        self.assertEqual(
            sp.factor(final_g + (t - 1) * (t + 1) / (2 * t**2)), 0
        )
        # In the strict domain f>0 and t>0, the displayed f formula forces
        # t>1; the displayed g is then strictly negative.

    def test_boundary_point_catches_strict_positivity(self) -> None:
        variables, equations = probe.schur_pair_equations(
            self.pattern["basis_masks"], self.pattern["pair_targets"]
        )
        point = dict(zip(variables, (1, 3, 1, 1, 2, 0, 0), strict=True))
        self.assertTrue(all(sp.expand(value.subs(point)) == 0 for value in equations.values()))
        self.assertFalse(all(point[variable] > 0 for variable in variables))

    def test_full_probe_concludes_exact_but_not_yet_production_union(self) -> None:
        report = probe.build_report()
        self.assertEqual(report["status"], "EXACT_CONTRADICTION")
        self.assertEqual(
            report["theorem_credit"], "EXPLORATORY_NOT_YET_PRODUCTION_UNION"
        )
        self.assertFalse(report["semantics"]["floating_point_enters_contradiction"])
        self.assertTrue(report["semantics"]["candidate_nonedges_optional"])

    def test_probe_imports_no_k7_production_module(self) -> None:
        source = inspect.getsource(probe)
        self.assertNotIn("import d6_k7_", source)
        self.assertNotIn("from d6_k7_", source)


if __name__ == "__main__":
    unittest.main()
