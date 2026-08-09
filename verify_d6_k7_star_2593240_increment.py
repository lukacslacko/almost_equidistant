#!/usr/bin/env python3
"""Independent checker for the exact graph-2593240 K7 star increment.

This checker imports neither the production builder nor the exploratory probe.
It reconstructs both K7 seed/cover quantifiers and every exact Schur entry from
the pinned artifacts.  Its algebra uses a one-stage lexicographic elimination,
different from the producer's two-stage grevlex membership route.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import subprocess
import sys
import time
from datetime import UTC, datetime
from itertools import combinations
from pathlib import Path
from typing import Iterable, Sequence

import sympy as sp


ROOT = Path(__file__).resolve().parent
TARGET_INDEX = 2_593_240
SEEDS = (
    (2, 4, 7, 10, 13, 15, 18),
    (6, 7, 8, 9, 12, 15, 18),
)
SEED_MASKS = tuple(sum(1 << vertex for vertex in seed) for seed in SEEDS)
CLIQUES = (
    (0, 1, 2, 7, 9, 10),
    (0, 1, 3, 7, 9, 10),
)
REMAINDERS = tuple(
    tuple(vertex for vertex in range(12) if vertex not in set(clique))
    for clique in CLIQUES
)
EXPECTED_PROPAGATED = (
    (100, 127, 5, 98, 25, 19, 26, 108, 11, 20, 38, 72),
    (127, 98, 21, 10, 28, 97, 25, 102, 13, 18, 35, 68),
)
EXPECTED_COVERS = ((0, 8), (0, 32))
EXPECTED_RAW_COVERS_SHA256 = (
    "fb44017cc67555b5954484a5bd012ba2c8b4ff80e5ad878eb517a14ac791d769",
    "f915ae089f9885b1d800c698306a8180f57ebe9d8818481094ca07996930a6bd",
)
EXPECTED_BASIS_MASKS = (
    (43, 30, 54, 58, 46, 11),
    (57, 29, 43, 53, 45, 11),
)
PAIR_ORDER = tuple(combinations(range(6), 2))
EXPECTED_PAIR_TARGETS = (
    tuple(map(int, "011111111110111")),
    tuple(map(int, "111110111111101")),
)
REMAINDER_ISOMORPHISM = (2, 1, 3, 0, 4, 5)
COORDINATE_ISOMORPHISM = (1, 0, 2, 3, 4, 5)
EXPECTED_ADJACENCY_SHA256 = (
    "56f6a6a8bacf5a67e8fd0a0688b86b610d2d8c0842060b6423f8d484f340ed82"
)
UPSTREAM = {
    "d6_current_residue_manifest_v6.json": (
        "c513bcc40c037af23340dce3cc91d593a101dd1d073f93feab82abd556efadc9"
    ),
    "d6_current_residue_manifest_v6_verification.json": (
        "3d11cbd20cb59b42a8c8fd5252fb3ef3dea649576982ea8017622fad3279a7bc"
    ),
    "d6_k7_one_two_star_increment_report.json": (
        "12b3d18b1ea81961f58d831d4c7c322fbceb2300e7533b64161e8011fbe9d1ec"
    ),
    "d6_k7_one_two_star_increment_verification.json": (
        "c6f40578685036cb1156cb0ea8bf06d6004a70d5f376aff5ac1e16a1c0ca09eb"
    ),
    "d6_k7_schur_3949382_increment_report.json": (
        "30efa40d1db8a3aef03af9e56aa915a04594aabe8b72e5a1c72489e2a8c952b8"
    ),
    "d6_k7_schur_3949382_increment_verification.json": (
        "18c9980baa74596523eba9b447af1ae5830a9558d80981e1fa3ff6c984bcb5bc"
    ),
}
FROZEN_EXPLORATORY = {
    "probe_d6_k7_2593240_star_boundary.py": (
        "7d88aafdda609a85255476cb81d72ee0685880c8b17e0823f3ffc4c690ad912f"
    ),
    "test_probe_d6_k7_2593240_star_boundary.py": (
        "30e5949baa4180edc6578d25054fd175c0bd81b248470158bb39e8bcc7b6c2a1"
    ),
    "d6_k7_2593240_star_boundary.md": (
        "cbaa5e0bcb1f313136fff96c17cb33b251465232f062df7efa1fa692c8b93d2e"
    ),
}
SOURCE_FILES = {
    *FROZEN_EXPLORATORY,
    "build_d6_k7_star_2593240_increment.py",
    "verify_d6_k7_star_2593240_increment.py",
    "test_d6_k7_star_2593240_increment.py",
    "d6_k7_star_2593240_increment.md",
}
EXPECTED_INPUT = (
    316173,
    2581209,
    2593240,
    3595554,
    3648882,
    3729907,
    3935560,
    3936177,
    3936310,
    3936435,
    3945490,
    3945555,
    3945557,
    3945564,
    3947605,
)
EXPECTED_INPUT_SHA256 = (
    "715c52011421fb6d8721d52339ed3437ef1cec344ddf3f03a3b0b618f35b00d7"
)
EXPECTED_OUTPUT_SHA256 = (
    "3f038f15c97483aaf529300a86da00459d04e4a8f49823b8849fa7a4a91dd7cb"
)
EXPECTED_SEMANTICS = {
    "all_actual_support_branches_for_both_seeds_are_quantified": True,
    "both_required_K7_seeds_checked": True,
    "candidate_nonedges_optional": True,
    "floating_point_enters_rejection": False,
    "one_infeasible_required_K7_seed_rejects_graph": True,
    "propagated_masks_are_support_supersets": True,
    "rank_one_schur_support_closure_checked": True,
    "required_edges_give_unit_K_entries": True,
    "strict_positive_diagonal_parameters_encode_distinct_points": True,
    "zero_K_entries_require_disjoint_propagated_support_supersets": True,
}
EXPECTED_STAR_SEMANTICS = {
    "allowed_unpinned_coordinates_may_be_zero": True,
    "candidate_nonedges_optional": True,
    "cap500000_campaign_used": False,
    "floating_point_enters_rejection": False,
    "graph_rejected_if_any_k7_seed_is_infeasible": True,
    "normalization_factor_t_nonzero_on_cover_complement_N": True,
    "only_required_edges_enter_star_equations": True,
    "positive_sqrt7_embedding_checked_exactly": True,
    "propagated_masks_are_support_supersets": True,
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_hash(value: object) -> str:
    raw = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")
    return hashlib.sha256(raw).hexdigest()


def atomic_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def bits(mask: int) -> tuple[int, ...]:
    return tuple(index for index in range(mask.bit_length()) if mask & (1 << index))


def validate_adjacency(adjacency: Sequence[int]) -> None:
    require(len(adjacency) == 19, "graph order")
    full = (1 << len(adjacency)) - 1
    for vertex, row in enumerate(adjacency):
        require(not row & ~full, "adjacency outside graph")
        require(not row & (1 << vertex), "adjacency loop")
        for other in range(vertex):
            require(
                bool(row & (1 << other))
                == bool(adjacency[other] & (1 << vertex)),
                "asymmetric adjacency",
            )


def is_clique(adjacency: Sequence[int], vertices: Sequence[int]) -> bool:
    return all(
        bool(adjacency[first] & (1 << second))
        for first, second in combinations(vertices, 2)
    )


def induced_graph(adjacency: Sequence[int], vertices: Sequence[int]) -> tuple[int, ...]:
    position = {vertex: index for index, vertex in enumerate(vertices)}
    return tuple(
        sum(
            1 << position[other]
            for other in vertices
            if adjacency[vertex] & (1 << other)
        )
        for vertex in vertices
    )


def reconstruct_seed(
    adjacency: Sequence[int], seed: Sequence[int]
) -> tuple[tuple[int, ...], tuple[int, ...], tuple[int, ...], int]:
    require(is_clique(adjacency, seed), "required K7 seed")
    seed_set = set(seed)
    seed_position = {vertex: position for position, vertex in enumerate(seed)}
    outside = tuple(vertex for vertex in range(19) if vertex not in seed_set)
    defects = tuple(
        sum(
            1 << seed_position[seed_vertex]
            for seed_vertex in seed
            if not adjacency[vertex] & (1 << seed_vertex)
        )
        for vertex in outside
    )
    ladj = [0] * len(outside)
    for first, second in combinations(range(len(outside)), 2):
        if (
            adjacency[outside[first]] & (1 << outside[second])
            and not defects[first] & defects[second]
        ):
            ladj[first] |= 1 << second
            ladj[second] |= 1 << first
    eligible = sum(
        1 << index for index, defect in enumerate(defects) if defect.bit_count() >= 3
    )
    return outside, defects, tuple(ladj), eligible


def is_vertex_cover(ladj: Sequence[int], zmask: int) -> bool:
    remaining = ((1 << len(ladj)) - 1) & ~zmask
    return all(not ladj[vertex] & remaining for vertex in bits(remaining))


def raw_covers(ladj: Sequence[int], eligible: int) -> tuple[int, ...]:
    return tuple(
        zmask
        for zmask in range(1 << len(ladj))
        if not zmask & ~eligible
        and zmask.bit_count() <= 7
        and is_vertex_cover(ladj, zmask)
    )


def load_upstream() -> tuple[dict, dict, dict, dict, dict]:
    for name, expected in UPSTREAM.items():
        require(sha256(ROOT / name) == expected, f"upstream hash: {name}")
    manifest = json.loads(
        (ROOT / "d6_current_residue_manifest_v6.json").read_text(encoding="utf-8")
    )
    manifest_check = json.loads(
        (ROOT / "d6_current_residue_manifest_v6_verification.json").read_text(
            encoding="utf-8"
        )
    )
    star = json.loads(
        (ROOT / "d6_k7_one_two_star_increment_report.json").read_text(
            encoding="utf-8"
        )
    )
    star_check = json.loads(
        (ROOT / "d6_k7_one_two_star_increment_verification.json").read_text(
            encoding="utf-8"
        )
    )
    preceding = json.loads(
        (ROOT / "d6_k7_schur_3949382_increment_report.json").read_text(
            encoding="utf-8"
        )
    )
    preceding_check = json.loads(
        (ROOT / "d6_k7_schur_3949382_increment_verification.json").read_text(
            encoding="utf-8"
        )
    )
    require(
        manifest.get("schema") == "d6-current-certified-residue-v6"
        and manifest.get("status") == "COMPLETE_MIXED_CERTIFICATE_UNION",
        "manifest boundary",
    )
    require(
        manifest_check.get("status") == "PASS"
        and manifest_check.get("manifest", {}).get("sha256")
        == UPSTREAM["d6_current_residue_manifest_v6.json"],
        "manifest checker boundary",
    )
    require(
        star.get("kind")
        == "d6_k7_one_free_neighbour_two_free_center_increment"
        and star.get("status") == "COMPLETE"
        and star.get("semantics") == EXPECTED_STAR_SEMANTICS,
        "star report boundary",
    )
    require(
        star_check.get("kind") == "d6_k7_one_two_star_increment_verification"
        and star_check.get("status") == "PASS"
        and star_check.get("report", {}).get("sha256")
        == UPSTREAM["d6_k7_one_two_star_increment_report.json"]
        and star_check.get("semantics") == EXPECTED_STAR_SEMANTICS
        and star_check.get("checked", {}).get("eligible_covers") == 19_932
        and star_check.get("checked", {}).get("current_passing_families") == 88,
        "star checker boundary",
    )
    require(
        preceding.get("kind") == "d6_k7_schur_3949382_exact_increment"
        and preceding.get("status") == "COMPLETE_EXACT_REJECTION"
        and preceding.get("summary", {}).get("ordered_survivor_indices")
        == list(EXPECTED_INPUT)
        and preceding.get("summary", {}).get("ordered_survivor_indices_sha256")
        == EXPECTED_INPUT_SHA256,
        "preceding exact increment",
    )
    require(
        preceding_check.get("kind")
        == "d6_k7_schur_3949382_increment_verification"
        and preceding_check.get("status") == "PASS"
        and preceding_check.get("report", {}).get("sha256")
        == UPSTREAM["d6_k7_schur_3949382_increment_report.json"]
        and preceding_check.get("conclusion", {}).get(
            "ordered_survivor_indices_sha256"
        )
        == EXPECTED_INPUT_SHA256,
        "preceding independent verification",
    )
    return manifest, star, star_check, preceding, preceding_check


def validate_report_schema(report: dict) -> None:
    output = [index for index in EXPECTED_INPUT if index != TARGET_INDEX]
    require(report.get("schema") == 1, "report schema")
    require(
        report.get("kind") == "d6_k7_rankone_star_2593240_exact_increment",
        "report kind",
    )
    require(report.get("status") == "COMPLETE_EXACT_REJECTION", "report status")
    require(report.get("semantics") == EXPECTED_SEMANTICS, "report semantics")
    require(report.get("upstream_sha256") == dict(sorted(UPSTREAM.items())), "upstream roots")
    require(
        report.get("frozen_exploratory_sha256")
        == dict(sorted(FROZEN_EXPLORATORY.items())),
        "exploratory roots",
    )
    require(
        report.get("input", {}).get("artifact")
        == "d6_k7_schur_3949382_increment_report.json"
        and report.get("input", {}).get("ordered_indices") == list(EXPECTED_INPUT)
        and report.get("input", {}).get("ordered_indices_sha256")
        == EXPECTED_INPUT_SHA256,
        "input residue",
    )
    summary = report.get("summary", {})
    require(
        summary.get("graphs_rejected") == 1
        and summary.get("rejected_indices") == [TARGET_INDEX]
        and summary.get("rejected_indices_sha256") == stable_hash([TARGET_INDEX])
        and summary.get("graphs_surviving") == 14
        and summary.get("ordered_survivor_indices") == output
        and summary.get("ordered_survivor_indices_sha256") == EXPECTED_OUTPUT_SHA256,
        "union-ready summary",
    )
    require(
        report.get("certificate_sha256") == stable_hash(report.get("certificate")),
        "certificate hash",
    )


def validate_source_boundary(report: dict) -> dict:
    sources = report.get("source_sha256")
    require(isinstance(sources, dict) and set(sources) == SOURCE_FILES, "source set")
    for name, expected in FROZEN_EXPLORATORY.items():
        require(sources.get(name) == expected, f"frozen source in report: {name}")
    launch = report.get("execution", {}).get("git", {})
    require(
        launch.get("branch") == "codex/dimension6"
        and launch.get("proof_and_checker_sources_equal_committed_blobs") is True,
        "launch branch/source claim",
    )
    commit = launch.get("commit")
    require(isinstance(commit, str) and len(commit) == 40, "launch commit")
    for name, expected in sources.items():
        require(sha256(ROOT / name) == expected, f"working source hash: {name}")
        blob = subprocess.run(
            ["git", "show", f"{commit}:{name}"],
            cwd=ROOT,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        ).stdout
        require(hashlib.sha256(blob).hexdigest() == expected, f"committed source: {name}")
    porcelain = "\n".join(launch.get("porcelain_lines", ()))
    require(
        hashlib.sha256(porcelain.encode("utf-8")).hexdigest()
        == launch.get("porcelain_sha256"),
        "launch porcelain binding",
    )
    return {
        "commit": commit,
        "branch": launch["branch"],
        "source_files": len(sources),
        "proof_and_checker_sources_equal_committed_blobs": True,
    }


def mapped_mask(mask: int) -> int:
    return sum(
        1 << COORDINATE_ISOMORPHISM[coordinate] for coordinate in bits(mask)
    )


def independent_quantifier_and_patterns(manifest: dict, star: dict) -> dict:
    graph = next(
        record
        for record in manifest["classes"]["K7"]["graphs"]
        if int(record["index"]) == TARGET_INDEX
    )
    adjacency = tuple(map(int, graph["adjacency"]))
    validate_adjacency(adjacency)
    require(stable_hash(list(adjacency)) == EXPECTED_ADJACENCY_SHA256, "adjacency hash")
    k7s = tuple(
        chosen
        for chosen in combinations(range(19), 7)
        if is_clique(adjacency, chosen)
    )
    require(k7s == SEEDS, "complete required-K7 list")

    graph_record = next(
        record for record in star["records"] if int(record["index"]) == TARGET_INDEX
    )
    require(graph_record.get("decision") == "SURVIVOR", "upstream target status")
    archived_by_seed = {
        int(record["seed_mask"]): record for record in graph_record["seeds"]
    }
    quantifier = []
    patterns = []
    reason_counts = {
        "required_unit_edges": 0,
        "disjoint_support_zeros": 0,
        "unresolved_optional_entries": 0,
    }
    for position, seed in enumerate(SEEDS):
        outside, _defects, ladj, eligible = reconstruct_seed(adjacency, seed)
        covers = raw_covers(ladj, eligible)
        require(len(covers) == 502, f"raw cover count at seed {position}")
        require(
            stable_hash(list(covers)) == EXPECTED_RAW_COVERS_SHA256[position],
            f"raw cover hash at seed {position}",
        )
        seed_record = archived_by_seed[SEED_MASKS[position]]
        require(tuple(seed_record["seed"]) == seed, f"archived seed {position}")
        archived = seed_record["current_covers"]
        require(
            tuple(int(row["zmask"]) for row in archived) == EXPECTED_COVERS[position],
            f"current covers at seed {position}",
        )
        by_z = {int(row["zmask"]): row for row in archived}
        nonzero = EXPECTED_COVERS[position][1]
        require(
            by_z[nonzero]["status"] == "INFEASIBLE"
            and int(by_z[nonzero]["current_passing_families"]) == 0,
            f"infeasible current cover at seed {position}",
        )
        zero = by_z[0]
        require(
            zero["status"] == "PASSING"
            and int(zero["current_passing_families"]) == 1
            and int(zero["star_passing_families"]) == 1,
            f"sole current family at seed {position}",
        )
        witness = zero["first_passing_witness"]
        require(witness.get("z_supports") == [], f"z=0 witness at seed {position}")
        propagated = tuple(map(int, witness["propagated_masks"]))
        require(
            propagated == EXPECTED_PROPAGATED[position],
            f"propagated supports at seed {position}",
        )

        graph_n = induced_graph(adjacency, outside)
        clique = CLIQUES[position]
        remainder = REMAINDERS[position]
        require(is_clique(graph_n, clique), f"required K6 at seed {position}")

        def exact_entry(first: int, second: int) -> int:
            if graph_n[first] & (1 << second):
                require(
                    propagated[first] & propagated[second],
                    "required edge with disjoint supports",
                )
                reason_counts["required_unit_edges"] += 1
                return 1
            require(
                not propagated[first] & propagated[second],
                f"optional unresolved entry {position}:{first},{second}",
            )
            reason_counts["disjoint_support_zeros"] += 1
            return 0

        basis_masks = []
        for vertex in remainder:
            mask = 0
            for coordinate, basis_vertex in enumerate(clique):
                mask |= exact_entry(vertex, basis_vertex) << coordinate
            basis_masks.append(mask)
        pair_targets = tuple(
            exact_entry(remainder[first], remainder[second])
            for first, second in PAIR_ORDER
        )
        require(
            tuple(basis_masks) == EXPECTED_BASIS_MASKS[position],
            f"basis masks at seed {position}",
        )
        require(
            pair_targets == EXPECTED_PAIR_TARGETS[position],
            f"pair targets at seed {position}",
        )
        quantifier.append(
            {
                "seed": list(seed),
                "outside": list(outside),
                "raw_eligible_covers": len(covers),
                "raw_eligible_covers_sha256": stable_hash(list(covers)),
                "current_covers": list(EXPECTED_COVERS[position]),
                "prior_layer_eliminated_raw_covers": len(covers) - len(archived),
                "current_infeasible_covers": [nonzero],
                "sole_new_branch": 0,
                "sole_new_branch_passing_families": 1,
            }
        )
        patterns.append(
            {
                "local_rank6_clique": list(clique),
                "global_rank6_clique": [outside[index] for index in clique],
                "local_remainder": list(remainder),
                "global_remainder": [outside[index] for index in remainder],
                "propagated_support_supersets": list(propagated),
                "basis_masks": basis_masks,
                "pair_order": [list(pair) for pair in PAIR_ORDER],
                "pair_targets": list(pair_targets),
                "pair_target_string": "".join(map(str, pair_targets)),
                "zero_entry_justification": (
                    "Every zero in the displayed basis/pair pattern has disjoint "
                    "propagated support supersets; it is not inferred from a "
                    "candidate nonedge alone."
                ),
            }
        )

    for old_remainder, new_remainder in enumerate(REMAINDER_ISOMORPHISM):
        require(
            mapped_mask(patterns[0]["basis_masks"][old_remainder])
            == patterns[1]["basis_masks"][new_remainder],
            "basis-mask isomorphism",
        )
    first_targets = dict(zip(PAIR_ORDER, patterns[0]["pair_targets"], strict=True))
    second_targets = dict(zip(PAIR_ORDER, patterns[1]["pair_targets"], strict=True))
    for (old_first, old_second), target in first_targets.items():
        new_pair = tuple(
            sorted(
                (
                    REMAINDER_ISOMORPHISM[old_first],
                    REMAINDER_ISOMORPHISM[old_second],
                )
            )
        )
        require(second_targets[new_pair] == target, "remainder-pair isomorphism")

    return {
        "target": {
            "index": TARGET_INDEX,
            "adjacency_sha256": EXPECTED_ADJACENCY_SHA256,
            "required_k7_seeds": [list(seed) for seed in SEEDS],
            "required_k7_seed_masks": list(SEED_MASKS),
        },
        "quantifier": quantifier,
        "patterns": patterns,
        "system_isomorphism": {
            "first_remainder_to_second": list(REMAINDER_ISOMORPHISM),
            "first_coordinate_to_second": list(COORDINATE_ISOMORPHISM),
            "checked_exactly": True,
        },
        "entry_reasons": reason_counts,
    }


def mask_sum(mask: int, variables: Sequence[sp.Symbol]) -> sp.Expr:
    return sum((variables[index] for index in bits(mask)), sp.Integer(0))


def polynomial_zero(expression: sp.Expr, variables: Iterable[sp.Symbol], label: str) -> None:
    require(
        sp.Poly(sp.together(expression), *variables, domain=sp.QQ).is_zero,
        label,
    )


def independent_algebra() -> dict:
    """One-stage lex elimination, independent of the archived grevlex route."""

    variables = sp.symbols("a b c d e f")
    a, b, c, d, e, f = variables
    total = 1 + sum(variables)
    masks = EXPECTED_BASIS_MASKS[0]
    sums = [mask_sum(mask, variables) for mask in masks]
    numerators: dict[tuple[int, int], sp.Expr] = {}
    for pair, target in zip(PAIR_ORDER, EXPECTED_PAIR_TARGETS[0], strict=True):
        first, second = pair
        intersection = mask_sum(masks[first] & masks[second], variables)
        numerators[pair] = sp.expand(
            target * total - total * intersection + sums[first] * sums[second]
        )
    tetrads = []
    tetrad_by_label = {}
    for vertices in combinations(range(6), 4):
        first, second, third, fourth = vertices
        common = numerators[(first, second)] * numerators[(third, fourth)]
        first_equation = sp.expand(
            common - numerators[(first, third)] * numerators[(second, fourth)]
        )
        second_equation = sp.expand(
            common - numerators[(first, fourth)] * numerators[(second, third)]
        )
        tetrad_by_label[(vertices, 0)] = first_equation
        tetrad_by_label[(vertices, 1)] = second_equation
        tetrads.extend((first_equation, second_equation))

    factor = -f * (a + 1) * (c - e) * total
    polynomial_zero(
        tetrad_by_label[((0, 3, 4, 5), 0)] - factor,
        variables,
        "e=c factor identity",
    )
    # Strict positivity gives e=c.  Divide the remaining tetrads by the
    # positive common factor T and compute a direct lex triangularization.
    substituted_total = total.subs(e, c)
    reduced = tuple(
        sp.cancel(equation.subs(e, c) / substituted_total) for equation in tetrads
    )
    lex_variables = (a, b, d, f, c)
    basis = sp.groebner(
        reduced, *lex_variables, order="lex", domain=sp.QQ
    )
    targets = (
        (c - 1) ** 2 * (c + 1) * (2 * c + 1) * (3 * c + 1),
        (c - f) * (2 * c + 1),
        6 * c**4 + 5 * c**3 - 5 * c**2 + 7 * c - 12 * d - 1,
        36 * b + 102 * c**4 + 13 * c**3 - 217 * c**2 - 85 * c + 36 * f - 29,
        108 * a + 102 * c**4 + 13 * c**3 - 217 * c**2 - 49 * c + 79,
    )
    for index, target in enumerate(targets):
        _quotients, remainder = basis.reduce(target)
        polynomial_zero(remainder, lex_variables, f"lex target {index}")

    point = {
        a: sp.Rational(2, 3),
        b: 5,
        c: 1,
        d: 1,
        e: 1,
        f: 1,
    }
    require(
        all(sp.expand(equation.subs(point)) == 0 for equation in tetrads),
        "positive tetrad control",
    )
    evaluated = {
        pair: sp.factor(expression.subs(point))
        for pair, expression in numerators.items()
    }
    expected_nonzero = {
        (0, 1): sp.Rational(-8, 3),
        (0, 2): sp.Integer(8),
        (0, 3): sp.Rational(-8, 3),
        (0, 4): sp.Rational(-8, 3),
        (0, 5): sp.Rational(-28, 3),
    }
    require(
        {pair: value for pair, value in evaluated.items() if value}
        == expected_nonzero,
        "star off-diagonal support",
    )
    require(
        evaluated[(0, 1)] != 0
        and evaluated[(0, 2)] != 0
        and evaluated[(1, 2)] == 0,
        "rank-one support closure",
    )
    return {
        "method": "independent_one_stage_lex_elimination_and_support_closure",
        "variables": [str(variable) for variable in variables],
        "pair_numerators": [str(numerators[pair]) for pair in PAIR_ORDER],
        "tetrads_reconstructed": len(tetrads),
        "first_factor_identity": str(sp.factor(factor)),
        "lex_variable_order": [str(variable) for variable in lex_variables],
        "lex_groebner_basis": [
            str(sp.factor(polynomial.as_expr())) for polynomial in basis.polys
        ],
        "triangular_targets": [str(sp.factor(target)) for target in targets],
        "strict_sign_chain": [
            "-f*(a+1)*T*(c-e)=0 and strict positivity give e=c",
            "the lex univariate factor and c>0 give c=1",
            "(c-f)*(2*c+1)=0 gives f=1",
            "the remaining triangular relations give d=1, b=5, a=2/3",
            "g_01 and g_02 are nonzero while g_12=0",
            "a rank-one factor h*h^T cannot have that off-diagonal support",
        ],
        "unique_positive_tetrad_point": ["2/3", "5", "1", "1", "1", "1"],
        "star_nonzero_entries": {
            f"{first},{second}": str(value)
            for (first, second), value in expected_nonzero.items()
        },
        "status": "EXACT_RANK_ONE_SUPPORT_CONTRADICTION",
    }


def validate_report_certificate(report: dict, rebuilt: dict, algebra: dict) -> None:
    certificate = report.get("certificate", {})
    require(certificate.get("target") == rebuilt["target"], "target certificate")
    require(
        certificate.get("quantifier") == rebuilt["quantifier"],
        "both-seed quantifier certificate",
    )
    require(
        certificate.get("exact_patterns") == rebuilt["patterns"],
        "both exact patterns",
    )
    require(
        certificate.get("system_isomorphism") == rebuilt["system_isomorphism"],
        "system isomorphism",
    )
    archived_algebra = certificate.get("algebra", {})
    require(
        archived_algebra.get("pair_numerators") == algebra["pair_numerators"],
        "archived pair numerators",
    )
    require(
        archived_algebra.get("tetrads_checked") == algebra["tetrads_reconstructed"],
        "archived tetrad count",
    )
    require(
        archived_algebra.get("unique_positive_tetrad_point")
        == algebra["unique_positive_tetrad_point"],
        "archived positive point",
    )
    require(
        archived_algebra.get("star_nonzero_entries")
        == algebra["star_nonzero_entries"],
        "archived star entries",
    )
    require(
        certificate.get("rank_argument")
        == {
            "normalized_gram_rank_upper_bound": 7,
            "required_clique_order": 6,
            "clique_principal_block_positive_definite": True,
            "schur_complement_rank_upper_bound": 1,
            "schur_complement_positive_semidefinite": True,
            "strict_domain": "a,b,c,d,e,f > 0",
        },
        "rank argument",
    )


def verify_report(
    report_path: Path,
    expected_report_sha256: str,
    *,
    enforce_source_boundary: bool = True,
) -> dict:
    report_path = report_path.resolve()
    report_hash = sha256(report_path)
    require(report_hash == expected_report_sha256, "report hash")
    report = json.loads(report_path.read_text(encoding="utf-8"))
    validate_report_schema(report)
    provenance = (
        validate_source_boundary(report)
        if enforce_source_boundary
        else {"status": "SKIPPED_FOR_TEST"}
    )
    manifest, star, star_check, preceding, preceding_check = load_upstream()
    require(
        tuple(map(int, preceding["summary"]["ordered_survivor_indices"]))
        == EXPECTED_INPUT,
        "independently loaded input list",
    )
    rebuilt = independent_quantifier_and_patterns(manifest, star)
    algebra = independent_algebra()
    validate_report_certificate(report, rebuilt, algebra)
    return {
        "schema": 1,
        "kind": "d6_k7_star_2593240_increment_verification",
        "status": "PASS",
        "report": {"path": report_path.name, "sha256": report_hash},
        "semantics": EXPECTED_SEMANTICS,
        "source_boundary": provenance,
        "upstream_sha256": dict(sorted(UPSTREAM.items())),
        "checked": {
            "target_index": TARGET_INDEX,
            "required_k7_seeds": len(SEEDS),
            "raw_eligible_covers": sum(
                row["raw_eligible_covers"] for row in rebuilt["quantifier"]
            ),
            "prior_eliminated_raw_covers": sum(
                row["prior_layer_eliminated_raw_covers"]
                for row in rebuilt["quantifier"]
            ),
            "upstream_current_covers": sum(
                len(row["current_covers"]) for row in rebuilt["quantifier"]
            ),
            "upstream_infeasible_covers": sum(
                len(row["current_infeasible_covers"])
                for row in rebuilt["quantifier"]
            ),
            "new_families": sum(
                row["sole_new_branch_passing_families"]
                for row in rebuilt["quantifier"]
            ),
            "basis_entries": len(SEEDS) * 6 * 6,
            "remainder_pair_entries": len(SEEDS) * len(PAIR_ORDER),
            "unresolved_optional_entries": rebuilt["entry_reasons"][
                "unresolved_optional_entries"
            ],
            "upstream_independently_checked_eligible_covers": star_check["checked"][
                "eligible_covers"
            ],
            "upstream_independently_checked_current_passing_families": star_check[
                "checked"
            ]["current_passing_families"],
            "preceding_increment_verification_status": preceding_check["status"],
            "input_graphs": len(EXPECTED_INPUT),
            "output_graphs": len(EXPECTED_INPUT) - 1,
        },
        "independent_algebra": algebra,
        "conclusion": {
            "rejected_indices": [TARGET_INDEX],
            "ordered_survivor_indices_sha256": EXPECTED_OUTPUT_SHA256,
            "floating_point_used": False,
        },
        "execution": {
            "command": " ".join([sys.executable, *sys.argv]),
            "finished_utc": datetime.now(UTC).isoformat(),
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python_version": sys.version,
            "verifier_sha256": sha256(Path(__file__).resolve()),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--report",
        type=Path,
        default=ROOT / "d6_k7_star_2593240_increment_report.json",
    )
    parser.add_argument("--report-sha256", required=True)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "d6_k7_star_2593240_increment_verification.json",
    )
    args = parser.parse_args()
    result = verify_report(args.report, args.report_sha256)
    atomic_json(args.output.resolve(), result)
    print(
        json.dumps(
            {
                "output": str(args.output.resolve()),
                "sha256": sha256(args.output.resolve()),
                "status": result["status"],
                "report_sha256": result["report"]["sha256"],
                "rejected_indices": result["conclusion"]["rejected_indices"],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
