#!/usr/bin/env python3
"""Exact pilot for hereditary PSD--Z support Hall on the K6 residue.

For every connected side-graph block to which the PSD Z-matrix theorem
applies, every *proper* principal submatrix is positive definite.  Hence any
proper subset of the corresponding defect vectors is linearly independent.
The full block retains the strongest already-proved fused rank lower bound.

For one generic Lorentz orientation this pilot chooses at most one subset
from each connected side-graph block, independently on both Lorentz sides,
and any subset of the singleton Z0 blocks.  These chosen spans are mutually
orthogonal, so their rank lower bounds add.  Every such choice must fit in
the union of its allowed seed-coordinate masks.  Choices inside one connected
block are alternatives, never additive blocks; this avoids double counting
overlapping subsets of the same span.

The search is exact.  A 64-state dynamic program compresses choices having
the same allowed-coordinate union by retaining the largest rank lower bound;
this is dominance, not a heuristic prune.  Candidate nonedges remain
unconstrained, and allowed defect coordinates may be zero.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Sequence

import d6_k6_fused_side_rank as frozen
import d6_k6_psd_zmatrix as parent
from d6_k6_bipartite_rank_reference import (
    ZeroForcingSolver,
    induced_required_graph,
    lorentz_components,
)
from d6_k6_lorentz_reference import (
    COORDINATES,
    K6LorentzInstance,
    build_instance,
    clique_masks,
    find_clique_mask,
    support_matching,
    validate_graph,
    vertices,
)
from d6_k6_normal_inertia import InertiaCache
from test_d6_k6_lorentz import lower_bound_18_graph


ROOT = Path(__file__).resolve().parent
SCHEMA = "d6-k6-psd-z-hereditary-pilot-v1"
EXPECTED_PARENT_REPORT_SHA256 = (
    "abb920fbd9eb848502cca37e32c2aa86bf220ab1773e6006ea1fbee4bc2f4e30"
)
EXPECTED_INPUT = 861
EXPECTED_INPUT_SHA256 = (
    "09ebce17d2b72fa6514fc6a8a938376626d373161d2e1b4dd8d31c5e00c18db5"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_hash(value: object) -> str:
    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


def atomic_json(path: Path, value: object) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


@dataclass(frozen=True)
class VectorChoice:
    """One selected subset from a single orthogonal span block."""

    rank_lower: int
    allowed_coordinates: int
    selected_vertices: tuple[int, ...]
    kind: str


@dataclass(frozen=True)
class HereditaryGroup:
    """Alternative subset choices inside one connected side block or Z0."""

    label: str
    options: tuple[VectorChoice, ...]
    raw_option_count: int


@dataclass(frozen=True)
class HereditaryHallResult:
    passed: bool
    groups: int
    raw_combinations: int
    compressed_options: int
    dp_transitions: int
    terminal_states: int
    first_failure: dict | None


def allowed_union(instance: K6LorentzInstance, local_vertices: Sequence[int]) -> int:
    answer = 0
    for local in local_vertices:
        answer |= instance.defects[local]
    return answer


def side_groups(
    instance: K6LorentzInstance,
    side_label: str,
    side,
) -> tuple[HereditaryGroup, ...]:
    """Build exact subset alternatives for each connected F component."""

    local_by_absolute = {
        absolute: local for local, absolute in enumerate(instance.outside)
    }
    answer = []
    for number, component in enumerate(side.components):
        local_vertices = tuple(
            local_by_absolute[absolute] for absolute in component["vertices"]
        )
        size = len(local_vertices)
        applicable = bool(component["psd_zmatrix_applicable"])
        full_rank = int(component["componentwise_fused_rank_lower"])
        if not 0 <= full_rank <= size:
            raise AssertionError("invalid full-component rank lower bound")
        if applicable and full_rank < size - 1:
            raise AssertionError("applicable PSD Z block lost its n-1 bound")

        raw = [VectorChoice(0, 0, (), "empty")]
        if applicable:
            for selected in range(1, 1 << size):
                chosen = tuple(
                    local_vertices[position]
                    for position in range(size)
                    if selected & (1 << position)
                )
                full = selected == (1 << size) - 1
                raw.append(VectorChoice(
                    full_rank if full else len(chosen),
                    allowed_union(instance, chosen),
                    tuple(instance.outside[local] for local in chosen),
                    "full" if full else "proper",
                ))
        else:
            raw.append(VectorChoice(
                full_rank,
                allowed_union(instance, local_vertices),
                tuple(instance.outside[local] for local in local_vertices),
                "full_nonapplicable",
            ))

        # Equal coordinate unions are future-equivalent.  The largest rank
        # dominates all smaller ranks for detecting a Hall violation.
        best: dict[int, VectorChoice] = {}
        for choice in raw:
            old = best.get(choice.allowed_coordinates)
            if old is None or choice.rank_lower > old.rank_lower:
                best[choice.allowed_coordinates] = choice
        options = tuple(
            best[mask] for mask in sorted(best, key=lambda mask: (mask.bit_count(), mask))
        )
        answer.append(HereditaryGroup(
            f"{side_label}:F{number}", options, len(raw)
        ))
    return tuple(answer)


def z0_groups(
    instance: K6LorentzInstance, z0: int
) -> tuple[HereditaryGroup, ...]:
    return tuple(
        HereditaryGroup(
            f"Z0:{instance.outside[local]}",
            (
                VectorChoice(0, 0, (), "empty"),
                VectorChoice(
                    1,
                    instance.defects[local],
                    (instance.outside[local],),
                    "singleton",
                ),
            ),
            2,
        )
        for local in vertices(z0)
    )


def check_hereditary_hall(
    groups: Sequence[HereditaryGroup],
) -> HereditaryHallResult:
    """Check all one-choice-per-span Hall inequalities by exact dominance DP."""

    # coordinate union -> (largest rank, selected nonempty group choices)
    states: dict[int, tuple[int, tuple[tuple[str, VectorChoice], ...]]] = {
        0: (0, ())
    }
    transitions = 0
    raw_combinations = 1
    compressed_options = 0
    for group in groups:
        raw_combinations *= group.raw_option_count
        compressed_options += len(group.options)
        next_states: dict[
            int, tuple[int, tuple[tuple[str, VectorChoice], ...]]
        ] = {}
        for coordinates, (rank, witness) in states.items():
            for choice in group.options:
                transitions += 1
                new_coordinates = coordinates | choice.allowed_coordinates
                new_rank = rank + choice.rank_lower
                new_witness = witness
                if choice.selected_vertices:
                    new_witness += ((group.label, choice),)
                if new_rank > new_coordinates.bit_count():
                    return HereditaryHallResult(
                        False,
                        len(groups),
                        raw_combinations,
                        compressed_options,
                        transitions,
                        len(next_states),
                        {
                            "rank_lower": new_rank,
                            "coordinate_capacity": new_coordinates.bit_count(),
                            "allowed_coordinates": vertices(new_coordinates),
                            "selected_groups": [
                                {
                                    "group": label,
                                    "kind": selected.kind,
                                    "rank_lower": selected.rank_lower,
                                    "selected_vertices": list(
                                        selected.selected_vertices
                                    ),
                                    "allowed_coordinates": vertices(
                                        selected.allowed_coordinates
                                    ),
                                }
                                for label, selected in new_witness
                            ],
                        },
                    )
                old = next_states.get(new_coordinates)
                if old is None or new_rank > old[0]:
                    next_states[new_coordinates] = (new_rank, new_witness)
        states = next_states
    return HereditaryHallResult(
        True,
        len(groups),
        raw_combinations,
        compressed_options,
        transitions,
        len(states),
        None,
    )


def check_orientation(
    instance: K6LorentzInstance,
    z0: int,
    side_a,
    side_b,
) -> tuple[HereditaryHallResult, tuple[HereditaryGroup, ...]]:
    groups = (
        side_groups(instance, "A", side_a)
        + side_groups(instance, "B", side_b)
        + z0_groups(instance, z0)
    )
    return check_hereditary_hall(groups), groups


def check_component(
    adj: Sequence[int],
    instance: K6LorentzInstance,
    z0: int,
    component,
    inertia: InertiaCache,
    zero_forcing: ZeroForcingSolver,
) -> dict:
    if not component.bipartite:
        raise ValueError("hereditary generic checker requires bipartite L")

    graph_a, absolute_a = induced_required_graph(adj, instance, component.side_a)
    graph_b, absolute_b = induced_required_graph(adj, instance, component.side_b)
    orientations = []
    for name, sign_a, sign_b in (
        ("A_positive", "positive", "negative"),
        ("A_negative", "negative", "positive"),
    ):
        side_a = parent.side_rank_lower(
            graph_a, absolute_a, sign_a, inertia, zero_forcing
        )
        side_b = parent.side_rank_lower(
            graph_b, absolute_b, sign_b, inertia, zero_forcing
        )
        old = parent.check_block_subsets(frozen.make_blocks(
            instance,
            z0,
            component.side_a,
            component.side_b,
            side_a.rank_lower,
            side_b.rank_lower,
        ))
        hereditary, groups = check_orientation(instance, z0, side_a, side_b)
        if not old.passed and hereditary.passed:
            raise AssertionError("hereditary Hall must subsume the parent system")
        orientations.append({
            "case": name,
            "parent_passed": old.passed,
            "hereditary_passed": hereditary.passed,
            "strict_new_failure": old.passed and not hereditary.passed,
            "groups": hereditary.groups,
            "raw_combinations": hereditary.raw_combinations,
            "compressed_options": hereditary.compressed_options,
            "dp_transitions": hereditary.dp_transitions,
            "terminal_states": hereditary.terminal_states,
            "first_failure": hereditary.first_failure,
            "applicable_groups": sum(
                len(group.options) > 2 for group in groups
            ),
        })

    light_selected = z0 | component.component
    light_dimension = light_selected.bit_count()
    light_passed = (
        light_dimension <= COORDINATES
        and support_matching(light_selected, instance.defects) is not None
    )
    parent_passed = any(row["parent_passed"] for row in orientations) or light_passed
    hereditary_passed = (
        any(row["hereditary_passed"] for row in orientations) or light_passed
    )
    if not parent_passed and hereditary_passed:
        raise AssertionError("hereditary component decision revived parent failure")
    return {
        "component": frozen.absolute_vertices(instance, component.component),
        "orientations": orientations,
        "lightlike_passed": light_passed,
        "parent_passed": parent_passed,
        "hereditary_passed": hereditary_passed,
        "strict_new_failure": parent_passed and not hereditary_passed,
    }


def check_hereditary_system(
    adj: Sequence[int],
    instance: K6LorentzInstance,
    z0: int,
    components,
    inertia: InertiaCache,
    zero_forcing: ZeroForcingSolver,
) -> tuple[bool, bool, dict, dict | None]:
    """Return parent/new pass booleans, counters, and first new failure."""

    parent_passed = True
    hereditary_passed = True
    first_failure = None
    counts = {
        "components_checked": 0,
        "components_parent_failed": 0,
        "components_strict_new_failed": 0,
        "orientations_checked": 0,
        "orientations_parent_failed": 0,
        "orientations_strict_new_failed": 0,
        "applicable_groups": 0,
        "raw_combinations": 0,
        "compressed_options": 0,
        "dp_transitions": 0,
    }
    for component in components:
        if not component.bipartite:
            continue
        result = check_component(
            adj, instance, z0, component, inertia, zero_forcing
        )
        counts["components_checked"] += 1
        counts["components_parent_failed"] += not result["parent_passed"]
        counts["components_strict_new_failed"] += result["strict_new_failure"]
        for orientation in result["orientations"]:
            counts["orientations_checked"] += 1
            counts["orientations_parent_failed"] += not orientation["parent_passed"]
            counts["orientations_strict_new_failed"] += orientation[
                "strict_new_failure"
            ]
            counts["applicable_groups"] += orientation["applicable_groups"]
            counts["raw_combinations"] += orientation["raw_combinations"]
            counts["compressed_options"] += orientation["compressed_options"]
            counts["dp_transitions"] += orientation["dp_transitions"]
        parent_passed &= result["parent_passed"]
        hereditary_passed &= result["hereditary_passed"]
        if not result["hereditary_passed"] and first_failure is None:
            first_failure = result
    if not parent_passed and hereditary_passed:
        raise AssertionError("hereditary system revived parent failure")
    return parent_passed, hereditary_passed, counts, first_failure


def add_counts(target: dict[str, int], source: dict[str, int]) -> None:
    for key, value in source.items():
        target[key] = target.get(key, 0) + value


def evaluate_seed(adj: Sequence[int], instance: K6LorentzInstance) -> dict:
    inertia = InertiaCache()
    zero_forcing = ZeroForcingSolver()
    choices = []
    totals: dict[str, int] = {}
    parent_common = 0
    hereditary_alone = 0
    hereditary_common = 0
    for z0 in parent.z0_subsets(instance.eligible_z0_mask):
        totals["z0_considered"] = totals.get("z0_considered", 0) + 1
        matchable = support_matching(z0, instance.defects) is not None
        if not matchable:
            choices.append({
                "Z0": frozen.absolute_vertices(instance, z0),
                "matchable": False,
            })
            continue
        totals["z0_matchable"] = totals.get("z0_matchable", 0) + 1
        components = lorentz_components(instance, z0)
        old_passed, new_passed, local, failure = check_hereditary_system(
            adj, instance, z0, components, inertia, zero_forcing
        )
        add_counts(totals, local)
        nonbipartite = frozen.check_nonbipartite_system(
            instance, z0, components
        )
        nonbipartite_passed = nonbipartite.passed
        old_common = old_passed and nonbipartite_passed
        new_common = new_passed and nonbipartite_passed
        if new_common and not old_common:
            raise AssertionError("hereditary conjunction revived parent conjunction")
        parent_common += old_common
        hereditary_alone += new_passed
        hereditary_common += new_common
        choices.append({
            "Z0": frozen.absolute_vertices(instance, z0),
            "matchable": True,
            "parent_psd_zmatrix_passed": old_passed,
            "hereditary_passed": new_passed,
            "nonbipartite_passed": nonbipartite_passed,
            "parent_common_passed": old_common,
            "hereditary_common_passed": new_common,
            "first_hereditary_failure": failure,
        })

        # A common hereditary/nonbipartite Z0 is simultaneously a witness for
        # the standalone hereditary system and for the coupled system.  Later
        # Z0 choices cannot change either existential seed decision.
        if new_common:
            return {
                "parent_feasible": True,
                "hereditary_standalone_feasible": True,
                "hereditary_combined_feasible": True,
                "parent_common_Z0": parent_common,
                "hereditary_passing_Z0": hereditary_alone,
                "hereditary_common_Z0": hereditary_common,
                "counts": totals,
                "choices": None,
            }

    return {
        "parent_feasible": bool(parent_common),
        "hereditary_standalone_feasible": bool(hereditary_alone),
        "hereditary_combined_feasible": bool(hereditary_common),
        "parent_common_Z0": parent_common,
        "hereditary_passing_Z0": hereditary_alone,
        "hereditary_common_Z0": hereditary_common,
        "counts": totals,
        "choices": choices if not hereditary_common else None,
    }


def evaluate_graph(adj: Sequence[int]) -> dict:
    validate_graph(adj, require_alpha_two=True)
    if find_clique_mask(adj, 7):
        raise ValueError("hereditary K6 pilot target must be K6-only")
    totals: dict[str, int] = {}
    seeds_checked = 0
    combined_impossible = 0
    standalone_impossible = 0
    first_combined = None
    first_standalone = None
    for seed_mask in clique_masks(adj, COORDINATES):
        seeds_checked += 1
        instance = build_instance(adj, vertices(seed_mask))
        decision = evaluate_seed(adj, instance)
        if not decision["parent_feasible"]:
            raise AssertionError("pinned PSD Z survivor failed its parent system")
        add_counts(totals, decision["counts"])
        if not decision["hereditary_standalone_feasible"]:
            standalone_impossible += 1
            if first_standalone is None:
                first_standalone = {
                    "seed": list(instance.seed),
                    "decision": decision,
                }
        if not decision["hereditary_combined_feasible"]:
            combined_impossible += 1
            if first_combined is None:
                first_combined = {
                    "seed": list(instance.seed),
                    "decision": decision,
                }
            # If the same seed also fails the hereditary system without the
            # nonbipartite conjunction, both graph-level classifications are
            # already final.  (For a joint-only failure we continue, because
            # a later seed might still give a standalone rejection.)
            if not decision["hereditary_standalone_feasible"]:
                break
    rejected = combined_impossible > 0
    standalone_rejected = standalone_impossible > 0
    return {
        "rejected": rejected,
        "hereditary_standalone_rejected": standalone_rejected,
        "joint_only_rejected": rejected and not standalone_rejected,
        "seeds_checked": seeds_checked,
        "combined_impossible_seeds": combined_impossible,
        "standalone_impossible_seeds": standalone_impossible,
        **totals,
        "first_combined_failure": first_combined,
        "first_standalone_failure": first_standalone,
    }


def evaluate_record(record: dict) -> dict:
    return {
        "index": record["index"],
        "decision": evaluate_graph(record["adjacency"]),
    }


def determinant(matrix: Sequence[Sequence[int]]) -> int:
    """Exact Bareiss determinant for the small integer controls."""

    work = [list(row) for row in matrix]
    n = len(work)
    if n == 0:
        return 1
    sign = 1
    previous = 1
    for pivot_index in range(n - 1):
        pivot_row = next(
            (row for row in range(pivot_index, n) if work[row][pivot_index]),
            None,
        )
        if pivot_row is None:
            return 0
        if pivot_row != pivot_index:
            work[pivot_index], work[pivot_row] = work[pivot_row], work[pivot_index]
            sign = -sign
        pivot = work[pivot_index][pivot_index]
        for row in range(pivot_index + 1, n):
            for column in range(pivot_index + 1, n):
                numerator = (
                    work[row][column] * pivot
                    - work[row][pivot_index] * work[pivot_index][column]
                )
                if numerator % previous:
                    raise AssertionError("Bareiss exact division failed")
                work[row][column] = numerator // previous
        previous = pivot
    return sign * work[-1][-1]


def principal(matrix: Sequence[Sequence[int]], selected: Sequence[int]) -> list[list[int]]:
    return [[matrix[row][column] for column in selected] for row in selected]


def synthetic_controls() -> dict:
    singular = (
        (1, -1, 0),
        (-1, 2, -1),
        (0, -1, 1),
    )
    positive_definite = (
        (2, -1, 0),
        (-1, 2, -1),
        (0, -1, 2),
    )
    proper_singular = []
    proper_pd = []
    for size in (1, 2):
        for chosen in combinations(range(3), size):
            proper_singular.append(determinant(principal(singular, chosen)))
            proper_pd.append(determinant(principal(positive_definite, chosen)))
    if determinant(singular) != 0 or min(proper_singular) <= 0:
        raise AssertionError("singular irreducible PSD Z control failed")
    if determinant(positive_definite) <= 0 or min(proper_pd) <= 0:
        raise AssertionError("positive-definite PSD Z control failed")

    # Signature switching the path Laplacian changes its negative edge entries
    # to positive ones and preserves every principal determinant.
    signature = (1, -1, 1)
    switched = tuple(tuple(
        signature[row] * singular[row][column] * signature[column]
        for column in range(3)
    ) for row in range(3))
    switched_proper = [
        determinant(principal(switched, chosen))
        for size in (1, 2)
        for chosen in combinations(range(3), size)
    ]
    if determinant(switched) != 0 or switched_proper != proper_singular:
        raise AssertionError("bipartite signature-switch control failed")

    # The full old block fits in two coordinates (rank 2 <= capacity 2), but
    # the proper pair {0,1} is independently rank 2 in one allowed coordinate.
    hall_group = HereditaryGroup(
        "synthetic:F0",
        (
            VectorChoice(0, 0, (), "empty"),
            VectorChoice(1, 1, (0,), "proper"),
            VectorChoice(1, 1, (1,), "proper"),
            VectorChoice(2, 1, (0, 1), "proper"),
            VectorChoice(1, 2, (2,), "proper"),
            VectorChoice(2, 3, (0, 1, 2), "full"),
        ),
        8,
    )
    hall = check_hereditary_hall((hall_group,))
    if hall.passed or hall.first_failure is None:
        raise AssertionError("hereditary support-Hall strict control failed")
    if hall.first_failure["rank_lower"] != 2:
        raise AssertionError("unexpected hereditary Hall control witness")

    return {
        "singular_path_laplacian_determinant": determinant(singular),
        "singular_path_proper_principal_determinants": proper_singular,
        "positive_definite_path_determinant": determinant(positive_definite),
        "positive_definite_path_proper_principal_determinants": proper_pd,
        "signature_switched_proper_determinants": switched_proper,
        "reducible_counterexample_proper_zero": determinant(((0,),)) == 0,
        "hereditary_hall_strict_control": hall.first_failure,
    }


def load_input() -> tuple[list[dict], list[int]]:
    if sha256(ROOT / "d6_k6_psd_zmatrix_report.json") != EXPECTED_PARENT_REPORT_SHA256:
        raise ValueError("parent PSD Z-matrix report hash changed")
    _, records, _ = parent.verify_inputs()
    report = json.loads(
        (ROOT / "d6_k6_psd_zmatrix_report.json").read_text(encoding="utf-8")
    )
    indices = [
        item["index"]
        for item in report["graph_results"]
        if not item["decision"]["rejected"]
    ]
    if len(indices) != EXPECTED_INPUT or stable_hash(indices) != EXPECTED_INPUT_SHA256:
        raise ValueError("pinned 861-graph survivor boundary changed")
    by_index = {record["index"]: record for record in records}
    return [by_index[index] for index in indices], indices


def run(output: Path, workers: int, limit: int | None) -> dict:
    controls = synthetic_controls()
    records, all_indices = load_input()
    if limit is not None:
        if limit < 1:
            raise ValueError("limit must be positive")
        records = records[:limit]
    indices = [record["index"] for record in records]
    started = time.perf_counter()
    if workers == 1:
        graph_results = [evaluate_record(record) for record in records]
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            graph_results = list(pool.map(evaluate_record, records, chunksize=1))
    wall = time.perf_counter() - started

    rejected = [row for row in graph_results if row["decision"]["rejected"]]
    standalone = [
        row for row in graph_results
        if row["decision"]["hereditary_standalone_rejected"]
    ]
    joint_only = [
        row for row in graph_results if row["decision"]["joint_only_rejected"]
    ]
    numeric_fields = sorted({
        key
        for row in graph_results
        for key, value in row["decision"].items()
        if isinstance(value, int) and not isinstance(value, bool)
    })
    totals = {
        key: sum(row["decision"].get(key, 0) for row in graph_results)
        for key in numeric_fields
    }

    positive = evaluate_graph(lower_bound_18_graph())
    if positive["rejected"] or positive["seeds_checked"] != 32:
        raise AssertionError("known realizable 18-point control failed")
    report = {
        "schema": SCHEMA,
        "status": "COMPLETE" if limit is None else "PILOT",
        "description": (
            "Exact hereditary principal-submatrix PSD-Z support-Hall pilot, "
            "coupled on the same Z0 to the frozen nonbipartite system."
        ),
        "parent_input_graphs": len(all_indices),
        "parent_input_indices_sha256": stable_hash(all_indices),
        "graphs_checked": len(records),
        "checked_indices_sha256": stable_hash(indices),
        "graphs_rejected_marginal": len(rejected),
        "graphs_rejected_standalone": len(standalone),
        "graphs_rejected_joint_only": len(joint_only),
        "rejected_indices": [row["index"] for row in rejected],
        "standalone_rejected_indices": [row["index"] for row in standalone],
        "joint_only_rejected_indices": [row["index"] for row in joint_only],
        **totals,
        "graph_results": graph_results,
        "positive_control": {
            "passed": not positive["rejected"],
            "K6_seeds": positive["seeds_checked"],
            "standalone_impossible_seeds": positive[
                "standalone_impossible_seeds"
            ],
            "combined_impossible_seeds": positive["combined_impossible_seeds"],
        },
        "synthetic_controls": controls,
        "semantics": {
            "candidate_nonedges": "unconstrained and may be unit",
            "allowed_defects": "upper bounds; allowed coordinates may be zero",
            "generic_component_choices": (
                "at most one subset per connected side span; ranks never add "
                "between two overlapping subsets of that span"
            ),
            "lightlike": "retained as a separate alternative",
            "arithmetic": "integer decisions only; no floating-point screen",
        },
        "sources": {
            "d6_k6_psd_zmatrix_report.json": EXPECTED_PARENT_REPORT_SHA256,
            Path(__file__).name: sha256(Path(__file__).resolve()),
        },
        "runtime": {
            "command": " ".join(sys.argv),
            "workers": workers,
            "wall_seconds": wall,
            "python": sys.version,
            "platform": platform.platform(),
        },
    }
    atomic_json(output, report)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "d6_k6_psd_z_hereditary_pilot_report.json",
    )
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--limit", type=int)
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("workers must be positive")
    report = run(args.output.resolve(), args.workers, args.limit)
    print(json.dumps({
        "status": report["status"],
        "graphs_checked": report["graphs_checked"],
        "graphs_rejected_marginal": report["graphs_rejected_marginal"],
        "graphs_rejected_standalone": report["graphs_rejected_standalone"],
        "graphs_rejected_joint_only": report["graphs_rejected_joint_only"],
        "rejected_indices": report["rejected_indices"],
        "output": str(args.output.resolve()),
    }, sort_keys=True))


if __name__ == "__main__":
    main()
