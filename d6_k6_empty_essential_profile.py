#!/usr/bin/env python3
"""Exact empty-essential profiler on the pinned 805-graph K6-only residue.

For every required K6 seed, this program asks whether the current exact
Lorentz/Hall relaxation has an all-nonempty branch.  Such a witness makes
virtual-K7 promotion irrelevant for rejecting that seed.  Only when every
passing branch needs one or two empty K6 defect vectors does the program
exhaust all common ``Z0`` choices and retain the identities of the possible
empty vertices.

Two state spaces are reported.  ``budget`` uses only the inherited global
cardinality bound of two empty vertices.  ``opposite_apex`` additionally uses
the exact fact that two distinct empty K6 defects are the opposite common K6
neighbours, whose squared distance is 7/3; hence their candidate pair cannot
be a required unit edge.  Candidate nonedges are never required to be
non-unit outside this branch.

All computations are exact integer and bit-mask logic.  Nonbipartite Lorentz
components remain ignored conservatively, as in the parent filter.  A
surviving state is only a necessary relaxation, not a realization.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import sys
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Iterable, Sequence

import d6_k6_arbitrary_subset_hall as hall_parent
import d6_k6_empty_support_budget as empty_parent
import d6_k6_psd_z_hereditary as hall
import d6_k6_psd_zmatrix as ranks
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
    vertices,
)
from d6_k6_normal_inertia import InertiaCache
from test_d6_k6_lorentz import lower_bound_18_graph


ROOT = Path(__file__).resolve().parent
PARENT_REPORT = ROOT / "d6_k6_empty_support_budget_report.json"
PARENT_VERIFICATION = ROOT / "d6_k6_empty_support_budget_verification.json"
EXPECTED_PARENT_REPORT_SHA256 = (
    "1860dbe69ae74b55203ea283cf83afff929a4b1f5cbfcd30e09f0138688eba9f"
)
EXPECTED_PARENT_SOURCE_SHA256 = (
    "284b8c3bda4d2a43581665ade3d52d9295b29d460418454e3a760897cabc5fac"
)
EXPECTED_INPUT = 805
EXPECTED_INPUT_SHA256 = (
    "53cad35eedf4f087d1e1dadbffcde54cbbbf04dbd5ee6ea6e10506329c958818"
)
REPORT_SCHEMA = "d6-k6-empty-essential-profile-v1"
DEFAULT_WORKERS = 11
GLOBAL_EMPTY_BUDGET = 2


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
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def worker_initializer() -> None:
    for name in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[name] = "1"


def load_residue() -> tuple[list[dict], list[int], dict[str, str]]:
    """Bind and reconstruct the independently checked ordered 805 residue."""

    if sha256(PARENT_REPORT) != EXPECTED_PARENT_REPORT_SHA256:
        raise ValueError("empty-support parent report hash changed")
    if sha256(ROOT / "d6_k6_empty_support_budget.py") != EXPECTED_PARENT_SOURCE_SHA256:
        raise ValueError("empty-support parent source hash changed")
    _, records, parent_indices, _ = empty_parent.verify_parent_input()
    report = json.loads(PARENT_REPORT.read_text(encoding="utf-8"))
    verification = json.loads(PARENT_VERIFICATION.read_text(encoding="utf-8"))
    residue = [int(index) for index in report.get("ordered_residue_indices", [])]
    if (
        report.get("schema") != "d6-k6-empty-support-budget-v1"
        or report.get("status") != "COMPLETE"
        or report.get("input_graphs") != 822
        or report.get("graphs_rejected") != 17
        or report.get("graphs_surviving") != EXPECTED_INPUT
        or report.get("ordered_input_indices") != parent_indices
        or len(residue) != EXPECTED_INPUT
        or stable_hash(residue) != EXPECTED_INPUT_SHA256
        or report.get("ordered_residue_indices_sha256") != EXPECTED_INPUT_SHA256
    ):
        raise ValueError("empty-support parent residue boundary changed")
    if (
        verification.get("schema")
        != "d6-k6-empty-support-budget-verification-v1"
        or verification.get("status") != "PASS"
        or verification.get("report_sha256") != EXPECTED_PARENT_REPORT_SHA256
        or verification.get("graphs_recomputed") != 822
        or verification.get("graphs_rejected") != 17
        or verification.get("graphs_surviving") != EXPECTED_INPUT
        or verification.get("ordered_residue_indices_sha256")
        != EXPECTED_INPUT_SHA256
    ):
        raise ValueError("empty-support parent independent verification changed")
    by_index = {int(record["index"]): record for record in records}
    selected = [by_index[index] for index in residue]
    for record in selected:
        if find_clique_mask(record["adjacency"], 7):
            raise ValueError("805-graph profile input unexpectedly contains a K7")
    return selected, residue, {
        PARENT_REPORT.name: sha256(PARENT_REPORT),
        PARENT_VERIFICATION.name: sha256(PARENT_VERIFICATION),
        "d6_k6_empty_support_budget.py": sha256(
            ROOT / "d6_k6_empty_support_budget.py"
        ),
    }


def _zero_capable_singletons(
    instance: K6LorentzInstance, groups: Sequence[hall.Group], side
) -> tuple[tuple[int, int, int], ...]:
    """Return ``(group index, absolute vertex, allowed mask)`` triples."""

    if len(groups) != len(side.components):
        raise ValueError("side group/component count mismatch")
    local = {absolute: i for i, absolute in enumerate(instance.outside)}
    answer = []
    for group_index, detail in enumerate(side.components):
        absolute = tuple(int(vertex) for vertex in detail["vertices"])
        rank_lower = int(detail["componentwise_fused_rank_lower"])
        if len(absolute) == 1 and rank_lower == 0:
            vertex = absolute[0]
            answer.append((group_index, vertex, instance.defects[local[vertex]]))
        elif len(absolute) == 1 and rank_lower != 1:
            raise AssertionError("singleton block rank lower is neither zero nor one")
    return tuple(answer)


def _force_singleton_status(
    groups: Sequence[hall.Group],
    zero_capable: Sequence[tuple[int, int, int]],
    chosen_empty_group: int | None,
) -> tuple[hall.Group, ...]:
    by_group = {group: (vertex, mask) for group, vertex, mask in zero_capable}
    if chosen_empty_group is not None and chosen_empty_group not in by_group:
        raise ValueError("chosen empty group is not zero-capable")
    answer = []
    for group_index, group in enumerate(groups):
        item = by_group.get(group_index)
        if item is None:
            answer.append(group)
            continue
        vertex, mask = item
        omit = hall.Choice(0, 0, (), "omit")
        if group_index == chosen_empty_group:
            answer.append(hall.Group(group.label, (omit,), True))
        else:
            nonempty = hall.Choice(
                1, mask, (vertex,), "forced_nonempty_negative_singleton"
            )
            answer.append(hall.Group(group.label, (omit, nonempty), True))
    return tuple(answer)


def orientation_empty_choices(
    instance: K6LorentzInstance,
    z0: int,
    side_a,
    side_b,
    subset_rank: hall_parent.PositiveSubsetRank,
    negative_side: str,
) -> tuple[int | None, ...]:
    """Enumerate every Hall-passing all-nonempty/singleton-empty choice."""

    groups_a, _ = hall_parent.make_side_groups(instance, "A", side_a, subset_rank)
    groups_b, _ = hall_parent.make_side_groups(instance, "B", side_b, subset_rank)
    if negative_side == "A":
        negative_groups, negative_rank = groups_a, side_a
    elif negative_side == "B":
        negative_groups, negative_rank = groups_b, side_b
    else:
        raise ValueError("negative_side must be A or B")
    zero_capable = _zero_capable_singletons(
        instance, negative_groups, negative_rank
    )
    z0_groups = tuple(
        hall.singleton_group(
            f"Z0:{instance.outside[local]}",
            instance.outside[local],
            instance.defects[local],
        )
        for local in vertices(z0)
    )
    by_group = {group: vertex for group, vertex, _ in zero_capable}
    passing: list[int | None] = []
    for chosen in (None, *(group for group, _, _ in zero_capable)):
        if negative_side == "A":
            test_a = _force_singleton_status(groups_a, zero_capable, chosen)
            test_b = groups_b
        else:
            test_a = groups_a
            test_b = _force_singleton_status(groups_b, zero_capable, chosen)
        if hall.check_groups(test_a + test_b + z0_groups).passed:
            passing.append(None if chosen is None else by_group[chosen])
    return tuple(passing)


def component_empty_choices(
    adjacency: Sequence[int],
    instance: K6LorentzInstance,
    z0: int,
    component,
    inertia: InertiaCache,
    zero_forcing: ZeroForcingSolver,
    subset_rank: hall_parent.PositiveSubsetRank,
) -> tuple[int | None, ...]:
    """Union every exact generic orientation and lightlike alternative."""

    graph_a, absolute_a = induced_required_graph(
        adjacency, instance, component.side_a
    )
    graph_b, absolute_b = induced_required_graph(
        adjacency, instance, component.side_b
    )
    choices: set[int | None] = set()
    for sign_a, sign_b, negative_side in (
        ("positive", "negative", "B"),
        ("negative", "positive", "A"),
    ):
        side_a = ranks.side_rank_lower(
            graph_a, absolute_a, sign_a, inertia, zero_forcing
        )
        side_b = ranks.side_rank_lower(
            graph_b, absolute_b, sign_b, inertia, zero_forcing
        )
        choices.update(
            orientation_empty_choices(
                instance, z0, side_a, side_b, subset_rank, negative_side
            )
        )
    light = z0 | component.component
    if (
        light.bit_count() <= COORDINATES
        and support_matching(light, instance.defects) is not None
    ):
        choices.add(None)
    ordered = tuple(sorted(choice for choice in choices if choice is not None))
    return (None, *ordered) if None in choices else ordered


def extend_states(
    states: Iterable[tuple[int, ...]],
    options: Sequence[int | None],
    adjacency: Sequence[int],
    require_opposite_pair: bool,
) -> set[tuple[int, ...]]:
    """Combine one component while enforcing the global two-empty budget."""

    following: set[tuple[int, ...]] = set()
    for state in states:
        for option in options:
            if option is None:
                following.add(tuple(state))
                continue
            if option in state:
                raise AssertionError("one vertex appeared in two components")
            chosen = tuple(sorted((*state, int(option))))
            if len(chosen) > GLOBAL_EMPTY_BUDGET:
                continue
            if (
                require_opposite_pair
                and len(chosen) == 2
                and adjacency[chosen[0]] & (1 << chosen[1])
            ):
                continue
            following.add(chosen)
    return following


def classify_seed(
    adjacency: Sequence[int], instance: K6LorentzInstance
) -> dict:
    """Find an all-nonempty witness or exhaust every possible empty state."""

    inertia = InertiaCache()
    zero_forcing = ZeroForcingSolver()
    subset_rank = hall_parent.PositiveSubsetRank(inertia, zero_forcing)
    budget_states: set[tuple[int, ...]] = set()
    opposite_states: set[tuple[int, ...]] = set()
    counts: Counter[str] = Counter()
    all_nonempty_witness = None

    for z0 in ranks.z0_subsets(instance.eligible_z0_mask):
        counts["z0_considered"] += 1
        if support_matching(z0, instance.defects) is None:
            counts["z0_unmatchable"] += 1
            continue
        counts["z0_matchable"] += 1
        local_budget: set[tuple[int, ...]] = {()}
        local_opposite: set[tuple[int, ...]] = {()}
        for component in lorentz_components(instance, z0):
            if not component.bipartite:
                counts["nonbipartite_components_ignored"] += 1
                continue
            counts["bipartite_components"] += 1
            counts["generic_orientations"] += 2
            choices = component_empty_choices(
                adjacency,
                instance,
                z0,
                component,
                inertia,
                zero_forcing,
                subset_rank,
            )
            counts["component_distinct_empty_options"] += len(choices)
            local_budget = extend_states(
                local_budget, choices, adjacency, require_opposite_pair=False
            )
            local_opposite = extend_states(
                local_opposite, choices, adjacency, require_opposite_pair=True
            )
            if not local_budget:
                counts["z0_component_impossible"] += 1
                break
        if local_budget:
            counts["z0_budget_feasible"] += 1
            budget_states.update(local_budget)
        if local_opposite:
            counts["z0_opposite_apex_feasible"] += 1
            opposite_states.update(local_opposite)
        if () in local_opposite:
            all_nonempty_witness = [
                instance.outside[local] for local in vertices(z0)
            ]
            counts["early_all_nonempty_witness"] += 1
            break

    if all_nonempty_witness is not None:
        return {
            "classification": "ALL_NONEMPTY_WITNESS",
            "budget_feasible": True,
            "opposite_apex_feasible": True,
            "empty_essential": False,
            "all_nonempty_witness_Z0": all_nonempty_witness,
            "exhaustive_empty_state_profile": False,
            "possible_empty_cardinalities": [0],
            "possible_empty_states": [[]],
            "promotion_apices": [],
            "promotion_pairs": [],
            "singleton_empty_apices": [],
            "mandatory_empty_vertices": [],
            "counts": dict(counts),
        }

    ordered_budget = sorted(budget_states, key=lambda state: (len(state), state))
    ordered_opposite = sorted(
        opposite_states, key=lambda state: (len(state), state)
    )
    if not ordered_budget:
        classification = "BUDGET_IMPOSSIBLE"
    elif not ordered_opposite:
        classification = "OPPOSITE_APEX_IMPOSSIBLE"
    else:
        classification = "EMPTY_ESSENTIAL"
    promotion_apices = sorted({vertex for state in ordered_opposite for vertex in state})
    singleton_apices = sorted(state[0] for state in ordered_opposite if len(state) == 1)
    pairs = [list(state) for state in ordered_opposite if len(state) == 2]
    mandatory = []
    if ordered_opposite:
        mandatory = sorted(set(ordered_opposite[0]).intersection(*map(set, ordered_opposite[1:])))
    return {
        "classification": classification,
        "budget_feasible": bool(ordered_budget),
        "opposite_apex_feasible": bool(ordered_opposite),
        "empty_essential": classification == "EMPTY_ESSENTIAL",
        "all_nonempty_witness_Z0": None,
        "exhaustive_empty_state_profile": True,
        "budget_possible_empty_cardinalities": sorted({len(state) for state in ordered_budget}),
        "budget_possible_empty_state_count": len(ordered_budget),
        "possible_empty_cardinalities": sorted({len(state) for state in ordered_opposite}),
        "possible_empty_states": [list(state) for state in ordered_opposite],
        "promotion_apices": promotion_apices,
        "promotion_pairs": pairs,
        "singleton_empty_apices": singleton_apices,
        "mandatory_empty_vertices": mandatory,
        "counts": dict(counts),
    }


def evaluate_record(record: dict) -> dict:
    adjacency = tuple(map(int, record["adjacency"]))
    seeds = []
    counts: Counter[str] = Counter()
    classifications: Counter[str] = Counter()
    for seed_mask in clique_masks(adjacency, COORDINATES):
        instance = build_instance(adjacency, vertices(seed_mask))
        decision = classify_seed(adjacency, instance)
        classifications[decision["classification"]] += 1
        counts.update(decision["counts"])
        # Keep the 25,243 routine witnesses compact.  The complete state space
        # is retained only for exceptional (essential or impossible) seeds.
        if decision["classification"] == "ALL_NONEMPTY_WITNESS":
            seeds.append(
                {
                    "seed_mask": seed_mask,
                    "classification": decision["classification"],
                    "all_nonempty_witness_Z0": decision[
                        "all_nonempty_witness_Z0"
                    ],
                }
            )
        else:
            compact = dict(decision)
            compact.pop("counts")
            seeds.append({"seed_mask": seed_mask, "seed": list(instance.seed), **compact})
    return {
        "index": int(record["index"]),
        "k6_seeds": len(seeds),
        "classification_counts": dict(classifications),
        "newly_impossible_seed_masks": [
            row["seed_mask"]
            for row in seeds
            if row["classification"]
            in {"BUDGET_IMPOSSIBLE", "OPPOSITE_APEX_IMPOSSIBLE"}
        ],
        "empty_essential_seed_masks": [
            row["seed_mask"]
            for row in seeds
            if row["classification"] == "EMPTY_ESSENTIAL"
        ],
        "counts": dict(counts),
        "seed_results": seeds,
    }


def positive_control() -> dict:
    adjacency = lower_bound_18_graph()
    rows = []
    for seed_mask in clique_masks(adjacency, COORDINATES):
        instance = build_instance(adjacency, vertices(seed_mask))
        decision = classify_seed(adjacency, instance)
        rows.append(
            {
                "seed_mask": seed_mask,
                "classification": decision["classification"],
                "all_nonempty_witness_Z0": decision["all_nonempty_witness_Z0"],
            }
        )
    if len(rows) != 32 or any(
        row["classification"] != "ALL_NONEMPTY_WITNESS" for row in rows
    ):
        raise AssertionError("known realizable 18-point control lost all-nonempty branch")
    return {"passed": True, "K6_seeds": len(rows), "seed_results": rows}


def aggregate(results: Sequence[dict]) -> dict:
    totals: Counter[str] = Counter()
    classifications: Counter[str] = Counter()
    effort: Counter[str] = Counter()
    essential_keys = []
    promotion_keys = []
    pair_keys = []
    impossible_keys = []
    for result in results:
        totals.update(
            graphs=1,
            k6_seeds=result["k6_seeds"],
            graphs_with_empty_essential_seed=bool(result["empty_essential_seed_masks"]),
            graphs_with_newly_impossible_seed=bool(result["newly_impossible_seed_masks"]),
        )
        classifications.update(result["classification_counts"])
        effort.update(result["counts"])
        for row in result["seed_results"]:
            key = [result["index"], row["seed_mask"]]
            if row["classification"] == "EMPTY_ESSENTIAL":
                essential_keys.append(key)
                promotion_keys.extend(
                    [*key, apex] for apex in row["promotion_apices"]
                )
                pair_keys.extend(
                    [*key, *pair] for pair in row["promotion_pairs"]
                )
            if row["classification"] in {
                "BUDGET_IMPOSSIBLE",
                "OPPOSITE_APEX_IMPOSSIBLE",
            }:
                impossible_keys.append(key)
    return {
        "totals": dict(totals),
        "seed_classification_counts": dict(classifications),
        "effort": dict(effort),
        "empty_essential_seed_keys": essential_keys,
        "empty_essential_seed_keys_sha256": stable_hash(essential_keys),
        "virtual_K7_promotion_keys": promotion_keys,
        "virtual_K7_promotion_keys_sha256": stable_hash(promotion_keys),
        "opposite_apex_pair_keys": pair_keys,
        "opposite_apex_pair_keys_sha256": stable_hash(pair_keys),
        "newly_impossible_seed_keys": impossible_keys,
        "newly_impossible_seed_keys_sha256": stable_hash(impossible_keys),
        "newly_rejected_indices": [
            result["index"] for result in results if result["newly_impossible_seed_masks"]
        ],
    }


def run(output: Path, workers: int) -> dict:
    records, indices, dependencies = load_residue()
    started = time.perf_counter()
    if workers == 1:
        results = [evaluate_record(record) for record in records]
    else:
        with ProcessPoolExecutor(
            max_workers=workers, initializer=worker_initializer
        ) as pool:
            results = list(pool.map(evaluate_record, records, chunksize=1))
    if [row["index"] for row in results] != indices:
        raise AssertionError("process pool changed ordered residue")
    summary = aggregate(results)
    report = {
        "schema": REPORT_SCHEMA,
        "status": "COMPLETE",
        "description": (
            "Exact targeted profile of all-nonempty versus empty-essential "
            "K6 Hall branches on the pinned 805-graph residue."
        ),
        "input_graphs": len(indices),
        "ordered_input_indices": indices,
        "input_indices_sha256": stable_hash(indices),
        "aggregate": summary,
        "graph_results": results,
        "positive_control": positive_control(),
        "sources": {
            **dependencies,
            Path(__file__).name: sha256(Path(__file__).resolve()),
        },
        "semantics": {
            "candidate_nonedges": (
                "unconstrained and may be unit; a candidate nonedge is used "
                "only as permission for the branch-forced non-unit opposite pair"
            ),
            "allowed_defects": "upper bounds; allowed coordinates may be zero",
            "all_nonempty_branch": "explicitly searched first and never discarded",
            "maximum_empty_vertices": GLOBAL_EMPTY_BUDGET,
            "two_empty_vertices": (
                "opposite common K6 neighbours with squared distance 7/3, "
                "so their pair cannot be a required unit edge"
            ),
            "empty_essential": (
                "every state surviving this necessary Hall/opposite-apex "
                "relaxation contains at least one empty defect"
            ),
            "promotion_apex": (
                "appears in at least one surviving state of an empty-essential seed"
            ),
            "nonbipartite_components": "ignored conservatively",
            "surviving_state": "necessary relaxation only, not a realization",
            "arithmetic": "exact integer and bit-mask logic only",
        },
        "runtime": {
            "command": list(sys.argv),
            "workers": workers,
            "wall_seconds": time.perf_counter() - started,
            "python": sys.version,
            "platform": platform.platform(),
        },
    }
    atomic_json(output, report)
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "d6_k6_empty_essential_profile.json",
    )
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    args = parser.parse_args()
    if args.workers < 1 or args.workers > DEFAULT_WORKERS:
        parser.error(f"workers must lie in 1..{DEFAULT_WORKERS}")
    report = run(args.output.resolve(), args.workers)
    print(json.dumps(report["aggregate"], sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
