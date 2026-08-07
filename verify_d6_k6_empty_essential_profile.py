#!/usr/bin/env python3
"""Independent replay of the 805-graph K6 empty-essential profile.

This checker does not import the production profiler or the empty-support
production kernel.  It reconstructs the ordered 805-graph residue from the
pinned parent artifacts, independently implements the Hall state traversal,
recomputes every K6 seed, and exactly compares the complete report.
"""

from __future__ import annotations

import argparse
import ast
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
import d6_k6_psd_z_hereditary as groups_parent
import d6_k6_psd_zmatrix as ranks
from d6_k6_bipartite_rank_reference import (
    ZeroForcingSolver,
    induced_required_graph,
    lorentz_components,
)
from d6_k6_lorentz_reference import (
    COORDINATES,
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
VERIFICATION_SCHEMA = "d6-k6-empty-essential-profile-verification-v1"
GLOBAL_EMPTY_BUDGET = 2
DEFAULT_WORKERS = 11


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


def assert_import_independence() -> None:
    tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    imports = {
        node.module
        for node in ast.walk(tree)
        if isinstance(node, ast.ImportFrom) and node.module
    }
    imports.update(
        alias.name
        for node in ast.walk(tree)
        if isinstance(node, ast.Import)
        for alias in node.names
    )
    forbidden = {
        "d6_k6_empty_essential_profile",
        "d6_k6_empty_support_budget",
        "probe_d6_k6_virtual_k7_empty_conjunction",
    }
    if imports & forbidden:
        raise AssertionError(
            f"independent verifier imports production code: {imports & forbidden}"
        )


def reconstruct_residue() -> tuple[list[dict], list[int], dict[str, str]]:
    """Rebuild 805 input without importing the empty-support producer."""

    if sha256(PARENT_REPORT) != EXPECTED_PARENT_REPORT_SHA256:
        raise ValueError("empty-support parent report hash changed")
    if sha256(ROOT / "d6_k6_empty_support_budget.py") != EXPECTED_PARENT_SOURCE_SHA256:
        raise ValueError("empty-support parent source hash changed")
    _, records, arbitrary_indices = hall_parent.verify_inputs()
    parent = json.loads(PARENT_REPORT.read_text(encoding="utf-8"))
    verification = json.loads(PARENT_VERIFICATION.read_text(encoding="utf-8"))
    arbitrary_report = json.loads(
        (ROOT / "d6_k6_arbitrary_subset_hall_report.json").read_text(
            encoding="utf-8"
        )
    )
    inherited = {int(index) for index in arbitrary_report["rejected_indices"]}
    expected_822 = [index for index in arbitrary_indices if index not in inherited]
    residue = [int(index) for index in parent.get("ordered_residue_indices", [])]
    if (
        parent.get("schema") != "d6-k6-empty-support-budget-v1"
        or parent.get("status") != "COMPLETE"
        or parent.get("ordered_input_indices") != expected_822
        or parent.get("graphs_rejected") != 17
        or parent.get("graphs_surviving") != EXPECTED_INPUT
        or len(residue) != EXPECTED_INPUT
        or stable_hash(residue) != EXPECTED_INPUT_SHA256
        or parent.get("ordered_residue_indices_sha256") != EXPECTED_INPUT_SHA256
    ):
        raise ValueError("independent 805-graph reconstruction failed")
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
        raise ValueError("parent independent verification boundary changed")
    by_index = {int(record["index"]): record for record in records}
    selected = [by_index[index] for index in residue]
    if any(find_clique_mask(record["adjacency"], 7) for record in selected):
        raise ValueError("independent input unexpectedly contains K7")
    return selected, residue, {
        PARENT_REPORT.name: sha256(PARENT_REPORT),
        PARENT_VERIFICATION.name: sha256(PARENT_VERIFICATION),
        "d6_k6_empty_support_budget.py": sha256(
            ROOT / "d6_k6_empty_support_budget.py"
        ),
    }


def independent_hall(groups: Sequence[groups_parent.Group]) -> bool:
    """Fresh 64-state Hall traversal over every selectable subblock."""

    states: dict[int, int] = {0: 0}
    for group in groups:
        following: dict[int, int] = {}
        for coordinates, rank in states.items():
            for choice in group.choices:
                union = coordinates | choice.coordinates
                total = rank + choice.rank
                if total > union.bit_count():
                    return False
                if total > following.get(union, -1):
                    following[union] = total
        states = following
    return True


def singleton_data(instance, side) -> dict[int, tuple[int, int]]:
    local = {absolute: i for i, absolute in enumerate(instance.outside)}
    answer = {}
    for group_index, detail in enumerate(side.components):
        absolute = tuple(int(vertex) for vertex in detail["vertices"])
        rank = int(detail["componentwise_fused_rank_lower"])
        if len(absolute) == 1 and rank == 0:
            vertex = absolute[0]
            answer[group_index] = (vertex, instance.defects[local[vertex]])
        elif len(absolute) == 1 and rank != 1:
            raise AssertionError("independent singleton rank lower is invalid")
    return answer


def impose_singleton_status(
    groups: Sequence[groups_parent.Group],
    zero: dict[int, tuple[int, int]],
    empty_group: int | None,
) -> tuple[groups_parent.Group, ...]:
    answer = []
    for group_index, group in enumerate(groups):
        if group_index not in zero:
            answer.append(group)
            continue
        vertex, mask = zero[group_index]
        omit = groups_parent.Choice(0, 0, (), "independent_omit")
        choices = (omit,)
        if group_index != empty_group:
            choices += (
                groups_parent.Choice(
                    1, mask, (vertex,), "independent_forced_nonempty"
                ),
            )
        answer.append(groups_parent.Group(group.label, choices, True))
    return tuple(answer)


def orientation_choices(
    instance, z0, side_a, side_b, subset_rank, negative_side: str
) -> tuple[int | None, ...]:
    groups_a, _ = hall_parent.make_side_groups(
        instance, "A", side_a, subset_rank
    )
    groups_b, _ = hall_parent.make_side_groups(
        instance, "B", side_b, subset_rank
    )
    zero = singleton_data(instance, side_a if negative_side == "A" else side_b)
    zgroups = tuple(
        groups_parent.singleton_group(
            f"Z0:{instance.outside[local]}",
            instance.outside[local],
            instance.defects[local],
        )
        for local in vertices(z0)
    )
    passing: list[int | None] = []
    for chosen in (None, *zero.keys()):
        if negative_side == "A":
            test_a = impose_singleton_status(groups_a, zero, chosen)
            test_b = groups_b
        else:
            test_a = groups_a
            test_b = impose_singleton_status(groups_b, zero, chosen)
        if independent_hall(test_a + test_b + zgroups):
            passing.append(None if chosen is None else zero[chosen][0])
    return tuple(passing)


def component_choices(
    adjacency, instance, z0, component, inertia, forcing, subset_rank
) -> tuple[int | None, ...]:
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
            graph_a, absolute_a, sign_a, inertia, forcing
        )
        side_b = ranks.side_rank_lower(
            graph_b, absolute_b, sign_b, inertia, forcing
        )
        choices.update(
            orientation_choices(
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


def independent_extend(
    states: Iterable[tuple[int, ...]],
    options: Sequence[int | None],
    adjacency: Sequence[int],
    require_nonedge_pair: bool,
) -> set[tuple[int, ...]]:
    following = set()
    for state in states:
        for option in options:
            if option is None:
                following.add(tuple(state))
                continue
            if option in state:
                raise AssertionError("duplicate empty vertex across components")
            chosen = tuple(sorted((*state, int(option))))
            if len(chosen) > GLOBAL_EMPTY_BUDGET:
                continue
            if (
                require_nonedge_pair
                and len(chosen) == 2
                and adjacency[chosen[0]] & (1 << chosen[1])
            ):
                continue
            following.add(chosen)
    return following


def classify_seed(adjacency, instance) -> dict:
    inertia = InertiaCache()
    forcing = ZeroForcingSolver()
    subset = hall_parent.PositiveSubsetRank(inertia, forcing)
    all_budget: set[tuple[int, ...]] = set()
    all_opposite: set[tuple[int, ...]] = set()
    counts: Counter[str] = Counter()

    for z0 in ranks.z0_subsets(instance.eligible_z0_mask):
        counts["z0_considered"] += 1
        if support_matching(z0, instance.defects) is None:
            counts["z0_unmatchable"] += 1
            continue
        counts["z0_matchable"] += 1
        budget: set[tuple[int, ...]] = {()}
        opposite: set[tuple[int, ...]] = {()}
        for component in lorentz_components(instance, z0):
            if not component.bipartite:
                counts["nonbipartite_components_ignored"] += 1
                continue
            counts["bipartite_components"] += 1
            counts["generic_orientations"] += 2
            options = component_choices(
                adjacency, instance, z0, component, inertia, forcing, subset
            )
            counts["component_distinct_empty_options"] += len(options)
            budget = independent_extend(
                budget, options, adjacency, require_nonedge_pair=False
            )
            opposite = independent_extend(
                opposite, options, adjacency, require_nonedge_pair=True
            )
            if not budget:
                counts["z0_component_impossible"] += 1
                break
        if budget:
            counts["z0_budget_feasible"] += 1
            all_budget.update(budget)
        if opposite:
            counts["z0_opposite_apex_feasible"] += 1
            all_opposite.update(opposite)
        if () in opposite:
            counts["early_all_nonempty_witness"] += 1
            return {
                "classification": "ALL_NONEMPTY_WITNESS",
                "all_nonempty_witness_Z0": [
                    instance.outside[local] for local in vertices(z0)
                ],
                "counts": dict(counts),
            }

    budget_rows = sorted(all_budget, key=lambda state: (len(state), state))
    opposite_rows = sorted(all_opposite, key=lambda state: (len(state), state))
    if not budget_rows:
        classification = "BUDGET_IMPOSSIBLE"
    elif not opposite_rows:
        classification = "OPPOSITE_APEX_IMPOSSIBLE"
    else:
        classification = "EMPTY_ESSENTIAL"
    promotion_apices = sorted(
        {vertex for state in opposite_rows for vertex in state}
    )
    singleton_apices = sorted(
        state[0] for state in opposite_rows if len(state) == 1
    )
    mandatory = []
    if opposite_rows:
        mandatory = sorted(
            set(opposite_rows[0]).intersection(*map(set, opposite_rows[1:]))
        )
    return {
        "classification": classification,
        "seed_detail": {
            "classification": classification,
            "budget_feasible": bool(budget_rows),
            "opposite_apex_feasible": bool(opposite_rows),
            "empty_essential": classification == "EMPTY_ESSENTIAL",
            "all_nonempty_witness_Z0": None,
            "exhaustive_empty_state_profile": True,
            "budget_possible_empty_cardinalities": sorted(
                {len(state) for state in budget_rows}
            ),
            "budget_possible_empty_state_count": len(budget_rows),
            "possible_empty_cardinalities": sorted(
                {len(state) for state in opposite_rows}
            ),
            "possible_empty_states": [list(state) for state in opposite_rows],
            "promotion_apices": promotion_apices,
            "promotion_pairs": [
                list(state) for state in opposite_rows if len(state) == 2
            ],
            "singleton_empty_apices": singleton_apices,
            "mandatory_empty_vertices": mandatory,
        },
        "counts": dict(counts),
    }


def evaluate_record(record: dict) -> dict:
    adjacency = tuple(map(int, record["adjacency"]))
    seed_rows = []
    classifications: Counter[str] = Counter()
    effort: Counter[str] = Counter()
    for seed_mask in clique_masks(adjacency, COORDINATES):
        instance = build_instance(adjacency, vertices(seed_mask))
        result = classify_seed(adjacency, instance)
        classification = result["classification"]
        classifications[classification] += 1
        effort.update(result["counts"])
        if classification == "ALL_NONEMPTY_WITNESS":
            seed_rows.append(
                {
                    "seed_mask": seed_mask,
                    "classification": classification,
                    "all_nonempty_witness_Z0": result[
                        "all_nonempty_witness_Z0"
                    ],
                }
            )
        else:
            seed_rows.append(
                {
                    "seed_mask": seed_mask,
                    "seed": list(instance.seed),
                    **result["seed_detail"],
                }
            )
    impossible = [
        row["seed_mask"]
        for row in seed_rows
        if row["classification"]
        in {"BUDGET_IMPOSSIBLE", "OPPOSITE_APEX_IMPOSSIBLE"}
    ]
    essential = [
        row["seed_mask"]
        for row in seed_rows
        if row["classification"] == "EMPTY_ESSENTIAL"
    ]
    return {
        "index": int(record["index"]),
        "k6_seeds": len(seed_rows),
        "classification_counts": dict(classifications),
        "newly_impossible_seed_masks": impossible,
        "empty_essential_seed_masks": essential,
        "counts": dict(effort),
        "seed_results": seed_rows,
    }


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
            graphs_with_empty_essential_seed=bool(
                result["empty_essential_seed_masks"]
            ),
            graphs_with_newly_impossible_seed=bool(
                result["newly_impossible_seed_masks"]
            ),
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
            result["index"]
            for result in results
            if result["newly_impossible_seed_masks"]
        ],
    }


def positive_control() -> dict:
    adjacency = lower_bound_18_graph()
    rows = []
    for seed_mask in clique_masks(adjacency, COORDINATES):
        instance = build_instance(adjacency, vertices(seed_mask))
        result = classify_seed(adjacency, instance)
        rows.append(
            {
                "seed_mask": seed_mask,
                "classification": result["classification"],
                "all_nonempty_witness_Z0": result.get(
                    "all_nonempty_witness_Z0"
                ),
            }
        )
    if len(rows) != 32 or any(
        row["classification"] != "ALL_NONEMPTY_WITNESS" for row in rows
    ):
        raise AssertionError("independent positive control failed")
    return {"passed": True, "K6_seeds": len(rows), "seed_results": rows}


class VerificationInfrastructureAbort(RuntimeError):
    pass


def verify(
    report_path: Path,
    expected_report_sha256: str,
    output: Path,
    workers: int,
) -> dict:
    if (
        len(expected_report_sha256) != 64
        or any(character not in "0123456789abcdef" for character in expected_report_sha256)
    ):
        raise ValueError("--expected-report-sha256 must be 64 lowercase hex digits")
    actual_hash = sha256(report_path)
    if actual_hash != expected_report_sha256:
        raise ValueError(
            f"explicit report SHA-256 boundary mismatch: actual={actual_hash}"
        )
    assert_import_independence()
    records, indices, dependencies = reconstruct_residue()
    report = json.loads(report_path.read_text(encoding="utf-8"))
    producer = ROOT / "d6_k6_empty_essential_profile.py"
    producer_hash = sha256(producer)
    if (
        report.get("schema") != REPORT_SCHEMA
        or report.get("status") != "COMPLETE"
        or report.get("input_graphs") != EXPECTED_INPUT
        or report.get("ordered_input_indices") != indices
        or report.get("input_indices_sha256") != EXPECTED_INPUT_SHA256
        or report.get("sources", {}).get(producer.name) != producer_hash
    ):
        raise ValueError("production profile boundary changed")
    for name, digest in dependencies.items():
        if report.get("sources", {}).get(name) != digest:
            raise ValueError(f"production report dependency mismatch: {name}")

    started = time.perf_counter()
    try:
        if workers == 1:
            recomputed = [evaluate_record(record) for record in records]
        else:
            with ProcessPoolExecutor(
                max_workers=workers, initializer=worker_initializer
            ) as pool:
                recomputed = list(
                    pool.map(evaluate_record, records, chunksize=1)
                )
    except BaseException as error:
        raise VerificationInfrastructureAbort(
            "fresh 805-graph verifier pool aborted; no conclusion inferred"
        ) from error
    if recomputed != report.get("graph_results"):
        for observed, archived in zip(recomputed, report.get("graph_results", [])):
            if observed != archived:
                raise AssertionError(
                    f"independent graph/seed profile mismatch at {observed['index']}"
                )
        raise AssertionError("independent graph result length/order mismatch")
    summary = aggregate(recomputed)
    if summary != report.get("aggregate"):
        raise AssertionError("independent aggregate profile mismatch")
    control = positive_control()
    if control != report.get("positive_control"):
        raise AssertionError("independent positive control rows differ")

    result = {
        "schema": VERIFICATION_SCHEMA,
        "status": "PASS",
        "method": (
            "fresh ordered 805-graph process pool; production profiler and "
            "empty-support producer not imported; independent Hall DP and "
            "empty-state traversal; exact comparison of all 25,354 K6 seeds"
        ),
        "report_sha256": actual_hash,
        "production_source_sha256": producer_hash,
        "verifier_source_sha256": sha256(Path(__file__).resolve()),
        "graphs_recomputed": len(recomputed),
        "K6_seeds_recomputed": summary["totals"]["k6_seeds"],
        "empty_essential_seeds": summary["seed_classification_counts"].get(
            "EMPTY_ESSENTIAL", 0
        ),
        "graphs_with_empty_essential_seed": summary["totals"][
            "graphs_with_empty_essential_seed"
        ],
        "newly_impossible_seeds": len(
            summary["newly_impossible_seed_keys"]
        ),
        "newly_rejected_indices": summary["newly_rejected_indices"],
        "virtual_K7_promotion_keys": len(
            summary["virtual_K7_promotion_keys"]
        ),
        "opposite_apex_pair_keys": len(
            summary["opposite_apex_pair_keys"]
        ),
        "positive_control": control,
        "dependency_hashes": dependencies,
        "import_independence": {
            "production_profiler_imported": False,
            "empty_support_producer_imported": False,
            "checked_by_AST": True,
        },
        "runtime": {
            "command": list(sys.argv),
            "workers": workers,
            "wall_seconds": time.perf_counter() - started,
            "python": sys.version,
            "platform": platform.platform(),
        },
    }
    atomic_json(output, result)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--report",
        type=Path,
        default=ROOT / "d6_k6_empty_essential_profile.json",
    )
    parser.add_argument("--expected-report-sha256", required=True)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "d6_k6_empty_essential_profile_verification.json",
    )
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    args = parser.parse_args()
    if args.workers < 1 or args.workers > DEFAULT_WORKERS:
        parser.error(f"workers must lie in 1..{DEFAULT_WORKERS}")
    try:
        result = verify(
            args.report.resolve(),
            args.expected_report_sha256,
            args.output.resolve(),
            args.workers,
        )
    except VerificationInfrastructureAbort as error:
        atomic_json(
            args.output.resolve(),
            {
                "schema": VERIFICATION_SCHEMA,
                "status": "INFRA_ABORT",
                "meaning": "no mathematical conclusion is inferred",
                "message": str(error),
                "report_sha256": sha256(args.report.resolve()),
                "verifier_source_sha256": sha256(Path(__file__).resolve()),
            },
        )
        print(str(error), file=sys.stderr)
        return 2
    print(
        json.dumps(
            {
                "status": result["status"],
                "graphs_recomputed": result["graphs_recomputed"],
                "K6_seeds_recomputed": result["K6_seeds_recomputed"],
                "empty_essential_seeds": result["empty_essential_seeds"],
                "virtual_K7_promotion_keys": result[
                    "virtual_K7_promotion_keys"
                ],
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
