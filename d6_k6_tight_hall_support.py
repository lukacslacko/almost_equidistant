#!/usr/bin/env python3
"""Exact tight-Hall support propagation on the 756-graph K6 residue.

This layer keeps the complete K6 quantifier over zero Lorentz factors,
generic sign orientations, the one possible empty defect on a negative side,
the lightlike alternative, and the global two-empty-point constraint.  Its
new rule is purely exact.

If a selection of mutually orthogonal defect subspans has total proved rank
equal to the size of the union of their *allowed* seed-coordinate supports,
then the actual subspans fill that coordinate space.  Every other orthogonal
side span is supported outside it.  The reduced supports are then checked for
nonzero Gram entries and for a same-Lorentz-line singleton collision.

Candidate nonedges remain unconstrained and every allowed coordinate may be
zero.  Nonbipartite Lorentz components are ignored conservatively.
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

import d6_k6_arbitrary_subset_hall as arbitrary
import d6_k6_psd_z_hereditary as hall
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
INPUT_MANIFEST = ROOT / "d6_current_residue_manifest_v5.json"
INPUT_VERIFICATION = ROOT / "d6_current_residue_manifest_v5_verification.json"
REPORT_SCHEMA = "d6-k6-tight-hall-support-v1"
CERTIFICATE_SCHEMA = "d6-k6-tight-hall-support-certificates-v1"
EXPECTED_INPUT = 756
EXPECTED_INPUT_SHA256 = "2cfedbc83f6ff01b6e386440fb7ea066b274d22ffbc52371c964e99d17cdcb90"
GLOBAL_EMPTY_BUDGET = 2
DEFAULT_WORKERS = 11

EXPECTED_DEPENDENCIES = {
    "d6_current_residue_manifest_v5.json": "164b2a12813845ef1f6d6ca5e3713ed1d2edac4be58563c1648a39ff8580c0b5",
    "d6_current_residue_manifest_v5_verification.json": "e871431b8b28fc04f0e58922ff3a6386054d15d37de5c9f9e67601ede47e4b46",
    "d6_k6_arbitrary_subset_hall.py": "cadc0048fcc7a77299bc2950e2eca57bfed87d54ed6ffd59886d155b2c9a88ff",
    "d6_k6_psd_z_hereditary.py": "7b869895cb9b4f7f3bbf1a1d1b11ec71507b6355c8ec0b7d7f8eaec2a51b180d",
    "d6_k6_psd_zmatrix.py": "24c2b3abfce5074560183d17f77778f96298bef6e38ece5bc937ccb2a4cba736",
    "d6_k6_bipartite_rank_reference.py": "e760fef74c423a89f1e1f9d08909fa1223e0389900cf3d45b0140ee296dea922",
    "d6_k6_lorentz_reference.py": "ef6a0877a5987119c92ec1bd9271a6a47774fc9e4992c10e1f28b6d5373592c2",
    "d6_k6_normal_inertia.py": "2d00ab40eb97aa79ffec3b8134b83c4fd5f703cbd691f094a25bf3ff9f3ee023",
    "test_d6_k6_lorentz.py": "f4fc8544a4311f742267808e6c785f82139accb9da99b496d206283b52831a47",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_hash(value: object) -> str:
    raw = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
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


def load_input() -> tuple[list[dict], list[int], dict[str, str]]:
    observed = {
        name: sha256(ROOT / name) for name in EXPECTED_DEPENDENCIES
    }
    if observed != EXPECTED_DEPENDENCIES:
        raise ValueError(f"tight-Hall dependency boundary changed: {observed}")
    manifest = json.loads(INPUT_MANIFEST.read_text(encoding="utf-8"))
    verification = json.loads(INPUT_VERIFICATION.read_text(encoding="utf-8"))
    block = manifest.get("classes", {}).get("K6_only", {})
    records = block.get("graphs", [])
    indices = [int(record["index"]) for record in records]
    if (
        manifest.get("schema") != "d6-current-exact-residue-v5"
        or manifest.get("status") != "COMPLETE_EXACT_FILTER_UNION"
        or block.get("count") != EXPECTED_INPUT
        or block.get("indices") != indices
        or block.get("indices_sha256") != EXPECTED_INPUT_SHA256
        or stable_hash(indices) != EXPECTED_INPUT_SHA256
        or len(records) != EXPECTED_INPUT
    ):
        raise ValueError("v5 K6-only input boundary changed")
    if (
        verification.get("schema") != "d6-current-exact-residue-v5-verification-v1"
        or verification.get("status") != "PASS"
        or verification.get("manifest", {}).get("sha256")
        != EXPECTED_DEPENDENCIES[INPUT_MANIFEST.name]
        or verification.get("counts", {}).get("K6_only") != EXPECTED_INPUT
    ):
        raise ValueError("v5 independent verification boundary changed")
    for record in records:
        adjacency = tuple(map(int, record["adjacency"]))
        if find_clique_mask(adjacency, 7):
            raise ValueError("K6-only production input contains a K7")
    return records, indices, observed


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
            raise AssertionError("singleton rank lower is neither zero nor one")
    return answer


def impose_singleton_status(
    groups: Sequence[hall.Group],
    zero: dict[int, tuple[int, int]],
    empty_group: int | None,
) -> tuple[hall.Group, ...]:
    answer = []
    for group_index, group in enumerate(groups):
        if group_index not in zero:
            answer.append(group)
            continue
        vertex, mask = zero[group_index]
        omit = hall.Choice(0, 0, (), "omit")
        choices = (omit,)
        if group_index != empty_group:
            choices += (hall.Choice(1, mask, (vertex,), "forced_nonempty"),)
        answer.append(hall.Group(group.label, choices, True))
    return tuple(answer)


def tight_states(
    groups: Sequence[hall.Group],
) -> dict[int, tuple[tuple[str, hall.Choice], ...]]:
    """DP all Hall-tight coordinate unions and one attaining selection."""

    states: dict[int, tuple[int, tuple[tuple[str, hall.Choice], ...]]] = {0: (0, ())}
    for group in groups:
        following: dict[int, tuple[int, tuple[tuple[str, hall.Choice], ...]]] = {}
        for coordinates, (rank, witness) in states.items():
            for choice in group.choices:
                union = coordinates | choice.coordinates
                total = rank + choice.rank
                if total > union.bit_count():
                    raise AssertionError("orientation advertised as Hall-passing")
                selected = witness
                if choice.selected:
                    selected += ((group.label, choice),)
                old = following.get(union)
                if old is None or total > old[0]:
                    following[union] = (total, selected)
        states = following
    return {
        union: witness
        for union, (rank, witness) in states.items()
        if rank == union.bit_count()
    }


def support_obstruction(
    graph: Sequence[int],
    absolute: Sequence[int],
    defects: dict[int, int],
    removed: int,
) -> dict | None:
    reduced = {vertex: defects[vertex] & ~removed for vertex in absolute}
    for row, vertex in enumerate(absolute):
        if reduced[vertex] == 0:
            return {
                "kind": "nonempty_diagonal_lost_all_support",
                "vertices": [vertex],
            }
        for column in range(row):
            other = absolute[column]
            if graph[row] & (1 << column) and not (
                reduced[vertex] & reduced[other]
            ):
                return {
                    "kind": "required_nonzero_gram_entry_lost_common_support",
                    "vertices": [other, vertex],
                }
    singleton: dict[int, list[int]] = {}
    for vertex in absolute:
        mask = reduced[vertex]
        if mask.bit_count() == 1:
            singleton.setdefault(mask, []).append(vertex)
    for mask, repeated in singleton.items():
        if len(repeated) >= 2:
            return {
                "kind": "same_line_singleton_collision",
                "vertices": repeated[:2],
                "coordinate": vertices(mask)[0],
            }
    return None


def orientation_state(
    adjacency: Sequence[int],
    instance,
    z0: int,
    component,
    sign_a: str,
    sign_b: str,
    chosen_empty: int | None,
    inertia: InertiaCache,
    forcing: ZeroForcingSolver,
    subset_rank: arbitrary.PositiveSubsetRank,
) -> tuple[str, dict | None]:
    graph_a, absolute_a = induced_required_graph(
        adjacency, instance, component.side_a
    )
    graph_b, absolute_b = induced_required_graph(
        adjacency, instance, component.side_b
    )
    side_a = ranks.side_rank_lower(
        graph_a, absolute_a, sign_a, inertia, forcing
    )
    side_b = ranks.side_rank_lower(
        graph_b, absolute_b, sign_b, inertia, forcing
    )
    groups_a, _ = arbitrary.make_side_groups(instance, "A", side_a, subset_rank)
    groups_b, _ = arbitrary.make_side_groups(instance, "B", side_b, subset_rank)
    negative_name, negative_groups, negative_side = (
        ("A", groups_a, side_a)
        if sign_a == "negative"
        else ("B", groups_b, side_b)
    )
    zero = singleton_data(instance, negative_side)
    empty_group = None
    if chosen_empty is not None:
        empty_group = next(
            (group for group, (vertex, _) in zero.items() if vertex == chosen_empty),
            None,
        )
        if empty_group is None:
            raise ValueError("chosen empty vertex is not eligible")
    if negative_name == "A":
        groups_a = impose_singleton_status(groups_a, zero, empty_group)
    else:
        groups_b = impose_singleton_status(groups_b, zero, empty_group)
    empty_label = (
        None if empty_group is None else f"{negative_name}:F{empty_group}"
    )
    zgroups = tuple(
        hall.singleton_group(
            f"Z0:{instance.outside[local]}",
            instance.outside[local],
            instance.defects[local],
        )
        for local in vertices(z0)
    )
    groups = groups_a + groups_b + zgroups
    parent = hall.check_groups(groups)
    if not parent.passed:
        return "PARENT_HALL_FAIL", parent.first_failure

    span_data = {}
    for prefix, side in (("A", side_a), ("B", side_b)):
        for number, detail in enumerate(side.components):
            label = f"{prefix}:F{number}"
            if label != empty_label:
                span_data[label] = (
                    tuple(detail["adjacency_rows"]),
                    tuple(detail["vertices"]),
                )
    defects = {
        absolute: instance.defects[local]
        for local, absolute in enumerate(instance.outside)
    }
    for target, (graph, absolute) in span_data.items():
        other = tuple(group for group in groups if group.label != target)
        for removed, tight in sorted(
            tight_states(other).items(),
            key=lambda item: (item[0].bit_count(), item[0]),
            reverse=True,
        ):
            obstruction = support_obstruction(
                graph, absolute, defects, removed
            )
            if obstruction is None:
                continue
            return "TIGHT_SUPPORT_FAIL", {
                "target_group": target,
                "target_vertices": list(absolute),
                "removed_coordinates": vertices(removed),
                "tight_rank": removed.bit_count(),
                "tight_groups": [
                    {
                        "group": label,
                        "rank": choice.rank,
                        "vertices": list(choice.selected),
                        "coordinates": vertices(choice.coordinates),
                    }
                    for label, choice in tight
                ],
                "obstruction": obstruction,
            }
    return "PASS", None


def component_options(
    adjacency: Sequence[int], instance, z0: int, component,
    inertia: InertiaCache, forcing: ZeroForcingSolver,
    subset_rank: arbitrary.PositiveSubsetRank,
) -> tuple[tuple[int | None, ...], list[dict], Counter[str]]:
    choices: set[int | None] = set()
    branches = []
    counts: Counter[str] = Counter()
    graph_a, absolute_a = induced_required_graph(
        adjacency, instance, component.side_a
    )
    graph_b, absolute_b = induced_required_graph(
        adjacency, instance, component.side_b
    )
    for sign_a, sign_b in (
        ("positive", "negative"),
        ("negative", "positive"),
    ):
        side_a = ranks.side_rank_lower(
            graph_a, absolute_a, sign_a, inertia, forcing
        )
        side_b = ranks.side_rank_lower(
            graph_b, absolute_b, sign_b, inertia, forcing
        )
        negative_side = side_a if sign_a == "negative" else side_b
        possible = tuple(vertex for vertex, _ in singleton_data(instance, negative_side).values())
        for chosen in (None, *possible):
            status, certificate = orientation_state(
                adjacency, instance, z0, component, sign_a, sign_b, chosen,
                inertia, forcing, subset_rank,
            )
            counts[f"generic_state_{status}"] += 1
            branches.append({
                "sign_A": sign_a,
                "chosen_empty": chosen,
                "status": status,
                "certificate": certificate,
            })
            if status == "PASS":
                choices.add(chosen)
    light = z0 | component.component
    lightlike = (
        light.bit_count() <= COORDINATES
        and support_matching(light, instance.defects) is not None
    )
    if lightlike:
        choices.add(None)
        counts["lightlike_pass"] += 1
    ordered = tuple(sorted(item for item in choices if item is not None))
    return (
        ((None, *ordered) if None in choices else ordered),
        branches,
        counts,
    )


def extend_states(
    states: Iterable[tuple[int, ...]],
    options: Sequence[int | None],
    adjacency: Sequence[int],
) -> set[tuple[int, ...]]:
    following = set()
    for state in states:
        for option in options:
            if option is None:
                following.add(tuple(state))
                continue
            chosen = tuple(sorted((*state, int(option))))
            if len(chosen) > GLOBAL_EMPTY_BUDGET:
                continue
            if len(chosen) == 2 and adjacency[chosen[0]] & (1 << chosen[1]):
                continue
            following.add(chosen)
    return following


def solve_seed(adjacency: Sequence[int], instance) -> dict:
    inertia = InertiaCache()
    forcing = ZeroForcingSolver()
    subset_rank = arbitrary.PositiveSubsetRank(inertia, forcing)
    counts: Counter[str] = Counter()
    failure_rows = []
    for z0 in ranks.z0_subsets(instance.eligible_z0_mask):
        counts["z0_considered"] += 1
        z0_vertices = [instance.outside[local] for local in vertices(z0)]
        if support_matching(z0, instance.defects) is None:
            counts["z0_unmatchable"] += 1
            failure_rows.append({"Z0": z0_vertices, "matchable": False})
            continue
        states: set[tuple[int, ...]] = {()}
        failed = None
        for component in lorentz_components(instance, z0):
            if not component.bipartite:
                counts["nonbipartite_components_ignored"] += 1
                continue
            counts["bipartite_components"] += 1
            options, branches, local = component_options(
                adjacency, instance, z0, component,
                inertia, forcing, subset_rank,
            )
            counts.update(local)
            if not options:
                states = set()
            else:
                states = extend_states(states, options, adjacency)
            if not states:
                failed = {
                    "component": [
                        instance.outside[local]
                        for local in vertices(component.component)
                    ],
                    "options": [
                        option for option in options if option is not None
                    ],
                    "all_nonempty_option": None in options,
                    "lightlike_passed": (
                        z0 | component.component
                    ).bit_count() <= COORDINATES and support_matching(
                        z0 | component.component, instance.defects
                    ) is not None,
                    "branches": branches,
                }
                break
        if states:
            counts["z0_feasible"] += 1
            return {
                "feasible": True,
                "witness": {
                    "Z0": z0_vertices,
                    "empty_state": list(min(states, key=lambda x: (len(x), x))),
                },
                "counts": dict(counts),
                "failure_rows": None,
            }
        failure_rows.append({
            "Z0": z0_vertices,
            "matchable": True,
            "failure": failed,
        })
    return {
        "feasible": False,
        "witness": None,
        "counts": dict(counts),
        "failure_rows": failure_rows,
    }


def evaluate_record(record: dict) -> dict:
    adjacency = tuple(map(int, record["adjacency"]))
    result_counts: Counter[str] = Counter()
    seeds_checked = 0
    for seed_mask in clique_masks(adjacency, COORDINATES):
        seeds_checked += 1
        instance = build_instance(adjacency, vertices(seed_mask))
        decision = solve_seed(adjacency, instance)
        result_counts.update(decision["counts"])
        if not decision["feasible"]:
            return {
                "index": int(record["index"]),
                "rejected": True,
                "seeds_checked": seeds_checked,
                "first_impossible_seed_mask": seed_mask,
                "first_impossible_seed": list(instance.seed),
                "counts": dict(result_counts),
                "certificate": {
                    "seed_mask": seed_mask,
                    "seed": list(instance.seed),
                    "rows": decision["failure_rows"],
                },
            }
    return {
        "index": int(record["index"]),
        "rejected": False,
        "seeds_checked": seeds_checked,
        "first_impossible_seed_mask": None,
        "first_impossible_seed": None,
        "counts": dict(result_counts),
        "certificate": None,
    }


def positive_control() -> dict:
    adjacency = lower_bound_18_graph()
    rows = []
    for seed_mask in clique_masks(adjacency, COORDINATES):
        instance = build_instance(adjacency, vertices(seed_mask))
        decision = solve_seed(adjacency, instance)
        rows.append({
            "seed_mask": seed_mask,
            "seed": list(instance.seed),
            "passed": decision["feasible"],
            "witness": decision["witness"],
        })
    if len(rows) != 32 or not all(row["passed"] for row in rows):
        raise AssertionError("known realizable 18-point control rejected")
    return {"passed": True, "K6_seeds": len(rows), "seed_results": rows}


def run(workers: int, output: Path, certificates: Path) -> dict:
    records, indices, dependencies = load_input()
    started = time.monotonic()
    if workers == 1:
        results = [evaluate_record(record) for record in records]
    else:
        with ProcessPoolExecutor(
            max_workers=workers, initializer=worker_initializer
        ) as pool:
            results = list(pool.map(evaluate_record, records, chunksize=1))
    rejected = [row["index"] for row in results if row["rejected"]]
    residue = [row["index"] for row in results if not row["rejected"]]
    counts: Counter[str] = Counter()
    for row in results:
        counts.update(row["counts"])
    archive = {
        "schema": CERTIFICATE_SCHEMA,
        "production_source_sha256": sha256(Path(__file__)),
        "ordered_input_indices_sha256": EXPECTED_INPUT_SHA256,
        "rejected_graphs": [
            {"index": row["index"], **row["certificate"]}
            for row in results if row["rejected"]
        ],
    }
    atomic_json(certificates, archive)
    report = {
        "schema": REPORT_SCHEMA,
        "status": "COMPLETE",
        "description": "exact tight-Hall equality support propagation",
        "semantics": {
            "candidate_nonedges": "unconstrained_and_may_be_unit",
            "allowed_supports": "upper_bounds_coordinates_may_be_zero",
            "distinct_points": "required",
            "nonbipartite_components": "ignored_conservatively",
            "arithmetic": "exact_integer_bitmask_graph_logic",
        },
        "production_source_sha256": sha256(Path(__file__)),
        "dependencies": dependencies,
        "ordered_input_indices": indices,
        "ordered_input_indices_sha256": EXPECTED_INPUT_SHA256,
        "input_graphs": len(indices),
        "graphs_rejected": len(rejected),
        "graphs_surviving": len(residue),
        "rejected_indices": rejected,
        "rejected_indices_sha256": stable_hash(rejected),
        "ordered_residue_indices": residue,
        "ordered_residue_indices_sha256": stable_hash(residue),
        "totals": {
            "seeds_checked": sum(row["seeds_checked"] for row in results),
            **dict(counts),
        },
        "certificate_archive": {
            "path": certificates.name,
            "sha256": sha256(certificates),
            "rejected_graphs": len(archive["rejected_graphs"]),
        },
        "positive_18_control": positive_control(),
        "runtime": {
            "workers": workers,
            "wall_seconds": time.monotonic() - started,
            "python": sys.version,
            "platform": platform.platform(),
        },
        "graph_results": [
            {key: value for key, value in row.items() if key != "certificate"}
            for row in results
        ],
        "nonclaims": [
            "a survivor is not a realization",
            "this incremental layer does not settle the K6 or dimension-six problem",
        ],
    }
    atomic_json(output, report)
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    parser.add_argument(
        "--output", type=Path, default=ROOT / "d6_k6_tight_hall_support_report.json"
    )
    parser.add_argument(
        "--certificates", type=Path,
        default=ROOT / "d6_k6_tight_hall_support_certificates.json",
    )
    args = parser.parse_args()
    report = run(args.workers, args.output, args.certificates)
    print(json.dumps({
        "status": report["status"],
        "input": report["input_graphs"],
        "rejected": report["graphs_rejected"],
        "surviving": report["graphs_surviving"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
