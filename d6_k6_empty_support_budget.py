#!/usr/bin/env python3
"""Exact empty-defect support budget on the pinned 822-graph K6 residue.

This layer strengthens the generic bordered-K6 Hall systems.  In one generic
Lorentz component, the positive-norm side has no zero defect vector.  On the
negative-norm side only a singleton required-edge block can have zero defect,
and at most one such singleton can be zero because all factors on that side
share one projective Lorentz line.  Across all components there are only the
two empty-defect points r_+ and r_-, so at most two components may require an
empty singleton.

The implementation is exact integer/bit-mask logic.  It deliberately ignores
non-bipartite Lorentz components here; doing so only weakens this new screen,
and the pinned parent boundary has already checked their earlier constraints.
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
import traceback
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import d6_k6_arbitrary_subset_hall as parent
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
    support_matching,
    vertices,
)
from d6_k6_normal_inertia import InertiaCache
from test_d6_k6_lorentz import lower_bound_18_graph


ROOT = Path(__file__).resolve().parent
EXPECTED_INPUT = 822
EXPECTED_INPUT_SHA256 = (
    "cf94aac41eba28904ee6254f9c50996e73d38dcc15316fb6a411b29ff5cd05d8"
)
PARENT_REPORT = "d6_k6_arbitrary_subset_hall_report.json"
GLOBAL_EMPTY_BUDGET = 2
REPORT_SCHEMA = "d6-k6-empty-support-budget-v1"
CERTIFICATE_SCHEMA = "d6-k6-empty-support-budget-certificates-v1"
CHECKPOINT_SCHEMA = "d6-k6-empty-support-budget-checkpoint-v1"
ABORT_SCHEMA = "d6-k6-empty-support-budget-infrastructure-abort-v1"
PILOT_SCHEMA = "d6-k6-empty-support-budget-pilot-v2"
DEFAULT_WORKERS = 11

# This layer starts from the *independently verified* 822 survivors of the
# arbitrary-subset Hall production result.  Pinning all four parent artifacts
# prevents a same-name replacement from silently changing this theorem input.
EXPECTED_DEPENDENCIES = {
    "d6_k6_arbitrary_subset_hall.py": (
        "cadc0048fcc7a77299bc2950e2eca57bfed87d54ed6ffd59886d155b2c9a88ff"
    ),
    "d6_k6_arbitrary_subset_hall_report.json": (
        "bcc4baecf28a95443a5b13e87a796cbba0a97cd42fa76573119bfd81ef3735fc"
    ),
    "d6_k6_arbitrary_subset_hall_certificates.json": (
        "b7fb04f769ddaf5d36630372f5134f2e4a32e572f620b26059667e9cf65f4a23"
    ),
    "d6_k6_arbitrary_subset_hall_verification.json": (
        "e6042edfff3683043903f9be544322379e85673207e35ee87349637ff4c73058"
    ),
    "d6_k6_psd_z_hereditary.py": (
        "7b869895cb9b4f7f3bbf1a1d1b11ec71507b6355c8ec0b7d7f8eaec2a51b180d"
    ),
    "d6_k6_psd_zmatrix.py": (
        "24c2b3abfce5074560183d17f77778f96298bef6e38ece5bc937ccb2a4cba736"
    ),
    "d6_k6_bipartite_rank_reference.py": (
        "e760fef74c423a89f1e1f9d08909fa1223e0389900cf3d45b0140ee296dea922"
    ),
    "d6_k6_lorentz_reference.py": (
        "ef6a0877a5987119c92ec1bd9271a6a47774fc9e4992c10e1f28b6d5373592c2"
    ),
    "d6_k6_normal_inertia.py": (
        "2d00ab40eb97aa79ffec3b8134b83c4fd5f703cbd691f094a25bf3ff9f3ee023"
    ),
    "test_d6_k6_lorentz.py": (
        "f4fc8544a4311f742267808e6c785f82139accb9da99b496d206283b52831a47"
    ),
}

SEMANTIC_CONFIG = {
    "global_empty_budget": GLOBAL_EMPTY_BUDGET,
    "candidate_nonedges": "unconstrained_and_may_be_unit",
    "allowed_defects": "upper_bounds_coordinates_may_be_zero",
    "nonbipartite_components": "ignored_conservatively",
    "generic_component_empty_rule": "at_most_one_negative_singleton",
    "arithmetic": "exact_integer_bitmask_only",
}


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
class EmptySingleton:
    group_index: int
    absolute_vertex: int
    allowed_coordinates: int


@dataclass(frozen=True)
class OrientationBudget:
    possible: bool
    minimum_empty: int | None
    negative_singletons: tuple[int, ...]
    chosen_empty: int | None


@dataclass(frozen=True)
class ComponentBudget:
    possible: bool
    minimum_empty: int | None
    alternatives: tuple[dict, ...]


def global_component_costs_pass(
    costs: Sequence[int | None], budget: int = GLOBAL_EMPTY_BUDGET
) -> bool:
    """Whether independent component alternatives fit the global empty budget."""

    return all(cost is not None for cost in costs) and sum(
        int(cost) for cost in costs if cost is not None
    ) <= budget


def zero_capable_singletons(
    instance: K6LorentzInstance,
    groups: Sequence[hall.Group],
    side,
) -> tuple[EmptySingleton, ...]:
    """Return negative-side singleton Gram blocks whose rank lower bound is 0."""

    if len(groups) != len(side.components):
        raise ValueError("side group/component count mismatch")
    local_by_absolute = {
        absolute: local for local, absolute in enumerate(instance.outside)
    }
    answer = []
    for group_index, (group, detail) in enumerate(zip(groups, side.components)):
        absolute = tuple(detail["vertices"])
        rank_lower = int(detail["componentwise_fused_rank_lower"])
        if len(absolute) == 1 and rank_lower == 0:
            vertex = absolute[0]
            answer.append(EmptySingleton(
                group_index,
                vertex,
                instance.defects[local_by_absolute[vertex]],
            ))
        elif len(absolute) == 1 and rank_lower != 1:
            raise AssertionError("a singleton Gram block must have rank lower 0 or 1")
    return tuple(answer)


def force_singleton_status(
    groups: Sequence[hall.Group],
    zero_capable: Sequence[EmptySingleton],
    chosen_empty_group: int | None,
) -> tuple[hall.Group, ...]:
    """Force all but the selected zero-capable singleton to be nonzero rank 1."""

    by_group = {item.group_index: item for item in zero_capable}
    if chosen_empty_group is not None and chosen_empty_group not in by_group:
        raise ValueError("chosen empty group is not zero-capable")
    answer = []
    for group_index, group in enumerate(groups):
        item = by_group.get(group_index)
        if item is None:
            answer.append(group)
            continue
        omitted = hall.Choice(0, 0, (), "omit")
        if group_index == chosen_empty_group:
            # A known zero vector adds neither rank nor actual support.
            answer.append(hall.Group(group.label, (omitted,), True))
        else:
            nonempty = hall.Choice(
                1,
                item.allowed_coordinates,
                (item.absolute_vertex,),
                "forced_nonempty_negative_singleton",
            )
            answer.append(hall.Group(group.label, (omitted, nonempty), True))
    return tuple(answer)


def orientation_minimum_empty(
    instance: K6LorentzInstance,
    z0: int,
    side_a,
    side_b,
    subset_rank: parent.PositiveSubsetRank,
    negative_side: str,
) -> OrientationBudget:
    """Minimum empty singletons needed by one fixed generic orientation."""

    groups_a, _ = parent.make_side_groups(instance, "A", side_a, subset_rank)
    groups_b, _ = parent.make_side_groups(instance, "B", side_b, subset_rank)
    if negative_side == "A":
        negative_groups, negative_rank = groups_a, side_a
    elif negative_side == "B":
        negative_groups, negative_rank = groups_b, side_b
    else:
        raise ValueError("negative_side must be A or B")
    zero_capable = zero_capable_singletons(
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
    # Cost zero first, then each possible unique empty point on this line.
    choices = (None,) + tuple(item.group_index for item in zero_capable)
    for chosen in choices:
        if negative_side == "A":
            test_a = force_singleton_status(groups_a, zero_capable, chosen)
            test_b = groups_b
        else:
            test_a = groups_a
            test_b = force_singleton_status(groups_b, zero_capable, chosen)
        result = hall.check_groups(test_a + test_b + z0_groups)
        if result.passed:
            absolute = None
            if chosen is not None:
                absolute = next(
                    item.absolute_vertex
                    for item in zero_capable
                    if item.group_index == chosen
                )
            return OrientationBudget(
                True,
                int(chosen is not None),
                tuple(item.absolute_vertex for item in zero_capable),
                absolute,
            )
    return OrientationBudget(
        False,
        None,
        tuple(item.absolute_vertex for item in zero_capable),
        None,
    )


def component_minimum_empty(
    adjacency: Sequence[int],
    instance: K6LorentzInstance,
    z0: int,
    component,
    inertia: InertiaCache,
    zero_forcing: ZeroForcingSolver,
    subset_rank: parent.PositiveSubsetRank,
) -> ComponentBudget:
    """Minimize empty-support cost over two generic orientations and lightlike."""

    graph_a, absolute_a = induced_required_graph(
        adjacency, instance, component.side_a
    )
    graph_b, absolute_b = induced_required_graph(
        adjacency, instance, component.side_b
    )
    alternatives = []
    for name, sign_a, sign_b, negative_side in (
        ("A_positive", "positive", "negative", "B"),
        ("A_negative", "negative", "positive", "A"),
    ):
        side_a = ranks.side_rank_lower(
            graph_a, absolute_a, sign_a, inertia, zero_forcing
        )
        side_b = ranks.side_rank_lower(
            graph_b, absolute_b, sign_b, inertia, zero_forcing
        )
        decision = orientation_minimum_empty(
            instance, z0, side_a, side_b, subset_rank, negative_side
        )
        alternatives.append({
            "kind": name,
            "possible": decision.possible,
            "minimum_empty": decision.minimum_empty,
            "negative_singletons": list(decision.negative_singletons),
            "chosen_empty": decision.chosen_empty,
        })

    light = z0 | component.component
    light_passed = (
        light.bit_count() <= COORDINATES
        and support_matching(light, instance.defects) is not None
    )
    alternatives.append({
        "kind": "lightlike",
        "possible": light_passed,
        "minimum_empty": 0 if light_passed else None,
        "orthonormal_vectors": light.bit_count(),
    })
    costs = [
        row["minimum_empty"]
        for row in alternatives
        if row["possible"]
    ]
    return ComponentBudget(
        bool(costs), min(costs) if costs else None, tuple(alternatives)
    )


def solve_seed(adjacency: Sequence[int], instance: K6LorentzInstance) -> dict:
    """Exhaust Z0 and apply the global two-empty-point budget conservatively."""

    inertia = InertiaCache()
    zero_forcing = ZeroForcingSolver()
    subset_rank = parent.PositiveSubsetRank(inertia, zero_forcing)
    rows = []
    counts = {"z0": 0, "components": 0, "orientations": 0}
    for z0 in ranks.z0_subsets(instance.eligible_z0_mask):
        counts["z0"] += 1
        if support_matching(z0, instance.defects) is None:
            rows.append({
                "Z0": parent.frozen.absolute_vertices(instance, z0),
                "matchable": False,
                "terminal": "Z0_ALLOWED_SUPPORT_UNMATCHABLE",
            })
            continue
        total_empty = 0
        component_costs: list[int | None] = []
        component_rows = []
        possible = True
        terminal = None
        for component in lorentz_components(instance, z0):
            if not component.bipartite:
                # Conservatively ignore this component in the new budget.
                continue
            counts["components"] += 1
            counts["orientations"] += 2
            decision = component_minimum_empty(
                adjacency, instance, z0, component,
                inertia, zero_forcing, subset_rank,
            )
            component_rows.append({
                "vertices": parent.frozen.absolute_vertices(
                    instance, component.component
                ),
                "possible": decision.possible,
                "minimum_empty": decision.minimum_empty,
                "alternatives": list(decision.alternatives),
            })
            if not decision.possible:
                component_costs.append(None)
                possible = False
                terminal = "COMPONENT_IMPOSSIBLE"
                break
            assert decision.minimum_empty is not None
            component_costs.append(decision.minimum_empty)
            total_empty += decision.minimum_empty
            if not global_component_costs_pass(component_costs):
                possible = False
                terminal = "GLOBAL_EMPTY_BUDGET_EXCEEDED"
                break
        if possible:
            return {
                "feasible": True,
                "counts": counts,
                "chosen_Z0": parent.frozen.absolute_vertices(instance, z0),
                "minimum_empty_used": total_empty,
                "certificate_rows": None,
            }
        rows.append({
            "Z0": parent.frozen.absolute_vertices(instance, z0),
            "matchable": True,
            "terminal": terminal,
            "minimum_empty_before_failure": total_empty,
            "global_empty_budget": GLOBAL_EMPTY_BUDGET,
            "components": component_rows,
        })
    return {
        "feasible": False,
        "counts": counts,
        "chosen_Z0": None,
        "minimum_empty_used": None,
        "certificate_rows": rows,
    }


def evaluate_graph(adjacency: Sequence[int]) -> dict:
    totals = {"z0": 0, "components": 0, "orientations": 0}
    seeds = 0
    for seed_mask in clique_masks(adjacency, COORDINATES):
        seeds += 1
        instance = build_instance(adjacency, vertices(seed_mask))
        decision = solve_seed(adjacency, instance)
        for name in totals:
            totals[name] += decision["counts"][name]
        if not decision["feasible"]:
            return {
                "rejected": True,
                "seeds_checked": seeds,
                **totals,
                "first_impossible_seed": list(instance.seed),
                "certificate": {
                    "seed": list(instance.seed),
                    "outside": list(instance.outside),
                    "allowed_defects": {
                        str(instance.outside[local]): [
                            instance.seed[coordinate]
                            for coordinate in vertices(mask)
                        ]
                        for local, mask in enumerate(instance.defects)
                    },
                    "choices": decision["certificate_rows"],
                },
            }
    return {
        "rejected": False,
        "seeds_checked": seeds,
        **totals,
        "first_impossible_seed": None,
        "certificate": None,
    }


def dependency_hashes() -> dict[str, str]:
    observed = {
        name: sha256(ROOT / name) for name in EXPECTED_DEPENDENCIES
    }
    if observed != EXPECTED_DEPENDENCIES:
        raise ValueError(
            "K6 empty-support dependency hash mismatch: "
            f"observed={observed}"
        )
    return observed


def verify_parent_input() -> tuple[dict[str, str], list[dict], list[int], list[int]]:
    """Reconstruct and bind the ordered 822-graph parent complement."""

    observed = dependency_hashes()
    _, records, indices = parent.verify_inputs()
    report = json.loads((ROOT / PARENT_REPORT).read_text(encoding="utf-8"))
    archive = json.loads(
        (ROOT / "d6_k6_arbitrary_subset_hall_certificates.json").read_text(
            encoding="utf-8"
        )
    )
    verification = json.loads(
        (ROOT / "d6_k6_arbitrary_subset_hall_verification.json").read_text(
            encoding="utf-8"
        )
    )
    if (
        report.get("schema") != "d6-k6-arbitrary-subset-hall-v1"
        or report.get("status") != "COMPLETE"
        or report.get("input_graphs") != 831
    ):
        raise ValueError("arbitrary-subset parent report is not complete")
    if (
        archive.get("schema")
        != "d6-k6-arbitrary-subset-hall-certificates-v1"
        or archive.get("status") != "COMPLETE"
    ):
        raise ValueError("arbitrary-subset parent archive is not complete")
    if (
        verification.get("schema")
        != "d6-k6-arbitrary-subset-hall-verification-v1"
        or verification.get("status") != "PASS"
    ):
        raise ValueError("arbitrary-subset parent verification is not PASS")
    rejected = list(report["rejected_indices"])
    if (
        len(rejected) != 9
        or rejected != archive.get("rejected_indices")
        or rejected != verification.get("rejected_indices")
    ):
        raise ValueError("parent rejection manifests disagree")
    if (
        report.get("certificate_archive", {}).get("sha256")
        != EXPECTED_DEPENDENCIES[
            "d6_k6_arbitrary_subset_hall_certificates.json"
        ]
    ):
        raise ValueError("parent report does not bind the pinned archive")
    rejected_set = set(rejected)
    survivors = [index for index in indices if index not in rejected_set]
    if (
        len(survivors) != EXPECTED_INPUT
        or stable_hash(survivors) != EXPECTED_INPUT_SHA256
    ):
        raise ValueError("822-graph empty-support input boundary changed")
    by_index = {record["index"]: record for record in records}
    return observed, [by_index[index] for index in survivors], survivors, rejected


def semantic_config_sha256() -> str:
    return stable_hash(SEMANTIC_CONFIG)


def checkpoint_payload(
    source: str,
    indices: Sequence[int],
    dependencies: dict[str, str],
    completed: Sequence[dict],
    status: str = "RUNNING",
) -> dict:
    return {
        "schema": CHECKPOINT_SCHEMA,
        "status": status,
        "production_source_sha256": source,
        "semantic_config_sha256": semantic_config_sha256(),
        "dependency_hashes": dependencies,
        "input_graphs": len(indices),
        "input_indices_sha256": stable_hash(indices),
        "completed": list(completed),
    }


def validate_checkpoint_payload(
    saved: object,
    source: str,
    indices: Sequence[int],
    dependencies: dict[str, str],
) -> list[dict]:
    if not isinstance(saved, dict):
        raise ValueError("checkpoint is not a JSON object")
    if (
        saved.get("schema") != CHECKPOINT_SCHEMA
        or saved.get("status") not in {"RUNNING", "COMPLETE"}
        or saved.get("production_source_sha256") != source
        or saved.get("semantic_config_sha256") != semantic_config_sha256()
        or saved.get("dependency_hashes") != dependencies
        or saved.get("input_graphs") != len(indices)
        or saved.get("input_indices_sha256") != stable_hash(indices)
    ):
        raise ValueError("empty-support checkpoint boundary mismatch")
    completed = saved.get("completed")
    if not isinstance(completed, list) or len(completed) > len(indices):
        raise ValueError("checkpoint completed payload is malformed")
    if [row.get("index") for row in completed] != list(indices[:len(completed)]):
        raise ValueError("checkpoint is not an ordered input prefix")
    for row in completed:
        if not isinstance(row.get("decision"), dict):
            raise ValueError("checkpoint decision row is malformed")
    return completed


def load_checkpoint(
    path: Path,
    source: str,
    indices: Sequence[int],
    dependencies: dict[str, str],
) -> list[dict]:
    return validate_checkpoint_payload(
        json.loads(path.read_text(encoding="utf-8")),
        source,
        indices,
        dependencies,
    )


def worker_initializer() -> None:
    # All kernels are exact combinatorics.  Avoid accidental nested BLAS pools
    # in imported controls while the outer process pool owns the parallelism.
    for name in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[name] = "1"


def evaluate_record(record: dict) -> dict:
    return {
        "index": record["index"],
        "decision": evaluate_graph(record["adjacency"]),
    }


def positive_seed_controls() -> list[dict]:
    adjacency = lower_bound_18_graph()
    rows = []
    for seed_mask in clique_masks(adjacency, COORDINATES):
        instance = build_instance(adjacency, vertices(seed_mask))
        decision = solve_seed(adjacency, instance)
        rows.append({
            "seed": list(instance.seed),
            "passed": bool(decision["feasible"]),
            "z0": decision["counts"]["z0"],
            "components": decision["counts"]["components"],
            "orientations": decision["counts"]["orientations"],
        })
    if len(rows) != 32 or not all(row["passed"] for row in rows):
        raise AssertionError(
            "known realizable 18-point control did not pass all 32 K6 seeds"
        )
    return rows


def synthetic_controls() -> dict:
    observed = {
        "two_components_cost_one": global_component_costs_pass([1, 1]),
        "three_components_cost_one": global_component_costs_pass([1, 1, 1]),
        "impossible_component": global_component_costs_pass([0, None]),
    }
    expected = {
        "two_components_cost_one": True,
        "three_components_cost_one": False,
        "impossible_component": False,
    }
    if observed != expected:
        raise AssertionError(f"synthetic empty-budget controls failed: {observed}")
    return observed


def git_provenance(committed_sources: dict[str, str]) -> dict:
    """Require each theorem-level source to equal its blob at HEAD."""

    def git_text(*args: str) -> str:
        result = subprocess.run(
            ("git", *args), cwd=ROOT, check=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        )
        return result.stdout.strip()

    head = git_text("rev-parse", "HEAD")
    branch = git_text("branch", "--show-current")
    for name, expected in committed_sources.items():
        result = subprocess.run(
            ("git", "show", f"HEAD:{name}"), cwd=ROOT, check=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
        observed = hashlib.sha256(result.stdout).hexdigest()
        if observed != expected:
            raise ValueError(
                f"{name} is not source-bound to HEAD {head}: "
                f"blob={observed}, working={expected}"
            )
    return {
        "commit": head,
        "branch": branch,
        "committed_source_sha256": committed_sources,
    }


def stratified_selection(records: Sequence[dict], per_stratum: int = 5) -> list[dict]:
    """Four K6-count strata, five stable-hash choices from each."""

    scored = []
    for record in records:
        count = sum(1 for _ in clique_masks(record["adjacency"], COORDINATES))
        priority = hashlib.sha256(str(record["index"]).encode("ascii")).digest()
        scored.append((count, priority, record))
    scored.sort(key=lambda row: (row[0], row[1]))
    answer = []
    for stratum in range(4):
        lo = len(scored) * stratum // 4
        hi = len(scored) * (stratum + 1) // 4
        block = sorted(scored[lo:hi], key=lambda row: row[1])
        answer.extend(row[2] for row in block[:per_stratum])
    return answer


def pilot(output: Path) -> dict:
    observed, records, indices, inherited = verify_parent_input()
    selected = stratified_selection(records)
    by_index = {record["index"]: record for record in records}
    # Preserve the independently located first witness even if hash sampling
    # does not select it.
    if 202_556 not in {record["index"] for record in selected}:
        selected.append(by_index[202_556])
    started = time.monotonic()
    decisions = [
        {"index": record["index"], "decision": evaluate_graph(record["adjacency"])}
        for record in selected
    ]
    positive = positive_seed_controls()
    report = {
        "schema": PILOT_SCHEMA,
        "status": "PILOT_ONLY",
        "input_graphs": len(indices),
        "input_indices_sha256": stable_hash(indices),
        "inherited_parent_rejections": inherited,
        "selection": [row["index"] for row in decisions],
        "rejected_indices": [
            row["index"] for row in decisions if row["decision"]["rejected"]
        ],
        "graph_results": decisions,
        "positive_control": {
            "passed": True,
            "K6_seeds": len(positive),
            "seed_results": positive,
        },
        "synthetic_controls": synthetic_controls(),
        "runtime_seconds": time.monotonic() - started,
        "sources": {
            **observed,
            "d6_k6_empty_support_budget.py": sha256(Path(__file__)),
        },
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
        },
        "semantics": {
            "allowed_defects": "upper bounds; coordinates may be zero",
            "nonbipartite_components": "ignored conservatively by this layer",
            "survivor": "filter non-rejection only",
        },
    }
    atomic_json(output, report)
    return report


class InfrastructureAbort(RuntimeError):
    pass


def infrastructure_abort_payload(
    source: str,
    dependencies: dict[str, str],
    indices: Sequence[int],
    completed: Sequence[dict],
    checkpoint: Path,
    error: BaseException,
) -> dict:
    return {
        "schema": ABORT_SCHEMA,
        "status": "INFRA_ABORT",
        "meaning": "no mathematical rejection is inferred from this abort",
        "production_source_sha256": source,
        "semantic_config_sha256": semantic_config_sha256(),
        "dependency_hashes": dependencies,
        "input_indices_sha256": stable_hash(indices),
        "completed_ordered_prefix": len(completed),
        "completed_indices_sha256": stable_hash(
            [row["index"] for row in completed]
        ),
        "checkpoint": {
            "path": checkpoint.name,
            "sha256": sha256(checkpoint),
        },
        "exception": {
            "type": type(error).__name__,
            "message": str(error),
            "traceback": traceback.format_exc(),
        },
    }


def production_commands(
    output: Path,
    certificates: Path,
    checkpoint: Path,
    abort: Path,
    workers: int,
    checkpoint_every: int,
) -> dict[str, list[str]]:
    base = [
        "python3", Path(__file__).name, "--full",
        "--output", output.name,
        "--certificates", certificates.name,
        "--checkpoint", checkpoint.name,
        "--abort", abort.name,
        "--workers", str(workers),
        "--checkpoint-every", str(checkpoint_every),
    ]
    return {
        "initial": base,
        "resume": [*base, "--resume"],
        "independent_verification": [
            "python3", "verify_d6_k6_empty_support_budget.py",
            "--report", output.name,
            "--certificates", certificates.name,
            "--expected-report-sha256", "<sha256-of-report>",
            "--workers", str(workers),
            "--output", "d6_k6_empty_support_budget_verification.json",
        ],
    }


def run_production(
    output: Path,
    certificates: Path,
    checkpoint: Path,
    abort: Path,
    workers: int,
    checkpoint_every: int,
    resume: bool,
) -> dict:
    dependencies, records, indices, inherited = verify_parent_input()
    source = sha256(Path(__file__).resolve())
    verifier_source = sha256(ROOT / "verify_d6_k6_empty_support_budget.py")
    git = git_provenance({
        Path(__file__).name: source,
        "verify_d6_k6_empty_support_budget.py": verifier_source,
    })
    controls = synthetic_controls()
    positive = positive_seed_controls()

    if resume:
        if not checkpoint.exists():
            raise ValueError("--resume requested but checkpoint does not exist")
        completed = load_checkpoint(
            checkpoint, source, indices, dependencies
        )
    else:
        if checkpoint.exists():
            raise FileExistsError(
                f"checkpoint already exists: {checkpoint}; use --resume or a new path"
            )
        completed: list[dict] = []
        atomic_json(
            checkpoint,
            checkpoint_payload(source, indices, dependencies, completed),
        )

    pending = records[len(completed):]
    try:
        if workers == 1:
            iterator = map(evaluate_record, pending)
            for row in iterator:
                expected = indices[len(completed)]
                if row["index"] != expected:
                    raise RuntimeError(
                        f"worker ordering failure: got {row['index']}, expected {expected}"
                    )
                completed.append(row)
                if len(completed) % checkpoint_every == 0:
                    atomic_json(
                        checkpoint,
                        checkpoint_payload(
                            source, indices, dependencies, completed
                        ),
                    )
        else:
            with ProcessPoolExecutor(
                max_workers=workers, initializer=worker_initializer
            ) as pool:
                for row in pool.map(
                    evaluate_record, pending, chunksize=1
                ):
                    expected = indices[len(completed)]
                    if row["index"] != expected:
                        raise RuntimeError(
                            "worker ordering failure: "
                            f"got {row['index']}, expected {expected}"
                        )
                    completed.append(row)
                    if len(completed) % checkpoint_every == 0:
                        atomic_json(
                            checkpoint,
                            checkpoint_payload(
                                source, indices, dependencies, completed
                            ),
                        )
    except BaseException as error:
        atomic_json(
            checkpoint,
            checkpoint_payload(source, indices, dependencies, completed),
        )
        atomic_json(
            abort,
            infrastructure_abort_payload(
                source, dependencies, indices, completed, checkpoint, error
            ),
        )
        raise InfrastructureAbort(
            f"production infrastructure aborted after {len(completed)}/"
            f"{len(indices)} graphs; resume from {checkpoint}"
        ) from error

    if len(completed) != len(indices):
        raise InfrastructureAbort(
            f"ordered production ended at {len(completed)}/{len(indices)}"
        )

    rejected_rows = [
        row for row in completed if row["decision"]["rejected"]
    ]
    rejected_indices = [row["index"] for row in rejected_rows]
    residue_indices = [
        row["index"] for row in completed if not row["decision"]["rejected"]
    ]
    archive = {
        "schema": CERTIFICATE_SCHEMA,
        "status": "COMPLETE",
        "production_source_sha256": source,
        "semantic_config_sha256": semantic_config_sha256(),
        "dependency_hashes": dependencies,
        "input_graphs": len(indices),
        "input_indices_sha256": stable_hash(indices),
        "rejected_indices": rejected_indices,
        "certificates": [
            {
                "index": row["index"],
                "certificate": row["decision"]["certificate"],
            }
            for row in rejected_rows
        ],
    }
    atomic_json(certificates, archive)

    graph_results = []
    for row in completed:
        decision = dict(row["decision"])
        decision.pop("certificate")
        graph_results.append({"index": row["index"], "decision": decision})

    atomic_json(
        checkpoint,
        checkpoint_payload(
            source, indices, dependencies, completed, status="COMPLETE"
        ),
    )
    totals = {
        name: sum(row["decision"][name] for row in graph_results)
        for name in ("seeds_checked", "z0", "components", "orientations")
    }
    report = {
        "schema": REPORT_SCHEMA,
        "status": "COMPLETE",
        "description": (
            "Exact K6 empty-defect distinctness budget layered on the pinned "
            "822-graph arbitrary-subset Hall residue."
        ),
        "input_graphs": len(indices),
        "input_indices_sha256": stable_hash(indices),
        "ordered_input_indices": list(indices),
        "inherited_parent_rejections": inherited,
        "graphs_rejected": len(rejected_indices),
        "graphs_surviving": len(residue_indices),
        "rejected_indices": rejected_indices,
        "ordered_residue_indices": residue_indices,
        "ordered_residue_indices_sha256": stable_hash(residue_indices),
        **totals,
        "graph_results": graph_results,
        "certificate_archive": {
            "path": certificates.name,
            "sha256": sha256(certificates),
            "bytes": certificates.stat().st_size,
            "certificates": len(rejected_indices),
        },
        "checkpoint": {
            "path": checkpoint.name,
            "sha256": sha256(checkpoint),
            "completed_ordered_prefix": len(completed),
            "status": "COMPLETE",
        },
        "abort_artifact_on_failure": abort.name,
        "positive_control": {
            "passed": True,
            "K6_seeds": len(positive),
            "seed_results": positive,
        },
        "synthetic_controls": controls,
        "sources": {
            **dependencies,
            Path(__file__).name: source,
            "verify_d6_k6_empty_support_budget.py": verifier_source,
        },
        "git": git,
        "semantic_config": SEMANTIC_CONFIG,
        "semantic_config_sha256": semantic_config_sha256(),
        "semantics": {
            "candidate_nonedges": "unconstrained and may also be unit",
            "allowed_defects": "upper bounds; allowed coordinates may be zero",
            "distinctness": (
                "at most one empty negative singleton per generic Lorentz "
                "line and at most two empty K6 common neighbours globally"
            ),
            "nonbipartite_components": "ignored conservatively by this layer",
            "survivor": "filter non-rejection only",
            "arithmetic": "exact integer and bit-mask logic only",
        },
        "commands": production_commands(
            output, certificates, checkpoint, abort,
            workers, checkpoint_every,
        ),
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "workers": workers,
            "checkpoint_every": checkpoint_every,
        },
        "determinism": (
            "ordered ProcessPoolExecutor.map; decisions and artifacts contain "
            "no timestamps or measured runtimes"
        ),
    }
    atomic_json(output, report)
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument(
        "--pilot", action="store_true", help="run the bounded audit pilot"
    )
    mode.add_argument(
        "--full", action="store_true", help="run all 822 pinned parent survivors"
    )
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--certificates", type=Path,
        default=ROOT / "d6_k6_empty_support_budget_certificates.json",
    )
    parser.add_argument(
        "--checkpoint", type=Path,
        default=ROOT / "d6_k6_empty_support_budget_checkpoint.json",
    )
    parser.add_argument(
        "--abort", type=Path,
        default=ROOT / "d6_k6_empty_support_budget_abort.json",
    )
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    parser.add_argument("--checkpoint-every", type=int, default=32)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    if args.workers < 1 or args.workers > DEFAULT_WORKERS:
        parser.error(f"workers must lie in 1..{DEFAULT_WORKERS}")
    if args.checkpoint_every < 1:
        parser.error("checkpoint-every must be positive")
    if args.pilot:
        if args.resume:
            parser.error("--resume is only valid with --full")
        output = (
            args.output
            if args.output is not None
            else ROOT / "d6_k6_empty_support_budget_pilot_report.json"
        )
        report = pilot(output.resolve())
        print(json.dumps({
            "status": report["status"],
            "selection": len(report["selection"]),
            "rejected_indices": report["rejected_indices"],
            "runtime_seconds": report["runtime_seconds"],
        }, sort_keys=True))
        return 0

    output = (
        args.output
        if args.output is not None
        else ROOT / "d6_k6_empty_support_budget_report.json"
    )
    try:
        report = run_production(
            output.resolve(), args.certificates.resolve(),
            args.checkpoint.resolve(), args.abort.resolve(),
            args.workers, args.checkpoint_every, args.resume,
        )
    except InfrastructureAbort as error:
        print(str(error), file=sys.stderr)
        return 2
    print(json.dumps({
        "status": report["status"],
        "input_graphs": report["input_graphs"],
        "graphs_rejected": report["graphs_rejected"],
        "graphs_surviving": report["graphs_surviving"],
        "rejected_indices": report["rejected_indices"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
