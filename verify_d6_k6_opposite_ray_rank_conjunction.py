#!/usr/bin/env python3
"""Independent verifier for the K6 opposite-light-ray rank conjunction.

The checker imports neither the new producer nor either discovery probe.  It
uses the frozen independent tight-same-Z0 kernel, independently reconstructs
the projection support graph and its component-rank bound, enumerates ray
colorings directly, and adds a separately ordered actual-support DFS with the
tight same-bin support rule.
"""

from __future__ import annotations

import argparse
import ast
import gzip
import hashlib
import json
import os
import platform
import shlex
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from itertools import combinations
from pathlib import Path

import verify_d6_k6_tight_same_z0 as prior_verifier


base = prior_verifier.base
zf_reference = prior_verifier.zf_reference
ROOT = Path(__file__).resolve().parent
REPORT_SCHEMA = "d6-k6-opposite-ray-rank-conjunction-v1"
CERTIFICATE_SCHEMA = "d6-k6-opposite-ray-rank-conjunction-certificates-v1"
VERIFICATION_SCHEMA = "d6-k6-opposite-ray-rank-conjunction-verification-v1"
EXPECTED_INPUT = 623
EXPECTED_INPUT_SHA256 = (
    "b31aeac00b91d0d843ea51909c64f2c4ca3da49e5d2d33792312a45aa9f25e79"
)
EXPECTED_REJECTIONS = 372
EXPECTED_REJECTED_SHA256 = (
    "968195801d0311b0e43349e6deb8b08f61e5167b76f26edd8905369a6e574c63"
)
EXPECTED_RESIDUE_SHA256 = (
    "a23da7bd9ae21eb7c078a08a7bce92f22d397d5c71f1e4f51d9cf8e59bbfa907"
)
DEFAULT_WORKERS = max(1, (os.cpu_count() or 2) - 1)
PACKAGE_SOURCES = (
    "d6_k6_opposite_ray_rank_conjunction.py",
    "verify_d6_k6_opposite_ray_rank_conjunction.py",
    "test_d6_k6_opposite_ray_rank_conjunction.py",
    "d6_k6_opposite_ray_rank_conjunction.md",
)
EXPECTED_PRODUCTION_DEPENDENCIES = {
    "d6_k6_tight_same_z0.py": (
        "b761ded129805bd957c5ecd2d6cbcdbaadb6972d750d18260035ca175c5b949e"
    ),
    "d6_k6_tight_same_z0_report.json": (
        "2f2fb63e91518717f3a340c4e9401f720a3393ee6816535eadff1b89d681be04"
    ),
    "d6_k6_tight_same_z0_certificates.json": (
        "ad408071d49d95724c6b9a45afca947b937ea1ce8b79fb2b67cd810a91756a11"
    ),
    "d6_k6_tight_same_z0_verification.json": (
        "d87ecd688301e2e811e8eb3f3fe52111a050765643de6a5dbb1f1c06582b311e"
    ),
    "d6_k6_bipartite_rank_reference.py": (
        "e760fef74c423a89f1e1f9d08909fa1223e0389900cf3d45b0140ee296dea922"
    ),
    "d6_k6_support_reference.py": (
        "d6481137ef49d88154882285660dd755eb7cb652ab87280f05c3792de9742577"
    ),
}
EXPECTED_PARENT_ARTIFACTS = {
    "d6_k6_tight_same_z0.py": (
        "b761ded129805bd957c5ecd2d6cbcdbaadb6972d750d18260035ca175c5b949e"
    ),
    "d6_k6_tight_same_z0_report.json": (
        "2f2fb63e91518717f3a340c4e9401f720a3393ee6816535eadff1b89d681be04"
    ),
    "d6_k6_tight_same_z0_certificates.json": (
        "ad408071d49d95724c6b9a45afca947b937ea1ce8b79fb2b67cd810a91756a11"
    ),
    "d6_k6_tight_same_z0_verification.json": (
        "d87ecd688301e2e811e8eb3f3fe52111a050765643de6a5dbb1f1c06582b311e"
    ),
    "verify_d6_k6_tight_same_z0.py": (
        "b9a601d760aab00ec548354160e382e434a22b02d7dafd9e165be851e9686c51"
    ),
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


def load_gzip_json(path: Path) -> tuple[dict, str, int]:
    with gzip.open(path, "rb") as stream:
        raw = stream.read()
    return json.loads(raw), hashlib.sha256(raw).hexdigest(), len(raw)


def worker_initializer() -> None:
    for name in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[name] = "1"
    prior_verifier.prior_verifier.activate_kernel()


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
        "d6_k6_opposite_ray_rank_conjunction",
        "probe_d6_k6_opposite_ray_rank_conjunction",
        "probe_d6_k6_opposite_ray_rank",
        "d6_k6_support_reference",
        "d6_k6_bipartite_rank_reference",
    }
    if imports & forbidden:
        raise AssertionError(
            f"independent checker imports production code: {imports & forbidden}"
        )


def validate_porcelain_binding(git: dict) -> dict:
    """Require a correctly hashed launch status containing only untracked dirt."""

    lines = git.get("porcelain_lines")
    if not isinstance(lines, list) or not all(
        isinstance(line, str) for line in lines
    ):
        raise ValueError("launch porcelain lines are malformed")
    porcelain = "\n".join(lines)
    if (
        hashlib.sha256(porcelain.encode("utf-8")).hexdigest()
        != git.get("porcelain_sha256")
        or git.get("dirty") != bool(lines)
        or not all(line.startswith("?? ") for line in lines)
    ):
        raise ValueError("tracked source was dirty or porcelain binding differs")
    return {
        "dirty": bool(lines),
        "porcelain_sha256": git["porcelain_sha256"],
        "untracked_entries": len(lines),
    }


def validate_git_provenance(report: dict) -> dict:
    """Check that the launch commit contains every exact package source."""

    expected_sources = {name: sha256(ROOT / name) for name in PACKAGE_SOURCES}
    if report.get("source_sha256") != expected_sources:
        raise ValueError("working package source boundary differs")
    git = report.get("execution", {}).get("git", {})
    if (
        not git.get("available")
        or git.get("branch") != "codex/dimension6"
        or not git.get("tracked_clean")
        or git.get("committed_source_sha256") != expected_sources
    ):
        raise ValueError("launch git provenance unavailable or wrong branch")
    commit = git.get("commit")
    if not isinstance(commit, str) or len(commit) != 40:
        raise ValueError("recorded launch commit is malformed")
    for name, expected in expected_sources.items():
        try:
            blob = subprocess.run(
                ["git", "show", f"{commit}:{name}"],
                cwd=ROOT,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            ).stdout
        except (OSError, subprocess.SubprocessError) as error:
            raise ValueError(
                f"cannot read committed package source: {name}"
            ) from error
        if hashlib.sha256(blob).hexdigest() != expected:
            raise ValueError(f"committed package source differs: {name}")
    porcelain = validate_porcelain_binding(git)
    return {
        "commit": commit,
        "branch": git["branch"],
        "tracked_sources_clean_at_launch": True,
        **porcelain,
    }


def load_input() -> tuple[list[dict], list[int]]:
    observed = {
        name: sha256(ROOT / name) for name in EXPECTED_PARENT_ARTIFACTS
    }
    if observed != EXPECTED_PARENT_ARTIFACTS:
        raise ValueError(f"frozen parent/independent artifacts changed: {observed}")
    report = json.loads(
        (ROOT / "d6_k6_tight_same_z0_report.json").read_text(encoding="utf-8")
    )
    verification = json.loads(
        (ROOT / "d6_k6_tight_same_z0_verification.json").read_text(
            encoding="utf-8"
        )
    )
    if (
        report.get("schema") != "d6-k6-tight-same-z0-v1"
        or report.get("status") != "COMPLETE"
        or report.get("graphs_rejected") != 2
        or report.get("graphs_surviving") != EXPECTED_INPUT
        or report.get("ordered_residue_indices_sha256")
        != EXPECTED_INPUT_SHA256
        or verification.get("schema")
        != "d6-k6-tight-same-z0-verification-v1"
        or verification.get("status") != "PASS"
        or verification.get("graphs_rejected") != 2
        or verification.get("graphs_surviving") != EXPECTED_INPUT
        or verification.get("ordered_residue_indices_sha256")
        != EXPECTED_INPUT_SHA256
    ):
        raise ValueError("frozen tight same-Z0 boundary changed")
    parent_records, _ = prior_verifier.load_input()
    by_index = {int(record["index"]): record for record in parent_records}
    indices = [int(index) for index in report["ordered_residue_indices"]]
    if (
        len(indices) != EXPECTED_INPUT
        or stable_hash(indices) != EXPECTED_INPUT_SHA256
        or any(index not in by_index for index in indices)
    ):
        raise ValueError("independent 623 input reconstruction changed")
    return [by_index[index] for index in indices], indices


def independent_projection_check(adjacency, instance, z0, source, basis) -> dict:
    """Independent edge-list/union-find reconstruction of one direction."""

    source_light = source & ~z0
    basis_light = basis & ~z0
    source_locals = tuple(base.bits(source_light))
    basis_locals = tuple(base.bits(basis_light))
    basis_absolute = [instance.outside[local] for local in basis_locals]
    source_absolute = [instance.outside[local] for local in source_locals]
    cross_neighbors = []
    for absolute in source_absolute:
        cross_neighbors.append({
            target for target in basis_absolute
            if adjacency[absolute] & (1 << target)
        })
    parent = list(range(len(source_locals)))

    def root(item: int) -> int:
        while parent[item] != item:
            parent[item] = parent[parent[item]]
            item = parent[item]
        return item

    def union(left: int, right: int) -> None:
        left_root = root(left)
        right_root = root(right)
        if left_root == right_root:
            return
        if left_root < right_root:
            parent[right_root] = left_root
        else:
            parent[left_root] = right_root

    for right in range(len(source_locals)):
        for left in range(right):
            if cross_neighbors[left] & cross_neighbors[right]:
                union(left, right)
    groups: dict[int, list[int]] = {}
    for local in range(len(source_locals)):
        groups.setdefault(root(local), []).append(local)
    components = [
        group for _, group in sorted(groups.items()) if len(group) >= 2
    ]
    zero_cross = [
        position for position, neighbors in enumerate(cross_neighbors)
        if not neighbors
    ]
    component_lower = sum(len(component) - 1 for component in components)
    lower = len(zero_cross) + component_lower
    capacity = base.COORDINATES - basis.bit_count()
    return {
        "passed": lower <= capacity,
        "basis_size_including_Z0": basis.bit_count(),
        "rank_capacity": capacity,
        "rank_lower": lower,
        "zero_cross_degree_rank": len(zero_cross),
        "component_rank_lower": component_lower,
        "source": source_absolute,
        "basis": [instance.outside[local] for local in base.bits(basis)],
        "basis_light": basis_absolute,
        "zero_cross_degree_source": [
            source_absolute[position] for position in zero_cross
        ],
        "common_neighbour_components": [
            [source_absolute[position] for position in component]
            for component in components
        ],
    }


def assigned_hall_passes(selected, assignments, outside_size) -> bool:
    assigned = selected & sum(1 << local for local in assignments)
    masks = [0] * outside_size
    for local, support in assignments.items():
        masks[local] = support
    return base.hall_matchable(assigned, masks)


def tight_supports_pass(selected, assignments) -> bool:
    """Combination-ordered independent transcription of tight spans."""

    assigned = [local for local in base.bits(selected) if local in assignments]
    for size in range(1, len(assigned)):
        for chosen in combinations(assigned, size):
            chosen_set = set(chosen)
            union = 0
            outside_union = 0
            for local in assigned:
                if local in chosen_set:
                    union |= assignments[local]
                else:
                    outside_union |= assignments[local]
            if union.bit_count() == size and union & outside_union:
                return False
    return True


def independent_cross_edge_local_pass(
    adjacency, instance, z0, bins, assignments,
) -> bool:
    left_assigned = [
        local for local in base.bits(bins[0] & ~z0)
        if local in assignments
    ]
    right_assigned = [
        local for local in base.bits(bins[1] & ~z0)
        if local in assignments
    ]
    required_pairs = [
        (left, right)
        for left in left_assigned
        for right in right_assigned
        if adjacency[instance.outside[right]]
        & (1 << instance.outside[left])
    ]
    for left, right in required_pairs:
        left_coordinates = set(base.bits(assignments[left]))
        right_coordinates = set(base.bits(assignments[right]))
        if left_coordinates.isdisjoint(right_coordinates):
            return False
        if len(left_coordinates) == len(right_coordinates) == 1:
            return False
    return True


def mask_family_hall(masks) -> bool:
    for size in range(1, len(masks) + 1):
        for chosen in combinations(range(len(masks)), size):
            union = 0
            for number in chosen:
                union |= masks[number]
            if union.bit_count() < size:
                return False
    return True


def independent_cross_components(adjacency, instance, left, right):
    left_locals = base.bits(left)
    right_locals = base.bits(right)
    left_neighbors = {
        local: {
            other for other in right_locals
            if adjacency[instance.outside[local]]
            & (1 << instance.outside[other])
        }
        for local in left_locals
    }
    right_neighbors = {
        local: {
            other for other in left_locals
            if adjacency[instance.outside[local]]
            & (1 << instance.outside[other])
        }
        for local in right_locals
    }
    unseen = {
        (0, local) for local in left_locals if left_neighbors[local]
    } | {
        (1, local) for local in right_locals if right_neighbors[local]
    }
    answer = []
    while unseen:
        root = min(unseen)
        reached = {root}
        frontier = [root]
        while frontier:
            side, local = frontier.pop()
            neighbors = (
                {(1, item) for item in left_neighbors[local]}
                if side == 0
                else {(0, item) for item in right_neighbors[local]}
            )
            for item in sorted(neighbors & unseen):
                if item not in reached:
                    reached.add(item)
                    frontier.append(item)
        unseen -= reached
        answer.append((
            sorted(local for side, local in reached if side == 0),
            sorted(local for side, local in reached if side == 1),
        ))
    return answer


def independent_private_coordinates(locals_, assignments):
    answer = 0
    for coordinate in range(base.COORDINATES):
        occurrences = [
            local for local in locals_
            if assignments[local] & (1 << coordinate)
        ]
        if len(occurrences) == 1:
            answer |= 1 << coordinate
    return answer


def independent_isolated_edge_line_energy_pass(z0_masks, edge_masks):
    tagged = [(mask, 0) for mask in z0_masks]
    tagged.extend((mask, 1) for mask in edge_masks)
    for size in range(1, len(tagged) + 1):
        for chosen in combinations(range(len(tagged)), size):
            if not any(tagged[number][1] for number in chosen):
                continue
            coordinates = set()
            for number in chosen:
                coordinates.update(base.bits(tagged[number][0]))
            if len(coordinates) == size:
                return False
    return True


def independent_cross_intersection_pass(
    adjacency, instance, z0, bins, assignments,
):
    left = bins[0] & ~z0
    right = bins[1] & ~z0
    components = independent_cross_components(
        adjacency, instance, left, right
    )
    data = []
    used_left = set()
    used_right = set()
    for left_component, right_component in components:
        left_union = 0
        for local in left_component:
            left_union |= assignments[local]
        right_union = 0
        for local in right_component:
            right_union |= assignments[local]
        intersection = left_union & right_union
        eligible = bool(intersection)
        if eligible:
            eligible = not (
                independent_private_coordinates(
                    left_component, assignments
                ) & ~right_union
                or independent_private_coordinates(
                    right_component, assignments
                ) & ~left_union
            )
        data.append((
            len(left_component) + len(right_component),
            left_union | right_union,
            intersection,
            eligible,
            len(left_component) == len(right_component) == 1,
        ))
        used_left.update(left_component)
        used_right.update(right_component)
    singletons = [
        (1, assignments[local]) for local in base.bits(left)
        if local not in used_left
    ] + [
        (1, assignments[local]) for local in base.bits(right)
        if local not in used_right
    ] + [
        (1, assignments[local]) for local in base.bits(z0)
    ]
    z0_supports = [assignments[local] for local in base.bits(z0)]
    for selected in range(1 << len(data)):
        if any(
            selected & (1 << number) and not component[3]
            for number, component in enumerate(data)
        ):
            continue
        shared = [
            component[2]
            for number, component in enumerate(data)
            if selected & (1 << number)
        ]
        if not mask_family_hall([*z0_supports, *shared]):
            continue
        edge_lines = [
            component[2]
            for number, component in enumerate(data)
            if selected & (1 << number) and component[4]
        ]
        if not independent_isolated_edge_line_energy_pass(
            z0_supports, edge_lines
        ):
            continue
        blocks = [
            (
                component[0] - int(bool(selected & (1 << number))),
                component[1],
            )
            for number, component in enumerate(data)
        ] + singletons
        passed = True
        for size in range(1, len(blocks) + 1):
            for chosen in combinations(range(len(blocks)), size):
                demand = sum(blocks[number][0] for number in chosen)
                union = 0
                for number in chosen:
                    union |= blocks[number][1]
                if demand > union.bit_count():
                    passed = False
                    break
            if not passed:
                break
        if passed:
            return True
    return False


def independent_actual_support_exists(adjacency, instance, z0, bins) -> bool:
    if bins[0] & bins[1] != z0:
        raise ValueError("independent light-ray bins intersect outside Z0")
    involved = bins[0] | bins[1]
    domains = {
        local: base.actual_support_domains(
            instance.defects[local], bool(z0 & (1 << local))
        )
        for local in base.bits(involved)
    }
    peers = {}
    for local in base.bits(involved):
        peer_mask = 0
        for bin_mask in bins:
            if bin_mask & (1 << local):
                peer_mask |= bin_mask & ~(1 << local)
        peers[local] = peer_mask
    assignments = {}

    def choices(local):
        return [
            support for support in domains[local]
            if all(
                (support & other_support).bit_count() != 1
                for other, other_support in assignments.items()
                if peers[local] & (1 << other)
            )
        ]

    def visit(unassigned):
        if not unassigned:
            return (
                all(
                    tight_supports_pass(bin_mask, assignments)
                    for bin_mask in bins
                )
                and independent_cross_intersection_pass(
                    adjacency, instance, z0, bins, assignments
                )
            )
        ranked = []
        for local in base.bits(unassigned):
            local_choices = choices(local)
            ranked.append((len(local_choices), -local, local_choices, local))
        _, _, local_choices, local = min(ranked)
        if not local_choices:
            return False
        bit = 1 << local
        for support in local_choices:
            assignments[local] = support
            if (
                all(
                    assigned_hall_passes(
                        bin_mask, assignments, len(instance.outside)
                    )
                    and tight_supports_pass(bin_mask, assignments)
                    for bin_mask in bins
                )
                and independent_cross_edge_local_pass(
                    adjacency, instance, z0, bins, assignments
                )
                and visit(unassigned ^ bit)
            ):
                return True
            assignments.pop(local, None)
        return False

    return visit(involved)


def light_ray_decision(adjacency, instance, z0, local_components):
    odd = [
        component.vertices
        for component in local_components
        if not component.bipartite
    ]
    raw_colorings = 1 if not odd else 1 << (len(odd) - 1)
    failures = []
    for encoded in range(raw_colorings):
        if odd:
            bins = [z0 | odd[0], z0]
            colors = [0]
            for position, component in enumerate(odd[1:]):
                color = (encoded >> position) & 1
                colors.append(color)
                bins[color] |= component
        else:
            bins = [z0, z0]
            colors = []
        pair = (bins[0], bins[1])
        if not all(base.hall_matchable(bin_mask, instance.defects) for bin_mask in pair):
            continue
        checks = (
            independent_projection_check(
                adjacency, instance, z0, pair[0], pair[1]
            ),
            independent_projection_check(
                adjacency, instance, z0, pair[1], pair[0]
            ),
        )
        rendered_bins = [
            [instance.outside[local] for local in base.bits(bin_mask)]
            for bin_mask in pair
        ]
        if (
            pair[0].bit_count() == base.COORDINATES
            and pair[1].bit_count() == base.COORDINATES
            and bool(pair[0] & ~z0)
        ):
            failures.append({
                "component_colors": colors,
                "bins": rendered_bins,
                "reason": "both_saturated_light_ray_arithmetic",
                "Z0_size": z0.bit_count(),
                "light_vertices_per_bin": base.COORDINATES - z0.bit_count(),
                "projection_checks": list(checks),
            })
            continue
        if not all(check["passed"] for check in checks):
            failures.append({
                "component_colors": colors,
                "bins": rendered_bins,
                "reason": "opposite_ray_projection_rank",
                "projection_checks": list(checks),
            })
            continue
        if independent_actual_support_exists(
            adjacency, instance, z0, pair
        ):
            return True, None
        failures.append({
            "component_colors": colors,
            "bins": rendered_bins,
            "reason": "joint_actual_support",
            "projection_checks": list(checks),
        })
    return False, {
        "raw_colorings": raw_colorings,
        "matchable_colorings": len(failures),
        "failures": failures,
    }


def solve_seed(adjacency, instance):
    inertia = base.SympyInertiaCache()
    forcing = zf_reference.IndependentZeroForcing()
    cache = {}
    failures = []
    for z0 in base.z0_subsets(instance.eligible_z0):
        zvertices = [instance.outside[local] for local in base.bits(z0)]
        if not base.hall_matchable(z0, instance.defects):
            failures.append({"Z0": zvertices, "reason": "Z0_matching"})
            continue
        states, local_components = prior_verifier.bipartite_states(
            adjacency, instance, z0, inertia, forcing, cache
        )
        if not states:
            failures.append({
                "Z0": zvertices,
                "reason": "tight_bipartite_system",
            })
            continue
        passed, coloring_certificate = light_ray_decision(
            adjacency, instance, z0, local_components
        )
        if passed:
            return {"feasible": True, "failures": None}
        failures.append({
            "Z0": zvertices,
            "reason": "opposite_ray_rank_and_actual_support_same_coloring",
            "colorings": coloring_certificate,
        })
    return {"feasible": False, "failures": failures}


def evaluate_record(record):
    prior_verifier.prior_verifier.activate_kernel()
    adjacency = tuple(map(int, record["adjacency"]))
    base.validate_graph(adjacency)
    seeds_checked = 0
    for seed_mask in base.clique_masks(adjacency, base.COORDINATES):
        seeds_checked += 1
        instance = base.build_instance(adjacency, seed_mask)
        decision = solve_seed(adjacency, instance)
        if not decision["feasible"]:
            return {
                "index": int(record["index"]),
                "rejected": True,
                "seeds_checked": seeds_checked,
                "seed_mask": seed_mask,
                "seed": list(instance.seed),
                "failures": decision["failures"],
            }
    return {
        "index": int(record["index"]),
        "rejected": False,
        "seeds_checked": seeds_checked,
        "seed_mask": None,
        "seed": None,
        "failures": None,
    }


def positive_control() -> dict:
    adjacency = base.lower_bound_18_graph()
    passed = 0
    for seed_mask in base.clique_masks(adjacency, base.COORDINATES):
        instance = base.build_instance(adjacency, seed_mask)
        if not solve_seed(adjacency, instance)["feasible"]:
            raise AssertionError("independent opposite-ray positive control failed")
        passed += 1
    if passed != 32:
        raise AssertionError("independent positive-control seed count changed")
    return {"passed": True, "K6_seeds": passed}


def verify(report_path, certificates_path, output, workers):
    assert_import_independence()
    records, indices = load_input()
    report = json.loads(report_path.read_text(encoding="utf-8"))
    archive, raw_sha256, raw_bytes = load_gzip_json(certificates_path)
    source_path = ROOT / "d6_k6_opposite_ray_rank_conjunction.py"
    package_sources = {name: sha256(ROOT / name) for name in PACKAGE_SOURCES}
    if (
        report.get("schema") != REPORT_SCHEMA
        or report.get("status") != "COMPLETE"
        or report.get("production_source_sha256") != sha256(source_path)
        or report.get("source_sha256") != package_sources
        or report.get("dependencies") != EXPECTED_PRODUCTION_DEPENDENCIES
        or report.get("ordered_input_indices") != indices
        or report.get("ordered_input_indices_sha256") != EXPECTED_INPUT_SHA256
        or report.get("input_graphs") != EXPECTED_INPUT
        or report.get("graphs_rejected") != EXPECTED_REJECTIONS
        or report.get("graphs_surviving") != EXPECTED_INPUT - EXPECTED_REJECTIONS
        or report.get("rejected_indices_sha256") != EXPECTED_REJECTED_SHA256
        or report.get("ordered_residue_indices_sha256")
        != EXPECTED_RESIDUE_SHA256
        or report.get("certificate_archive", {}).get("sha256")
        != sha256(certificates_path)
        or report.get("certificate_archive", {}).get("uncompressed_sha256")
        != raw_sha256
        or report.get("certificate_archive", {}).get("uncompressed_bytes")
        != raw_bytes
        or report.get("positive_18_control")
        != {"passed": True, "K6_seeds": 32}
        or archive.get("schema") != CERTIFICATE_SCHEMA
        or archive.get("production_source_sha256") != sha256(source_path)
        or archive.get("source_sha256") != package_sources
        or archive.get("ordered_input_indices_sha256") != EXPECTED_INPUT_SHA256
    ):
        raise ValueError("opposite-ray production boundary mismatch")
    provenance = validate_git_provenance(report)
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
    if (
        len(rejected) != EXPECTED_REJECTIONS
        or stable_hash(rejected) != EXPECTED_REJECTED_SHA256
        or stable_hash(residue) != EXPECTED_RESIDUE_SHA256
        or rejected != report.get("rejected_indices")
        or stable_hash(rejected) != report.get("rejected_indices_sha256")
        or residue != report.get("ordered_residue_indices")
        or stable_hash(residue) != report.get("ordered_residue_indices_sha256")
        or len(rejected) != report.get("graphs_rejected")
        or len(residue) != report.get("graphs_surviving")
    ):
        raise AssertionError("independent opposite-ray decisions differ")
    saved_rows = {int(row["index"]): row for row in report["graph_results"]}
    certificates = {int(row["index"]): row for row in archive["rejected_graphs"]}
    if set(saved_rows) != set(indices) or len(saved_rows) != len(indices):
        raise AssertionError("production graph-result coverage differs")
    if (
        set(certificates) != set(rejected)
        or len(certificates) != len(archive["rejected_graphs"])
    ):
        raise AssertionError("certificate archive index coverage differs")
    for row in results:
        saved = saved_rows[row["index"]]
        if (
            row["rejected"] != saved["rejected"]
            or row["seeds_checked"] != saved["seeds_checked"]
            or row["seed_mask"] != saved["first_impossible_seed_mask"]
            or row["seed"] != saved["first_impossible_seed"]
        ):
            raise AssertionError(f"graph detail differs at {row['index']}")
        if row["rejected"]:
            certificate = certificates[row["index"]]
            if (
                certificate["seed_mask"] != row["seed_mask"]
                or certificate["seed"] != row["seed"]
                or certificate["Z0_failures"] != row["failures"]
            ):
                raise AssertionError(f"certificate differs at {row['index']}")
    control = positive_control()
    result = {
        "schema": VERIFICATION_SCHEMA,
        "status": "PASS",
        "report": {"path": report_path.name, "sha256": sha256(report_path)},
        "certificates": {
            "path": certificates_path.name,
            "sha256": sha256(certificates_path),
            "uncompressed_sha256": raw_sha256,
        },
        "sources": {
            Path(__file__).name: sha256(Path(__file__)),
            source_path.name: sha256(source_path),
            **{
                name: sha256(ROOT / name)
                for name in EXPECTED_PARENT_ARTIFACTS
            },
        },
        "checks": {
            "import_independence": True,
            "independent_projection_support_graph_and_union_find": True,
            "independent_PSD_Z_component_rank_bound": True,
            "independent_both_saturated_integer_type_equation": True,
            "independent_tight_actual_support_DFS": True,
            "independent_opposite_edge_support_overlap_and_values": True,
            "independent_cross_component_intersection_Hall": True,
            "independent_isolated_edge_all_ones_energy": True,
            "same_Z0_and_same_ray_coloring": True,
            "all_graph_decisions": True,
            "all_rejected_Z0_and_coloring_partitions": True,
            "positive_18_control": True,
        },
        "provenance": provenance,
        "graphs_recomputed": len(results),
        "graphs_rejected": len(rejected),
        "graphs_surviving": len(residue),
        "rejected_indices": rejected,
        "rejected_indices_sha256": stable_hash(rejected),
        "ordered_residue_indices_sha256": stable_hash(residue),
        "positive_18_control": control,
        "runtime": {
            "command": shlex.join([sys.executable, *sys.argv]),
            "workers": workers,
            "wall_seconds": time.monotonic() - started,
            "python": sys.version,
            "platform": platform.platform(),
        },
    }
    atomic_json(output, result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--report", type=Path,
        default=ROOT / "d6_k6_opposite_ray_rank_conjunction_report.json",
    )
    parser.add_argument(
        "--certificates", type=Path,
        default=(
            ROOT / "d6_k6_opposite_ray_rank_conjunction_certificates.json.gz"
        ),
    )
    parser.add_argument(
        "--output", type=Path,
        default=ROOT / "d6_k6_opposite_ray_rank_conjunction_verification.json",
    )
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("workers must be positive")
    result = verify(args.report, args.certificates, args.output, args.workers)
    print(json.dumps({
        "status": result["status"],
        "rejected": result["graphs_rejected"],
        "surviving": result["graphs_surviving"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
