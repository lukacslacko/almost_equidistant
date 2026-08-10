#!/usr/bin/env python3
"""Exact opposite-light-ray projection rank on the certified K6 residue.

For one fixed K6 seed, Z0, and coloring of forced nonbipartite Lorentz
components onto the two light rays, this layer conjoins three necessary
systems on their correct shared quantifiers:

* the frozen strongest bipartite-component system uses the same Z0;
* opposite-ray projection residuals obey an exact PSD-Z rank bound;
* the frozen actual-support CSP uses the same light-ray coloring.

Candidate nonedges remain unconstrained and every allowed defect coordinate
may vanish.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import platform
import shlex
import subprocess
import sys
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from datetime import UTC, datetime
from itertools import combinations
from pathlib import Path
from typing import Sequence

import d6_k6_tight_same_z0 as parent


ROOT = Path(__file__).resolve().parent
PARENT_REPORT = ROOT / "d6_k6_tight_same_z0_report.json"
PARENT_CERTIFICATES = ROOT / "d6_k6_tight_same_z0_certificates.json"
PARENT_VERIFICATION = ROOT / "d6_k6_tight_same_z0_verification.json"
REPORT_SCHEMA = "d6-k6-opposite-ray-rank-conjunction-v1"
CERTIFICATE_SCHEMA = "d6-k6-opposite-ray-rank-conjunction-certificates-v1"
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

EXPECTED_DEPENDENCIES = {
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


@dataclass
class SupportResult:
    feasible: bool
    dfs_nodes: int
    supports: dict[int, int] | None


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


def atomic_gzip_json(path: Path, value: object) -> tuple[str, int]:
    """Write deterministic gzip JSON and return raw SHA-256 and byte count."""

    raw = (
        json.dumps(value, indent=2, sort_keys=True).encode("utf-8") + b"\n"
    )
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    with temporary.open("wb") as raw_stream:
        with gzip.GzipFile(
            filename="", mode="wb", compresslevel=9, mtime=0,
            fileobj=raw_stream,
        ) as stream:
            stream.write(raw)
        raw_stream.flush()
        os.fsync(raw_stream.fileno())
    os.replace(temporary, path)
    return hashlib.sha256(raw).hexdigest(), len(raw)


def worker_initializer() -> None:
    for name in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[name] = "1"
    parent.prior.activate_kernel()


def git_provenance() -> dict:
    """Require every package source to match the recorded commit blob."""

    def git(*arguments: str, binary: bool = False):
        return subprocess.run(
            ["git", *arguments],
            cwd=ROOT,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=not binary,
        ).stdout

    try:
        commit = git("rev-parse", "HEAD").strip()
        branch = git("branch", "--show-current").strip()
        status = git(
            "status", "--porcelain=v1", "--untracked-files=all"
        ).strip()
    except (OSError, subprocess.SubprocessError) as error:
        raise RuntimeError(f"git provenance unavailable: {error}") from error
    lines = status.splitlines() if status else []
    if branch != "codex/dimension6":
        raise ValueError(f"wrong theorem branch at launch: {branch!r}")
    if not all(line.startswith("?? ") for line in lines):
        raise ValueError("tracked or staged worktree dirt blocks theorem launch")
    source_hashes = {name: sha256(ROOT / name) for name in PACKAGE_SOURCES}
    for name, expected in source_hashes.items():
        try:
            blob = git("show", f"{commit}:{name}", binary=True)
        except (OSError, subprocess.SubprocessError) as error:
            raise RuntimeError(f"source is not committed at launch: {name}") from error
        if hashlib.sha256(blob).hexdigest() != expected:
            raise ValueError(f"committed source differs at launch: {name}")
    return {
        "available": True,
        "commit": commit,
        "branch": branch,
        "tracked_clean": True,
        "dirty": bool(lines),
        "porcelain_sha256": hashlib.sha256(status.encode("utf-8")).hexdigest(),
        "porcelain_lines": lines,
        "committed_source_sha256": source_hashes,
    }


def load_input() -> tuple[list[dict], list[int], dict[str, str]]:
    observed = {
        name: sha256(ROOT / name) for name in EXPECTED_DEPENDENCIES
    }
    if observed != EXPECTED_DEPENDENCIES:
        raise ValueError(f"opposite-ray dependency boundary changed: {observed}")
    report = json.loads(PARENT_REPORT.read_text(encoding="utf-8"))
    verification = json.loads(PARENT_VERIFICATION.read_text(encoding="utf-8"))
    if (
        report.get("schema") != parent.REPORT_SCHEMA
        or report.get("status") != "COMPLETE"
        or report.get("graphs_rejected") != 2
        or report.get("graphs_surviving") != EXPECTED_INPUT
        or report.get("ordered_residue_indices_sha256")
        != EXPECTED_INPUT_SHA256
        or report.get("certificate_archive", {}).get("sha256")
        != sha256(PARENT_CERTIFICATES)
        or verification.get("schema")
        != "d6-k6-tight-same-z0-verification-v1"
        or verification.get("status") != "PASS"
        or verification.get("graphs_rejected") != 2
        or verification.get("graphs_surviving") != EXPECTED_INPUT
        or verification.get("ordered_residue_indices_sha256")
        != EXPECTED_INPUT_SHA256
    ):
        raise ValueError("frozen tight same-Z0 result boundary changed")
    parent_records, _, _ = parent.load_input()
    by_index = {int(record["index"]): record for record in parent_records}
    indices = [int(index) for index in report["ordered_residue_indices"]]
    if (
        len(indices) != EXPECTED_INPUT
        or stable_hash(indices) != EXPECTED_INPUT_SHA256
        or any(index not in by_index for index in indices)
    ):
        raise ValueError("certified 623-graph residue changed")
    return [by_index[index] for index in indices], indices, observed


def common_neighbour_rows(
    adjacency: Sequence[int], instance, source: int, target: int,
) -> tuple[tuple[int, ...], tuple[int, ...], tuple[int, ...]]:
    """Common-required-neighbour graph on local ``source`` vertices."""

    source_locals = tuple(parent.kernel.vertices(source))
    target_absolute = sum(
        1 << instance.outside[local]
        for local in parent.kernel.vertices(target)
    )
    rows = [0] * len(source_locals)
    zero_cross_degree = []
    for position, source_local in enumerate(source_locals):
        source_absolute = instance.outside[source_local]
        if not adjacency[source_absolute] & target_absolute:
            zero_cross_degree.append(position)
    for right, right_local in enumerate(source_locals):
        right_absolute = instance.outside[right_local]
        for left in range(right):
            left_absolute = instance.outside[source_locals[left]]
            if (
                adjacency[right_absolute]
                & adjacency[left_absolute]
                & target_absolute
            ):
                rows[right] |= 1 << left
                rows[left] |= 1 << right
    return source_locals, tuple(rows), tuple(zero_cross_degree)


def psd_z_component_rank_lower(rows: Sequence[int]) -> tuple[int, list[list[int]]]:
    """Return ``sum(|H|-1)`` over nontrivial support components."""

    unseen = (1 << len(rows)) - 1
    lower = 0
    components = []
    while unseen:
        root = unseen & -unseen
        reached = root
        frontier = root
        while frontier:
            bit = frontier & -frontier
            frontier ^= bit
            local = bit.bit_length() - 1
            new = rows[local] & unseen & ~reached
            reached |= new
            frontier |= new
        unseen &= ~reached
        if reached.bit_count() >= 2:
            lower += reached.bit_count() - 1
            components.append(parent.kernel.vertices(reached))
    return lower, components


def directed_projection_check(
    adjacency: Sequence[int], instance, z0: int, source: int, basis: int,
) -> dict:
    """Project one nonzero light-ray bin off the opposite orthonormal bin."""

    source_light = source & ~z0
    basis_light = basis & ~z0
    source_locals, rows, zero_cross_degree = common_neighbour_rows(
        adjacency, instance, source_light, basis_light
    )
    component_lower, components = psd_z_component_rank_lower(rows)
    lower = len(zero_cross_degree) + component_lower
    capacity = parent.kernel.COORDINATES - basis.bit_count()
    return {
        "passed": lower <= capacity,
        "basis_size_including_Z0": basis.bit_count(),
        "rank_capacity": capacity,
        "rank_lower": lower,
        "zero_cross_degree_rank": len(zero_cross_degree),
        "component_rank_lower": component_lower,
        "source": [instance.outside[local] for local in source_locals],
        "basis": [
            instance.outside[local] for local in parent.kernel.vertices(basis)
        ],
        "basis_light": [
            instance.outside[local]
            for local in parent.kernel.vertices(basis_light)
        ],
        "zero_cross_degree_source": [
            instance.outside[source_locals[position]]
            for position in zero_cross_degree
        ],
        "common_neighbour_components": [
            [instance.outside[source_locals[local]] for local in component]
            for component in components
        ],
    }


def coloring_projection_checks(
    adjacency: Sequence[int], instance, z0: int, bins: tuple[int, int],
) -> tuple[dict, dict]:
    return (
        directed_projection_check(adjacency, instance, z0, bins[0], bins[1]),
        directed_projection_check(adjacency, instance, z0, bins[1], bins[0]),
    )


def both_saturated_with_light(z0: int, bins: tuple[int, int]) -> bool:
    """Both ON bins are full bases and have a nonzero light-ray part."""

    return (
        bins[0].bit_count() == parent.kernel.COORDINATES
        and bins[1].bit_count() == parent.kernel.COORDINATES
        and bool(bins[0] & ~z0)
    )


def support_domains(allowed: int, zero_factor: bool) -> tuple[int, ...]:
    """Exhaust actual nonzero subsets of an allowed defect mask."""

    minimum = 3 if zero_factor else 1
    choices = []
    subset = allowed
    while subset:
        if subset.bit_count() >= minimum:
            choices.append(subset)
        subset = (subset - 1) & allowed
    choices.sort(key=lambda mask: (mask.bit_count(), mask))
    return tuple(choices)


def assigned_hall_passes(
    selected: int, assignments: dict[int, int], outside_size: int,
) -> bool:
    assigned = selected & sum(1 << local for local in assignments)
    masks = [0] * outside_size
    for local, support in assignments.items():
        masks[local] = support
    return parent.kernel.support_matching(assigned, masks) is not None


def tight_orthogonal_supports_pass(
    selected: int, assignments: dict[int, int],
) -> bool:
    """A tight actual-support subfamily fills its coordinate union."""

    assigned = [
        local for local in parent.kernel.vertices(selected)
        if local in assignments
    ]
    for encoded in range(1, 1 << len(assigned)):
        if encoded == (1 << len(assigned)) - 1:
            continue
        union = 0
        count = 0
        outside_union = 0
        for position, local in enumerate(assigned):
            if encoded & (1 << position):
                union |= assignments[local]
                count += 1
            else:
                outside_union |= assignments[local]
        if union.bit_count() == count and union & outside_union:
            return False
    return True


def opposite_required_edge_supports_pass(
    adjacency: Sequence[int], instance, z0: int,
    bins: tuple[int, int], assignments: dict[int, int],
) -> bool:
    """Check local nonzero-dot support rules across opposite light rays."""

    left = bins[0] & ~z0
    right = bins[1] & ~z0
    for left_local in parent.kernel.vertices(left):
        if left_local not in assignments:
            continue
        left_absolute = instance.outside[left_local]
        left_support = assignments[left_local]
        for right_local in parent.kernel.vertices(right):
            if (
                right_local not in assignments
                or not (
                    adjacency[left_absolute]
                    & (1 << instance.outside[right_local])
                )
            ):
                continue
            right_support = assignments[right_local]
            if not left_support & right_support:
                return False
            if left_support.bit_count() == right_support.bit_count() == 1:
                return False
    return True


def masks_matchable(masks: Sequence[int]) -> bool:
    """Exact augmenting-path matching for a short family of masks."""

    matched = [-1] * parent.kernel.COORDINATES

    def augment(row: int, seen: int) -> bool:
        for coordinate in parent.kernel.vertices(masks[row]):
            bit = 1 << coordinate
            if seen & bit:
                continue
            if matched[coordinate] < 0 or augment(matched[coordinate], seen | bit):
                matched[coordinate] = row
                return True
        return False

    return all(augment(row, 0) for row in range(len(masks)))


def cross_edge_components(
    adjacency: Sequence[int], instance, left: int, right: int,
) -> list[tuple[list[int], list[int]]]:
    """Required-edge components meeting both opposite light-ray bins."""

    left_vertices = parent.kernel.vertices(left)
    right_vertices = parent.kernel.vertices(right)
    combined = [(0, local) for local in left_vertices] + [
        (1, local) for local in right_vertices
    ]
    rows = [0] * len(combined)
    for left_position, left_local in enumerate(left_vertices):
        left_absolute = instance.outside[left_local]
        for right_offset, right_local in enumerate(right_vertices):
            right_absolute = instance.outside[right_local]
            if adjacency[left_absolute] & (1 << right_absolute):
                right_position = len(left_vertices) + right_offset
                rows[left_position] |= 1 << right_position
                rows[right_position] |= 1 << left_position
    unseen = sum(1 << position for position, row in enumerate(rows) if row)
    components = []
    while unseen:
        root = unseen & -unseen
        reached = root
        frontier = root
        while frontier:
            bit = frontier & -frontier
            frontier ^= bit
            position = bit.bit_length() - 1
            new = rows[position] & unseen & ~reached
            reached |= new
            frontier |= new
        unseen &= ~reached
        component_left = []
        component_right = []
        for position in parent.kernel.vertices(reached):
            side, local = combined[position]
            (component_left if side == 0 else component_right).append(local)
        components.append((component_left, component_right))
    return components


def private_support_coordinates(
    locals_: Sequence[int], assignments: dict[int, int],
) -> int:
    private = 0
    for coordinate in range(parent.kernel.COORDINATES):
        if sum(
            bool(assignments[local] & (1 << coordinate))
            for local in locals_
        ) == 1:
            private |= 1 << coordinate
    return private


def isolated_edge_line_tight_supports_pass(
    z0_supports: Sequence[int], isolated_edge_supports: Sequence[int],
) -> bool:
    """Exclude a coordinate-tight family containing a saturated edge line."""

    tagged = [
        (support, False) for support in z0_supports
    ] + [
        (support, True) for support in isolated_edge_supports
    ]
    for selected in range(1, 1 << len(tagged)):
        if not any(
            selected & (1 << number) and is_edge
            for number, (_, is_edge) in enumerate(tagged)
        ):
            continue
        support_union = 0
        line_count = 0
        for number, (support, _) in enumerate(tagged):
            if selected & (1 << number):
                support_union |= support
                line_count += 1
        if support_union.bit_count() == line_count:
            return False
    return True


def cross_component_intersection_supports_pass(
    adjacency: Sequence[int], instance, z0: int,
    bins: tuple[int, int], assignments: dict[int, int],
) -> bool:
    """Exhaust the possible 0/1 cross-span intersections componentwise."""

    left = bins[0] & ~z0
    right = bins[1] & ~z0
    component_data = []
    used_left = set()
    used_right = set()
    for left_component, right_component in cross_edge_components(
        adjacency, instance, left, right
    ):
        left_union = 0
        for local in left_component:
            left_union |= assignments[local]
        right_union = 0
        for local in right_component:
            right_union |= assignments[local]
        intersection = left_union & right_union
        eligible = bool(intersection)
        if eligible:
            private_left = private_support_coordinates(
                left_component, assignments
            )
            private_right = private_support_coordinates(
                right_component, assignments
            )
            eligible = not (
                private_left & ~right_union or private_right & ~left_union
            )
        component_data.append({
            "count": len(left_component) + len(right_component),
            "support": left_union | right_union,
            "intersection": intersection,
            "eligible": eligible,
            "isolated_edge": (
                len(left_component) == len(right_component) == 1
            ),
        })
        used_left.update(left_component)
        used_right.update(right_component)
    singleton_blocks = []
    for local in parent.kernel.vertices(left):
        if local not in used_left:
            singleton_blocks.append((1, assignments[local]))
    for local in parent.kernel.vertices(right):
        if local not in used_right:
            singleton_blocks.append((1, assignments[local]))
    z0_locals = parent.kernel.vertices(z0)
    for local in z0_locals:
        singleton_blocks.append((1, assignments[local]))
    z0_supports = [assignments[local] for local in z0_locals]

    for selected in range(1 << len(component_data)):
        if any(
            selected & (1 << number) and not data["eligible"]
            for number, data in enumerate(component_data)
        ):
            continue
        intersections = [
            data["intersection"]
            for number, data in enumerate(component_data)
            if selected & (1 << number)
        ]
        if not masks_matchable([*z0_supports, *intersections]):
            continue
        isolated_edge_supports = [
            data["intersection"]
            for number, data in enumerate(component_data)
            if (
                selected & (1 << number)
                and data["isolated_edge"]
            )
        ]
        if not isolated_edge_line_tight_supports_pass(
            z0_supports, isolated_edge_supports
        ):
            continue
        blocks = [
            (
                data["count"] - bool(selected & (1 << number)),
                data["support"],
            )
            for number, data in enumerate(component_data)
        ]
        blocks.extend(singleton_blocks)
        if all(
            sum(
                rank for rank, support in blocks
                if not support & ~container
            ) <= container.bit_count()
            for container in range(1 << parent.kernel.COORDINATES)
        ):
            return True
    return False


def tight_actual_supports_for_bins(
    adjacency: Sequence[int], instance, z0: int, bins: tuple[int, int],
) -> SupportResult:
    """Actual-support CSP with same-bin tight-span propagation."""

    if bins[0] & bins[1] != z0:
        raise ValueError("light-ray bins intersect outside Z0")
    involved = bins[0] | bins[1]
    domains = {
        local: support_domains(
            instance.defects[local], bool(z0 & (1 << local))
        )
        for local in parent.kernel.vertices(involved)
    }
    peers = {}
    for local in parent.kernel.vertices(involved):
        peer_mask = 0
        for bin_mask in bins:
            if bin_mask & (1 << local):
                peer_mask |= bin_mask & ~(1 << local)
        peers[local] = peer_mask
    assignments: dict[int, int] = {}
    nodes = 0

    def compatible_choices(local: int) -> list[int]:
        return [
            support for support in domains[local]
            if all(
                (support & other_support).bit_count() != 1
                for other, other_support in assignments.items()
                if peers[local] & (1 << other)
            )
        ]

    def visit(unassigned: int) -> bool:
        nonlocal nodes
        nodes += 1
        if not unassigned:
            return (
                all(
                    tight_orthogonal_supports_pass(bin_mask, assignments)
                    for bin_mask in bins
                )
                and cross_component_intersection_supports_pass(
                    adjacency, instance, z0, bins, assignments
                )
            )
        best = None
        best_choices = None
        best_key = None
        for local in parent.kernel.vertices(unassigned):
            choices = compatible_choices(local)
            key = (
                len(choices),
                -(peers[local] & unassigned).bit_count(),
                local,
            )
            if best_key is None or key < best_key:
                best = local
                best_choices = choices
                best_key = key
            if not choices:
                return False
        assert best is not None and best_choices is not None
        bit = 1 << best
        for support in best_choices:
            assignments[best] = support
            if (
                all(
                    assigned_hall_passes(
                        bin_mask, assignments, len(instance.outside)
                    )
                    and tight_orthogonal_supports_pass(bin_mask, assignments)
                    for bin_mask in bins
                )
                and opposite_required_edge_supports_pass(
                    adjacency, instance, z0, bins, assignments
                )
                and visit(unassigned ^ bit)
            ):
                return True
            assignments.pop(best, None)
        return False

    feasible = visit(involved)
    return SupportResult(
        feasible, nodes, assignments.copy() if feasible else None
    )


def light_ray_decision(
    adjacency: Sequence[int], instance, z0: int, components,
) -> tuple[dict | None, Counter[str], dict]:
    """Conjoin projection rank and actual supports on one ray coloring."""

    odd = [component.component for component in components if not component.bipartite]
    raw_colorings = 1 if not odd else 1 << (len(odd) - 1)
    counts: Counter[str] = Counter()
    failures = []
    for colors, bins in parent._pure_colorings(instance, z0, odd):
        counts["matchable_colorings"] += 1
        checks = coloring_projection_checks(adjacency, instance, z0, bins)
        for check in checks:
            counts[
                f"directed_basis_size_{check['basis_size_including_Z0']}"
            ] += 1
        if both_saturated_with_light(z0, bins):
            counts["both_saturated_arithmetic_failed_colorings"] += 1
            failures.append({
                "component_colors": list(colors),
                "bins": [
                    [
                        instance.outside[local]
                        for local in parent.kernel.vertices(bin_mask)
                    ]
                    for bin_mask in bins
                ],
                "reason": "both_saturated_light_ray_arithmetic",
                "Z0_size": z0.bit_count(),
                "light_vertices_per_bin": (
                    parent.kernel.COORDINATES - z0.bit_count()
                ),
                "projection_checks": list(checks),
            })
            continue
        failed_checks = [check for check in checks if not check["passed"]]
        if failed_checks:
            counts["projection_rank_failed_colorings"] += 1
            for check in failed_checks:
                counts[
                    f"projection_failure_basis_size_{check['basis_size_including_Z0']}"
                ] += 1
            failures.append({
                "component_colors": list(colors),
                "bins": [
                    [
                        instance.outside[local]
                        for local in parent.kernel.vertices(bin_mask)
                    ]
                    for bin_mask in bins
                ],
                "reason": "opposite_ray_projection_rank",
                "projection_checks": list(checks),
            })
            continue
        counts["projection_rank_passed_colorings"] += 1
        support = tight_actual_supports_for_bins(
            adjacency, instance, z0, bins
        )
        counts["actual_support_searches"] += 1
        counts["actual_support_dfs_nodes"] += support.dfs_nodes
        if support.feasible:
            counts["joint_coloring_passed"] += 1
            return ({
                "component_colors": list(colors),
                "bins": [
                    [
                        instance.outside[local]
                        for local in parent.kernel.vertices(bin_mask)
                    ]
                    for bin_mask in bins
                ],
                "projection_checks": list(checks),
            }, counts, {
                "raw_colorings": raw_colorings,
                "matchable_colorings": counts["matchable_colorings"],
                "failures": failures,
            })
        counts["actual_support_failed_colorings"] += 1
        failures.append({
            "component_colors": list(colors),
            "bins": [
                [
                    instance.outside[local]
                    for local in parent.kernel.vertices(bin_mask)
                ]
                for bin_mask in bins
            ],
            "reason": "joint_actual_support",
            "projection_checks": list(checks),
        })
    return None, counts, {
        "raw_colorings": raw_colorings,
        "matchable_colorings": counts["matchable_colorings"],
        "failures": failures,
    }


def solve_seed(adjacency: Sequence[int], instance) -> dict:
    inertia = parent.kernel.InertiaCache()
    forcing = parent.kernel.ZeroForcingSolver()
    subset_rank = parent.kernel.arbitrary.PositiveSubsetRank(inertia, forcing)
    counts: Counter[str] = Counter()
    failures = []
    for z0 in parent.kernel.ranks.z0_subsets(instance.eligible_z0_mask):
        counts["z0_considered"] += 1
        zvertices = [
            instance.outside[local] for local in parent.kernel.vertices(z0)
        ]
        if parent.kernel.support_matching(z0, instance.defects) is None:
            counts["z0_unmatchable"] += 1
            failures.append({"Z0": zvertices, "reason": "Z0_matching"})
            continue
        states, components, local = parent.bipartite_states(
            adjacency, instance, z0, inertia, forcing, subset_rank
        )
        counts.update(local)
        if not states:
            counts["z0_bipartite_failed"] += 1
            failures.append({
                "Z0": zvertices,
                "reason": "tight_bipartite_system",
            })
            continue
        counts["z0_bipartite_passed"] += 1
        witness, local, coloring_certificate = light_ray_decision(
            adjacency, instance, z0, components
        )
        counts.update(local)
        if witness is None:
            counts["z0_light_ray_conjunction_failed"] += 1
            failures.append({
                "Z0": zvertices,
                "reason": "opposite_ray_rank_and_actual_support_same_coloring",
                "colorings": coloring_certificate,
            })
            continue
        counts["z0_full_conjunction_passed"] += 1
        return {
            "feasible": True,
            "counts": dict(counts),
            "failures": None,
        }
    return {
        "feasible": False,
        "counts": dict(counts),
        "failures": failures,
    }


def evaluate_record(record: dict) -> dict:
    parent.prior.activate_kernel()
    adjacency = tuple(map(int, record["adjacency"]))
    total: Counter[str] = Counter()
    seeds_checked = 0
    for seed_mask in parent.kernel.clique_masks(
        adjacency, parent.kernel.COORDINATES
    ):
        seeds_checked += 1
        instance = parent.kernel.build_instance(
            adjacency, parent.kernel.vertices(seed_mask)
        )
        decision = solve_seed(adjacency, instance)
        total.update(decision["counts"])
        if decision["feasible"]:
            continue
        return {
            "index": int(record["index"]),
            "rejected": True,
            "seeds_checked": seeds_checked,
            "first_impossible_seed_mask": seed_mask,
            "first_impossible_seed": list(instance.seed),
            "counts": dict(total),
            "certificate": {
                "seed_mask": seed_mask,
                "seed": list(instance.seed),
                "Z0_failures": decision["failures"],
            },
        }
    return {
        "index": int(record["index"]),
        "rejected": False,
        "seeds_checked": seeds_checked,
        "first_impossible_seed_mask": None,
        "first_impossible_seed": None,
        "counts": dict(total),
        "certificate": None,
    }


def positive_control() -> dict:
    parent.prior.activate_kernel()
    adjacency = parent.kernel.lower_bound_18_graph()
    passed = 0
    for seed_mask in parent.kernel.clique_masks(
        adjacency, parent.kernel.COORDINATES
    ):
        instance = parent.kernel.build_instance(
            adjacency, parent.kernel.vertices(seed_mask)
        )
        if not solve_seed(adjacency, instance)["feasible"]:
            raise AssertionError("opposite-ray rank rejects positive 18 control")
        passed += 1
    if passed != 32:
        raise AssertionError("positive-control K6 seed count changed")
    return {"passed": True, "K6_seeds": passed}


def run(workers: int, output: Path, certificates: Path) -> dict:
    provenance = git_provenance()
    records, indices, dependencies = load_input()
    started_utc = datetime.now(UTC).isoformat()
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
    ):
        raise AssertionError("production result differs from source-bound census")
    totals: Counter[str] = Counter()
    for row in results:
        totals.update(row["counts"])
    archive = {
        "schema": CERTIFICATE_SCHEMA,
        "production_source_sha256": sha256(Path(__file__)),
        "source_sha256": {
            name: sha256(ROOT / name) for name in PACKAGE_SOURCES
        },
        "ordered_input_indices_sha256": EXPECTED_INPUT_SHA256,
        "rejected_graphs": [
            {"index": row["index"], **row["certificate"]}
            for row in results if row["rejected"]
        ],
    }
    raw_certificate_sha256, raw_certificate_bytes = atomic_gzip_json(
        certificates, archive
    )
    report = {
        "schema": REPORT_SCHEMA,
        "status": "COMPLETE",
        "description": (
            "same-Z0, same-coloring opposite-light-ray projection rank, "
            "actual-support, cross-component, and all-ones conjunction"
        ),
        "semantics": {
            "candidate_nonedges": "unconstrained_and_may_be_unit",
            "allowed_supports": "upper_bounds_actual_subsets_exhausted",
            "shared_quantifiers": "same_Z0_and_same_two_ray_coloring",
            "arithmetic": "exact_integer_graph_bitmask_Hall_inertia_and_DFS",
            "cross_component_intersections": (
                "all_componentwise_zero_or_one_span_intersection_choices"
            ),
            "opposite_required_edges": (
                "actual_support_overlap_and_no_double_singleton"
            ),
            "isolated_edge_line_energy": (
                "all_tagged_coordinate_tight_subfamilies_exhausted"
            ),
        },
        "production_source_sha256": sha256(Path(__file__)),
        "source_sha256": {
            name: sha256(ROOT / name) for name in PACKAGE_SOURCES
        },
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
            **dict(totals),
        },
        "certificate_archive": {
            "path": certificates.name,
            "sha256": sha256(certificates),
            "uncompressed_sha256": raw_certificate_sha256,
            "uncompressed_bytes": raw_certificate_bytes,
            "rejected_graphs": len(archive["rejected_graphs"]),
        },
        "positive_18_control": positive_control(),
        "execution": {
            "command": shlex.join([sys.executable, *sys.argv]),
            "workers": workers,
            "logical_cpus": os.cpu_count(),
            "started_utc": started_utc,
            "finished_utc": datetime.now(UTC).isoformat(),
            "wall_seconds": time.monotonic() - started,
            "python": sys.version,
            "platform": platform.platform(),
            "machine": platform.machine(),
            "git": provenance,
        },
        "graph_results": [
            {key: value for key, value in row.items() if key != "certificate"}
            for row in results
        ],
        "nonclaims": [
            "a survivor is not a realization",
            "different bipartite empty-state and ray-support witnesses may remain",
            "passing support masks do not assert coefficient realizability",
            "this layer does not settle the K6 or dimension-six problem",
        ],
    }
    atomic_json(output, report)
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    parser.add_argument(
        "--output", type=Path,
        default=ROOT / "d6_k6_opposite_ray_rank_conjunction_report.json",
    )
    parser.add_argument(
        "--certificates", type=Path,
        default=(
            ROOT / "d6_k6_opposite_ray_rank_conjunction_certificates.json.gz"
        ),
    )
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("workers must be positive")
    report = run(args.workers, args.output, args.certificates)
    print(json.dumps({
        "status": report["status"],
        "input": report["input_graphs"],
        "rejected": report["graphs_rejected"],
        "surviving": report["graphs_surviving"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
