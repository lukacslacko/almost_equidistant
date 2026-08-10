#!/usr/bin/env python3
"""Exact saturated-singleton-basis extension on the certified K6 residue.

This production layer starts from the frozen 251-graph residue of
``d6_k6_opposite_ray_rank_conjunction.py``.  Inside that package's complete
actual-support DFS, after all older leaf tests pass, it applies one further
necessary condition whenever a nonzero light-ray bin consists of six
singleton supports.  Those six vectors are the coordinate basis.  Every
other outside vector is then reconstructed up to one of two exact quadratic
roots, leaving a finite two-sign CSP.

The test is deliberately inside the actual-support quantifier: failure
rejects only that support leaf, and the DFS continues through every other
allowed actual-support assignment.  Candidate nonedges remain unconstrained,
and an allowed defect coordinate is never required to be nonzero.
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
from pathlib import Path
from typing import Sequence

import d6_k6_opposite_ray_rank_conjunction as parent


ROOT = Path(__file__).resolve().parent
PARENT_REPORT = ROOT / "d6_k6_opposite_ray_rank_conjunction_report.json"
PARENT_CERTIFICATES = (
    ROOT / "d6_k6_opposite_ray_rank_conjunction_certificates.json.gz"
)
PARENT_VERIFICATION = (
    ROOT / "d6_k6_opposite_ray_rank_conjunction_verification.json"
)
REPORT_SCHEMA = "d6-k6-saturated-singleton-basis-v1"
CERTIFICATE_SCHEMA = "d6-k6-saturated-singleton-basis-certificates-v1"
EXPECTED_INPUT = 251
EXPECTED_INPUT_SHA256 = (
    "a23da7bd9ae21eb7c078a08a7bce92f22d397d5c71f1e4f51d9cf8e59bbfa907"
)
EXPECTED_REJECTIONS = 2
EXPECTED_REJECTED_INDICES = (3_138_618, 3_673_988)
EXPECTED_REJECTED_SHA256 = (
    "45ce66ccf4db5fe1a565a576cd16fc719e51b20c6ea9fee858dc3682463e9e84"
)
EXPECTED_RESIDUE_SHA256 = (
    "f274aab3f2d408fce34c7cd7eb4518621fa9d634a8943e3ee1a38d76df9f3d31"
)
EXPECTED_SINGLETON_BRANCH_COUNTS = {
    "saturated_singleton_branches_checked": 291,
    "saturated_singleton_branches_failed": 291,
}
DEFAULT_WORKERS = max(1, (os.cpu_count() or 2) - 1)
PACKAGE_SOURCES = (
    "d6_k6_saturated_singleton_basis.py",
    "verify_d6_k6_saturated_singleton_basis.py",
    "test_d6_k6_saturated_singleton_basis.py",
    "d6_k6_saturated_singleton_basis.md",
)
EXPECTED_DEPENDENCIES = {
    "d6_k6_opposite_ray_rank_conjunction.py": (
        "1306cdcfc3784a4e1cd0e8370edf30c94e7f64ec1564d45b16b777b851d50eaf"
    ),
    "d6_k6_opposite_ray_rank_conjunction_report.json": (
        "3415f48f498758a5d0daac9403012e35186eb7b620fe9c0324e4674f7a11fcc8"
    ),
    "d6_k6_opposite_ray_rank_conjunction_certificates.json.gz": (
        "99a86290da7221f0d2a07ee15ed93432b19dcd58a9ba2bcebaa2956167ac1e37"
    ),
    "d6_k6_opposite_ray_rank_conjunction_verification.json": (
        "e9e64d93faad85af321011b7f48842d13c416f3e61ec802455972d13d48b334e"
    ),
}


@dataclass
class SupportResult:
    feasible: bool
    dfs_nodes: int
    supports: dict[int, int] | None
    singleton_counts: dict[str, int]
    singleton_witnesses: list[dict]


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
    """Write deterministic gzip JSON and return raw hash and byte count."""

    raw = json.dumps(value, indent=2, sort_keys=True).encode("utf-8") + b"\n"
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
    parent.parent.prior.activate_kernel()


def git_provenance() -> dict:
    """Require every package source to equal its launch-commit blob."""

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
        status = git("status", "--porcelain=v1", "--untracked-files=all").strip()
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
    """Load only the frozen, independently verified 251-graph residue."""

    observed = {name: sha256(ROOT / name) for name in EXPECTED_DEPENDENCIES}
    if observed != EXPECTED_DEPENDENCIES:
        raise ValueError(f"opposite-ray dependency boundary changed: {observed}")
    report = json.loads(PARENT_REPORT.read_text(encoding="utf-8"))
    verification = json.loads(PARENT_VERIFICATION.read_text(encoding="utf-8"))
    if (
        report.get("schema") != parent.REPORT_SCHEMA
        or report.get("status") != "COMPLETE"
        or report.get("input_graphs") != 623
        or report.get("graphs_rejected") != 372
        or report.get("graphs_surviving") != EXPECTED_INPUT
        or report.get("ordered_residue_indices_sha256") != EXPECTED_INPUT_SHA256
        or report.get("certificate_archive", {}).get("sha256")
        != sha256(PARENT_CERTIFICATES)
        or verification.get("schema")
        != "d6-k6-opposite-ray-rank-conjunction-verification-v1"
        or verification.get("status") != "PASS"
        or verification.get("graphs_rejected") != 372
        or verification.get("graphs_surviving") != EXPECTED_INPUT
        or verification.get("ordered_residue_indices_sha256")
        != EXPECTED_INPUT_SHA256
    ):
        raise ValueError("frozen opposite-ray result boundary changed")
    parent_records, _, _ = parent.load_input()
    by_index = {int(record["index"]): record for record in parent_records}
    indices = [int(index) for index in report["ordered_residue_indices"]]
    if (
        len(indices) != EXPECTED_INPUT
        or stable_hash(indices) != EXPECTED_INPUT_SHA256
        or any(index not in by_index for index in indices)
    ):
        raise ValueError("certified 251-graph residue changed")
    return [by_index[index] for index in indices], indices, observed


def sign_csp(
    neighbourhoods: Sequence[int], remaining_adjacency: Sequence[int],
) -> tuple[bool, dict]:
    """Exhaust the two exact roots for all reconstructed remaining points."""

    count = len(neighbourhoods)
    if len(remaining_adjacency) != count:
        raise ValueError("neighbourhood and adjacency sizes differ")
    for signs in range(1 << count):
        valid = True
        for right in range(count):
            degree_right = neighbourhoods[right].bit_count()
            sign_right = 1 if signs & (1 << right) else -1
            for left in range(right):
                degree_left = neighbourhoods[left].bit_count()
                sign_left = 1 if signs & (1 << left) else -1
                required = bool(remaining_adjacency[right] & (1 << left))
                if required:
                    if degree_left == degree_right == 0:
                        valid = False
                        break
                    if degree_left == 0 or degree_right == 0:
                        if max(degree_left, degree_right) != 4:
                            valid = False
                            break
                        if sign_left != sign_right:
                            valid = False
                            break
                    else:
                        common = (
                            neighbourhoods[left] & neighbourhoods[right]
                        ).bit_count()
                        if degree_left + degree_right - 2 * common != 4:
                            valid = False
                            break
                        if sign_left == sign_right:
                            valid = False
                            break
                if (
                    neighbourhoods[left] == neighbourhoods[right]
                    and sign_left == sign_right
                ):
                    # Equal root and equal reconstructed u give equal points.
                    valid = False
                    break
            if not valid:
                break
        if valid:
            return True, {
                "signs": [
                    1 if signs & (1 << number) else -1
                    for number in range(count)
                ]
            }
    return False, {
        "neighbourhoods": [int(mask) for mask in neighbourhoods],
        "remaining_adjacency": [int(row) for row in remaining_adjacency],
        "sign_assignments_exhausted": 1 << count,
    }


def saturated_singleton_extension(
    adjacency: Sequence[int], instance, basis: int,
    assignments: dict[int, int],
) -> tuple[bool, dict]:
    """Check the exact extension CSP for one six-singleton ray basis."""

    basis_locals = tuple(parent.parent.kernel.vertices(basis))
    if len(basis_locals) != parent.parent.kernel.COORDINATES:
        raise ValueError("basis must have six local vertices")
    coordinate_of: dict[int, int] = {}
    used_coordinates = set()
    for local in basis_locals:
        support = assignments[local]
        if support.bit_count() != 1:
            return True, {"applicable": False}
        coordinate = support.bit_length() - 1
        if coordinate in used_coordinates:
            return False, {"reason": "duplicate_singleton_basis_coordinate"}
        coordinate_of[local] = coordinate
        used_coordinates.add(coordinate)
    if len(used_coordinates) != parent.parent.kernel.COORDINATES:
        return False, {"reason": "singleton_bin_is_not_coordinate_basis"}

    remaining = tuple(
        local for local in range(len(instance.outside))
        if local not in basis_locals
    )
    neighbourhoods: list[int] = []
    for local in remaining:
        absolute = instance.outside[local]
        reconstructed = 0
        for basis_local in basis_locals:
            if adjacency[absolute] & (1 << instance.outside[basis_local]):
                reconstructed |= 1 << coordinate_of[basis_local]
        if reconstructed and reconstructed & ~instance.defects[local]:
            return False, {
                "reason": "reconstructed_support_not_allowed_by_seed_defects",
                "vertex": absolute,
                "reconstructed_support": parent.parent.kernel.vertices(
                    reconstructed
                ),
                "allowed_support": parent.parent.kernel.vertices(
                    instance.defects[local]
                ),
            }
        if local in assignments and assignments[local] != reconstructed:
            return False, {
                "reason": "assigned_light_support_differs_from_reconstruction",
                "vertex": absolute,
                "reconstructed_support": parent.parent.kernel.vertices(
                    reconstructed
                ),
                "assigned_support": parent.parent.kernel.vertices(
                    assignments[local]
                ),
            }
        neighbourhoods.append(reconstructed)

    position = {local: number for number, local in enumerate(remaining)}
    remaining_rows = []
    for local in remaining:
        absolute = instance.outside[local]
        remaining_rows.append(sum(
            1 << position[other]
            for other in remaining
            if adjacency[absolute] & (1 << instance.outside[other])
        ))
    feasible, detail = sign_csp(neighbourhoods, remaining_rows)
    detail.update({
        "applicable": True,
        "basis_vertices_by_coordinate": {
            str(coordinate): instance.outside[local]
            for local, coordinate in coordinate_of.items()
        },
        "remaining_vertices": [instance.outside[local] for local in remaining],
        "remaining_basis_neighbourhoods": {
            str(instance.outside[local]): parent.parent.kernel.vertices(
                neighbourhoods[number]
            )
            for number, local in enumerate(remaining)
        },
    })
    return feasible, detail


def singleton_basis_leaf_passes(
    adjacency: Sequence[int], instance, z0: int,
    bins: tuple[int, int], assignments: dict[int, int],
    counts: Counter[str], witnesses: list[dict],
) -> bool:
    """Apply every applicable basis test at one complete support leaf."""

    if z0:
        return True
    for basis in bins:
        if basis.bit_count() != parent.parent.kernel.COORDINATES:
            continue
        if any(
            assignments[local].bit_count() != 1
            for local in parent.parent.kernel.vertices(basis)
        ):
            continue
        counts["saturated_singleton_branches_checked"] += 1
        feasible, detail = saturated_singleton_extension(
            adjacency, instance, basis, assignments
        )
        if feasible:
            counts["saturated_singleton_branches_passed"] += 1
            continue
        counts["saturated_singleton_branches_failed"] += 1
        if len(witnesses) < 8:
            witnesses.append(detail)
        return False
    return True


def tight_actual_supports_for_bins(
    adjacency: Sequence[int], instance, z0: int, bins: tuple[int, int],
) -> SupportResult:
    """Parent actual-support DFS plus the leaf-local singleton-basis rule."""

    if bins[0] & bins[1] != z0:
        raise ValueError("light-ray bins intersect outside Z0")
    involved = bins[0] | bins[1]
    domains = {
        local: parent.support_domains(
            instance.defects[local], bool(z0 & (1 << local))
        )
        for local in parent.parent.kernel.vertices(involved)
    }
    peers = {}
    for local in parent.parent.kernel.vertices(involved):
        peer_mask = 0
        for bin_mask in bins:
            if bin_mask & (1 << local):
                peer_mask |= bin_mask & ~(1 << local)
        peers[local] = peer_mask
    assignments: dict[int, int] = {}
    nodes = 0
    singleton_counts: Counter[str] = Counter()
    singleton_witnesses: list[dict] = []

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
                    parent.tight_orthogonal_supports_pass(
                        bin_mask, assignments
                    )
                    for bin_mask in bins
                )
                and parent.cross_component_intersection_supports_pass(
                    adjacency, instance, z0, bins, assignments
                )
                and singleton_basis_leaf_passes(
                    adjacency, instance, z0, bins, assignments,
                    singleton_counts, singleton_witnesses,
                )
            )
        best = None
        best_choices = None
        best_key = None
        for local in parent.parent.kernel.vertices(unassigned):
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
                    parent.assigned_hall_passes(
                        bin_mask, assignments, len(instance.outside)
                    )
                    and parent.tight_orthogonal_supports_pass(
                        bin_mask, assignments
                    )
                    for bin_mask in bins
                )
                and parent.opposite_required_edge_supports_pass(
                    adjacency, instance, z0, bins, assignments
                )
                and visit(unassigned ^ bit)
            ):
                return True
            assignments.pop(best, None)
        return False

    feasible = visit(involved)
    return SupportResult(
        feasible=feasible,
        dfs_nodes=nodes,
        supports=assignments.copy() if feasible else None,
        singleton_counts=dict(singleton_counts),
        singleton_witnesses=singleton_witnesses,
    )


def light_ray_decision(
    adjacency: Sequence[int], instance, z0: int, components,
) -> tuple[dict | None, Counter[str], dict, Counter[str], list[dict]]:
    """Conjoin the new support-leaf rule with the frozen coloring checks."""

    odd = [component.component for component in components if not component.bipartite]
    raw_colorings = 1 if not odd else 1 << (len(odd) - 1)
    counts: Counter[str] = Counter()
    singleton_counts: Counter[str] = Counter()
    singleton_witnesses: list[dict] = []
    failures = []
    for colors, bins in parent.parent._pure_colorings(instance, z0, odd):
        counts["matchable_colorings"] += 1
        checks = parent.coloring_projection_checks(adjacency, instance, z0, bins)
        for check in checks:
            counts[
                f"directed_basis_size_{check['basis_size_including_Z0']}"
            ] += 1
        rendered_bins = [
            [
                instance.outside[local]
                for local in parent.parent.kernel.vertices(bin_mask)
            ]
            for bin_mask in bins
        ]
        if parent.both_saturated_with_light(z0, bins):
            counts["both_saturated_arithmetic_failed_colorings"] += 1
            failures.append({
                "component_colors": list(colors),
                "bins": rendered_bins,
                "reason": "both_saturated_light_ray_arithmetic",
                "Z0_size": z0.bit_count(),
                "light_vertices_per_bin": (
                    parent.parent.kernel.COORDINATES - z0.bit_count()
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
                "bins": rendered_bins,
                "reason": "opposite_ray_projection_rank",
                "projection_checks": list(checks),
            })
            continue
        counts["projection_rank_passed_colorings"] += 1
        support = tight_actual_supports_for_bins(adjacency, instance, z0, bins)
        counts["actual_support_searches"] += 1
        counts["actual_support_dfs_nodes"] += support.dfs_nodes
        singleton_counts.update(support.singleton_counts)
        if len(singleton_witnesses) < 8:
            singleton_witnesses.extend(
                support.singleton_witnesses[: 8 - len(singleton_witnesses)]
            )
        if support.feasible:
            counts["joint_coloring_passed"] += 1
            return ({
                "component_colors": list(colors),
                "bins": rendered_bins,
                "projection_checks": list(checks),
            }, counts, {
                "raw_colorings": raw_colorings,
                "matchable_colorings": counts["matchable_colorings"],
                "failures": failures,
            }, singleton_counts, singleton_witnesses)
        counts["actual_support_failed_colorings"] += 1
        failures.append({
            "component_colors": list(colors),
            "bins": rendered_bins,
            "reason": "joint_actual_support",
            "projection_checks": list(checks),
        })
    return None, counts, {
        "raw_colorings": raw_colorings,
        "matchable_colorings": counts["matchable_colorings"],
        "failures": failures,
    }, singleton_counts, singleton_witnesses


def solve_seed(adjacency: Sequence[int], instance) -> dict:
    inertia = parent.parent.kernel.InertiaCache()
    forcing = parent.parent.kernel.ZeroForcingSolver()
    subset_rank = parent.parent.kernel.arbitrary.PositiveSubsetRank(
        inertia, forcing
    )
    counts: Counter[str] = Counter()
    singleton_counts: Counter[str] = Counter()
    singleton_witnesses: list[dict] = []
    failures = []
    for z0 in parent.parent.kernel.ranks.z0_subsets(instance.eligible_z0_mask):
        counts["z0_considered"] += 1
        zvertices = [
            instance.outside[local]
            for local in parent.parent.kernel.vertices(z0)
        ]
        if parent.parent.kernel.support_matching(z0, instance.defects) is None:
            counts["z0_unmatchable"] += 1
            failures.append({"Z0": zvertices, "reason": "Z0_matching"})
            continue
        states, components, local = parent.parent.bipartite_states(
            adjacency, instance, z0, inertia, forcing, subset_rank
        )
        counts.update(local)
        if not states:
            counts["z0_bipartite_failed"] += 1
            failures.append({"Z0": zvertices, "reason": "tight_bipartite_system"})
            continue
        counts["z0_bipartite_passed"] += 1
        witness, local, coloring_certificate, local_singleton, local_witnesses = (
            light_ray_decision(adjacency, instance, z0, components)
        )
        counts.update(local)
        singleton_counts.update(local_singleton)
        if len(singleton_witnesses) < 8:
            singleton_witnesses.extend(
                local_witnesses[: 8 - len(singleton_witnesses)]
            )
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
            "singleton_counts": dict(singleton_counts),
            "singleton_witnesses": singleton_witnesses,
        }
    return {
        "feasible": False,
        "counts": dict(counts),
        "failures": failures,
        "singleton_counts": dict(singleton_counts),
        "singleton_witnesses": singleton_witnesses,
    }


def evaluate_record(record: dict) -> dict:
    parent.parent.prior.activate_kernel()
    adjacency = tuple(map(int, record["adjacency"]))
    total: Counter[str] = Counter()
    singleton_total: Counter[str] = Counter()
    singleton_witnesses: list[dict] = []
    seeds_checked = 0
    for seed_mask in parent.parent.kernel.clique_masks(
        adjacency, parent.parent.kernel.COORDINATES
    ):
        seeds_checked += 1
        instance = parent.parent.kernel.build_instance(
            adjacency, parent.parent.kernel.vertices(seed_mask)
        )
        decision = solve_seed(adjacency, instance)
        total.update(decision["counts"])
        singleton_total.update(decision["singleton_counts"])
        if len(singleton_witnesses) < 8:
            singleton_witnesses.extend(
                decision["singleton_witnesses"][: 8 - len(singleton_witnesses)]
            )
        if decision["feasible"]:
            continue
        return {
            "index": int(record["index"]),
            "rejected": True,
            "seeds_checked": seeds_checked,
            "first_impossible_seed_mask": seed_mask,
            "first_impossible_seed": list(instance.seed),
            "counts": dict(total),
            "singleton_basis_counts": dict(singleton_total),
            "singleton_basis_witnesses": singleton_witnesses,
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
        "singleton_basis_counts": dict(singleton_total),
        "singleton_basis_witnesses": singleton_witnesses,
        "certificate": None,
    }


def synthetic_controls() -> dict:
    neighbourhoods = (0b000011, 0b001100, 0b110000)
    triangle = (0b110, 0b101, 0b011)
    feasible, _ = sign_csp(neighbourhoods, triangle)
    if feasible:
        raise AssertionError("odd opposite-sign cycle control passed")
    path = (0b010, 0b101, 0b010)
    feasible, witness = sign_csp(neighbourhoods, path)
    if not feasible:
        raise AssertionError("bipartite opposite-sign path control failed")
    feasible, _ = sign_csp((0, 0b001111), (0b10, 0b01))
    if not feasible:
        raise AssertionError("empty/degree-four required-edge control failed")
    feasible, _ = sign_csp((0, 0b000111), (0b10, 0b01))
    if feasible:
        raise AssertionError("empty/degree-three required-edge control passed")
    return {
        "odd_opposite_sign_cycle_fails": True,
        "opposite_sign_path_passes": witness,
        "empty_degree_four_same_sign_passes": True,
        "empty_degree_three_edge_fails": True,
    }


def positive_control() -> dict:
    parent.parent.prior.activate_kernel()
    adjacency = parent.parent.kernel.lower_bound_18_graph()
    passed = 0
    for seed_mask in parent.parent.kernel.clique_masks(
        adjacency, parent.parent.kernel.COORDINATES
    ):
        instance = parent.parent.kernel.build_instance(
            adjacency, parent.parent.kernel.vertices(seed_mask)
        )
        if not solve_seed(adjacency, instance)["feasible"]:
            raise AssertionError("singleton-basis rule rejects positive 18 control")
        passed += 1
    if passed != 32:
        raise AssertionError("positive-control K6 seed count changed")
    return {"passed": True, "K6_seeds": passed}


def run(workers: int, output: Path, certificates: Path) -> dict:
    provenance = git_provenance()
    records, indices, dependencies = load_input()
    controls = synthetic_controls()
    started_utc = datetime.now(UTC).isoformat()
    started = time.monotonic()
    if workers == 1:
        results = [evaluate_record(record) for record in records]
    else:
        with ProcessPoolExecutor(
            max_workers=workers, initializer=worker_initializer
        ) as pool:
            results = list(pool.map(evaluate_record, records, chunksize=1))
    rejected = [int(row["index"]) for row in results if row["rejected"]]
    residue = [int(row["index"]) for row in results if not row["rejected"]]
    singleton_totals: Counter[str] = Counter()
    parent_totals: Counter[str] = Counter()
    for row in results:
        parent_totals.update(row["counts"])
        singleton_totals.update(row["singleton_basis_counts"])
    if (
        rejected != list(EXPECTED_REJECTED_INDICES)
        or stable_hash(rejected) != EXPECTED_REJECTED_SHA256
        or stable_hash(residue) != EXPECTED_RESIDUE_SHA256
        or dict(singleton_totals) != EXPECTED_SINGLETON_BRANCH_COUNTS
    ):
        raise AssertionError("production result differs from source-bound census")

    package_hashes = {name: sha256(ROOT / name) for name in PACKAGE_SOURCES}
    archive = {
        "schema": CERTIFICATE_SCHEMA,
        "production_source_sha256": sha256(Path(__file__)),
        "source_sha256": package_hashes,
        "ordered_input_indices_sha256": EXPECTED_INPUT_SHA256,
        "rejected_graphs": [
            {"index": row["index"], **row["certificate"]}
            for row in results if row["rejected"]
        ],
    }
    raw_certificate_sha256, raw_certificate_bytes = atomic_gzip_json(
        certificates, archive
    )
    positive = positive_control()
    report = {
        "schema": REPORT_SCHEMA,
        "status": "COMPLETE",
        "description": (
            "exact saturated-singleton K6 ray-basis reconstruction and "
            "two-root extension CSP inside the exhaustive actual-support DFS"
        ),
        "semantics": {
            "candidate_nonedges": "unconstrained_and_may_be_unit",
            "allowed_supports": "upper_bounds_actual_subsets_exhausted",
            "actual_support_quantifier": (
                "basis_rule_prunes_only_one_complete_support_leaf"
            ),
            "basis_trigger": "Z0_empty_and_one_ray_bin_has_six_singletons",
            "remaining_coefficients": (
                "exact_roots_t=(1+epsilon*sqrt(7))/3"
            ),
            "arithmetic": "exact_integer_bitmasks_and_exhaustive_two_sign_CSP",
            "distinctness": "equal_neighbourhood_equal_root_collisions_forbidden",
        },
        "production_source_sha256": sha256(Path(__file__)),
        "source_sha256": package_hashes,
        "dependencies": dependencies,
        "discovery_boundary": {
            "probe_source_sha256": (
                "b829259adc7270ae7f8a235335ef43e027fbea4f99ca4a2314ef10e09954cbed"
            ),
            "pilot_report_sha256": (
                "822c8051c371eb3b09bae99e90b62ba5c5ba1f443b1104804c400714223dc09d"
            ),
        },
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
            **dict(parent_totals),
            **dict(singleton_totals),
        },
        "certificate_archive": {
            "path": certificates.name,
            "sha256": sha256(certificates),
            "uncompressed_sha256": raw_certificate_sha256,
            "uncompressed_bytes": raw_certificate_bytes,
            "rejected_graphs": len(archive["rejected_graphs"]),
        },
        "positive_18_control": positive,
        "synthetic_controls": controls,
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
            "passing graph nonedges are not asserted to be nonunit",
            "the two-root CSP is invoked only at saturated singleton leaves",
            "this layer does not settle the K6 residue or dimension six",
        ],
    }
    atomic_json(output, report)
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    parser.add_argument(
        "--output", type=Path,
        default=ROOT / "d6_k6_saturated_singleton_basis_report.json",
    )
    parser.add_argument(
        "--certificates", type=Path,
        default=ROOT / "d6_k6_saturated_singleton_basis_certificates.json.gz",
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
