#!/usr/bin/env python3
"""Independent verifier for the saturated-singleton K6 basis layer.

The checker imports neither the new producer nor its discovery probe.  It
starts from the frozen independent verifier for the 251-graph opposite-ray
residue, independently replays that verifier's actual-support DFS, and adds
the new rule only at complete actual-support leaves.

For a saturated six-vertex light-ray bin with six distinct singleton actual
supports, the remaining seven vertices have reconstructed basis
neighbourhoods.  The producer enumerates their two algebraic signs directly;
this verifier instead translates required-edge and collision conditions to
XOR equations and solves them with parity union-find.  Candidate nonedges are
never required to be nonunit, and every allowed defect coordinate may vanish.
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
from typing import Sequence

import verify_d6_k6_opposite_ray_rank_conjunction as parent_verifier


base = parent_verifier.base
zf_reference = parent_verifier.zf_reference
ROOT = Path(__file__).resolve().parent
REPORT_SCHEMA = "d6-k6-saturated-singleton-basis-v1"
CERTIFICATE_SCHEMA = "d6-k6-saturated-singleton-basis-certificates-v1"
VERIFICATION_SCHEMA = "d6-k6-saturated-singleton-basis-verification-v1"
EXPECTED_INPUT = 251
EXPECTED_INPUT_SHA256 = (
    "a23da7bd9ae21eb7c078a08a7bce92f22d397d5c71f1e4f51d9cf8e59bbfa907"
)
EXPECTED_REJECTIONS = 2
EXPECTED_REJECTED_INDICES = [3138618, 3673988]
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
EXPECTED_PRODUCTION_DEPENDENCIES = {
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
EXPECTED_PARENT_ARTIFACTS = {
    **EXPECTED_PRODUCTION_DEPENDENCIES,
    "verify_d6_k6_opposite_ray_rank_conjunction.py": (
        "870c59dd1801b054814bf250bc07effc22e2b5a1f4c3c299af7d9b7b376b94d1"
    ),
}
EXPECTED_PRODUCTION_SYNTHETIC_CONTROLS = {
    "odd_opposite_sign_cycle_fails": True,
    "opposite_sign_path_passes": {"signs": [-1, 1, -1]},
    "empty_degree_four_same_sign_passes": True,
    "empty_degree_three_edge_fails": True,
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_hash(value: object) -> str:
    raw = json.dumps(value, sort_keys=True, separators=(",", ":")).encode(
        "ascii"
    )
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
    parent_verifier.worker_initializer()


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
        "d6_k6_saturated_singleton_basis",
        "probe_d6_k6_saturated_singleton_basis",
        "probe_d6_k6_opposite_ray_rank_conjunction",
        "d6_k6_opposite_ray_rank_conjunction",
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
        (
            ROOT / "d6_k6_opposite_ray_rank_conjunction_report.json"
        ).read_text(encoding="utf-8")
    )
    verification = json.loads(
        (
            ROOT / "d6_k6_opposite_ray_rank_conjunction_verification.json"
        ).read_text(encoding="utf-8")
    )
    if (
        report.get("schema") != "d6-k6-opposite-ray-rank-conjunction-v1"
        or report.get("status") != "COMPLETE"
        or report.get("graphs_rejected") != 372
        or report.get("graphs_surviving") != EXPECTED_INPUT
        or report.get("ordered_residue_indices_sha256")
        != EXPECTED_INPUT_SHA256
        or report.get("certificate_archive", {}).get("sha256")
        != sha256(
            ROOT / "d6_k6_opposite_ray_rank_conjunction_certificates.json.gz"
        )
        or verification.get("schema")
        != "d6-k6-opposite-ray-rank-conjunction-verification-v1"
        or verification.get("status") != "PASS"
        or verification.get("graphs_rejected") != 372
        or verification.get("graphs_surviving") != EXPECTED_INPUT
        or verification.get("ordered_residue_indices_sha256")
        != EXPECTED_INPUT_SHA256
    ):
        raise ValueError("frozen opposite-ray result boundary changed")
    parent_records, _ = parent_verifier.load_input()
    by_index = {int(record["index"]): record for record in parent_records}
    indices = [int(index) for index in report["ordered_residue_indices"]]
    if (
        len(indices) != EXPECTED_INPUT
        or stable_hash(indices) != EXPECTED_INPUT_SHA256
        or any(index not in by_index for index in indices)
    ):
        raise ValueError("independent 251-graph input reconstruction changed")
    return [by_index[index] for index in indices], indices


class ParityUnionFind:
    """Union-find for equations ``sign[left] XOR sign[right] = parity``."""

    def __init__(self, size: int) -> None:
        self.parent = list(range(size))
        self.rank = [0] * size
        self.to_parent = [0] * size

    def find(self, item: int) -> tuple[int, int]:
        parent = self.parent[item]
        if parent == item:
            return item, 0
        root, above = self.find(parent)
        self.to_parent[item] ^= above
        self.parent[item] = root
        return root, self.to_parent[item]

    def constrain(self, left: int, right: int, parity: int) -> bool:
        left_root, left_value = self.find(left)
        right_root, right_value = self.find(right)
        if left_root == right_root:
            return left_value ^ right_value == parity
        if self.rank[left_root] > self.rank[right_root]:
            left_root, right_root = right_root, left_root
        self.parent[left_root] = right_root
        self.to_parent[left_root] = left_value ^ right_value ^ parity
        if self.rank[left_root] == self.rank[right_root]:
            self.rank[right_root] += 1
        return True

    def witness(self) -> list[int]:
        """Choose root signs zero and return one bit for each vertex."""

        return [self.find(item)[1] for item in range(len(self.parent))]


def parity_constraints_feasible(
    vertex_count: int, constraints: Sequence[tuple[int, int, int]],
) -> tuple[bool, list[int] | None]:
    system = ParityUnionFind(vertex_count)
    for left, right, parity in constraints:
        if parity not in (0, 1):
            raise ValueError("parity must be zero or one")
        if not system.constrain(left, right, parity):
            return False, None
    return True, system.witness()


def reconstructed_sign_system(
    neighbourhoods: Sequence[int], remaining_adjacency: Sequence[int],
) -> tuple[bool, dict]:
    """Check the exact two-root system by XOR propagation, not enumeration."""

    count = len(neighbourhoods)
    if len(remaining_adjacency) != count:
        raise ValueError("neighbourhood and adjacency lengths differ")
    constraints: list[tuple[int, int, int]] = []
    for right in range(count):
        degree_right = neighbourhoods[right].bit_count()
        for left in range(right):
            degree_left = neighbourhoods[left].bit_count()
            required = bool(remaining_adjacency[right] & (1 << left))
            if required:
                if degree_left == degree_right == 0:
                    return False, {
                        "reason": "required_edge_between_empty_apices",
                        "pair": [left, right],
                    }
                if degree_left == 0 or degree_right == 0:
                    nonempty_degree = max(degree_left, degree_right)
                    if nonempty_degree != 4:
                        return False, {
                            "reason": "empty_apex_edge_has_wrong_basis_degree",
                            "pair": [left, right],
                            "nonempty_degree": nonempty_degree,
                        }
                    # Empty-apex delta equals the nonempty quadratic sign.
                    constraints.append((left, right, 0))
                else:
                    common = (
                        neighbourhoods[left] & neighbourhoods[right]
                    ).bit_count()
                    if degree_left + degree_right - 2 * common != 4:
                        return False, {
                            "reason": "nonempty_edge_neighbourhood_equation",
                            "pair": [left, right],
                            "degree_left": degree_left,
                            "degree_right": degree_right,
                            "common": common,
                        }
                    # The two roots must have opposite sqrt(7) signs.
                    constraints.append((left, right, 1))
            if neighbourhoods[left] == neighbourhoods[right]:
                # Equal reconstructed support and equal root would be a
                # collision, including two copies of the same empty apex.
                constraints.append((left, right, 1))
    feasible, signs = parity_constraints_feasible(count, constraints)
    if not feasible:
        return False, {
            "reason": "inconsistent_sign_parities",
            "constraints": [list(row) for row in constraints],
        }
    assert signs is not None
    return True, {
        "sign_bits": signs,
        "constraints": [list(row) for row in constraints],
    }


def saturated_singleton_extension(
    adjacency: Sequence[int], instance, basis: int,
    assignments: dict[int, int],
) -> tuple[bool, dict]:
    """Reconstruct all seven other points from a singleton coordinate basis."""

    basis_locals = tuple(base.bits(basis))
    if len(basis_locals) != base.COORDINATES:
        raise ValueError("basis must have six local vertices")
    coordinate_of: dict[int, int] = {}
    seen_coordinates = set()
    for local in basis_locals:
        support = assignments[local]
        if support.bit_count() != 1:
            return True, {"applicable": False}
        coordinate = support.bit_length() - 1
        if coordinate in seen_coordinates:
            return False, {
                "applicable": True,
                "reason": "duplicate_singleton_basis_coordinate",
            }
        coordinate_of[local] = coordinate
        seen_coordinates.add(coordinate)
    if len(seen_coordinates) != base.COORDINATES:
        return False, {
            "applicable": True,
            "reason": "singleton_bin_is_not_coordinate_basis",
        }

    basis_set = set(basis_locals)
    remaining = tuple(
        local
        for local in range(len(instance.outside))
        if local not in basis_set
    )
    neighbourhoods = []
    for local in remaining:
        absolute = instance.outside[local]
        support = 0
        for basis_local in basis_locals:
            if adjacency[absolute] & (1 << instance.outside[basis_local]):
                support |= 1 << coordinate_of[basis_local]
        if support & ~instance.defects[local]:
            return False, {
                "applicable": True,
                "reason": "reconstructed_support_not_allowed_by_seed_defects",
                "vertex": absolute,
                "reconstructed_support": list(base.bits(support)),
                "allowed_support": list(base.bits(instance.defects[local])),
            }
        if local in assignments and assignments[local] != support:
            return False, {
                "applicable": True,
                "reason": "assigned_light_support_differs_from_reconstruction",
                "vertex": absolute,
                "reconstructed_support": list(base.bits(support)),
                "assigned_support": list(base.bits(assignments[local])),
            }
        neighbourhoods.append(support)

    position = {local: number for number, local in enumerate(remaining)}
    remaining_rows = []
    for local in remaining:
        absolute = instance.outside[local]
        row = 0
        for other in remaining:
            if adjacency[absolute] & (1 << instance.outside[other]):
                row |= 1 << position[other]
        remaining_rows.append(row)
    feasible, detail = reconstructed_sign_system(
        neighbourhoods, remaining_rows
    )
    detail.update({
        "applicable": True,
        "basis_vertices_by_coordinate": {
            str(coordinate): instance.outside[local]
            for local, coordinate in coordinate_of.items()
        },
        "remaining_vertices": [instance.outside[local] for local in remaining],
        "remaining_basis_neighbourhoods": {
            str(instance.outside[local]): list(base.bits(neighbourhoods[number]))
            for number, local in enumerate(remaining)
        },
    })
    return feasible, detail


def saturated_singleton_leaves_pass(
    adjacency: Sequence[int], instance, z0: int, bins: tuple[int, int],
    assignments: dict[int, int],
) -> bool:
    """Apply the new theorem only to completed qualifying support branches."""

    if z0:
        return True
    for basis in bins:
        if basis.bit_count() != base.COORDINATES:
            continue
        if any(
            assignments[local].bit_count() != 1
            for local in base.bits(basis)
        ):
            continue
        feasible, _ = saturated_singleton_extension(
            adjacency, instance, basis, assignments
        )
        if not feasible:
            return False
    return True


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


def independent_actual_support_exists(adjacency, instance, z0, bins) -> bool:
    """Replay every actual-support alternative before rejecting a coloring."""

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
    assignments: dict[int, int] = {}

    def choices(local: int) -> list[int]:
        return [
            support
            for support in domains[local]
            if all(
                (support & other_support).bit_count() != 1
                for other, other_support in assignments.items()
                if peers[local] & (1 << other)
            )
        ]

    def visit(unassigned: int) -> bool:
        if not unassigned:
            if not all(
                tight_supports_pass(bin_mask, assignments)
                for bin_mask in bins
            ):
                return False
            if not parent_verifier.independent_cross_intersection_pass(
                adjacency, instance, z0, bins, assignments
            ):
                return False
            return saturated_singleton_leaves_pass(
                adjacency, instance, z0, bins, assignments
            )
        ranked = []
        for local in base.bits(unassigned):
            local_choices = choices(local)
            # Deliberately retain the independent verifier's ordering, which
            # differs from the producer's peer-degree tie breaker.
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
                and parent_verifier.independent_cross_edge_local_pass(
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
        if not all(
            base.hall_matchable(bin_mask, instance.defects)
            for bin_mask in pair
        ):
            continue
        checks = (
            parent_verifier.independent_projection_check(
                adjacency, instance, z0, pair[0], pair[1]
            ),
            parent_verifier.independent_projection_check(
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
        if independent_actual_support_exists(adjacency, instance, z0, pair):
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
        states, local_components = parent_verifier.prior_verifier.bipartite_states(
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
    parent_verifier.prior_verifier.prior_verifier.activate_kernel()
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
            raise AssertionError(
                "independent saturated-singleton positive control failed"
            )
        passed += 1
    if passed != 32:
        raise AssertionError("independent positive-control seed count changed")
    return {"passed": True, "K6_seeds": passed}


def synthetic_controls() -> dict:
    neighbourhoods = (0b000011, 0b001100, 0b110000)
    triangle = (0b110, 0b101, 0b011)
    feasible, _ = reconstructed_sign_system(neighbourhoods, triangle)
    if feasible:
        raise AssertionError("odd opposite-sign cycle control passed")
    path = (0b010, 0b101, 0b010)
    feasible, witness = reconstructed_sign_system(neighbourhoods, path)
    if not feasible:
        raise AssertionError("opposite-sign path control failed")
    feasible, _ = reconstructed_sign_system(
        (0, 0b001111), (0b10, 0b01)
    )
    if not feasible:
        raise AssertionError("empty/degree-four required-edge control failed")
    feasible, _ = reconstructed_sign_system(
        (0, 0b000111), (0b10, 0b01)
    )
    if feasible:
        raise AssertionError("empty/degree-three required-edge control passed")
    feasible, _ = reconstructed_sign_system((0, 0), (0, 0))
    if not feasible:
        raise AssertionError("two empty apices with opposite signs should pass")
    feasible, _ = reconstructed_sign_system((0, 0), (0b10, 0b01))
    if feasible:
        raise AssertionError("required edge between empty apices passed")
    return {
        "odd_opposite_sign_cycle_fails": True,
        "opposite_sign_path_passes": witness["sign_bits"],
        "empty_degree_four_same_sign_passes": True,
        "empty_degree_three_edge_fails": True,
        "duplicate_empty_apices_require_opposite_signs": True,
        "required_empty_apex_edge_fails": True,
    }


def verify(report_path, certificates_path, output, workers):
    assert_import_independence()
    records, indices = load_input()
    report = json.loads(report_path.read_text(encoding="utf-8"))
    archive, raw_sha256, raw_bytes = load_gzip_json(certificates_path)
    source_path = ROOT / "d6_k6_saturated_singleton_basis.py"
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
        or report.get("rejected_indices") != EXPECTED_REJECTED_INDICES
        or report.get("rejected_indices_sha256") != EXPECTED_REJECTED_SHA256
        or report.get("ordered_residue_indices_sha256")
        != EXPECTED_RESIDUE_SHA256
        or {
            key: report.get("totals", {}).get(key)
            for key in EXPECTED_SINGLETON_BRANCH_COUNTS
        }
        != EXPECTED_SINGLETON_BRANCH_COUNTS
        or report.get("certificate_archive", {}).get("sha256")
        != sha256(certificates_path)
        or report.get("certificate_archive", {}).get("uncompressed_sha256")
        != raw_sha256
        or report.get("certificate_archive", {}).get("uncompressed_bytes")
        != raw_bytes
        or report.get("certificate_archive", {}).get("rejected_graphs")
        != EXPECTED_REJECTIONS
        or report.get("positive_18_control")
        != {"passed": True, "K6_seeds": 32}
        or report.get("synthetic_controls")
        != EXPECTED_PRODUCTION_SYNTHETIC_CONTROLS
        or archive.get("schema") != CERTIFICATE_SCHEMA
        or archive.get("production_source_sha256") != sha256(source_path)
        or archive.get("source_sha256") != package_sources
        or archive.get("ordered_input_indices_sha256") != EXPECTED_INPUT_SHA256
    ):
        raise ValueError("saturated-singleton production boundary mismatch")
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
        rejected != EXPECTED_REJECTED_INDICES
        or len(rejected) != EXPECTED_REJECTIONS
        or stable_hash(rejected) != EXPECTED_REJECTED_SHA256
        or stable_hash(residue) != EXPECTED_RESIDUE_SHA256
        or rejected != report.get("rejected_indices")
        or residue != report.get("ordered_residue_indices")
        or len(rejected) != report.get("graphs_rejected")
        or len(residue) != report.get("graphs_surviving")
    ):
        raise AssertionError("independent saturated-singleton decisions differ")
    saved_rows = {int(row["index"]): row for row in report["graph_results"]}
    certificates = {
        int(row["index"]): row for row in archive["rejected_graphs"]
    }
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
    algebra_controls = synthetic_controls()
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
            "frozen_251_graph_parent_residue": True,
            "independent_parent_actual_support_DFS_replay": True,
            "actual_support_alternatives_exhausted": True,
            "optional_allowed_support_zeros_preserved": True,
            "saturated_rule_only_at_complete_support_leaves": True,
            "independent_basis_neighbourhood_reconstruction": True,
            "independent_parity_union_find_sign_CSP": True,
            "required_nonempty_edge_equation": True,
            "required_empty_apex_edge_equation": True,
            "distinct_point_collision_constraints": True,
            "all_graph_decisions": True,
            "all_rejected_Z0_and_coloring_partitions": True,
            "positive_18_control": True,
            "synthetic_algebra_controls": True,
        },
        "provenance": provenance,
        "graphs_recomputed": len(results),
        "graphs_rejected": len(rejected),
        "graphs_surviving": len(residue),
        "rejected_indices": rejected,
        "rejected_indices_sha256": stable_hash(rejected),
        "ordered_residue_indices_sha256": stable_hash(residue),
        "positive_18_control": control,
        "synthetic_controls": algebra_controls,
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
        "--report",
        type=Path,
        default=ROOT / "d6_k6_saturated_singleton_basis_report.json",
    )
    parser.add_argument(
        "--certificates",
        type=Path,
        default=(
            ROOT / "d6_k6_saturated_singleton_basis_certificates.json.gz"
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "d6_k6_saturated_singleton_basis_verification.json",
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
