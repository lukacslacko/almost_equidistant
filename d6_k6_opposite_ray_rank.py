#!/usr/bin/env python3
"""Exact opposite-light-ray projection-rank filter on the certified K6-623.

Fix a required K6 and one possible actual zero Lorentz-factor set ``Z0``.
Every connected nonbipartite component of the disjoint-defect Lorentz graph
after deleting ``Z0`` must lie on one of the two light rays.  In each ray bin
the defect vectors, together with ``Z0``, are orthonormal in R^6.

For opposite nonzero ray sets A and B, let X be their cross Gram block and
project A off the orthonormal bin ``B union Z0``.  The residual Gram matrix

    R = I_A - X X^T

is PSD of rank at most ``6-|B union Z0|``.  After diagonal congruence its
strict negative off-diagonal graph joins two A vertices exactly when they
share a required neighbour in B.  Every nontrivial connected component H
therefore contributes at least ``|H|-1`` to rank.  A source vertex with no
required B-neighbour has a zero row in X and contributes one more.  The
resulting necessary inequality is checked in both directions, for every
symmetry-reduced assignment of all nonbipartite components, and for every
eligible ``Z0``.

Only required unit edges, alpha(G)<=2, and allowed support masks are used.
Candidate nonedges remain unconstrained and allowed support entries may be
zero.  The computation uses only integer bit masks.
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
from datetime import UTC, datetime
from itertools import combinations
from pathlib import Path
from typing import Iterator, Sequence

from d6_k6_lorentz_reference import (
    COORDINATES,
    build_instance,
    clique_masks,
    nonbipartite_components,
    support_matching,
    validate_graph,
    vertices,
)
from test_d6_k6_lorentz import lower_bound_18_graph


ROOT = Path(__file__).resolve().parent
INPUT_MANIFEST = ROOT / "d6_current_residue_manifest_v6.json"
PARENT_REPORT = ROOT / "d6_k6_tight_same_z0_report.json"
PARENT_VERIFICATION = ROOT / "d6_k6_tight_same_z0_verification.json"

REPORT_SCHEMA = "d6-k6-opposite-light-ray-rank-v1"
CERTIFICATE_SCHEMA = "d6-k6-opposite-light-ray-rank-certificates-v1"
EXPECTED_INPUT = 623
EXPECTED_INPUT_SHA256 = (
    "b31aeac00b91d0d843ea51909c64f2c4ca3da49e5d2d33792312a45aa9f25e79"
)
EXPECTED_REJECTIONS = 171
EXPECTED_REJECTED_SHA256 = (
    "390bca6d3661ba4ab497d2ba23400d217d759702b071ecbf0438a9d03301bdf6"
)
EXPECTED_RESIDUE_SHA256 = (
    "cb62e002449e9be802ae37d66214904168de9528777d582cf6ddf28b350a5656"
)
DEFAULT_WORKERS = max(1, (os.cpu_count() or 2) - 1)
PACKAGE_SOURCES = (
    "d6_k6_opposite_ray_rank.py",
    "verify_d6_k6_opposite_ray_rank.py",
    "test_d6_k6_opposite_ray_rank.py",
    "d6_k6_opposite_ray_rank.md",
)

EXPECTED_DEPENDENCIES = {
    "build_d6_current_residue_manifest_v6.py": (
        "32345fc272efec6e04cab76a158411f65b4ea7cc6970ad4d1442dc8ca8471918"
    ),
    "verify_d6_current_residue_manifest_v6.py": (
        "664aa48f3bfed1acc983e2b4ccd1b5546e86ad8e5a57f89e73d20f20bd2c63fa"
    ),
    "d6_current_residue_manifest_v6.json": (
        "c513bcc40c037af23340dce3cc91d593a101dd1d073f93feab82abd556efadc9"
    ),
    "d6_current_residue_manifest_v6_verification.json": (
        "3d11cbd20cb59b42a8c8fd5252fb3ef3dea649576982ea8017622fad3279a7bc"
    ),
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
    "d6_k6_lorentz_reference.py": (
        "ef6a0877a5987119c92ec1bd9271a6a47774fc9e4992c10e1f28b6d5373592c2"
    ),
    "test_d6_k6_lorentz.py": (
        "f4fc8544a4311f742267808e6c785f82139accb9da99b496d206283b52831a47"
    ),
}


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


def rendered_json(value: object) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode("utf-8")


def atomic_gzip_json(path: Path, value: object) -> str:
    """Write deterministic gzip JSON and return the raw-payload SHA-256."""

    payload = rendered_json(value)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    with temporary.open("wb") as raw:
        with gzip.GzipFile(
            filename="", mode="wb", compresslevel=9, fileobj=raw, mtime=0
        ) as stream:
            stream.write(payload)
        raw.flush()
        os.fsync(raw.fileno())
    os.replace(temporary, path)
    return hashlib.sha256(payload).hexdigest()


def worker_initializer() -> None:
    for name in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[name] = "1"


def git_provenance() -> dict:
    """Bind a run to committed sources on the campaign branch.

    Untracked files are permitted because result artifacts and unrelated
    discovery probes may coexist in the shared worktree.  Any tracked or
    staged dirt blocks a theorem-level launch.
    """

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
        "dirty": bool(lines),
        "porcelain_sha256": hashlib.sha256(status.encode("utf-8")).hexdigest(),
        "porcelain_lines": lines,
        "committed_source_sha256": source_hashes,
    }


def eligible_z0_subsets(eligible: int) -> Iterator[int]:
    """Enumerate every eligible zero-factor set of size at most six."""

    eligible_vertices = vertices(eligible)
    for size in range(min(COORDINATES, len(eligible_vertices)) + 1):
        for chosen in combinations(eligible_vertices, size):
            yield sum(1 << local for local in chosen)


def absolute_vertices(instance, mask: int) -> list[int]:
    return [instance.outside[local] for local in vertices(mask)]


def common_neighbour_graph(
    adjacency: Sequence[int], instance, source: int, target_light: int
) -> tuple[int, ...]:
    """Graph on source: two vertices meet a required target-light vertex."""

    source_vertices = vertices(source)
    target_absolute = sum(
        1 << instance.outside[local] for local in vertices(target_light)
    )
    rows = [0] * len(source_vertices)
    for right in range(len(source_vertices)):
        absolute_right = instance.outside[source_vertices[right]]
        for left in range(right):
            absolute_left = instance.outside[source_vertices[left]]
            if adjacency[absolute_right] & adjacency[absolute_left] & target_absolute:
                rows[right] |= 1 << left
                rows[left] |= 1 << right
    return tuple(rows)


def psd_z_component_rank_lower(rows: Sequence[int]) -> tuple[int, list[list[int]]]:
    """Return sum |H|-1 over nontrivial connected support components."""

    unseen = (1 << len(rows)) - 1
    rank_lower = 0
    components: list[list[int]] = []
    while unseen:
        root = unseen & -unseen
        component = 0
        frontier = root
        while frontier:
            bit = frontier & -frontier
            frontier ^= bit
            local = bit.bit_length() - 1
            component |= bit
            frontier |= rows[local] & unseen & ~component
        unseen &= ~component
        if component.bit_count() >= 2:
            rank_lower += component.bit_count() - 1
            components.append(vertices(component))
    return rank_lower, components


def directed_ray_check(
    adjacency: Sequence[int], instance, z0: int, source: int, target: int
) -> dict:
    """Project the nonzero source-ray vectors off target union Z0."""

    source_light = source & ~z0
    target_light = target & ~z0
    source_vertices = vertices(source_light)
    rows = common_neighbour_graph(
        adjacency, instance, source_light, target_light
    )
    component_lower, local_components = psd_z_component_rank_lower(rows)
    target_absolute_mask = sum(
        1 << instance.outside[local] for local in vertices(target_light)
    )
    zero_cross_degree = [
        instance.outside[local]
        for local in source_vertices
        if not (adjacency[instance.outside[local]] & target_absolute_mask)
    ]
    lower = component_lower + len(zero_cross_degree)
    components = [
        [instance.outside[source_vertices[local]] for local in component]
        for component in local_components
    ]
    capacity = COORDINATES - target.bit_count()
    return {
        "passed": lower <= capacity,
        "source": absolute_vertices(instance, source_light),
        "orthonormal_target": absolute_vertices(instance, target),
        "target_light": absolute_vertices(instance, target_light),
        "common_neighbour_components": components,
        "component_rank_lower": component_lower,
        "zero_cross_degree_source": zero_cross_degree,
        "rank_lower": lower,
        "rank_capacity": capacity,
    }


def ray_assignment_decision(
    adjacency: Sequence[int], instance, z0: int
) -> tuple[dict, Counter[str]]:
    """Exhaust all nonbipartite-component assignments modulo ray exchange."""

    components = nonbipartite_components(instance, z0)
    counts: Counter[str] = Counter()
    counts["nonbipartite_components"] += len(components)
    expected = 1 if not components else 1 << (len(components) - 1)
    failures: list[dict] = []
    for encoded in range(expected):
        counts["ray_colorings_considered"] += 1
        if components:
            bins = [z0 | components[0], z0]
            colors = [0]
            for number, component in enumerate(components[1:]):
                color = (encoded >> number) & 1
                colors.append(color)
                bins[color] |= component
        else:
            bins = [z0, z0]
            colors = []
        matches = [
            support_matching(bin_mask, instance.defects) for bin_mask in bins
        ]
        row = {
            "colors": colors,
            "bins": [absolute_vertices(instance, bin_mask) for bin_mask in bins],
            "bin_matchable": [match is not None for match in matches],
        }
        if any(match is None for match in matches):
            counts["ray_coloring_bin_hall_fail"] += 1
            row["reason"] = "allowed_support_Hall"
            failures.append(row)
            continue
        checks = [
            directed_ray_check(adjacency, instance, z0, bins[0], bins[1]),
            directed_ray_check(adjacency, instance, z0, bins[1], bins[0]),
        ]
        counts["directed_ray_checks"] += 2
        counts["directed_ray_rank_fail"] += sum(
            not check["passed"] for check in checks
        )
        row["checks"] = checks
        if all(check["passed"] for check in checks):
            row["reason"] = "PASS"
            return ({
                "passed": True,
                "odd_components": [
                    absolute_vertices(instance, component)
                    for component in components
                ],
                "expected_colorings_mod_ray_swap": expected,
                "colorings_considered": int(counts["ray_colorings_considered"]),
                "witness": row,
            }, counts)
        row["reason"] = "opposite_ray_psd_Z_rank"
        failures.append(row)
    return ({
        "passed": False,
        "odd_components": [
            absolute_vertices(instance, component) for component in components
        ],
        "expected_colorings_mod_ray_swap": expected,
        "colorings_considered": int(counts["ray_colorings_considered"]),
        "coloring_failures": failures,
    }, counts)


def solve_seed(adjacency: Sequence[int], instance) -> dict:
    """Exhaust every eligible actual Z0; return the first passing witness."""

    counts: Counter[str] = Counter()
    failures: list[dict] = []
    for z0 in eligible_z0_subsets(instance.eligible_z0_mask):
        counts["z0_considered"] += 1
        z0_absolute = absolute_vertices(instance, z0)
        matching = support_matching(z0, instance.defects)
        if matching is None:
            counts["z0_unmatchable"] += 1
            failures.append({
                "Z0": z0_absolute,
                "Z0_matchable": False,
                "reason": "Z0_allowed_support_Hall",
            })
            continue
        ray, local = ray_assignment_decision(adjacency, instance, z0)
        counts.update(local)
        if ray["passed"]:
            counts["z0_feasible"] += 1
            return {
                "feasible": True,
                "witness": {"Z0": z0_absolute, "ray": ray},
                "counts": dict(counts),
                "failures": None,
            }
        counts["z0_ray_impossible"] += 1
        failures.append({
            "Z0": z0_absolute,
            "Z0_matchable": True,
            "reason": "all_ray_colorings_fail",
            "ray": ray,
        })
    return {
        "feasible": False,
        "witness": None,
        "counts": dict(counts),
        "failures": failures,
    }


def evaluate_record(record: dict) -> dict:
    adjacency = tuple(map(int, record["adjacency"]))
    validate_graph(adjacency, require_alpha_two=True)
    counts: Counter[str] = Counter()
    seeds_checked = 0
    for seed_mask in clique_masks(adjacency, COORDINATES):
        seeds_checked += 1
        instance = build_instance(adjacency, vertices(seed_mask))
        decision = solve_seed(adjacency, instance)
        counts.update(decision["counts"])
        if not decision["feasible"]:
            return {
                "index": int(record["index"]),
                "rejected": True,
                "seeds_checked": seeds_checked,
                "first_impossible_seed": list(instance.seed),
                "counts": dict(counts),
                "certificate": {
                    "seed": list(instance.seed),
                    "allowed_defects": {
                        str(instance.outside[local]): vertices(mask)
                        for local, mask in enumerate(instance.defects)
                    },
                    "eligible_Z0_vertices": absolute_vertices(
                        instance, instance.eligible_z0_mask
                    ),
                    "Z0_failures": decision["failures"],
                },
            }
    if not seeds_checked:
        raise ValueError(f"input graph {record['index']} has no required K6")
    return {
        "index": int(record["index"]),
        "rejected": False,
        "seeds_checked": seeds_checked,
        "first_impossible_seed": None,
        "counts": dict(counts),
        "certificate": None,
    }


def load_input() -> tuple[list[dict], list[int], dict[str, str]]:
    dependencies = {name: sha256(ROOT / name) for name in EXPECTED_DEPENDENCIES}
    if dependencies != EXPECTED_DEPENDENCIES:
        raise ValueError(f"opposite-ray dependency boundary changed: {dependencies}")
    manifest = json.loads(INPUT_MANIFEST.read_text(encoding="utf-8"))
    manifest_verification = json.loads(
        (ROOT / "d6_current_residue_manifest_v6_verification.json").read_text(
            encoding="utf-8"
        )
    )
    parent = json.loads(PARENT_REPORT.read_text(encoding="utf-8"))
    parent_verification = json.loads(PARENT_VERIFICATION.read_text(encoding="utf-8"))
    if (
        manifest.get("schema") != "d6-current-certified-residue-v6"
        or manifest.get("status") != "COMPLETE_MIXED_CERTIFICATE_UNION"
        or manifest.get("classes", {}).get("K6_only", {}).get("count") != 625
        or manifest_verification.get("schema")
        != "d6-current-certified-residue-v6-verification-v1"
        or manifest_verification.get("status") != "PASS"
        or parent.get("schema") != "d6-k6-tight-same-z0-v1"
        or parent.get("status") != "COMPLETE"
        or parent.get("input_graphs") != 625
        or parent.get("graphs_rejected") != 2
        or parent.get("graphs_surviving") != EXPECTED_INPUT
        or parent.get("ordered_residue_indices_sha256")
        != EXPECTED_INPUT_SHA256
        or parent_verification.get("schema")
        != "d6-k6-tight-same-z0-verification-v1"
        or parent_verification.get("status") != "PASS"
        or parent_verification.get("graphs_surviving") != EXPECTED_INPUT
        or parent_verification.get("ordered_residue_indices_sha256")
        != EXPECTED_INPUT_SHA256
    ):
        raise ValueError("certified 623-graph parent boundary changed")
    all_records = manifest["classes"]["K6_only"]["graphs"]
    by_index = {int(record["index"]): record for record in all_records}
    indices = [int(index) for index in parent["ordered_residue_indices"]]
    if (
        len(indices) != EXPECTED_INPUT
        or stable_hash(indices) != EXPECTED_INPUT_SHA256
        or len(by_index) != 625
        or any(index not in by_index for index in indices)
    ):
        raise ValueError("failed to reconstruct exact ordered K6-623 input")
    return [by_index[index] for index in indices], indices, dependencies


def positive_control() -> dict:
    adjacency = lower_bound_18_graph()
    record = {"index": -18, "adjacency": adjacency}
    result = evaluate_record(record)
    if result["rejected"] or result["seeds_checked"] != 32:
        raise AssertionError("opposite-ray theorem rejects standard 18-point set")
    return {
        "passed": True,
        "points": 18,
        "K6_seeds": result["seeds_checked"],
        "alpha_at_most_two": True,
    }


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
    counts: Counter[str] = Counter()
    for row in results:
        counts.update(row["counts"])
    archive = {
        "schema": CERTIFICATE_SCHEMA,
        "production_source_sha256": sha256(Path(__file__)),
        "ordered_input_indices_sha256": EXPECTED_INPUT_SHA256,
        "rejected_indices_sha256": EXPECTED_REJECTED_SHA256,
        "rejected_graphs": [
            {"index": row["index"], **row["certificate"]}
            for row in results if row["rejected"]
        ],
    }
    raw_certificate_sha256 = atomic_gzip_json(certificates, archive)
    report = {
        "schema": REPORT_SCHEMA,
        "status": "COMPLETE",
        "description": "exact opposite-light-ray PSD Z-matrix rank filter",
        "semantics": {
            "candidate_nonedges": "unconstrained_and_may_be_unit",
            "allowed_supports": "upper_bounds_coordinates_may_be_zero",
            "distinct_points": "required",
            "graph_precondition": "alpha_at_most_two",
            "arithmetic": "exact_integer_bitmasks_only",
        },
        "quantifiers": {
            "graph": "one impossible required K6 seed rejects",
            "seed": "every eligible actual Z0 is enumerated",
            "Z0": (
                "every nonbipartite-component assignment modulo global ray "
                "exchange is enumerated"
            ),
            "ray_assignment": "both directed projection-rank bounds must pass",
        },
        "rank_bound": (
            "zero_cross_degree_source_count + "
            "sum_nontrivial_common_neighbour_components(size-1) "
            "<= 6-orthonormal_target_size"
        ),
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
            **dict(counts),
        },
        "certificate_archive": {
            "path": certificates.name,
            "compressed_sha256": sha256(certificates),
            "raw_payload_sha256": raw_certificate_sha256,
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
            "the filter does not use the stronger identical-coloring actual-support conjunction",
            "this layer does not settle the K6-only or dimension-six problem",
        ],
    }
    atomic_json(output, report)
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    parser.add_argument(
        "--output", type=Path,
        default=ROOT / "d6_k6_opposite_ray_rank_report.json",
    )
    parser.add_argument(
        "--certificates", type=Path,
        default=ROOT / "d6_k6_opposite_ray_rank_certificates.json.gz",
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
