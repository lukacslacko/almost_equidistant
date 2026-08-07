#!/usr/bin/env python3
"""Independent verifier for the opposite-light-ray projection-rank layer.

The checker imports neither the producer nor any of its graph/Lorentz helper
modules.  It independently reconstructs the certified K6-623 input, K6
seeds, allowed masks, Lorentz graph, all eligible Z0 sets, nonbipartite
components, Hall decisions (by Hall's subset criterion), ray colorings, and
the common-neighbour component rank lower bound.  It then compares all 623
decisions and every exhaustive rejected-seed certificate byte-for-structure.
"""

from __future__ import annotations

import argparse
import ast
import gzip
import hashlib
import json
import os
import platform
import subprocess
import sys
import time
from collections import Counter, deque
from concurrent.futures import ProcessPoolExecutor
from itertools import combinations
from pathlib import Path
from typing import Iterator, Sequence


ROOT = Path(__file__).resolve().parent
REPORT_SCHEMA = "d6-k6-opposite-light-ray-rank-v1"
CERTIFICATE_SCHEMA = "d6-k6-opposite-light-ray-rank-certificates-v1"
VERIFICATION_SCHEMA = "d6-k6-opposite-light-ray-rank-verification-v1"
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
DIMENSION = 6
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


def raw_gzip_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with gzip.open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path, value: object) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


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
        "d6_k6_opposite_ray_rank",
        "probe_d6_k6_opposite_ray_rank",
        "d6_k6_lorentz_reference",
        "d6_k6_bipartite_rank_reference",
    }
    if imports & forbidden:
        raise AssertionError(
            f"independent checker imports production helpers: {imports & forbidden}"
        )


def validate_porcelain_binding(git: dict) -> dict:
    """Require a correctly hashed launch status with only untracked dirt."""

    lines = git.get("porcelain_lines")
    if not isinstance(lines, list) or not all(isinstance(line, str) for line in lines):
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
    """Check the recorded commit contains every exact package source blob."""

    source_hashes = report.get("source_sha256")
    expected_sources = {name: sha256(ROOT / name) for name in PACKAGE_SOURCES}
    if source_hashes != expected_sources:
        raise ValueError("working package source boundary differs")
    git = report.get("execution", {}).get("git", {})
    if (
        not git.get("available")
        or git.get("branch") != "codex/dimension6"
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
            raise ValueError(f"cannot read committed package source: {name}") from error
        if hashlib.sha256(blob).hexdigest() != expected:
            raise ValueError(f"committed package source differs: {name}")
    porcelain = validate_porcelain_binding(git)
    return {
        "commit": commit,
        "branch": git["branch"],
        "tracked_sources_clean_at_launch": True,
        **porcelain,
    }


def bits(mask: int) -> list[int]:
    answer: list[int] = []
    while mask:
        bit = mask & -mask
        answer.append(bit.bit_length() - 1)
        mask ^= bit
    return answer


def add_edge(adjacency: list[int], left: int, right: int) -> None:
    adjacency[left] |= 1 << right
    adjacency[right] |= 1 << left


def validate_graph(adjacency: Sequence[int]) -> None:
    n = len(adjacency)
    full = (1 << n) - 1
    for left, row in enumerate(adjacency):
        if row & ~full or row & (1 << left):
            raise ValueError("invalid simple graph row")
        for right in range(left):
            if bool(row & (1 << right)) != bool(
                adjacency[right] & (1 << left)
            ):
                raise ValueError("asymmetric graph")
    # Deliberately direct independent-triple check.
    for a, b, c in combinations(range(n), 3):
        if not (
            adjacency[a] & (1 << b)
            or adjacency[a] & (1 << c)
            or adjacency[b] & (1 << c)
        ):
            raise ValueError("opposite-ray theorem requires alpha(G)<=2")


def k6_seeds(adjacency: Sequence[int]) -> Iterator[tuple[int, ...]]:
    """Lexicographic clique DFS, transcribed independently."""

    n = len(adjacency)

    def visit(prefix: tuple[int, ...], candidates: tuple[int, ...]):
        need = DIMENSION - len(prefix)
        if not need:
            yield prefix
            return
        for offset, vertex in enumerate(candidates):
            if len(candidates) - offset < need:
                break
            later = tuple(
                other
                for other in candidates[offset + 1 :]
                if adjacency[vertex] & (1 << other)
            )
            yield from visit(prefix + (vertex,), later)

    yield from visit((), tuple(range(n)))


def build_seed_instance(adjacency: Sequence[int], seed: tuple[int, ...]) -> dict:
    seed_set = set(seed)
    outside = tuple(vertex for vertex in range(len(adjacency)) if vertex not in seed_set)
    defects = []
    for vertex in outside:
        mask = 0
        for coordinate, seed_vertex in enumerate(seed):
            if not (adjacency[vertex] & (1 << seed_vertex)):
                mask |= 1 << coordinate
        defects.append(mask)
    lorentz = [0] * len(outside)
    for right in range(len(outside)):
        for left in range(right):
            if (
                adjacency[outside[left]] & (1 << outside[right])
                and not (defects[left] & defects[right])
            ):
                lorentz[left] |= 1 << right
                lorentz[right] |= 1 << left
    eligible = sum(
        1 << local
        for local, defect in enumerate(defects)
        if defect.bit_count() >= 3
    )
    return {
        "seed": seed,
        "outside": outside,
        "defects": tuple(defects),
        "lorentz": tuple(lorentz),
        "eligible": eligible,
    }


def z0_subsets(eligible: int) -> Iterator[int]:
    eligible_vertices = bits(eligible)
    for size in range(min(DIMENSION, len(eligible_vertices)) + 1):
        for chosen in combinations(eligible_vertices, size):
            yield sum(1 << local for local in chosen)


def hall_matchable(selected: int, defects: Sequence[int]) -> bool:
    """Independent Hall-subset test, rather than the producer's matching DFS."""

    selected_vertices = bits(selected)
    if len(selected_vertices) > DIMENSION:
        return False
    for size in range(1, len(selected_vertices) + 1):
        for subset in combinations(selected_vertices, size):
            union = 0
            for local in subset:
                union |= defects[local]
            if union.bit_count() < size:
                return False
    return True


def odd_components(instance: dict, z0: int) -> list[int]:
    outside_count = len(instance["outside"])
    remaining = set(bits(((1 << outside_count) - 1) & ~z0))
    answer: list[int] = []
    while remaining:
        root = min(remaining)
        queue = deque([root])
        color = {root: 0}
        component: set[int] = set()
        nonbipartite = False
        while queue:
            vertex = queue.popleft()
            if vertex in component:
                continue
            component.add(vertex)
            for neighbour in bits(instance["lorentz"][vertex]):
                if neighbour not in remaining:
                    continue
                wanted = 1 - color[vertex]
                if neighbour in color:
                    if color[neighbour] != wanted:
                        nonbipartite = True
                else:
                    color[neighbour] = wanted
                    queue.append(neighbour)
        remaining -= component
        if nonbipartite:
            answer.append(sum(1 << local for local in component))
    return answer


def absolute_vertices(instance: dict, mask: int) -> list[int]:
    return [instance["outside"][local] for local in bits(mask)]


def common_graph(
    adjacency: Sequence[int], instance: dict, source: int, target_light: int
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    source_vertices = tuple(bits(source))
    target_vertices = tuple(bits(target_light))
    rows = [0] * len(source_vertices)
    for left, right in combinations(range(len(source_vertices)), 2):
        source_left = instance["outside"][source_vertices[left]]
        source_right = instance["outside"][source_vertices[right]]
        if any(
            adjacency[source_left] & (1 << instance["outside"][target])
            and adjacency[source_right] & (1 << instance["outside"][target])
            for target in target_vertices
        ):
            rows[left] |= 1 << right
            rows[right] |= 1 << left
    return source_vertices, tuple(rows)


def connected_nontrivial(rows: Sequence[int]) -> list[list[int]]:
    unseen = set(range(len(rows)))
    answer: list[list[int]] = []
    while unseen:
        root = min(unseen)
        queue = [root]
        component: set[int] = set()
        while queue:
            vertex = queue.pop()
            if vertex in component:
                continue
            component.add(vertex)
            queue.extend(neighbour for neighbour in bits(rows[vertex]) if neighbour in unseen)
        unseen -= component
        if len(component) >= 2:
            answer.append(sorted(component))
    return answer


def directed_check(
    adjacency: Sequence[int], instance: dict, z0: int, source: int, target: int
) -> dict:
    source_light = source & ~z0
    target_light = target & ~z0
    source_vertices, rows = common_graph(
        adjacency, instance, source_light, target_light
    )
    local_components = connected_nontrivial(rows)
    component_lower = sum(len(component) - 1 for component in local_components)
    target_absolute = tuple(
        instance["outside"][local] for local in bits(target_light)
    )
    zero_cross = [
        instance["outside"][local]
        for local in source_vertices
        if not any(
            adjacency[instance["outside"][local]] & (1 << target_vertex)
            for target_vertex in target_absolute
        )
    ]
    lower = component_lower + len(zero_cross)
    capacity = DIMENSION - target.bit_count()
    return {
        "passed": lower <= capacity,
        "source": absolute_vertices(instance, source_light),
        "orthonormal_target": absolute_vertices(instance, target),
        "target_light": absolute_vertices(instance, target_light),
        "common_neighbour_components": [
            [instance["outside"][source_vertices[local]] for local in component]
            for component in local_components
        ],
        "component_rank_lower": component_lower,
        "zero_cross_degree_source": zero_cross,
        "rank_lower": lower,
        "rank_capacity": capacity,
    }


def ray_decision(
    adjacency: Sequence[int], instance: dict, z0: int
) -> tuple[dict, Counter[str]]:
    components = odd_components(instance, z0)
    counts: Counter[str] = Counter()
    counts["nonbipartite_components"] += len(components)
    expected = 1 if not components else 1 << (len(components) - 1)
    failures = []
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
        matchable = [
            hall_matchable(bin_mask, instance["defects"]) for bin_mask in bins
        ]
        row = {
            "colors": colors,
            "bins": [absolute_vertices(instance, bin_mask) for bin_mask in bins],
            "bin_matchable": matchable,
        }
        if not all(matchable):
            counts["ray_coloring_bin_hall_fail"] += 1
            row["reason"] = "allowed_support_Hall"
            failures.append(row)
            continue
        checks = [
            directed_check(adjacency, instance, z0, bins[0], bins[1]),
            directed_check(adjacency, instance, z0, bins[1], bins[0]),
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


def solve_seed(adjacency: Sequence[int], instance: dict) -> dict:
    counts: Counter[str] = Counter()
    failures = []
    for z0 in z0_subsets(instance["eligible"]):
        counts["z0_considered"] += 1
        z0_absolute = absolute_vertices(instance, z0)
        if not hall_matchable(z0, instance["defects"]):
            counts["z0_unmatchable"] += 1
            failures.append({
                "Z0": z0_absolute,
                "Z0_matchable": False,
                "reason": "Z0_allowed_support_Hall",
            })
            continue
        ray, local = ray_decision(adjacency, instance, z0)
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
    validate_graph(adjacency)
    counts: Counter[str] = Counter()
    seeds_checked = 0
    for seed in k6_seeds(adjacency):
        seeds_checked += 1
        instance = build_seed_instance(adjacency, seed)
        decision = solve_seed(adjacency, instance)
        counts.update(decision["counts"])
        if not decision["feasible"]:
            return {
                "index": int(record["index"]),
                "rejected": True,
                "seeds_checked": seeds_checked,
                "first_impossible_seed": list(seed),
                "counts": dict(counts),
                "certificate": {
                    "seed": list(seed),
                    "allowed_defects": {
                        str(instance["outside"][local]): bits(mask)
                        for local, mask in enumerate(instance["defects"])
                    },
                    "eligible_Z0_vertices": absolute_vertices(
                        instance, instance["eligible"]
                    ),
                    "Z0_failures": decision["failures"],
                },
            }
    if not seeds_checked:
        raise ValueError("certified input unexpectedly has no K6")
    return {
        "index": int(record["index"]),
        "rejected": False,
        "seeds_checked": seeds_checked,
        "first_impossible_seed": None,
        "counts": dict(counts),
        "certificate": None,
    }


def load_input() -> tuple[list[dict], list[int]]:
    observed = {name: sha256(ROOT / name) for name in EXPECTED_DEPENDENCIES}
    if observed != EXPECTED_DEPENDENCIES:
        raise ValueError(f"frozen 623 boundary changed: {observed}")
    manifest = json.loads(
        (ROOT / "d6_current_residue_manifest_v6.json").read_text(encoding="utf-8")
    )
    manifest_verification = json.loads(
        (ROOT / "d6_current_residue_manifest_v6_verification.json").read_text(
            encoding="utf-8"
        )
    )
    parent = json.loads(
        (ROOT / "d6_k6_tight_same_z0_report.json").read_text(encoding="utf-8")
    )
    parent_verification = json.loads(
        (ROOT / "d6_k6_tight_same_z0_verification.json").read_text(
            encoding="utf-8"
        )
    )
    if (
        manifest.get("schema") != "d6-current-certified-residue-v6"
        or manifest.get("status") != "COMPLETE_MIXED_CERTIFICATE_UNION"
        or manifest.get("classes", {}).get("K6_only", {}).get("count") != 625
        or manifest_verification.get("status") != "PASS"
        or parent.get("schema") != "d6-k6-tight-same-z0-v1"
        or parent.get("status") != "COMPLETE"
        or parent.get("graphs_surviving") != EXPECTED_INPUT
        or parent.get("ordered_residue_indices_sha256") != EXPECTED_INPUT_SHA256
        or parent_verification.get("schema")
        != "d6-k6-tight-same-z0-verification-v1"
        or parent_verification.get("status") != "PASS"
        or parent_verification.get("graphs_surviving") != EXPECTED_INPUT
        or parent_verification.get("ordered_residue_indices_sha256")
        != EXPECTED_INPUT_SHA256
    ):
        raise ValueError("parent contents fail independent boundary checks")
    all_records = manifest["classes"]["K6_only"]["graphs"]
    by_index = {int(record["index"]): record for record in all_records}
    indices = [int(index) for index in parent["ordered_residue_indices"]]
    if (
        len(indices) != EXPECTED_INPUT
        or stable_hash(indices) != EXPECTED_INPUT_SHA256
        or len(by_index) != 625
        or any(index not in by_index for index in indices)
    ):
        raise ValueError("independent K6-623 reconstruction failed")
    return [by_index[index] for index in indices], indices


def standard_18_graph() -> list[int]:
    """Independent exact construction of the half-cube plus two poles."""

    words = [word for word in range(32) if word.bit_count() % 2]
    adjacency = [0] * 18
    for right, word in enumerate(words):
        for left, other in enumerate(words[:right]):
            if (word ^ other).bit_count() == 2:
                add_edge(adjacency, left, right)
    for pole in (16, 17):
        for base in range(16):
            add_edge(adjacency, base, pole)
    validate_graph(adjacency)
    return adjacency


def positive_control() -> dict:
    result = evaluate_record({"index": -18, "adjacency": standard_18_graph()})
    if result["rejected"] or result["seeds_checked"] != 32:
        raise AssertionError("independent theorem rejects standard 18-point set")
    return {
        "passed": True,
        "points": 18,
        "K6_seeds": result["seeds_checked"],
        "alpha_at_most_two": True,
    }


def audit_exhaustive_certificates(archive: dict) -> dict:
    z0_rows = 0
    coloring_rows = 0
    hall_failures = 0
    rank_failures = 0
    for graph in archive["rejected_graphs"]:
        eligible_count = len(graph["eligible_Z0_vertices"])
        expected_z0 = sum(
            1 for size in range(min(DIMENSION, eligible_count) + 1)
            for _ in combinations(range(eligible_count), size)
        )
        rows = graph["Z0_failures"]
        if len(rows) != expected_z0:
            raise AssertionError("certificate omits an eligible Z0")
        z0_rows += len(rows)
        for row in rows:
            if not row["Z0_matchable"]:
                hall_failures += 1
                continue
            ray = row["ray"]
            failures = ray["coloring_failures"]
            if (
                ray["passed"]
                or len(failures) != ray["expected_colorings_mod_ray_swap"]
                or len(failures) != ray["colorings_considered"]
            ):
                raise AssertionError("certificate omits a ray coloring")
            coloring_rows += len(failures)
            for failure in failures:
                if failure["reason"] == "allowed_support_Hall":
                    hall_failures += 1
                elif failure["reason"] == "opposite_ray_psd_Z_rank":
                    if all(check["passed"] for check in failure["checks"]):
                        raise AssertionError("rank failure has no failed direction")
                    rank_failures += 1
                else:
                    raise AssertionError("unknown rejected-coloring reason")
    return {
        "rejected_graphs": len(archive["rejected_graphs"]),
        "Z0_rows": z0_rows,
        "ray_coloring_rows": coloring_rows,
        "Hall_failure_rows": hall_failures,
        "rank_failure_rows": rank_failures,
    }


def verify(
    report_path: Path, certificates_path: Path, output: Path, workers: int
) -> dict:
    assert_import_independence()
    records, indices = load_input()
    report = json.loads(report_path.read_text(encoding="utf-8"))
    with gzip.open(certificates_path, "rt", encoding="utf-8") as stream:
        archive = json.load(stream)
    source_path = ROOT / "d6_k6_opposite_ray_rank.py"
    source_sha = sha256(source_path)
    if (
        report.get("schema") != REPORT_SCHEMA
        or report.get("status") != "COMPLETE"
        or report.get("production_source_sha256") != source_sha
        or report.get("source_sha256", {}).get(source_path.name) != source_sha
        or report.get("dependencies") != EXPECTED_DEPENDENCIES
        or report.get("ordered_input_indices") != indices
        or report.get("ordered_input_indices_sha256") != EXPECTED_INPUT_SHA256
        or report.get("graphs_rejected") != EXPECTED_REJECTIONS
        or report.get("rejected_indices_sha256") != EXPECTED_REJECTED_SHA256
        or report.get("ordered_residue_indices_sha256") != EXPECTED_RESIDUE_SHA256
        or report.get("certificate_archive", {}).get("compressed_sha256")
        != sha256(certificates_path)
        or report.get("certificate_archive", {}).get("raw_payload_sha256")
        != raw_gzip_sha256(certificates_path)
        or archive.get("schema") != CERTIFICATE_SCHEMA
        or archive.get("production_source_sha256") != source_sha
        or archive.get("ordered_input_indices_sha256") != EXPECTED_INPUT_SHA256
        or archive.get("rejected_indices_sha256") != EXPECTED_REJECTED_SHA256
    ):
        raise ValueError("production report/certificate boundary mismatch")
    provenance = validate_git_provenance(report)
    exhaustive = audit_exhaustive_certificates(archive)
    started = time.monotonic()
    if workers == 1:
        results = [evaluate_record(record) for record in records]
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            results = list(pool.map(evaluate_record, records, chunksize=1))
    rejected = [row["index"] for row in results if row["rejected"]]
    residue = [row["index"] for row in results if not row["rejected"]]
    graph_results = [
        {key: value for key, value in row.items() if key != "certificate"}
        for row in results
    ]
    recomputed_archive_graphs = [
        {"index": row["index"], **row["certificate"]}
        for row in results if row["rejected"]
    ]
    if (
        len(rejected) != EXPECTED_REJECTIONS
        or stable_hash(rejected) != EXPECTED_REJECTED_SHA256
        or stable_hash(residue) != EXPECTED_RESIDUE_SHA256
        or rejected != report.get("rejected_indices")
        or residue != report.get("ordered_residue_indices")
        or graph_results != report.get("graph_results")
        or recomputed_archive_graphs != archive.get("rejected_graphs")
    ):
        raise AssertionError("independent full-corpus replay disagrees")
    control = positive_control()
    if report.get("positive_18_control") != control:
        raise AssertionError("producer and verifier positive controls differ")
    verification = {
        "schema": VERIFICATION_SCHEMA,
        "status": "PASS",
        "checks": {
            "producer_imported": False,
            "input_reconstructed_independently": True,
            "Hall_algorithm": "all_subsets_union_cardinality",
            "all_graph_decisions_recomputed": True,
            "all_rejected_certificates_recomputed": True,
            "all_Z0_and_ray_colorings_archived": True,
            "integer_bitmask_arithmetic_only": True,
        },
        "report": {
            "path": report_path.name,
            "sha256": sha256(report_path),
        },
        "certificates": {
            "path": certificates_path.name,
            "compressed_sha256": sha256(certificates_path),
            "raw_payload_sha256": raw_gzip_sha256(certificates_path),
        },
        "sources": {
            "producer_sha256": source_sha,
            "verifier_sha256": sha256(Path(__file__)),
        },
        "provenance": provenance,
        "graphs_recomputed": len(results),
        "graphs_rejected": len(rejected),
        "graphs_surviving": len(residue),
        "rejected_indices": rejected,
        "rejected_indices_sha256": stable_hash(rejected),
        "ordered_residue_indices_sha256": stable_hash(residue),
        "exhaustive_certificate_audit": exhaustive,
        "positive_18_control": control,
        "runtime": {
            "workers": workers,
            "wall_seconds": time.monotonic() - started,
            "python": sys.version,
            "platform": platform.platform(),
        },
    }
    atomic_json(output, verification)
    return verification


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--report", type=Path,
        default=ROOT / "d6_k6_opposite_ray_rank_report.json",
    )
    parser.add_argument(
        "--certificates", type=Path,
        default=ROOT / "d6_k6_opposite_ray_rank_certificates.json.gz",
    )
    parser.add_argument(
        "--output", type=Path,
        default=ROOT / "d6_k6_opposite_ray_rank_verification.json",
    )
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    args = parser.parse_args()
    result = verify(args.report, args.certificates, args.output, args.workers)
    print(json.dumps({
        "status": result["status"],
        "graphs": result["graphs_recomputed"],
        "rejected": result["graphs_rejected"],
        "surviving": result["graphs_surviving"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
