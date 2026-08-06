#!/usr/bin/env python3
"""Independent exact verifier for the full K6 normal-inertia campaign.

This checker never imports the production evaluator or runner.  It rebuilds
K6 seeds, defect masks, zero-factor choices, Lorentz components, support
matchings, and every graph decision independently.  Most importantly, it
computes the inertia of each integer ``I+Adj`` matrix with SymPy exact
characteristic polynomials and Sturm root counts, rather than the production
evaluator's rational symmetric elimination.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import platform
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Iterator, Sequence

import sympy as sp

from test_d6_k6_lorentz import lower_bound_18_graph


ROOT = Path(__file__).resolve().parent
EXPECTED_SELECTED = 1_097
EXPECTED_INPUT_SHA256 = (
    "dfee979e6ca2759e80cf803c3054c8954d98543875d467644b91ee434c22e845"
)
EXPECTED_PRIOR_REPORT_SHA256 = (
    "ede04a71a7f9b6901bbc1e4c198a117797bfd90d4d1366912ca9efd1ca559bdc"
)
EXPECTED_DEPENDENCY_SHA256 = {
    "d6_k6_lorentz_reference.py": (
        "ef6a0877a5987119c92ec1bd9271a6a47774fc9e4992c10e1f28b6d5373592c2"
    ),
    "d6_k6_support_reference.py": (
        "d6481137ef49d88154882285660dd755eb7cb652ab87280f05c3792de9742577"
    ),
    "d6_k6_bipartite_rank_reference.py": (
        "e760fef74c423a89f1e1f9d08909fa1223e0389900cf3d45b0140ee296dea922"
    ),
    "d6_k6_normal_inertia.py": (
        "2d00ab40eb97aa79ffec3b8134b83c4fd5f703cbd691f094a25bf3ff9f3ee023"
    ),
    "d6_k6_normal_inertia_test.py": (
        "4c107490befec2919c0a7eabe8f79f61eb8b11e56e40800cc630b5fdd3db7e12"
    ),
    "d6_k6_normal_coordinates.md": (
        "dd755349da45bbaca2f6c85dcdd5d56e4c08042f5d137cb349204e81b9686ca7"
    ),
    "d6_k6_normal_inertia_sample.json": (
        "3bc04d4541e11597302631fb30466af19d4c8b18333d4dffabca2f059b06128d"
    ),
    "d6_k6_normal_inertia_full_runner.py": (
        "a3ff3859686b503dd7c5066dcc77d07d15b57ea8eb169df28efc182ea5db5c70"
    ),
    "test_d6_k6_lorentz.py": (
        "f4fc8544a4311f742267808e6c785f82139accb9da99b496d206283b52831a47"
    ),
}
DECISION_SCHEMA = "d6-k6-normal-inertia-full-decision-v1"
REPORT_SCHEMA = "d6-k6-normal-inertia-full-report-v1"
VERIFICATION_SCHEMA = "d6-k6-normal-inertia-full-verification-v1"
COORDINATES = 6
VARIABLE = sp.symbols("lambda")


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
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def bits(mask: int) -> list[int]:
    answer = []
    while mask:
        bit = mask & -mask
        mask ^= bit
        answer.append(bit.bit_length() - 1)
    return answer


def validate_graph(adj: Sequence[int]) -> None:
    size = len(adj)
    full = (1 << size) - 1
    for vertex, row in enumerate(adj):
        if row & ~full or row & (1 << vertex):
            raise ValueError("invalid adjacency mask")
        for other in range(vertex):
            if bool(row & (1 << other)) != bool(adj[other] & (1 << vertex)):
                raise ValueError("asymmetric adjacency")
    for first in range(size):
        nonneighbors = full & ~adj[first] & ~((1 << (first + 1)) - 1)
        remaining = nonneighbors
        while remaining:
            bit = remaining & -remaining
            remaining ^= bit
            second = bit.bit_length() - 1
            if remaining & ~adj[second]:
                raise ValueError("independent triple violates verifier precondition")


def clique_masks(adj: Sequence[int], size: int) -> Iterator[int]:
    def visit(candidates: int, need: int, selected: int) -> Iterator[int]:
        if need == 0:
            yield selected
            return
        while candidates.bit_count() >= need:
            bit = candidates & -candidates
            candidates ^= bit
            vertex = bit.bit_length() - 1
            yield from visit(candidates & adj[vertex], need - 1, selected | bit)

    yield from visit((1 << len(adj)) - 1, size, 0)


@dataclass(frozen=True)
class Instance:
    seed: tuple[int, ...]
    outside: tuple[int, ...]
    defects: tuple[int, ...]
    l_adj: tuple[int, ...]
    eligible_z0: int


def build_instance(adj: Sequence[int], seed_mask: int) -> Instance:
    seed = tuple(bits(seed_mask))
    if len(seed) != COORDINATES:
        raise ValueError("verifier seed does not have six vertices")
    outside = tuple(vertex for vertex in range(len(adj)) if not seed_mask & (1 << vertex))
    defects = []
    for vertex in outside:
        mask = 0
        for coordinate, seed_vertex in enumerate(seed):
            if not adj[vertex] & (1 << seed_vertex):
                mask |= 1 << coordinate
        defects.append(mask)
    l_adj = [0] * len(outside)
    for left, right in combinations(range(len(outside)), 2):
        if (
            adj[outside[left]] & (1 << outside[right])
            and not defects[left] & defects[right]
        ):
            l_adj[left] |= 1 << right
            l_adj[right] |= 1 << left
    eligible = sum(
        1 << vertex
        for vertex, mask in enumerate(defects)
        if mask.bit_count() >= 3
    )
    return Instance(seed, outside, tuple(defects), tuple(l_adj), eligible)


def z0_subsets(eligible: int) -> Iterator[int]:
    candidates = bits(eligible)
    for size in range(min(COORDINATES, len(candidates)) + 1):
        for selected in combinations(candidates, size):
            yield sum(1 << vertex for vertex in selected)


def support_matchable(selected: int, defects: Sequence[int]) -> bool:
    vertices = bits(selected)
    if len(vertices) > COORDINATES:
        return False
    vertices.sort(key=lambda vertex: (defects[vertex].bit_count(), vertex))

    def visit(position: int, used: int) -> bool:
        if position == len(vertices):
            return True
        available = defects[vertices[position]] & ~used
        while available:
            bit = available & -available
            available ^= bit
            if visit(position + 1, used | bit):
                return True
        return False

    return visit(0, 0)


@dataclass(frozen=True)
class Component:
    vertices: int
    bipartite: bool
    side_a: int
    side_b: int


def components(instance: Instance, deleted: int) -> tuple[Component, ...]:
    remaining = ((1 << len(instance.outside)) - 1) & ~deleted
    answer = []
    while remaining:
        root = remaining & -remaining
        reached = root
        frontier = root
        side_b = 0
        bipartite = True
        while frontier:
            bit = frontier & -frontier
            frontier ^= bit
            vertex = bit.bit_length() - 1
            neighbors = instance.l_adj[vertex] & remaining
            if side_b & bit:
                if neighbors & side_b:
                    bipartite = False
            elif neighbors & (reached & ~side_b):
                bipartite = False
            new = neighbors & ~reached
            if not side_b & bit:
                side_b |= new
            reached |= new
            frontier |= new
        side_b &= reached
        answer.append(Component(reached, bipartite, reached & ~side_b, side_b))
        remaining &= ~reached
    return tuple(answer)


def side_pattern(adj: Sequence[int], instance: Instance, selected: int) -> tuple[int, ...]:
    absolute = tuple(instance.outside[index] for index in bits(selected))
    return tuple(
        (1 << row)
        | sum(
            1 << column
            for column, other in enumerate(absolute)
            if adj[vertex] & (1 << other)
        )
        for row, vertex in enumerate(absolute)
    )


class SympyInertiaCache:
    def __init__(self) -> None:
        self.values: dict[tuple[int, ...], tuple[int, int, int]] = {}
        self.hits = 0

    def solve(self, rows: tuple[int, ...]) -> tuple[int, int, int]:
        if rows in self.values:
            self.hits += 1
            return self.values[rows]
        size = len(rows)
        matrix = sp.Matrix(
            size,
            size,
            lambda row, column: (rows[row] >> column) & 1,
        )
        polynomial = matrix.charpoly(VARIABLE).as_poly()
        positive = 0
        negative = 0
        zero = 0
        _, factors = polynomial.sqf_list()
        for factor, multiplicity in factors:
            current = factor
            if current.eval(0) == 0:
                zero += multiplicity
                current = current.exquo(sp.Poly(VARIABLE, VARIABLE))
            if current.degree() <= 0:
                continue
            negative += multiplicity * int(current.count_roots(-sp.oo, 0))
            positive += multiplicity * int(current.count_roots(0, sp.oo))
        if positive + negative + zero != size:
            raise AssertionError(
                f"inertia roots do not sum to size for {rows}: "
                f"{positive},{negative},{zero}"
            )
        answer = positive, negative, zero
        self.values[rows] = answer
        return answer


@dataclass(frozen=True)
class IndependentDecision:
    rejected: bool
    seeds_checked: int
    impossible_seeds: int
    z0_subsets_considered: int
    z0_support_matchable: int
    bipartite_components_checked: int
    inertia_failures: int
    inertia_cache_entries: int
    inertia_cache_hits: int
    first_impossible_seed: tuple[int, ...] | None


def solve_seed(
    adj: Sequence[int], instance: Instance, cache: SympyInertiaCache
) -> tuple[bool, tuple[int, int, int, int]]:
    considered = 0
    matchable = 0
    component_checks = 0
    failures = 0
    for z0 in z0_subsets(instance.eligible_z0):
        considered += 1
        if not support_matchable(z0, instance.defects):
            continue
        matchable += 1
        passed = True
        local_components = components(instance, z0)
        for component in local_components:
            if not component.bipartite:
                continue
            component_checks += 1
            pattern_a = side_pattern(adj, instance, component.side_a)
            pattern_b = side_pattern(adj, instance, component.side_b)
            positive_a, negative_a, _ = cache.solve(pattern_a)
            positive_b, negative_b, _ = cache.solve(pattern_b)
            size_a = component.side_a.bit_count()
            size_b = component.side_b.bit_count()
            first_orientation = (
                z0.bit_count()
                + size_a
                - negative_a
                + size_b
                - positive_b
            )
            second_orientation = (
                z0.bit_count()
                + size_a
                - positive_a
                + size_b
                - negative_b
            )
            if min(first_orientation, second_orientation) > COORDINATES:
                passed = False
        if passed:
            return True, (considered, matchable, component_checks, failures)
        failures += 1
    return False, (considered, matchable, component_checks, failures)


def evaluate_graph(adj: Sequence[int]) -> IndependentDecision:
    validate_graph(adj)
    if next(clique_masks(adj, 7), 0):
        raise ValueError("verification target unexpectedly contains a K7")
    cache = SympyInertiaCache()
    checked = 0
    impossible = 0
    considered = 0
    matchable = 0
    component_checks = 0
    failures = 0
    first_seed = None
    for seed_mask in clique_masks(adj, COORDINATES):
        checked += 1
        instance = build_instance(adj, seed_mask)
        feasible, counts = solve_seed(adj, instance, cache)
        considered += counts[0]
        matchable += counts[1]
        component_checks += counts[2]
        failures += counts[3]
        if not feasible:
            impossible += 1
            if first_seed is None:
                first_seed = instance.seed
    return IndependentDecision(
        impossible > 0,
        checked,
        impossible,
        considered,
        matchable,
        component_checks,
        failures,
        len(cache.values),
        cache.hits,
        first_seed,
    )


SCALAR_FIELDS = (
    "rejected",
    "seeds_checked",
    "impossible_seeds",
    "z0_subsets_considered",
    "z0_support_matchable",
    "bipartite_components_checked",
    "inertia_failures",
    "inertia_cache_entries",
    "inertia_cache_hits",
)


def verify_one(payload: tuple[dict, dict]) -> dict:
    record, production = payload
    independent = evaluate_graph(record["adjacency"])
    for name in SCALAR_FIELDS:
        observed = getattr(independent, name)
        expected = production["decision"][name]
        if observed != expected:
            raise AssertionError(
                f"graph {record['index']} field {name}: independent={observed}, "
                f"production={expected}"
            )
    witness = production["decision"].get("first_witness")
    production_seed = tuple(witness["seed"]) if witness is not None else None
    if independent.first_impossible_seed != production_seed:
        raise AssertionError(
            f"graph {record['index']} first impossible seed mismatch: "
            f"{independent.first_impossible_seed} != {production_seed}"
        )
    if production["decision"].get("applicable") is not True:
        raise AssertionError(f"graph {record['index']} production applicability false")
    return {
        "index": record["index"],
        **{name: getattr(independent, name) for name in SCALAR_FIELDS},
        "first_impossible_seed": independent.first_impossible_seed,
    }


def verify_dependencies(input_path: Path, prior_path: Path) -> dict[str, str]:
    observed = {
        input_path.name: sha256(input_path),
        prior_path.name: sha256(prior_path),
        **{name: sha256(ROOT / name) for name in EXPECTED_DEPENDENCY_SHA256},
    }
    expected = {
        input_path.name: EXPECTED_INPUT_SHA256,
        prior_path.name: EXPECTED_PRIOR_REPORT_SHA256,
        **EXPECTED_DEPENDENCY_SHA256,
    }
    if observed != expected:
        raise ValueError(f"dependency hash mismatch: observed={observed}, expected={expected}")
    return observed


def select_residue(input_path: Path, prior_path: Path) -> list[dict]:
    payload = json.loads(input_path.read_text(encoding="utf-8"))
    prior = json.loads(prior_path.read_text(encoding="utf-8"))
    rejected = {
        item["index"]
        for item in prior["graph_results"]
        if item["decision"]["rejected"]
    }
    if rejected != {461_363}:
        raise ValueError("prior report rejection set changed")
    selected = [record for record in payload["graphs"] if record["index"] not in rejected]
    if len(selected) != EXPECTED_SELECTED:
        raise ValueError("independent selector did not recover 1,097 graphs")
    return selected


def read_archive(path: Path) -> tuple[list[dict], str]:
    digest = hashlib.sha256()
    records = []
    with gzip.open(path, "rb") as stream:
        for line_number, line in enumerate(stream, 1):
            digest.update(line)
            try:
                record = json.loads(line)
            except json.JSONDecodeError as error:
                raise ValueError(f"invalid archive JSONL line {line_number}") from error
            records.append(record)
    return records, digest.hexdigest()


def aggregate(rows: Sequence[dict]) -> dict:
    rejected = [row["index"] for row in rows if row["decision"]["rejected"]]
    answer = {
        "graphs_rejected": len(rejected),
        "graphs_surviving": len(rows) - len(rejected),
        "rejected_indices": rejected,
    }
    for name in SCALAR_FIELDS[1:]:
        answer[name] = sum(int(row["decision"][name]) for row in rows)
    return answer


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input", type=Path, default=ROOT / "d6_k6_bipartite_rank_input.json"
    )
    parser.add_argument(
        "--prior-report", type=Path, default=ROOT / "d6_k6_bipartite_rank_report.json"
    )
    parser.add_argument(
        "--report", type=Path, default=ROOT / "d6_k6_normal_inertia_full_report.json"
    )
    parser.add_argument(
        "--archive",
        type=Path,
        default=ROOT / "d6_k6_normal_inertia_full_decisions.jsonl.gz",
    )
    parser.add_argument(
        "--checkpoint",
        type=Path,
        default=ROOT / "d6_k6_normal_inertia_full_checkpoint.json",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "d6_k6_normal_inertia_full_verification_report.json",
    )
    parser.add_argument("--workers", type=int, default=11)
    parser.add_argument("--chunksize", type=int, default=1)
    args = parser.parse_args()
    if min(args.workers, args.chunksize) < 1:
        parser.error("workers and chunksize must be positive")

    input_path = args.input.resolve()
    prior_path = args.prior_report.resolve()
    report_path = args.report.resolve()
    archive_path = args.archive.resolve()
    checkpoint_path = args.checkpoint.resolve()
    output_path = args.output.resolve()
    dependencies = verify_dependencies(input_path, prior_path)
    report = json.loads(report_path.read_text(encoding="utf-8"))
    if report.get("schema") != REPORT_SCHEMA or report.get("status") != "COMPLETE":
        raise ValueError("production report is not a COMPLETE expected-schema report")
    if sha256(archive_path) != report["artifacts"]["archive_sha256"]:
        raise ValueError("archive hash does not match production report")
    if sha256(checkpoint_path) != report["artifacts"]["checkpoint_sha256"]:
        raise ValueError("checkpoint hash does not match production report")
    checkpoint = json.loads(checkpoint_path.read_text(encoding="utf-8"))
    if checkpoint.get("status") != "COMPLETE" or checkpoint.get("committed_rows") != 1_097:
        raise ValueError("production checkpoint is not COMPLETE at 1,097 rows")

    production_rows, decompressed_hash = read_archive(archive_path)
    if decompressed_hash != report["artifacts"]["decisions_sha256"]:
        raise ValueError("decompressed decision hash does not match report")
    selected = select_residue(input_path, prior_path)
    if len(production_rows) != len(selected):
        raise ValueError("decision archive has the wrong number of rows")
    for position, (row, record) in enumerate(zip(production_rows, selected)):
        if row.get("schema") != DECISION_SCHEMA or row.get("index") != record["index"]:
            raise ValueError(f"decision archive prefix mismatch at row {position}")
    selection_hash = stable_hash(
        [{"index": row["index"], "adjacency": row["adjacency"]} for row in selected]
    )
    if selection_hash != report["selection"]["sha256"]:
        raise ValueError("independent selection hash does not match report")
    production_aggregate = aggregate(production_rows)
    for name, value in production_aggregate.items():
        if report.get(name) != value:
            raise ValueError(f"report aggregate mismatch for {name}")

    started = time.time()
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        verified = list(
            executor.map(
                verify_one,
                zip(selected, production_rows),
                chunksize=args.chunksize,
            )
        )
    independent_rejected = [row["index"] for row in verified if row["rejected"]]
    if independent_rejected != production_aggregate["rejected_indices"]:
        raise AssertionError("independent full rejection set differs")

    positive = evaluate_graph(lower_bound_18_graph())
    if positive.rejected or positive.seeds_checked != 32:
        raise AssertionError("known realizable 18-point positive control failed")

    verification = {
        "schema": VERIFICATION_SCHEMA,
        "status": "PASS",
        "method": (
            "independent graph/CSP reconstruction and exact SymPy "
            "characteristic-polynomial Sturm inertia"
        ),
        "graphs_recomputed": len(verified),
        "graphs_rejected": len(independent_rejected),
        "graphs_surviving": len(verified) - len(independent_rejected),
        "rejected_indices": independent_rejected,
        "all_scalar_decision_fields_matched": list(SCALAR_FIELDS),
        "all_first_impossible_seeds_matched": True,
        "positive_control": {
            "name": "known realizable 18-point construction",
            "K6_seeds": positive.seeds_checked,
            "passed": not positive.rejected,
        },
        "inputs": {
            **dependencies,
            report_path.name: sha256(report_path),
            archive_path.name: sha256(archive_path),
            checkpoint_path.name: sha256(checkpoint_path),
            Path(__file__).name: sha256(Path(__file__).resolve()),
        },
        "selection_sha256": selection_hash,
        "decompressed_decisions_sha256": decompressed_hash,
        "runtime": {
            "workers": args.workers,
            "chunksize": args.chunksize,
            "wall_seconds": time.time() - started,
            "python": sys.version,
            "sympy": sp.__version__,
            "platform": platform.platform(),
            "command": " ".join(sys.argv),
        },
        "arithmetic": (
            "exact integer characteristic polynomials, square-free "
            "factorization, and exact Sturm root counts"
        ),
    }
    atomic_json(output_path, verification)
    print(
        json.dumps(
            {
                "status": "PASS",
                "graphs_recomputed": len(verified),
                "graphs_rejected": len(independent_rejected),
                "output": str(output_path),
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
