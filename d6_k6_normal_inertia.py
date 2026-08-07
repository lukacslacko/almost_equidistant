#!/usr/bin/env python3
"""Exact normal-coordinate inertia screen for the dimension-six K6 residue.

Fix a required unit K6 and remove a possible zero Lorentz-factor set Z0.
For a connected bipartite component C=A union B of the disjoint-defect unit
edge graph, the nonzero Lorentz factors on either side are scalar multiples
of two orthogonal projective directions.  In the generic case their Lorentz
norms have opposite signs.  If ``H_A = I + Adj(G[A])``, the Gram matrix of
the actual defect vectors on A is

    I + kappa * D * H_A * D,

where D is nonsingular diagonal and kappa is the Lorentz norm of the chosen
direction divided by six.  Its nullity is at most the negative inertia of
H_A when kappa is positive, and at most the positive inertia when kappa is
negative.  The signs swap on B.  The A and B spans, and the |Z0|-dimensional
zero-factor span, are mutually orthogonal.  This gives the exact necessary
condition implemented below.

The lightlike case makes every vector of C orthonormal and therefore also
implies the (weaker) displayed inertia condition.  All inertia calculations
use rational symmetric elimination.  Floating point participates only in
wall-clock reporting.

This module is an isolated prototype.  It deliberately does not modify the
earlier production K6 references or their artifacts.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import sys
import time
from dataclasses import asdict, dataclass
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from typing import Iterator, Sequence

from d6_k6_bipartite_rank_reference import lorentz_components
from d6_k6_lorentz_reference import (
    COORDINATES,
    K6LorentzInstance,
    build_instance,
    clique_masks,
    find_clique_mask,
    support_matching,
    validate_graph,
    vertices,
)


ROOT = Path(__file__).resolve().parent


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def exact_inertia_rows(rows: Sequence[Sequence[int]]) -> tuple[int, int, int]:
    """Return (positive, negative, zero) inertia by exact congruence.

    A nonzero diagonal gives a one-by-one symmetric pivot.  If every
    remaining diagonal vanishes but an off-diagonal entry does not, the
    corresponding ``[[0,a],[a,0]]`` pivot has inertia (1,1).  Schur
    complements are over ``Fraction``, so Sylvester's law of inertia makes
    every returned count exact.
    """

    matrix = [[Fraction(value) for value in row] for row in rows]
    size = len(matrix)
    if any(len(row) != size for row in matrix):
        raise ValueError("inertia input must be square")
    if any(matrix[i][j] != matrix[j][i] for i in range(size) for j in range(i)):
        raise ValueError("inertia input must be symmetric")

    positive = 0
    negative = 0
    zero = 0
    while matrix:
        size = len(matrix)
        diagonal = next(
            (index for index in range(size) if matrix[index][index] != 0),
            None,
        )
        if diagonal is not None:
            order = [diagonal] + [
                index for index in range(size) if index != diagonal
            ]
            matrix = [
                [matrix[row][column] for column in order] for row in order
            ]
            pivot = matrix[0][0]
            if pivot > 0:
                positive += 1
            else:
                negative += 1
            matrix = [
                [
                    matrix[row][column]
                    - matrix[row][0] * matrix[0][column] / pivot
                    for column in range(1, size)
                ]
                for row in range(1, size)
            ]
            continue

        off_diagonal = next(
            (
                (row, column)
                for row in range(size)
                for column in range(row)
                if matrix[row][column] != 0
            ),
            None,
        )
        if off_diagonal is None:
            zero += size
            break

        first, second = off_diagonal
        order = [first, second] + [
            index for index in range(size) if index not in (first, second)
        ]
        matrix = [[matrix[row][column] for column in order] for row in order]
        pivot = matrix[0][1]
        positive += 1
        negative += 1
        matrix = [
            [
                matrix[row][column]
                - (
                    matrix[row][0] * matrix[1][column]
                    + matrix[row][1] * matrix[0][column]
                )
                / pivot
                for column in range(2, size)
            ]
            for row in range(2, size)
        ]

    if positive + negative + zero != len(rows):
        raise AssertionError("inertia counts do not sum to the matrix size")
    return positive, negative, zero


def _z0_subsets(eligible: int) -> Iterator[int]:
    choices = vertices(eligible)
    for size in range(min(COORDINATES, len(choices)) + 1):
        for selected in combinations(choices, size):
            yield sum(1 << vertex for vertex in selected)


class InertiaCache:
    """Cache exact inertia by a labeled ``I+Adj`` bit-row tuple."""

    def __init__(self) -> None:
        self.values: dict[tuple[int, ...], tuple[int, int, int]] = {}
        self.hits = 0

    def solve(self, rows: tuple[int, ...]) -> tuple[int, int, int]:
        cached = self.values.get(rows)
        if cached is not None:
            self.hits += 1
            return cached
        size = len(rows)
        matrix = tuple(
            tuple((rows[row] >> column) & 1 for column in range(size))
            for row in range(size)
        )
        answer = exact_inertia_rows(matrix)
        self.values[rows] = answer
        return answer


def side_unit_pattern(
    adj: Sequence[int], instance: K6LorentzInstance, selected: int
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """Return absolute vertices and bit rows of ``I+Adj(G[selected])``."""

    absolute = tuple(instance.outside[index] for index in vertices(selected))
    rows = tuple(
        (1 << row)
        | sum(
            1 << column
            for column, other in enumerate(absolute)
            if adj[vertex] & (1 << other)
        )
        for row, vertex in enumerate(absolute)
    )
    return absolute, rows


@dataclass(frozen=True)
class NormalInertiaCheck:
    component: tuple[int, ...]
    side_a: tuple[int, ...]
    side_b: tuple[int, ...]
    inertia_a: tuple[int, int, int]
    inertia_b: tuple[int, int, int]
    dimension_if_a_positive: int
    dimension_if_a_negative: int
    necessary_dimension: int
    passed: bool

    def jsonable(self) -> dict:
        return asdict(self)


def normal_dimension_lower_bound(
    size_a: int,
    inertia_a: tuple[int, int, int],
    size_b: int,
    inertia_b: tuple[int, int, int],
    zero_dimension: int = 0,
) -> tuple[int, int, int]:
    """Return the two Lorentz orientations and their necessary minimum."""

    if sum(inertia_a) != size_a or sum(inertia_b) != size_b:
        raise ValueError("inertia must sum to its side size")
    if zero_dimension < 0:
        raise ValueError("zero-factor dimension cannot be negative")
    positive_a, negative_a, _ = inertia_a
    positive_b, negative_b, _ = inertia_b
    if_a_positive = (
        zero_dimension
        + size_a
        - negative_a
        + size_b
        - positive_b
    )
    if_a_negative = (
        zero_dimension
        + size_a
        - positive_a
        + size_b
        - negative_b
    )
    return if_a_positive, if_a_negative, min(if_a_positive, if_a_negative)


def check_normal_inertia(
    adj: Sequence[int],
    instance: K6LorentzInstance,
    z0: int,
    cache: InertiaCache | None = None,
) -> tuple[bool, tuple[NormalInertiaCheck, ...]]:
    """Check every bipartite Lorentz component for one possible ``Z0``."""

    if cache is None:
        cache = InertiaCache()
    checks: list[NormalInertiaCheck] = []
    passed = True
    zero_dimension = z0.bit_count()
    for component in lorentz_components(instance, z0):
        if not component.bipartite:
            continue
        absolute_a, pattern_a = side_unit_pattern(
            adj, instance, component.side_a
        )
        absolute_b, pattern_b = side_unit_pattern(
            adj, instance, component.side_b
        )
        inertia_a = cache.solve(pattern_a)
        inertia_b = cache.solve(pattern_b)
        positive_a, negative_a, zero_a = inertia_a
        positive_b, negative_b, zero_b = inertia_b
        size_a = len(absolute_a)
        size_b = len(absolute_b)

        # If kappa_A>0, nullity on A is at most n_-(H_A), while
        # kappa_B<0 makes nullity on B at most n_+(H_B).  The other
        # orientation swaps the two signs.
        if_a_positive, if_a_negative, necessary = normal_dimension_lower_bound(
            size_a,
            inertia_a,
            size_b,
            inertia_b,
            zero_dimension,
        )
        item = NormalInertiaCheck(
            tuple(
                instance.outside[index]
                for index in vertices(component.component)
            ),
            absolute_a,
            absolute_b,
            (positive_a, negative_a, zero_a),
            (positive_b, negative_b, zero_b),
            if_a_positive,
            if_a_negative,
            necessary,
            necessary <= COORDINATES,
        )
        checks.append(item)
        if not item.passed:
            passed = False
    return passed, tuple(checks)


@dataclass(frozen=True)
class NormalSeedDecision:
    feasible: bool
    z0_subsets_considered: int
    z0_support_matchable: int
    z0_inertia_passed: int
    bipartite_components_checked: int
    inertia_failures: int
    chosen_z0: tuple[int, ...] | None
    witness: dict | None

    def jsonable(self) -> dict:
        return asdict(self)


def solve_normal_seed(
    adj: Sequence[int],
    instance: K6LorentzInstance,
    cache: InertiaCache | None = None,
) -> NormalSeedDecision:
    """Exhaust every eligible zero-factor set for one K6 seed."""

    if cache is None:
        cache = InertiaCache()
    considered = 0
    matchable = 0
    inertia_passed = 0
    components_checked = 0
    failures = 0
    failed_examples: list[dict] = []

    for z0 in _z0_subsets(instance.eligible_z0_mask):
        considered += 1
        if support_matching(z0, instance.defects) is None:
            continue
        matchable += 1
        passed, checks = check_normal_inertia(adj, instance, z0, cache)
        components_checked += len(checks)
        if passed:
            inertia_passed += 1
            return NormalSeedDecision(
                True,
                considered,
                matchable,
                inertia_passed,
                components_checked,
                failures,
                tuple(instance.outside[index] for index in vertices(z0)),
                None,
            )
        failures += 1
        if len(failed_examples) < 16:
            failed_examples.append(
                {
                    "Z0": [
                        instance.outside[index] for index in vertices(z0)
                    ],
                    "failed_components": [
                        check.jsonable() for check in checks if not check.passed
                    ],
                }
            )

    witness = instance.jsonable()
    witness.update(
        {
            "failure_kind": "no_Z0_passing_normal_inertia",
            "counts": {
                "Z0_subsets_considered": considered,
                "Z0_support_matchable": matchable,
                "Z0_inertia_passed": inertia_passed,
                "bipartite_components_checked": components_checked,
                "inertia_failures": failures,
            },
            "failed_examples": failed_examples,
        }
    )
    return NormalSeedDecision(
        False,
        considered,
        matchable,
        inertia_passed,
        components_checked,
        failures,
        None,
        witness,
    )


@dataclass(frozen=True)
class NormalGraphDecision:
    applicable: bool
    rejected: bool
    seeds_checked: int
    impossible_seeds: int
    z0_subsets_considered: int
    z0_support_matchable: int
    bipartite_components_checked: int
    inertia_failures: int
    inertia_cache_entries: int
    inertia_cache_hits: int
    first_witness: dict | None
    note: str | None = None

    def jsonable(self) -> dict:
        return asdict(self)


def evaluate_normal_graph(
    adj: Sequence[int], scan_all: bool = True
) -> NormalGraphDecision:
    """Reject iff one K6 seed exhausts every admissible ``Z0``."""

    validate_graph(adj, require_alpha_two=True)
    if find_clique_mask(adj, 7):
        return NormalGraphDecision(
            False, False, 0, 0, 0, 0, 0, 0, 0, 0, None,
            note="contains a required K7; prototype scope is K6-only",
        )
    cache = InertiaCache()
    checked = 0
    impossible = 0
    considered = 0
    matchable = 0
    components = 0
    failures = 0
    first_witness = None
    for seed_mask in clique_masks(adj, COORDINATES):
        checked += 1
        instance = build_instance(adj, vertices(seed_mask))
        decision = solve_normal_seed(adj, instance, cache)
        considered += decision.z0_subsets_considered
        matchable += decision.z0_support_matchable
        components += decision.bipartite_components_checked
        failures += decision.inertia_failures
        if not decision.feasible:
            impossible += 1
            if first_witness is None:
                first_witness = decision.witness
            if not scan_all:
                break
    return NormalGraphDecision(
        checked > 0,
        impossible > 0,
        checked,
        impossible,
        considered,
        matchable,
        components,
        failures,
        len(cache.values),
        cache.hits,
        first_witness,
        note=None if checked else "graph has no required K6 seed",
    )


def _atomic_json(path: Path, value: dict) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def _sample_records(records: list[dict], size: int) -> list[dict]:
    if size < 1 or size > len(records):
        raise ValueError("sample size must lie between one and residue size")
    return sorted(
        records,
        key=lambda record: hashlib.sha256(
            str(record["index"]).encode("ascii")
        ).digest(),
    )[:size]


def run_sample(
    input_path: Path,
    prior_report_path: Path,
    output_path: Path,
    sample_size: int,
) -> dict:
    """Profile a deterministic hash sample of the frozen K6 residue."""

    with input_path.open(encoding="utf-8") as stream:
        payload = json.load(stream)
    with prior_report_path.open(encoding="utf-8") as stream:
        prior = json.load(stream)
    survivors = {
        item["index"]
        for item in prior["graph_results"]
        if not item["decision"]["rejected"]
    }
    records = [
        record for record in payload["graphs"] if record["index"] in survivors
    ]
    if len(records) != prior["graphs_surviving"] or len(records) != 1_097:
        raise ValueError("prior report and K6 residue input do not select 1,097 graphs")
    selected = _sample_records(records, sample_size)
    started = time.time()
    results = []
    totals = {
        "seeds_checked": 0,
        "impossible_seeds": 0,
        "Z0_subsets_considered": 0,
        "Z0_support_matchable": 0,
        "bipartite_components_checked": 0,
        "inertia_failures": 0,
    }
    rejected = []
    for record in selected:
        decision = evaluate_normal_graph(record["adjacency"], scan_all=True)
        results.append({"index": record["index"], "decision": decision.jsonable()})
        if decision.rejected:
            rejected.append(record["index"])
        totals["seeds_checked"] += decision.seeds_checked
        totals["impossible_seeds"] += decision.impossible_seeds
        totals["Z0_subsets_considered"] += decision.z0_subsets_considered
        totals["Z0_support_matchable"] += decision.z0_support_matchable
        totals["bipartite_components_checked"] += decision.bipartite_components_checked
        totals["inertia_failures"] += decision.inertia_failures

    canonical_selection = json.dumps(
        [record["index"] for record in selected], separators=(",", ":")
    ).encode("ascii")
    report = {
        "schema": "d6-k6-normal-inertia-sample-v1",
        "status": "COMPLETE",
        "filter": "K6_bipartite_component_I_plus_adjacency_exact_inertia",
        "scope": "deterministic sample of the 1,097 exact prior K6 survivors",
        "sample_rule": "smallest SHA256(decimal corpus index), bytewise",
        "sample_size": sample_size,
        "sample_indices": [record["index"] for record in selected],
        "sample_indices_sha256": hashlib.sha256(canonical_selection).hexdigest(),
        "graphs_rejected": len(rejected),
        "graphs_surviving": sample_size - len(rejected),
        "rejected_indices": rejected,
        **totals,
        "graph_results": results,
        "inputs": {
            str(input_path): sha256_file(input_path),
            str(prior_report_path): sha256_file(prior_report_path),
            Path(__file__).name: sha256_file(Path(__file__).resolve()),
        },
        "runtime": {
            "wall_seconds": time.time() - started,
            "python": sys.version,
            "platform": platform.platform(),
            "command": " ".join(sys.argv),
            "workers": 1,
        },
        "arithmetic": "exact integers and fractions; wall time only is floating point",
    }
    _atomic_json(output_path, report)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=ROOT / "d6_k6_bipartite_rank_input.json",
    )
    parser.add_argument(
        "--prior-report",
        type=Path,
        default=ROOT / "d6_k6_bipartite_rank_report.json",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "d6_k6_normal_inertia_sample.json",
    )
    parser.add_argument("--sample-size", type=int, default=64)
    args = parser.parse_args()
    report = run_sample(
        args.input.resolve(),
        args.prior_report.resolve(),
        args.output.resolve(),
        args.sample_size,
    )
    print(json.dumps({
        "status": report["status"],
        "sample_size": report["sample_size"],
        "graphs_rejected": report["graphs_rejected"],
        "rejected_indices": report["rejected_indices"],
        "output": str(args.output.resolve()),
    }, sort_keys=True))


if __name__ == "__main__":
    main()
