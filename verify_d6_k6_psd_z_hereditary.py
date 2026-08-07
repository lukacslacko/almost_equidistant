#!/usr/bin/env python3
"""Independent verifier for hereditary PSD--Z support Hall.

This checker imports neither the production evaluator, its probe, nor a new
shared kernel.  It reconstructs Lorentz and side components, uses the frozen
independent SymPy/Sturm inertia and simultaneous zero-forcing implementations,
derives PSD--Z applicability by sign and bipartiteness, enumerates exactly one
subset per connected side span, and checks every production certificate.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
import os
import platform
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import verify_d6_k6_all_same_z0 as zf_reference
import verify_d6_k6_same_z0 as base
from test_d6_k6_lorentz import lower_bound_18_graph


ROOT = Path(__file__).resolve().parent
COORDINATES = 6
REPORT_SCHEMA = "d6-k6-psd-z-hereditary-v1"
CERTIFICATE_SCHEMA = "d6-k6-psd-z-hereditary-certificates-v1"
VERIFICATION_SCHEMA = "d6-k6-psd-z-hereditary-verification-v1"
EXPECTED_INPUT = 861
EXPECTED_INPUT_SHA256 = (
    "09ebce17d2b72fa6514fc6a8a938376626d373161d2e1b4dd8d31c5e00c18db5"
)
EXPECTED = {
    "verify_d6_k6_same_z0.py": (
        "4db05ff06910a11395e4bbdd2fd13161d6c19f58f81cc99851b8c37558d7000f"
    ),
    "verify_d6_k6_all_same_z0.py": (
        "9ad84985cffbc6f53562730f1aa536a6d7046709c21af0dfd1bfbab907a32162"
    ),
    "d6_k6_psd_zmatrix_report.json": (
        "abb920fbd9eb848502cca37e32c2aa86bf220ab1773e6006ea1fbee4bc2f4e30"
    ),
    "d6_k6_bipartite_rank_input.json": (
        "dfee979e6ca2759e80cf803c3054c8954d98543875d467644b91ee434c22e845"
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


def atomic_json(path: Path, value: object) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def connected_parts(rows: Sequence[int]) -> tuple[tuple[int, ...], ...]:
    unseen = set(range(len(rows)))
    answer = []
    while unseen:
        root = min(unseen)
        unseen.remove(root)
        reached = {root}
        frontier = [root]
        while frontier:
            vertex = frontier.pop()
            new = [other for other in sorted(unseen) if rows[vertex] & (1 << other)]
            for other in new:
                unseen.remove(other)
                reached.add(other)
                frontier.append(other)
        answer.append(tuple(sorted(reached)))
    return tuple(answer)


def induced_rows(rows: Sequence[int], selected: Sequence[int]) -> tuple[int, ...]:
    return tuple(
        sum(1 << column for column, other in enumerate(selected) if rows[vertex] & (1 << other))
        for vertex in selected
    )


def induced_adjacency(
    adj: Sequence[int], instance: base.Instance, selected: int
) -> tuple[tuple[int, ...], tuple[int, ...], tuple[int, ...]]:
    locals_ = tuple(base.bits(selected))
    absolute = tuple(instance.outside[local] for local in locals_)
    rows = tuple(
        sum(1 << column for column, other in enumerate(absolute) if adj[vertex] & (1 << other))
        for vertex in absolute
    )
    return rows, absolute, locals_


def is_bipartite(rows: Sequence[int]) -> bool:
    colors: dict[int, int] = {}
    for root in range(len(rows)):
        if root in colors:
            continue
        colors[root] = 0
        stack = [root]
        while stack:
            vertex = stack.pop()
            for other in base.bits(rows[vertex]):
                if other in colors:
                    if colors[other] == colors[vertex]:
                        return False
                else:
                    colors[other] = 1 - colors[vertex]
                    stack.append(other)
    return True


@dataclass(frozen=True)
class Option:
    rank: int
    coordinates: int
    selected: tuple[int, ...]
    kind: str


@dataclass(frozen=True)
class Span:
    label: str
    options: tuple[Option, ...]
    applicable: bool


def span_options(
    label: str,
    absolute: Sequence[int],
    masks: Sequence[int],
    full_rank: int,
    applicable: bool,
) -> Span:
    size = len(absolute)
    options = [Option(0, 0, (), "empty")]
    subsets = range(1, 1 << size) if applicable else ((1 << size) - 1,)
    for encoded in subsets:
        selected_positions = [i for i in range(size) if encoded & (1 << i)]
        coordinates = 0
        for i in selected_positions:
            coordinates |= masks[i]
        full = encoded == (1 << size) - 1
        options.append(Option(
            full_rank if full else len(selected_positions),
            coordinates,
            tuple(absolute[i] for i in selected_positions),
            "full" if full else "proper",
        ))
    return Span(label, tuple(options), applicable)


def exact_hall(spans: Sequence[Span]) -> tuple[bool, dict | None]:
    """Brute-force Cartesian product: one alternative from each span."""

    witness: list[tuple[str, Option]] = []

    def visit(position: int, rank: int, coordinates: int) -> dict | None:
        if position == len(spans):
            return None
        span = spans[position]
        for option in span.options:
            union = coordinates | option.coordinates
            total = rank + option.rank
            if option.selected:
                witness.append((span.label, option))
            if total > union.bit_count():
                failure = {
                    "rank_lower": total,
                    "coordinate_capacity": union.bit_count(),
                    "allowed_coordinates": base.bits(union),
                    "selected_groups": [
                        {
                            "group": label,
                            "kind": item.kind,
                            "rank_lower": item.rank,
                            "selected_vertices": list(item.selected),
                            "allowed_coordinates": base.bits(item.coordinates),
                        }
                        for label, item in witness
                    ],
                }
            else:
                failure = visit(position + 1, total, union)
            if option.selected:
                witness.pop()
            if failure is not None:
                return failure
        return None

    failure = visit(0, 0, 0)
    return failure is None, failure


def side_spans(
    adj: Sequence[int],
    instance: base.Instance,
    selected: int,
    label: str,
    sign: str,
    inertia: base.SympyInertiaCache,
    zero_forcing: zf_reference.IndependentZeroForcing,
) -> tuple[Span, ...]:
    rows, absolute, locals_ = induced_adjacency(adj, instance, selected)
    answer = []
    for number, part in enumerate(connected_parts(rows)):
        local_rows = induced_rows(rows, part)
        pattern = tuple(row | (1 << i) for i, row in enumerate(local_rows))
        positive, negative, _ = inertia.solve(pattern)
        size = len(part)
        zf_rank = size - zero_forcing.number(local_rows)
        inertia_rank = size - (negative if sign == "positive" else positive)
        bipartite = is_bipartite(local_rows)
        applicable = sign == "negative" or bipartite
        full_rank = max(zf_rank, inertia_rank, size - 1 if applicable else 0)
        part_absolute = tuple(absolute[i] for i in part)
        masks = tuple(instance.defects[locals_[i]] for i in part)
        answer.append(span_options(
            f"{label}:F{number}", part_absolute, masks, full_rank, applicable
        ))
    return tuple(answer)


def orientation_spans(
    adj: Sequence[int],
    instance: base.Instance,
    z0: int,
    component: base.Component,
    sign_a: str,
    sign_b: str,
    inertia: base.SympyInertiaCache,
    zero_forcing: zf_reference.IndependentZeroForcing,
) -> tuple[Span, ...]:
    spans = side_spans(adj, instance, component.side_a, "A", sign_a, inertia, zero_forcing)
    spans += side_spans(adj, instance, component.side_b, "B", sign_b, inertia, zero_forcing)
    spans += tuple(
        Span(
            f"Z0:{instance.outside[local]}",
            (Option(0, 0, (), "empty"), Option(1, instance.defects[local], (instance.outside[local],), "Z0")),
            True,
        )
        for local in base.bits(z0)
    )
    return spans


def check_component(
    adj: Sequence[int],
    instance: base.Instance,
    z0: int,
    component: base.Component,
    inertia: base.SympyInertiaCache,
    zero_forcing: zf_reference.IndependentZeroForcing,
) -> dict:
    orientations = []
    for name, sign_a, sign_b in (
        ("A_positive", "positive", "negative"),
        ("A_negative", "negative", "positive"),
    ):
        spans = orientation_spans(
            adj, instance, z0, component, sign_a, sign_b, inertia, zero_forcing
        )
        passed, failure = exact_hall(spans)
        orientations.append({"case": name, "passed": passed, "failure": failure, "spans": spans})
    light = z0 | component.vertices
    light_passed = light.bit_count() <= COORDINATES and base.hall_matchable(light, instance.defects)
    passed = any(row["passed"] for row in orientations) or light_passed
    return {
        "passed": passed,
        "orientations": orientations,
        "lightlike_passed": light_passed,
        "component": [instance.outside[local] for local in base.bits(component.vertices)],
    }


def solve_seed(adj: Sequence[int], instance: base.Instance) -> dict:
    inertia = base.SympyInertiaCache()
    zero_forcing = zf_reference.IndependentZeroForcing()
    choices = []
    considered = 0
    for z0 in base.z0_subsets(instance.eligible_z0):
        considered += 1
        matchable = base.hall_matchable(z0, instance.defects)
        if not matchable:
            choices.append({"Z0": [instance.outside[i] for i in base.bits(z0)], "matchable": False, "failure": None})
            continue
        failure = None
        for component in base.components(instance, z0):
            if not component.bipartite:
                continue
            result = check_component(adj, instance, z0, component, inertia, zero_forcing)
            if not result["passed"]:
                failure = result
                break
        if failure is None:
            return {"feasible": True, "z0_considered": considered, "choices": None}
        choices.append({
            "Z0": [instance.outside[i] for i in base.bits(z0)],
            "matchable": True,
            "failure": failure,
        })
    return {"feasible": False, "z0_considered": considered, "choices": choices}


def evaluate_graph(adj: Sequence[int]) -> dict:
    base.validate_graph(adj)
    if next(base.clique_masks(adj, 7), 0):
        raise ValueError("independent hereditary target contains K7")
    seeds_checked = 0
    z0_considered = 0
    for seed_mask in base.clique_masks(adj, COORDINATES):
        seeds_checked += 1
        instance = base.build_instance(adj, seed_mask)
        result = solve_seed(adj, instance)
        z0_considered += result["z0_considered"]
        if not result["feasible"]:
            return {
                "rejected": True,
                "seeds_checked": seeds_checked,
                "z0_considered": z0_considered,
                "first_impossible_seed": list(instance.seed),
                "instance": instance,
                "choices": result["choices"],
            }
    return {
        "rejected": False,
        "seeds_checked": seeds_checked,
        "z0_considered": z0_considered,
        "first_impossible_seed": None,
        "instance": None,
        "choices": None,
    }


def validate_archived_certificate(
    adj: Sequence[int], result: dict, archived: dict
) -> None:
    instance = result["instance"]
    if archived["seed"] != list(instance.seed):
        raise AssertionError("archived hereditary seed differs")
    expected_z0 = [row["Z0"] for row in result["choices"]]
    archived_z0 = [row["Z0"] for row in archived["choices"]]
    if archived_z0 != expected_z0:
        raise AssertionError("archived hereditary Z0 enumeration is incomplete")
    local_by_absolute = {absolute: local for local, absolute in enumerate(instance.outside)}
    for independent_row, row in zip(result["choices"], archived["choices"]):
        if row["matchable"] != independent_row["matchable"]:
            raise AssertionError("archived hereditary matchability differs")
        if not row["matchable"]:
            continue
        z0 = sum(1 << local_by_absolute[x] for x in row["Z0"])
        archived_failure = row["first_failed_component"]
        candidates = {
            tuple(instance.outside[i] for i in base.bits(component.vertices)): component
            for component in base.components(instance, z0)
            if component.bipartite
        }
        key = tuple(archived_failure["component"])
        if key not in candidates:
            raise AssertionError("archived failed component is not a bipartite L component")
        component = candidates[key]
        fresh = check_component(
            adj, instance, z0, component,
            base.SympyInertiaCache(), zf_reference.IndependentZeroForcing(),
        )
        if fresh["passed"]:
            raise AssertionError("archived component independently passes")
        light = z0 | component.vertices
        expected_light = light.bit_count() <= COORDINATES and base.hall_matchable(light, instance.defects)
        if archived_failure["lightlike"]["passed"] != expected_light:
            raise AssertionError("archived lightlike decision differs")
        fresh_by_case = {item["case"]: item for item in fresh["orientations"]}
        for archived_orientation in archived_failure["orientations"]:
            case = archived_orientation["case"]
            current = fresh_by_case[case]
            if archived_orientation["passed"] or current["passed"]:
                raise AssertionError("failed certificate orientation unexpectedly passes")
            witness = archived_orientation["first_failure"]
            catalog = {span.label: span for span in current["spans"]}
            rank = 0
            coordinates = 0
            used = set()
            for selected in witness["selected_groups"]:
                label = selected["group"]
                if label in used or label not in catalog:
                    raise AssertionError("certificate double-counts or invents a span")
                used.add(label)
                wanted = (
                    selected["rank_lower"],
                    tuple(selected["selected_vertices"]),
                    tuple(selected["allowed_coordinates"]),
                )
                valid = [
                    option for option in catalog[label].options
                    if (
                        option.rank,
                        option.selected,
                        tuple(base.bits(option.coordinates)),
                    ) == wanted
                ]
                if not valid:
                    raise AssertionError("certificate selected an invalid span subset")
                rank += valid[0].rank
                coordinates |= valid[0].coordinates
            if rank <= coordinates.bit_count():
                raise AssertionError("archived hereditary witness is not a Hall violation")
            if rank != witness["rank_lower"] or coordinates.bit_count() != witness["coordinate_capacity"]:
                raise AssertionError("archived hereditary witness totals differ")


def kernel_controls() -> dict:
    overlap = span_options("overlap", (0, 1, 2), (1, 2, 4), 2, True)
    one_passed, _ = exact_hall((overlap,))
    two_passed, _ = exact_hall((overlap, overlap))
    strict = span_options("strict", (0, 1, 2), (1, 1, 2), 2, True)
    strict_passed, strict_failure = exact_hall((strict,))
    if not one_passed or two_passed or strict_passed:
        raise AssertionError("independent hereditary span controls failed")
    # Fresh PF applicability controls: P3 is bipartite; C3 is not.
    if not is_bipartite((2, 5, 2)) or is_bipartite((6, 5, 3)):
        raise AssertionError("independent bipartite controls failed")
    return {
        "one_span_overlap_passed": one_passed,
        "two_distinct_spans_failed": not two_passed,
        "proper_subset_strict_failure": strict_failure,
        "P3_bipartite": True,
        "C3_bipartite": False,
    }


def verify(output: Path) -> dict:
    observed = {name: sha256(ROOT / name) for name in EXPECTED}
    if observed != EXPECTED:
        raise ValueError(f"independent hereditary dependency mismatch: {observed}")
    report = json.loads((ROOT / "d6_k6_psd_z_hereditary_report.json").read_text())
    archive = json.loads((ROOT / "d6_k6_psd_z_hereditary_certificates.json").read_text())
    if report.get("schema") != REPORT_SCHEMA or report.get("status") != "COMPLETE":
        raise ValueError("hereditary report schema/status differs")
    if archive.get("schema") != CERTIFICATE_SCHEMA or archive.get("status") != "COMPLETE":
        raise ValueError("hereditary archive schema/status differs")
    if report["certificate_archive"]["sha256"] != sha256(ROOT / report["certificate_archive"]["path"]):
        raise AssertionError("hereditary archive hash differs")
    parent_report = json.loads((ROOT / "d6_k6_psd_zmatrix_report.json").read_text())
    indices = [row["index"] for row in parent_report["graph_results"] if not row["decision"]["rejected"]]
    if len(indices) != EXPECTED_INPUT or base.stable_hash(indices) != EXPECTED_INPUT_SHA256:
        raise ValueError("independent hereditary input boundary differs")
    payload = json.loads((ROOT / "d6_k6_bipartite_rank_input.json").read_text())
    records = {row["index"]: row["adjacency"] for row in payload["graphs"]}
    production = {row["index"]: row["decision"] for row in report["graph_results"]}
    certificates = {row["index"]: row["certificate"] for row in archive["certificates"]}
    controls = kernel_controls()
    started = time.perf_counter()
    rejected = []
    summaries = []
    for index in indices:
        result = evaluate_graph(records[index])
        expected = production[index]
        for field in ("rejected", "seeds_checked", "first_impossible_seed"):
            if result[field] != expected[field]:
                raise AssertionError(f"graph {index} independent {field} differs")
        if not result["rejected"]:
            continue
        rejected.append(index)
        validate_archived_certificate(records[index], result, certificates[index])
        summaries.append({
            "index": index,
            "first_impossible_seed": result["first_impossible_seed"],
            "Z0_choices": len(result["choices"]),
        })
    if rejected != report["rejected_indices"] or rejected != archive["rejected_indices"]:
        raise AssertionError("independent hereditary rejection set differs")
    positive = evaluate_graph(lower_bound_18_graph())
    if positive["rejected"] or positive["seeds_checked"] != 32:
        raise AssertionError("known realizable 18-point control failed")
    tree = ast.parse(Path(__file__).read_text())
    imports = {
        node.module for node in ast.walk(tree)
        if isinstance(node, ast.ImportFrom) and node.module
    }
    imports.update(
        alias.name for node in ast.walk(tree) if isinstance(node, ast.Import)
        for alias in node.names
    )
    forbidden = {
        "d6_k6_psd_z_hereditary", "probe_d6_k6_psd_z_hereditary"
    }
    if imports & forbidden:
        raise AssertionError("independent verifier imported a forbidden kernel")
    result = {
        "schema": VERIFICATION_SCHEMA,
        "status": "PASS",
        "method": (
            "fresh Lorentz/side components and bipartitions; independent "
            "SymPy/Sturm inertia and simultaneous zero forcing; brute-force "
            "one-subset-per-span Hall; exhaustive Z0 and certificate replay"
        ),
        "graphs_recomputed": len(indices),
        "graphs_rejected": len(rejected),
        "graphs_surviving": len(indices) - len(rejected),
        "rejected_indices": rejected,
        "certificate_summaries": summaries,
        "kernel_controls": controls,
        "positive_control": {"passed": True, "K6_seeds": positive["seeds_checked"]},
        "inputs": {**observed, Path(__file__).name: sha256(Path(__file__).resolve())},
        "runtime": {
            "command": " ".join(sys.argv),
            "workers": 1,
            "wall_seconds": time.perf_counter() - started,
            "python": sys.version,
            "platform": platform.platform(),
        },
    }
    atomic_json(output, result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "d6_k6_psd_z_hereditary_verification.json")
    args = parser.parse_args()
    result = verify(args.output.resolve())
    print(json.dumps({
        "status": result["status"],
        "graphs_recomputed": result["graphs_recomputed"],
        "graphs_rejected": result["graphs_rejected"],
        "output": str(args.output.resolve()),
    }, sort_keys=True))


if __name__ == "__main__":
    main()
