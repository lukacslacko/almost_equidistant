#!/usr/bin/env python3
"""Independent verifier for exact K6 arbitrary-subset rank Hall.

This checker imports neither the new production evaluator nor the discovery
probe.  It extends the frozen independent hereditary checker with a fresh
positive-sign subset-rank reconstruction.  SymPy/Sturm inertia and the
independent simultaneous zero-forcing solver are retained; Hall alternatives
are brute-forced with exactly one subset per original connected span.
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
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Sequence

import verify_d6_k6_psd_z_hereditary as prior
from test_d6_k6_lorentz import lower_bound_18_graph


base = prior.base
zf_reference = prior.zf_reference
ROOT = Path(__file__).resolve().parent
COORDINATES = 6
REPORT_SCHEMA = "d6-k6-arbitrary-subset-hall-v1"
CERTIFICATE_SCHEMA = "d6-k6-arbitrary-subset-hall-certificates-v1"
VERIFICATION_SCHEMA = "d6-k6-arbitrary-subset-hall-verification-v1"
EXPECTED_INPUT = 831
EXPECTED_INPUT_SHA256 = (
    "2bcdad095c6bd3038a4bb1d117f2faadfa98113d814633351ec9e79ff553a44b"
)
EXPECTED = {
    "verify_d6_k6_psd_z_hereditary.py": (
        "297e749923ae60b84cda12170d19509da237b665965d01af6214465f4e4df17f"
    ),
    "verify_d6_k6_same_z0.py": (
        "4db05ff06910a11395e4bbdd2fd13161d6c19f58f81cc99851b8c37558d7000f"
    ),
    "verify_d6_k6_all_same_z0.py": (
        "9ad84985cffbc6f53562730f1aa536a6d7046709c21af0dfd1bfbab907a32162"
    ),
    "d6_k6_psd_z_hereditary_report.json": (
        "202a844d6505d3c68c9a0bea8a5d82a983a7e44b96cb9f3c711d80c47f6c00c1"
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


def components_in_mask(rows: Sequence[int], selected: int) -> tuple[int, ...]:
    remaining = selected
    answer = []
    while remaining:
        root = remaining & -remaining
        reached = {root.bit_length() - 1}
        frontier = list(reached)
        while frontier:
            vertex = frontier.pop()
            for other in base.bits(rows[vertex] & selected):
                if other not in reached:
                    reached.add(other)
                    frontier.append(other)
        mask = sum(1 << vertex for vertex in reached)
        answer.append(mask)
        remaining &= ~mask
    return tuple(answer)


def induced_from_mask(rows: Sequence[int], selected: int) -> tuple[int, ...]:
    chosen = base.bits(selected)
    return tuple(
        sum(1 << column for column, other in enumerate(chosen) if rows[vertex] & (1 << other))
        for vertex in chosen
    )


def positive_subset_rank(
    rows: Sequence[int],
    selected: int,
    inertia: base.SympyInertiaCache,
    zero_forcing: zf_reference.IndependentZeroForcing,
) -> int:
    total = 0
    for component in components_in_mask(rows, selected):
        local = induced_from_mask(rows, component)
        size = len(local)
        pattern = tuple(row | (1 << i) for i, row in enumerate(local))
        _, negative, _ = inertia.solve(pattern)
        zf_rank = size - zero_forcing.number(local)
        inertia_rank = size - negative
        psd_z_rank = size - 1 if prior.is_bipartite(local) else 0
        total += max(zf_rank, inertia_rank, psd_z_rank)
    return total


def side_spans(
    adj: Sequence[int],
    instance: base.Instance,
    selected: int,
    label: str,
    sign: str,
    inertia: base.SympyInertiaCache,
    zero_forcing: zf_reference.IndependentZeroForcing,
    cache: dict[tuple[int, str, str], tuple[prior.Span, ...]],
) -> tuple[prior.Span, ...]:
    key = selected, label, sign
    if key in cache:
        return cache[key]
    rows, absolute, locals_ = prior.induced_adjacency(adj, instance, selected)
    answer = []
    for number, part in enumerate(prior.connected_parts(rows)):
        local_rows = prior.induced_rows(rows, part)
        pattern = tuple(row | (1 << i) for i, row in enumerate(local_rows))
        positive, negative, _ = inertia.solve(pattern)
        size = len(part)
        zf_rank = size - zero_forcing.number(local_rows)
        inertia_rank = size - (negative if sign == "positive" else positive)
        bipartite = prior.is_bipartite(local_rows)
        applicable = sign == "negative" or bipartite
        full_rank = max(zf_rank, inertia_rank, size - 1 if applicable else 0)
        part_absolute = tuple(absolute[i] for i in part)
        masks = tuple(instance.defects[locals_[i]] for i in part)
        if applicable:
            span = prior.span_options(
                f"{label}:F{number}", part_absolute, masks, full_rank, True
            )
        else:
            raw = [prior.Option(0, 0, (), "empty")]
            for subset in range(1, 1 << size):
                subset_rank = positive_subset_rank(
                    local_rows, subset, inertia, zero_forcing
                )
                full = subset == (1 << size) - 1
                if full and subset_rank != full_rank:
                    raise AssertionError("independent full arbitrary rank differs")
                coordinates = 0
                chosen = []
                for i in base.bits(subset):
                    coordinates |= masks[i]
                    chosen.append(part_absolute[i])
                raw.append(prior.Option(
                    subset_rank,
                    coordinates,
                    tuple(chosen),
                    "full" if full else "arbitrary_positive_nonbip",
                ))
            span = prior.Span(f"{label}:F{number}", tuple(raw), False)
        answer.append(span)
    result = tuple(answer)
    cache[key] = result
    return result


def orientation_spans(
    adj: Sequence[int],
    instance: base.Instance,
    z0: int,
    component: base.Component,
    sign_a: str,
    sign_b: str,
    inertia: base.SympyInertiaCache,
    zero_forcing: zf_reference.IndependentZeroForcing,
    cache: dict[tuple[int, str, str], tuple[prior.Span, ...]],
) -> tuple[prior.Span, ...]:
    spans = side_spans(
        adj, instance, component.side_a, "A", sign_a, inertia, zero_forcing,
        cache,
    )
    spans += side_spans(
        adj, instance, component.side_b, "B", sign_b, inertia, zero_forcing,
        cache,
    )
    spans += tuple(
        prior.Span(
            f"Z0:{instance.outside[local]}",
            (
                prior.Option(0, 0, (), "empty"),
                prior.Option(1, instance.defects[local], (instance.outside[local],), "Z0"),
            ),
            True,
        )
        for local in base.bits(z0)
    )
    return spans


def independent_hall(
    spans: Sequence[prior.Span],
) -> tuple[bool, dict | None]:
    """Exact coordinate-container enumeration, independent of production DP.

    For each of the 64 coordinate sets C, choose independently in every
    original span the largest-rank subset alternative supported inside C.
    A violating combination exists iff one of these 64 maxima exceeds |C|.
    Every subset alternative is inspected; alternatives from one span are
    maximized, never added.
    """

    for container in range(1 << COORDINATES):
        rank = 0
        selected = []
        union = 0
        for span in spans:
            eligible = [
                option for option in span.options
                if not option.coordinates & ~container
            ]
            choice = max(eligible, key=lambda option: option.rank)
            rank += choice.rank
            union |= choice.coordinates
            if choice.selected:
                selected.append((span.label, choice))
        if rank > container.bit_count():
            return False, {
                "rank_lower": rank,
                "coordinate_capacity": union.bit_count(),
                "allowed_coordinates": base.bits(union),
                "selected_groups": [
                    {
                        "group": label,
                        "kind": option.kind,
                        "rank_lower": option.rank,
                        "selected_vertices": list(option.selected),
                        "allowed_coordinates": base.bits(option.coordinates),
                    }
                    for label, option in selected
                ],
            }
    return True, None


def check_component(
    adj: Sequence[int],
    instance: base.Instance,
    z0: int,
    component: base.Component,
    inertia: base.SympyInertiaCache,
    zero_forcing: zf_reference.IndependentZeroForcing,
    cache: dict[tuple[int, str, str], tuple[prior.Span, ...]],
) -> dict:
    orientations = []
    for name, sign_a, sign_b in (
        ("A_positive", "positive", "negative"),
        ("A_negative", "negative", "positive"),
    ):
        spans = orientation_spans(
            adj, instance, z0, component, sign_a, sign_b, inertia, zero_forcing,
            cache,
        )
        passed, failure = independent_hall(spans)
        orientations.append({
            "case": name,
            "passed": passed,
            "failure": failure,
            "spans": spans,
        })
    light = z0 | component.vertices
    light_passed = (
        light.bit_count() <= COORDINATES
        and base.hall_matchable(light, instance.defects)
    )
    return {
        "passed": any(item["passed"] for item in orientations) or light_passed,
        "orientations": orientations,
        "lightlike_passed": light_passed,
        "component": [instance.outside[i] for i in base.bits(component.vertices)],
    }


def solve_seed(adj: Sequence[int], instance: base.Instance) -> dict:
    inertia = base.SympyInertiaCache()
    zero_forcing = zf_reference.IndependentZeroForcing()
    cache: dict[tuple[int, str, str], tuple[prior.Span, ...]] = {}
    rows = []
    considered = 0
    for z0 in base.z0_subsets(instance.eligible_z0):
        considered += 1
        matchable = base.hall_matchable(z0, instance.defects)
        if not matchable:
            rows.append({
                "Z0": [instance.outside[i] for i in base.bits(z0)],
                "matchable": False,
                "failure": None,
            })
            continue
        failure = None
        for component in base.components(instance, z0):
            if not component.bipartite:
                continue
            result = check_component(
                adj, instance, z0, component, inertia, zero_forcing, cache
            )
            if not result["passed"]:
                failure = result
                break
        if failure is None:
            return {"feasible": True, "z0_considered": considered, "rows": None}
        rows.append({
            "Z0": [instance.outside[i] for i in base.bits(z0)],
            "matchable": True,
            "failure": failure,
        })
    return {"feasible": False, "z0_considered": considered, "rows": rows}


def evaluate_graph(adj: Sequence[int]) -> dict:
    base.validate_graph(adj)
    if next(base.clique_masks(adj, 7), 0):
        raise ValueError("independent arbitrary-subset target contains K7")
    seeds_checked = z0_considered = 0
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
                "rows": result["rows"],
            }
    return {
        "rejected": False,
        "seeds_checked": seeds_checked,
        "z0_considered": z0_considered,
        "first_impossible_seed": None,
        "instance": None,
        "rows": None,
    }


def evaluate_task(task: tuple[int, Sequence[int]]) -> dict:
    index, adjacency = task
    return {"index": index, "result": evaluate_graph(adjacency)}


def positive_seed_task(task: tuple[Sequence[int], int]) -> dict:
    adjacency, seed_mask = task
    instance = base.build_instance(adjacency, seed_mask)
    return {
        "seed": list(instance.seed),
        "passed": solve_seed(adjacency, instance)["feasible"],
    }


def validate_certificate(adj: Sequence[int], result: dict, archived: dict) -> None:
    instance = result["instance"]
    if archived["seed"] != list(instance.seed):
        raise AssertionError("archived arbitrary-subset seed differs")
    if archived["outside"] != list(instance.outside):
        raise AssertionError("archived arbitrary-subset outside order differs")
    expected_eligible = [
        instance.outside[i] for i in base.bits(instance.eligible_z0)
    ]
    if archived["eligible_Z0"] != expected_eligible:
        raise AssertionError("archived arbitrary-subset eligible Z0 differs")
    expected_defects = {
        str(instance.outside[local]): [
            instance.seed[coordinate]
            for coordinate in base.bits(mask)
        ]
        for local, mask in enumerate(instance.defects)
    }
    if archived["allowed_defects"] != expected_defects:
        raise AssertionError("archived arbitrary-subset allowed masks differ")
    if [row["Z0"] for row in archived["choices"]] != [row["Z0"] for row in result["rows"]]:
        raise AssertionError("archived arbitrary-subset Z0 list is incomplete")
    local_by_absolute = {absolute: local for local, absolute in enumerate(instance.outside)}
    inertia = base.SympyInertiaCache()
    zero_forcing = zf_reference.IndependentZeroForcing()
    cache: dict[tuple[int, str, str], tuple[prior.Span, ...]] = {}
    for independent_row, row in zip(result["rows"], archived["choices"]):
        if row["matchable"] != independent_row["matchable"]:
            raise AssertionError("archived arbitrary-subset matchability differs")
        if not row["matchable"]:
            continue
        z0 = sum(1 << local_by_absolute[x] for x in row["Z0"])
        archived_failure = row["first_failed_component"]
        components = {
            tuple(instance.outside[i] for i in base.bits(component.vertices)): component
            for component in base.components(instance, z0)
            if component.bipartite
        }
        key = tuple(archived_failure["component"])
        if key not in components:
            raise AssertionError("archived failed component is not present")
        component = components[key]
        fresh = check_component(
            adj, instance, z0, component,
            inertia, zero_forcing, cache,
        )
        if fresh["passed"]:
            raise AssertionError("archived arbitrary-subset component passes")
        light = z0 | component.vertices
        expected_light = (
            light.bit_count() <= COORDINATES
            and base.hall_matchable(light, instance.defects)
        )
        if archived_failure["lightlike"]["passed"] != expected_light:
            raise AssertionError("archived arbitrary-subset lightlike result differs")
        if [item["case"] for item in archived_failure["orientations"]] != [
            "A_positive", "A_negative"
        ]:
            raise AssertionError("archive omits a generic orientation")
        fresh_cases = {item["case"]: item for item in fresh["orientations"]}
        for saved in archived_failure["orientations"]:
            current = fresh_cases[saved["case"]]
            if saved["passed"] or current["passed"]:
                raise AssertionError("archived failed orientation passes")
            witness = saved["first_failure"]
            catalog = {span.label: span for span in current["spans"]}
            used = set()
            rank = coordinates = 0
            for selected in witness["selected_groups"]:
                label = selected["group"]
                if label in used or label not in catalog:
                    raise AssertionError("certificate double-counts/invents a span")
                used.add(label)
                wanted = (
                    selected["kind"], selected["rank_lower"],
                    tuple(selected["selected_vertices"]),
                    tuple(selected["allowed_coordinates"]),
                )
                valid = [
                    option for option in catalog[label].options
                    if (
                        option.kind, option.rank, option.selected,
                        tuple(base.bits(option.coordinates)),
                    ) == wanted
                ]
                if not valid:
                    raise AssertionError("certificate subset/rank is invalid")
                rank += valid[0].rank
                coordinates |= valid[0].coordinates
            if rank <= coordinates.bit_count():
                raise AssertionError("certificate is not a Hall violation")
            if rank != witness["rank_lower"] or base.bits(coordinates) != witness["allowed_coordinates"]:
                raise AssertionError("certificate Hall totals differ")


def kernel_controls() -> dict:
    inertia = base.SympyInertiaCache()
    zero_forcing = zf_reference.IndependentZeroForcing()
    c5 = (18, 5, 10, 20, 9)
    if positive_subset_rank(c5, 7, inertia, zero_forcing) != 2:
        raise AssertionError("independent C5/P3 rank control failed")
    masks = (1, 2, 1, 4, 8)
    raw = [prior.Option(0, 0, (), "empty")]
    for selected in range(1, 32):
        coordinates = 0
        for i in base.bits(selected):
            coordinates |= masks[i]
        raw.append(prior.Option(
            positive_subset_rank(c5, selected, inertia, zero_forcing),
            coordinates,
            tuple(base.bits(selected)),
            "full" if selected == 31 else "arbitrary_positive_nonbip",
        ))
    passed, failure = independent_hall((prior.Span("C5", tuple(raw), False),))
    if passed:
        raise AssertionError("independent arbitrary C5 Hall should fail")
    return {
        "C5_selected_P3_rank": 2,
        "arbitrary_C5_Hall_failed": True,
        "failure": failure,
        "one_vertex_extension_omitted": True,
    }


def verify(output: Path, workers: int) -> dict:
    observed = {name: sha256(ROOT / name) for name in EXPECTED}
    if observed != EXPECTED:
        raise ValueError(f"independent arbitrary-subset dependency mismatch: {observed}")
    report = json.loads((ROOT / "d6_k6_arbitrary_subset_hall_report.json").read_text())
    archive = json.loads((ROOT / "d6_k6_arbitrary_subset_hall_certificates.json").read_text())
    if report.get("schema") != REPORT_SCHEMA or report.get("status") != "COMPLETE":
        raise ValueError("arbitrary-subset report schema/status differs")
    if archive.get("schema") != CERTIFICATE_SCHEMA or archive.get("status") != "COMPLETE":
        raise ValueError("arbitrary-subset archive schema/status differs")
    if report["certificate_archive"]["sha256"] != sha256(ROOT / report["certificate_archive"]["path"]):
        raise AssertionError("arbitrary-subset archive hash differs")
    parent_report = json.loads((ROOT / "d6_k6_psd_z_hereditary_report.json").read_text())
    indices = [row["index"] for row in parent_report["graph_results"] if not row["decision"]["rejected"]]
    if len(indices) != EXPECTED_INPUT or base.stable_hash(indices) != EXPECTED_INPUT_SHA256:
        raise ValueError("independent 831 input boundary differs")
    if [row["index"] for row in report["graph_results"]] != indices:
        raise AssertionError("production graph rows are not the ordered input")
    if report.get("input_graphs") != len(indices) or report.get("input_indices_sha256") != EXPECTED_INPUT_SHA256:
        raise AssertionError("production arbitrary-subset input boundary differs")
    production_hash = sha256(ROOT / "d6_k6_arbitrary_subset_hall.py")
    if (
        report.get("sources", {}).get("d6_k6_arbitrary_subset_hall.py") != production_hash
        or archive.get("production_source_sha256") != production_hash
    ):
        raise AssertionError("production arbitrary-subset source boundary differs")
    payload = json.loads((ROOT / "d6_k6_bipartite_rank_input.json").read_text())
    adjacencies = {row["index"]: row["adjacency"] for row in payload["graphs"]}
    production = {row["index"]: row["decision"] for row in report["graph_results"]}
    certificates = {row["index"]: row["certificate"] for row in archive["certificates"]}
    if (
        report["rejected_indices"] != archive["rejected_indices"]
        or set(certificates) != set(report["rejected_indices"])
    ):
        raise AssertionError("arbitrary-subset report/archive indices differ")
    controls = kernel_controls()
    started = time.perf_counter()
    tasks = [(index, adjacencies[index]) for index in indices]
    if workers == 1:
        recomputed_rows = [evaluate_task(task) for task in tasks]
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            recomputed_rows = list(pool.map(evaluate_task, tasks, chunksize=1))
    recomputed = {row["index"]: row["result"] for row in recomputed_rows}
    rejected = []
    summaries = []
    for index in indices:
        result = recomputed[index]
        expected = production[index]
        for field in ("rejected", "seeds_checked", "first_impossible_seed"):
            if result[field] != expected[field]:
                raise AssertionError(f"graph {index} independent {field} differs")
        if not result["rejected"]:
            continue
        rejected.append(index)
        validate_certificate(adjacencies[index], result, certificates[index])
        summaries.append({
            "index": index,
            "first_impossible_seed": result["first_impossible_seed"],
            "Z0_choices": len(result["rows"]),
        })
    if rejected != report["rejected_indices"] or rejected != archive["rejected_indices"]:
        raise AssertionError("independent arbitrary-subset rejection set differs")
    positive_graph = lower_bound_18_graph()
    positive_seeds = list(base.clique_masks(positive_graph, COORDINATES))
    if workers == 1:
        positive_passes = [
            positive_seed_task((positive_graph, seed))
            for seed in positive_seeds
        ]
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            positive_passes = list(pool.map(
                positive_seed_task,
                [(positive_graph, seed) for seed in positive_seeds],
                chunksize=1,
            ))
    expected_positive_seed_rows = [
        list(base.bits(seed)) for seed in positive_seeds
    ]
    if (
        len(positive_seeds) != 32
        or [row["seed"] for row in positive_passes] != expected_positive_seed_rows
        or not all(row["passed"] for row in positive_passes)
    ):
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
    if imports & {"d6_k6_arbitrary_subset_hall", "probe_d6_k6_arbitrary_subset_hall"}:
        raise AssertionError("independent checker imported a forbidden kernel")
    result = {
        "schema": VERIFICATION_SCHEMA,
        "status": "PASS",
        "method": (
            "fresh arbitrary induced components and componentwise rank; "
            "independent SymPy/Sturm inertia and simultaneous zero forcing; "
            "64 coordinate containers with brute-force subset alternatives "
            "and one choice per original span; exhaustive Z0 replay"
        ),
        "graphs_recomputed": len(indices),
        "graphs_rejected": len(rejected),
        "graphs_surviving": len(indices) - len(rejected),
        "rejected_indices": rejected,
        "certificate_summaries": summaries,
        "kernel_controls": controls,
        "positive_control": {
            "passed": True,
            "K6_seeds": len(positive_seeds),
            "seed_results": positive_passes,
        },
        "inputs": {**observed, Path(__file__).name: sha256(Path(__file__).resolve())},
        "runtime": {
            "command": " ".join(sys.argv),
            "workers": workers,
            "wall_seconds": time.perf_counter() - started,
            "python": sys.version,
            "platform": platform.platform(),
        },
    }
    atomic_json(output, result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "d6_k6_arbitrary_subset_hall_verification.json")
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("workers must be positive")
    result = verify(args.output.resolve(), args.workers)
    print(json.dumps({
        "status": result["status"],
        "graphs_recomputed": result["graphs_recomputed"],
        "graphs_rejected": result["graphs_rejected"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
