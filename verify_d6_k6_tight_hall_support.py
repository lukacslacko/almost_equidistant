#!/usr/bin/env python3
"""Structurally independent verifier for tight-Hall K6 support propagation.

This checker imports neither the new producer nor either discovery probe.  It
uses the prior independent SymPy/Sturm inertia and simultaneous zero-forcing
reconstruction.  Unlike production's union DP, tight coordinate sets are
found by enumerating all 64 coordinate containers and maximizing separately
inside each original orthogonal span.
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
from typing import Iterable, Sequence

import verify_d6_k6_arbitrary_subset_hall as reference
from test_d6_k6_lorentz import lower_bound_18_graph


base = reference.base
prior = reference.prior
zf_reference = reference.zf_reference
ROOT = Path(__file__).resolve().parent
REPORT_SCHEMA = "d6-k6-tight-hall-support-v1"
CERTIFICATE_SCHEMA = "d6-k6-tight-hall-support-certificates-v1"
VERIFICATION_SCHEMA = "d6-k6-tight-hall-support-verification-v1"
EXPECTED_INPUT = 756
EXPECTED_INPUT_SHA256 = "2cfedbc83f6ff01b6e386440fb7ea066b274d22ffbc52371c964e99d17cdcb90"
EXPECTED_MANIFEST_SHA256 = "164b2a12813845ef1f6d6ca5e3713ed1d2edac4be58563c1648a39ff8580c0b5"
COORDINATES = 6
GLOBAL_EMPTY_BUDGET = 2


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


def worker_initializer() -> None:
    for name in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[name] = "1"


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
        "d6_k6_tight_hall_support",
        "probe_d6_k6_tight_support",
        "probe_d6_k6_rectangular_minor",
    }
    if imports & forbidden:
        raise AssertionError(f"independent checker imports new code: {imports & forbidden}")


def load_input() -> tuple[list[dict], list[int]]:
    manifest_path = ROOT / "d6_current_residue_manifest_v5.json"
    verification_path = ROOT / "d6_current_residue_manifest_v5_verification.json"
    if sha256(manifest_path) != EXPECTED_MANIFEST_SHA256:
        raise ValueError("v5 manifest changed")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    verification = json.loads(verification_path.read_text(encoding="utf-8"))
    block = manifest["classes"]["K6_only"]
    records = block["graphs"]
    indices = [int(record["index"]) for record in records]
    if (
        len(records) != EXPECTED_INPUT
        or block["indices"] != indices
        or block["indices_sha256"] != EXPECTED_INPUT_SHA256
        or stable_hash(indices) != EXPECTED_INPUT_SHA256
        or verification.get("status") != "PASS"
        or verification.get("manifest", {}).get("sha256") != EXPECTED_MANIFEST_SHA256
        or verification.get("counts", {}).get("K6_only") != EXPECTED_INPUT
    ):
        raise ValueError("independent v5 boundary check failed")
    return records, indices


def singleton_data(
    instance: base.Instance, spans: Sequence[prior.Span]
) -> dict[int, tuple[int, int]]:
    local = {absolute: i for i, absolute in enumerate(instance.outside)}
    answer = {}
    for span_index, span in enumerate(spans):
        full = [option for option in span.options if option.kind == "full"]
        if len(full) != 1:
            raise AssertionError("one full option expected per side span")
        option = full[0]
        if len(option.selected) == 1 and option.rank == 0:
            vertex = option.selected[0]
            answer[span_index] = (vertex, instance.defects[local[vertex]])
        elif len(option.selected) == 1 and option.rank != 1:
            raise AssertionError("independent singleton rank is invalid")
    return answer


def impose_singleton_status(
    spans: Sequence[prior.Span],
    zero: dict[int, tuple[int, int]],
    empty_span: int | None,
) -> tuple[prior.Span, ...]:
    answer = []
    for span_index, span in enumerate(spans):
        if span_index not in zero:
            answer.append(span)
            continue
        vertex, mask = zero[span_index]
        options = (prior.Option(0, 0, (), "omit"),)
        if span_index != empty_span:
            options += (prior.Option(1, mask, (vertex,), "forced_nonempty"),)
        answer.append(prior.Span(span.label, options, True))
    return tuple(answer)


def tight_containers(
    spans: Sequence[prior.Span],
) -> tuple[tuple[int, tuple[tuple[str, prior.Option], ...]], ...]:
    """Enumerate all tight unions by 64 independent container maxima."""

    answer = []
    for container in range(1 << COORDINATES):
        rank = 0
        union = 0
        witness = []
        for span in spans:
            eligible = [
                option
                for option in span.options
                if not option.coordinates & ~container
            ]
            choice = max(
                eligible,
                key=lambda option: (
                    option.rank,
                    option.coordinates.bit_count(),
                    option.coordinates,
                    option.selected,
                ),
            )
            rank += choice.rank
            union |= choice.coordinates
            if choice.selected:
                witness.append((span.label, choice))
        if rank > union.bit_count():
            raise AssertionError("parent Hall-passing orientation has a violation")
        if rank == container.bit_count():
            # Since union is contained in container and Hall gives rank<=|union|,
            # equality with |container| forces union=container.
            if union != container:
                raise AssertionError("tight container did not equal the actual union")
            answer.append((container, tuple(witness)))
    return tuple(answer)


def side_block_data(
    adjacency: Sequence[int],
    instance: base.Instance,
    selected: int,
    prefix: str,
) -> dict[str, tuple[tuple[int, ...], tuple[int, ...]]]:
    rows, absolute, _ = prior.induced_adjacency(adjacency, instance, selected)
    answer = {}
    for number, part in enumerate(prior.connected_parts(rows)):
        answer[f"{prefix}:F{number}"] = (
            prior.induced_rows(rows, part),
            tuple(absolute[i] for i in part),
        )
    return answer


def support_obstructed(
    graph: Sequence[int],
    absolute: Sequence[int],
    defects: dict[int, int],
    removed: int,
) -> bool:
    reduced = {vertex: defects[vertex] & ~removed for vertex in absolute}
    if any(mask == 0 for mask in reduced.values()):
        return True
    for row, vertex in enumerate(absolute):
        for column in range(row):
            other = absolute[column]
            if graph[row] & (1 << column) and not (
                reduced[vertex] & reduced[other]
            ):
                return True
    seen = set()
    for mask in reduced.values():
        if mask.bit_count() == 1:
            if mask in seen:
                return True
            seen.add(mask)
    return False


def orientation_state(
    adjacency: Sequence[int],
    instance: base.Instance,
    z0: int,
    component: base.Component,
    sign_a: str,
    sign_b: str,
    chosen_empty: int | None,
    inertia: base.SympyInertiaCache,
    forcing: zf_reference.IndependentZeroForcing,
    cache: dict,
) -> bool:
    spans_a = reference.side_spans(
        adjacency, instance, component.side_a, "A", sign_a,
        inertia, forcing, cache,
    )
    spans_b = reference.side_spans(
        adjacency, instance, component.side_b, "B", sign_b,
        inertia, forcing, cache,
    )
    negative_name, negative_spans = (
        ("A", spans_a) if sign_a == "negative" else ("B", spans_b)
    )
    zero = singleton_data(instance, negative_spans)
    empty_span = None
    if chosen_empty is not None:
        empty_span = next(
            (index for index, (vertex, _) in zero.items() if vertex == chosen_empty),
            None,
        )
        if empty_span is None:
            raise ValueError("independent chosen empty vertex is ineligible")
    if negative_name == "A":
        spans_a = impose_singleton_status(spans_a, zero, empty_span)
    else:
        spans_b = impose_singleton_status(spans_b, zero, empty_span)
    empty_label = None if empty_span is None else f"{negative_name}:F{empty_span}"
    zspans = tuple(
        prior.Span(
            f"Z0:{instance.outside[local]}",
            (
                prior.Option(0, 0, (), "empty"),
                prior.Option(
                    1,
                    instance.defects[local],
                    (instance.outside[local],),
                    "Z0",
                ),
            ),
            True,
        )
        for local in base.bits(z0)
    )
    spans = spans_a + spans_b + zspans
    parent_passed, _ = reference.independent_hall(spans)
    if not parent_passed:
        return False
    blocks = side_block_data(
        adjacency, instance, component.side_a, "A"
    )
    blocks.update(side_block_data(
        adjacency, instance, component.side_b, "B"
    ))
    if empty_label is not None:
        blocks.pop(empty_label)
    defects = {
        absolute: instance.defects[local]
        for local, absolute in enumerate(instance.outside)
    }
    for label, (graph, absolute) in blocks.items():
        other = tuple(span for span in spans if span.label != label)
        for removed, _ in tight_containers(other):
            if support_obstructed(graph, absolute, defects, removed):
                return False
    return True


def component_options(
    adjacency: Sequence[int],
    instance: base.Instance,
    z0: int,
    component: base.Component,
    inertia: base.SympyInertiaCache,
    forcing: zf_reference.IndependentZeroForcing,
    cache: dict,
) -> tuple[int | None, ...]:
    choices: set[int | None] = set()
    for sign_a, sign_b in (
        ("positive", "negative"),
        ("negative", "positive"),
    ):
        negative_selected = (
            component.side_a if sign_a == "negative" else component.side_b
        )
        negative_name = "A" if sign_a == "negative" else "B"
        negative_spans = reference.side_spans(
            adjacency, instance, negative_selected, negative_name, "negative",
            inertia, forcing, cache,
        )
        possible = tuple(vertex for vertex, _ in singleton_data(instance, negative_spans).values())
        for chosen in (None, *possible):
            if orientation_state(
                adjacency, instance, z0, component,
                sign_a, sign_b, chosen, inertia, forcing, cache,
            ):
                choices.add(chosen)
    light = z0 | component.vertices
    if (
        light.bit_count() <= COORDINATES
        and base.hall_matchable(light, instance.defects)
    ):
        choices.add(None)
    ordered = tuple(sorted(item for item in choices if item is not None))
    return (None, *ordered) if None in choices else ordered


def extend_states(
    states: Iterable[tuple[int, ...]],
    options: Sequence[int | None],
    adjacency: Sequence[int],
) -> set[tuple[int, ...]]:
    answer = set()
    for state in states:
        for option in options:
            if option is None:
                answer.add(tuple(state))
                continue
            chosen = tuple(sorted((*state, int(option))))
            if len(chosen) > GLOBAL_EMPTY_BUDGET:
                continue
            if len(chosen) == 2 and adjacency[chosen[0]] & (1 << chosen[1]):
                continue
            answer.add(chosen)
    return answer


def solve_seed(adjacency: Sequence[int], instance: base.Instance) -> dict:
    inertia = base.SympyInertiaCache()
    forcing = zf_reference.IndependentZeroForcing()
    cache = {}
    rows = []
    for z0 in base.z0_subsets(instance.eligible_z0):
        zvertices = [instance.outside[local] for local in base.bits(z0)]
        if not base.hall_matchable(z0, instance.defects):
            rows.append({"Z0": zvertices, "matchable": False})
            continue
        states: set[tuple[int, ...]] = {()}
        failed_component = None
        for component in base.components(instance, z0):
            if not component.bipartite:
                continue
            options = component_options(
                adjacency, instance, z0, component, inertia, forcing, cache
            )
            states = extend_states(states, options, adjacency)
            if not states:
                failed_component = [
                    instance.outside[local]
                    for local in base.bits(component.vertices)
                ]
                break
        if states:
            return {"feasible": True, "rows": None}
        rows.append({
            "Z0": zvertices,
            "matchable": True,
            "failed_component": failed_component,
        })
    return {"feasible": False, "rows": rows}


def evaluate_record(record: dict) -> dict:
    adjacency = tuple(map(int, record["adjacency"]))
    base.validate_graph(adjacency)
    seeds_checked = 0
    for seed_mask in base.clique_masks(adjacency, COORDINATES):
        seeds_checked += 1
        instance = base.build_instance(adjacency, seed_mask)
        result = solve_seed(adjacency, instance)
        if not result["feasible"]:
            return {
                "index": int(record["index"]),
                "rejected": True,
                "seeds_checked": seeds_checked,
                "seed_mask": seed_mask,
                "seed": list(instance.seed),
                "rows": result["rows"],
            }
    return {
        "index": int(record["index"]),
        "rejected": False,
        "seeds_checked": seeds_checked,
        "seed_mask": None,
        "seed": None,
        "rows": None,
    }


def positive_control() -> dict:
    adjacency = lower_bound_18_graph()
    passed = 0
    for seed_mask in base.clique_masks(adjacency, COORDINATES):
        instance = base.build_instance(adjacency, seed_mask)
        if not solve_seed(adjacency, instance)["feasible"]:
            raise AssertionError("independent checker rejects positive 18 control")
        passed += 1
    if passed != 32:
        raise AssertionError("positive control K6 seed count changed")
    return {"passed": True, "K6_seeds": passed}


def verify(
    report_path: Path,
    certificates_path: Path,
    output: Path,
    workers: int,
) -> dict:
    assert_import_independence()
    records, indices = load_input()
    report = json.loads(report_path.read_text(encoding="utf-8"))
    archive = json.loads(certificates_path.read_text(encoding="utf-8"))
    source_path = ROOT / "d6_k6_tight_hall_support.py"
    if (
        report.get("schema") != REPORT_SCHEMA
        or report.get("status") != "COMPLETE"
        or report.get("production_source_sha256") != sha256(source_path)
        or report.get("ordered_input_indices") != indices
        or report.get("ordered_input_indices_sha256") != EXPECTED_INPUT_SHA256
        or report.get("certificate_archive", {}).get("sha256")
        != sha256(certificates_path)
        or archive.get("schema") != CERTIFICATE_SCHEMA
        or archive.get("production_source_sha256") != sha256(source_path)
    ):
        raise ValueError("production artifact boundary mismatch")
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
        rejected != report.get("rejected_indices")
        or stable_hash(rejected) != report.get("rejected_indices_sha256")
        or residue != report.get("ordered_residue_indices")
        or stable_hash(residue) != report.get("ordered_residue_indices_sha256")
        or len(rejected) != report.get("graphs_rejected")
        or len(residue) != report.get("graphs_surviving")
    ):
        raise AssertionError("independent graph decisions differ")
    report_rows = {int(row["index"]): row for row in report["graph_results"]}
    archived = {int(row["index"]): row for row in archive["rejected_graphs"]}
    for row in results:
        saved = report_rows[row["index"]]
        if (
            row["rejected"] != saved["rejected"]
            or row["seeds_checked"] != saved["seeds_checked"]
            or row["seed_mask"] != saved["first_impossible_seed_mask"]
            or row["seed"] != saved["first_impossible_seed"]
        ):
            raise AssertionError(f"graph detail differs at {row['index']}")
        if row["rejected"]:
            certificate = archived[row["index"]]
            if (
                certificate["seed_mask"] != row["seed_mask"]
                or certificate["seed"] != row["seed"]
                or len(certificate["rows"]) != len(row["rows"])
                or [item["Z0"] for item in certificate["rows"]]
                != [item["Z0"] for item in row["rows"]]
            ):
                raise AssertionError(f"archive coverage differs at {row['index']}")
    result = {
        "schema": VERIFICATION_SCHEMA,
        "status": "PASS",
        "report": {"path": report_path.name, "sha256": sha256(report_path)},
        "certificates": {
            "path": certificates_path.name,
            "sha256": sha256(certificates_path),
        },
        "checks": {
            "import_independence": True,
            "source_and_input_pins": True,
            "container_maxima_instead_of_production_union_DP": True,
            "all_graph_decisions": True,
            "certificate_Z0_coverage": True,
            "positive_18_control": True,
        },
        "graphs_recomputed": len(results),
        "graphs_rejected": len(rejected),
        "graphs_surviving": len(residue),
        "rejected_indices_sha256": stable_hash(rejected),
        "ordered_residue_indices_sha256": stable_hash(residue),
        "positive_18_control": positive_control(),
        "runtime": {
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
        default=ROOT / "d6_k6_tight_hall_support_report.json",
    )
    parser.add_argument(
        "--certificates", type=Path,
        default=ROOT / "d6_k6_tight_hall_support_certificates.json",
    )
    parser.add_argument(
        "--output", type=Path,
        default=ROOT / "d6_k6_tight_hall_support_verification.json",
    )
    parser.add_argument("--workers", type=int, default=11)
    args = parser.parse_args()
    result = verify(args.report, args.certificates, args.output, args.workers)
    print(json.dumps({
        "status": result["status"],
        "rejected": result["graphs_rejected"],
        "surviving": result["graphs_surviving"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
