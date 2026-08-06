#!/usr/bin/env python3
"""Audit and aggregate the hardened dimension-six interval benchmarks.

This is a verifier, not a search driver.  It independently checks the three
64-graph benchmark reports and their atomic checkpoints, reconstructs the
sample from the hashed corpus/selection inputs, replays every claimed KILLED
slice and both controls with the pinned kernel, and resumes copies of the
checkpoint directories while forbidding graph recomputation.

The generated tracked summary preserves every winning order and slice record.
No ABORT or other unresolved result is promoted to a mathematical claim.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import random
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path
from typing import Sequence
from unittest.mock import patch

import cdriver6


ROOT = Path(__file__).resolve().parent
N = 19
CAPS = (1_000, 5_000, 20_000)
REPORT_TEMPLATE = ".runs/d6_interval_hardened_cap{cap}.json"
RESULT_STATUSES = ("KILLED", "ABORT", "UNRESOLVED", "INFRA_ERROR")
CORPUS_LINES = 3_971_787
CORPUS_SHA256 = "12bc7e87e6e67eb9ff2851982e206410784c6737a6397874c170cc1280b225b5"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def bytes_sha256(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def stable_hash(value: object) -> str:
    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")
    return bytes_sha256(encoded)


def resolve_recorded(path: str) -> Path:
    candidate = Path(path)
    return candidate if candidate.is_absolute() else ROOT / candidate


def load_json(path: Path) -> dict:
    value = json.loads(path.read_text(encoding="utf-8"))
    require(isinstance(value, dict), f"{path} is not a JSON object")
    return value


def validate_adjacency(adjacency: Sequence[int]) -> None:
    require(len(adjacency) == N, "wrong graph order")
    full = (1 << N) - 1
    for vertex, row in enumerate(adjacency):
        require(isinstance(row, int), "noninteger adjacency row")
        require(not (row & ~full), "adjacency bit outside graph")
        require(not (row & (1 << vertex)), "adjacency loop")
        for other in range(vertex):
            require(
                bool(row & (1 << other))
                == bool(adjacency[other] & (1 << vertex)),
                "asymmetric adjacency",
            )


def load_k7(path: Path, expected_hash: str) -> set[int]:
    require(sha256(path) == expected_hash, "K7 decision hash mismatch")
    survivors: set[int] = set()
    seen: set[int] = set()
    with gzip.open(path, "rt", encoding="ascii", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        require(reader.fieldnames is not None, "K7 decisions have no header")
        require("index" in reader.fieldnames, "K7 decisions omit index")
        require("decision" in reader.fieldnames, "unexpected K7 decision format")
        for row in reader:
            index = int(row["index"])
            require(index not in seen, "duplicate K7 decision")
            seen.add(index)
            require(row["decision"] in ("REJECTED", "SURVIVOR"), "bad K7 decision")
            if row["decision"] == "SURVIVOR":
                survivors.add(index)
    require(len(seen) == 113_136, "wrong K7 decision row count")
    require(len(survivors) == 17_764, "wrong K7 survivor count")
    return survivors


def load_k6(path: Path, expected_hash: str) -> set[int]:
    require(sha256(path) == expected_hash, "K6 report hash mismatch")
    report = load_json(path)
    results = report.get("graph_results")
    require(isinstance(results, list), "K6 report omits graph_results")
    survivors: set[int] = set()
    seen: set[int] = set()
    for result in results:
        index = int(result["index"])
        require(index not in seen, "duplicate K6 result")
        seen.add(index)
        if not bool(result.get("decision", {}).get("rejected")):
            survivors.add(index)
    require(len(results) == 1_106, "wrong K6 input row count")
    require(len(survivors) == 1_098, "wrong K6 survivor count")
    return survivors


def reconstruct_sample(configuration: dict) -> list[dict]:
    selection = configuration["selection"]
    layers = selection["selection_layers"]
    require(len(layers) == 2, "benchmark is not the two-layer base residue")
    require(layers[0]["kind"] == "K7_enhanced_rank_survivors", "wrong K7 layer")
    require(layers[1]["kind"] == "K6_support_survivors", "wrong K6 layer")
    k7 = load_k7(resolve_recorded(layers[0]["path"]), layers[0]["sha256"])
    k6 = load_k6(resolve_recorded(layers[1]["path"]), layers[1]["sha256"])
    require(not (k7 & k6), "K6 and K7 populations overlap")
    population = {index: "K7" for index in k7}
    population.update((index, "K6") for index in k6)

    corpus_path = resolve_recorded(configuration["inputs"]["corpus"])
    require(
        configuration["inputs"]["corpus_sha256"] == CORPUS_SHA256,
        "configuration has wrong corpus hash",
    )
    digest = hashlib.sha256()
    graphs: list[dict] = []
    line_count = 0
    with corpus_path.open("rb") as stream:
        for index, raw in enumerate(stream):
            digest.update(raw)
            line_count += 1
            if index not in population:
                continue
            fields = raw.split()
            require(
                len(fields) == N + 1 and fields[0] == str(N).encode("ascii"),
                f"malformed selected corpus row {index}",
            )
            adjacency = tuple(map(int, fields[1:]))
            validate_adjacency(adjacency)
            graphs.append({
                "index": index,
                "population": population[index],
                "adjacency": adjacency,
            })
    require(line_count == CORPUS_LINES, "wrong corpus line count")
    require(digest.hexdigest() == CORPUS_SHA256, "corpus hash mismatch")
    require(len(graphs) == len(population) == 18_862, "base residue coverage mismatch")

    sampling = configuration["sampling"]
    require(sampling == {
        "sample": 64,
        "sample_seed": 6003,
        "shards": 1,
        "shard": 0,
    }, "unexpected benchmark sampling parameters")
    random.Random(6003).shuffle(graphs)
    graphs = graphs[:64]
    graphs.sort(key=lambda graph: graph["index"])
    for ordinal, graph in enumerate(graphs):
        graph["ordinal"] = ordinal
    indices = [graph["index"] for graph in graphs]
    counts = Counter(graph["population"] for graph in graphs)
    require(selection["graphs"] == 64, "selection graph count mismatch")
    require(selection["indices_sha256"] == stable_hash(indices), "sample index hash mismatch")
    require(selection["population_counts"] == dict(counts), "sample population mismatch")
    return graphs


def validate_sources(configuration: dict) -> dict:
    sources = configuration["sources"]
    expected_files = {
        "runner": "run_d6_interval_residue.py",
        "cdriver6.py": "cdriver6.py",
        "ckernel6.c": "ckernel6.c",
        "ival.py": "ival.py",
        cdriver6.LIBNAME: cdriver6.LIBNAME,
    }
    for key, filename in expected_files.items():
        require(sha256(ROOT / filename) == sources[key], f"source hash mismatch: {key}")

    git = configuration["git"]
    require(git.get("available") is True, "git provenance unavailable")
    commit = git["commit"]
    for key in ("runner", "cdriver6.py", "ckernel6.c", "ival.py"):
        filename = expected_files[key]
        blob = subprocess.run(
            ["git", "show", f"{commit}:{filename}"],
            cwd=ROOT,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        ).stdout
        require(bytes_sha256(blob) == sources[key], f"commit blob mismatch: {filename}")
    porcelain = "\n".join(git["porcelain_lines"])
    require(bytes_sha256(porcelain.encode("utf-8")) == git["porcelain_sha256"],
            "recorded porcelain hash mismatch")
    require(git["dirty"] == bool(git["porcelain_lines"]), "dirty flag mismatch")
    require(
        all(line.startswith("?? ") for line in git["porcelain_lines"]),
        "benchmark launch had a modified tracked file",
    )
    return {
        "commit": commit,
        "branch": git["branch"],
        "tracked_sources_clean_at_launch": True,
        "porcelain_sha256": git["porcelain_sha256"],
    }


def ncircle(adjacency: Sequence[int], seed: Sequence[int], order: Sequence[int]) -> int:
    placed = set(seed)
    answer = 0
    for vertex in order:
        if sum(bool(adjacency[vertex] & (1 << other)) for other in placed) == 5:
            answer += 1
        placed.add(vertex)
    return answer


def expected_interval(part: int, slices: int) -> tuple[float, float]:
    return (
        2.0 * math.pi * part / slices,
        2.0 * math.pi * (part + 1) / slices,
    )


def validate_slice_records(
    records: object,
    slices: int,
    *,
    required_status: str | None = None,
    forbid_status: str | None = None,
) -> None:
    require(isinstance(records, list), "slice records are not a list")
    require(len(records) == slices, "wrong slice-record count")
    for part, record in enumerate(records):
        require(record["part"] == part, "slice part mismatch")
        lo, hi = expected_interval(part, slices)
        require(record["lo"] == lo and record["hi"] == hi, "slice endpoints mismatch")
        require(record["status"] in ("KILLED", "SURVIVORS", "ABORT"),
                "unknown kernel status")
        if required_status is not None:
            require(record["status"] == required_status, "required slice status failed")
        if forbid_status is not None:
            require(record["status"] != forbid_status, "forbidden slice status observed")
        require(isinstance(record["nodes"], int) and record["nodes"] >= 0,
                "invalid slice node count")
        require(
            isinstance(record["unresolved_cells"], int)
            and record["unresolved_cells"] >= 0,
            "invalid unresolved-cell count",
        )


def control_graph(negative: bool) -> tuple[int, ...]:
    adjacency = [0] * N
    clique_order = 8 if negative else 7
    for vertex in range(clique_order):
        for other in range(vertex):
            adjacency[vertex] |= 1 << other
            adjacency[other] |= 1 << vertex
    for vertex in range(clique_order, N):
        for other in range(5):
            adjacency[vertex] |= 1 << other
            adjacency[other] |= 1 << vertex
    validate_adjacency(adjacency)
    return tuple(adjacency)


def validate_controls(controls: dict, *, replay: bool) -> dict:
    require(controls["status"] == "PASS", "kernel controls did not pass")
    slices = controls["production_slices_exercised"]
    require(slices == 24, "controls did not exercise 24 slices")
    replayed = 0
    for name, negative in (("positive", False), ("negative", True)):
        control = controls[name]
        required = "KILLED" if negative else None
        forbidden = None if negative else "KILLED"
        validate_slice_records(
            control["records"], slices,
            required_status=required, forbid_status=forbidden,
        )
        observed = Counter(record["status"] for record in control["records"])
        require(control["status_counts"] == dict(observed), "control count mismatch")
        seed = control["seed"]
        order = control["placement_order"]
        require(sorted(seed + order) == list(range(N)), "control order is incomplete")
        if replay:
            adjacency = control_graph(negative)
            for record in control["records"]:
                got = cdriver6.decide6(
                    adjacency,
                    N,
                    seed=seed,
                    order=order,
                    th0=(record["lo"], record["hi"]),
                    max_nodes=1_000,
                )
                expected = (
                    record["status"],
                    record["nodes"],
                    record["unresolved_cells"],
                )
                require(got == expected, f"{name} control replay differs")
                replayed += 1
    return {
        "status": "PASS",
        "records_checked": 48,
        "records_replayed": replayed,
        "records_sha256": stable_hash({
            name: controls[name]["records"] for name in ("positive", "negative")
        }),
    }


def decision_text(results: Sequence[dict]) -> str:
    fields = (
        "ordinal", "index", "population", "status", "winning_order",
        "winning_circles", "orders_available", "orders_tried", "kernel_calls",
        "kernel_killed", "kernel_survivors", "kernel_aborts", "total_nodes",
        "unresolved_cells", "elapsed_seconds", "error_type", "error_message",
    )
    lines = ["\t".join(fields) + "\n"]
    for result in results:
        values: list[str] = []
        for field in fields:
            value = result.get(field)
            if field == "elapsed_seconds":
                value = f"{float(value):.9f}"
            if value is None:
                value = ""
            values.append(str(value).replace("\t", " ").replace("\n", " "))
        lines.append("\t".join(values) + "\n")
    return "".join(lines)


def validate_result_structure(
    result: dict, graph: dict, configuration: dict
) -> None:
    require(result["ordinal"] == graph["ordinal"], "result ordinal mismatch")
    require(result["index"] == graph["index"], "result index mismatch")
    require(result["population"] == graph["population"], "result population mismatch")
    require(result["status"] in RESULT_STATUSES, "bad result status")
    integer_fields = (
        "winning_order", "winning_circles", "orders_available", "orders_tried",
        "kernel_calls", "kernel_killed", "kernel_survivors", "kernel_aborts",
        "total_nodes", "unresolved_cells",
    )
    require(all(isinstance(result[name], int) for name in integer_fields),
            "noninteger result counter")
    require(all(result[name] >= 0 for name in integer_fields[2:]),
            "negative result counter")
    require(result["kernel_calls"] == sum(
        result[name] for name in ("kernel_killed", "kernel_survivors", "kernel_aborts")
    ), "kernel counters do not add")
    require(float(result["elapsed_seconds"]) >= 0, "negative elapsed time")

    search = configuration["search"]
    orders = cdriver6.gen_orders(
        graph["adjacency"], N, kmax=search["orders"]
    )
    require(result["orders_available"] == len(orders), "available order count mismatch")
    indexed = [
        (position, seed, order, ncircle(graph["adjacency"], seed, order))
        for position, (seed, order) in enumerate(orders)
    ]
    candidates = [entry for entry in indexed if search["include_bulk_order"] or entry[0] != 0]
    if search["zero_circle_only"]:
        candidates = [entry for entry in candidates if entry[3] == 0]
    if result["status"] == "KILLED":
        require(1 <= result["orders_tried"] <= len(candidates), "bad winning try count")
        position, seed, order, circles = candidates[result["orders_tried"] - 1]
        require(result["winning_order"] == position, "winning order position mismatch")
        require(result["winning_seed"] == list(seed), "winning seed mismatch")
        require(result["winning_placement_order"] == list(order), "placement order mismatch")
        require(result["winning_circles"] == circles, "winning circle count mismatch")
        validate_slice_records(
            result["winning_slice_records"], search["slices"],
            required_status="KILLED",
        )
    else:
        require(result["winning_order"] == -1, "non-KILLED has a winning order")
        require(result["winning_seed"] is None, "non-KILLED has winning seed")
        require(result["winning_placement_order"] is None,
                "non-KILLED has placement order")
        require(result["winning_slice_records"] is None,
                "non-KILLED has slice witness")
        require(result["orders_tried"] == len(candidates), "non-KILLED skipped an order")
    if result["status"] == "ABORT":
        require(result["kernel_aborts"] > 0, "ABORT has no kernel abort")
    if result["status"] == "UNRESOLVED":
        require(result["kernel_aborts"] == 0, "UNRESOLVED contains an abort")
    if result["status"] == "INFRA_ERROR":
        require(bool(result["error_type"]), "INFRA_ERROR omits exception type")


def validate_checkpoint_directory(
    report: dict, report_path: Path, graphs: list[dict]
) -> dict:
    checkpoint_dir = resolve_recorded(report["checkpoint_directory"])
    configuration = report["configuration"]
    config_hash = report["config_sha256"]
    campaign = load_json(checkpoint_dir / "campaign.json")
    require(campaign["schema"] == 1, "campaign schema mismatch")
    require(campaign["configuration"] == configuration, "campaign configuration mismatch")
    require(campaign["config_sha256"] == config_hash, "campaign config hash mismatch")

    entries: list[list[str]] = []
    for graph, result in zip(graphs, report["results"], strict=True):
        path = checkpoint_dir / "results" / (
            f"result_{graph['ordinal']:06d}_{graph['index']:07d}.json"
        )
        require(path.exists(), f"missing checkpoint {path.name}")
        value = load_json(path)
        require(value["schema"] == 1, "checkpoint schema mismatch")
        require(value["config_sha256"] == config_hash, "checkpoint config mismatch")
        require(value["search_slices"] == configuration["search"]["slices"],
                "checkpoint slice mismatch")
        require(value["result"] == result, "checkpoint/report result mismatch")
        entries.append([path.name, sha256(path)])
    require(
        report["checkpoint_result_index_sha256"] == stable_hash(entries),
        "checkpoint result-index hash mismatch",
    )

    decisions_path = resolve_recorded(report["decisions"]["path"])
    decisions_bytes = decision_text(report["results"]).encode("utf-8")
    require(decisions_path.read_bytes() == decisions_bytes, "decision TSV differs")
    require(sha256(decisions_path) == report["decisions"]["sha256"],
            "decision TSV hash mismatch")
    require(report["decisions"]["rows_excluding_header"] == len(graphs),
            "decision TSV row count mismatch")

    progress = load_json(checkpoint_dir / "progress.json")
    require(progress["config_sha256"] == config_hash, "progress config mismatch")
    require(progress["completed_graphs"] == len(graphs), "progress incomplete")
    require(progress["status_counts"] == report["status_counts"],
            "progress status mismatch")
    run_state = load_json(checkpoint_dir / "run_state.json")
    require(run_state["status"] == "COMPLETE", "campaign run state not COMPLETE")
    require(run_state["config_sha256"] == config_hash, "run-state config mismatch")
    require(run_state["report"] == str(Path(report["checkpoint_directory"]).parent /
                                      report_path.name), "run-state report path mismatch")
    require(run_state["report_sha256"] == sha256(report_path), "run-state report hash mismatch")
    require(run_state["status_counts"] == report["status_counts"],
            "run-state status mismatch")
    session_path = checkpoint_dir / "sessions" / run_state["session"]
    session = load_json(session_path)
    for key, value in report["session"].items():
        require(session[key] == value, f"session field mismatch: {key}")
    require(session["final_status"] == "COMPLETE", "session not complete")
    return {
        "directory": str(Path(report["checkpoint_directory"])),
        "campaign_sha256": sha256(checkpoint_dir / "campaign.json"),
        "result_index_sha256": stable_hash(entries),
        "progress_sha256": sha256(checkpoint_dir / "progress.json"),
        "run_state_sha256": sha256(checkpoint_dir / "run_state.json"),
        "session_sha256": sha256(session_path),
        "decision_TSV_sha256": sha256(decisions_path),
        "result_files": len(entries),
    }


def validate_report(
    report: dict,
    report_path: Path,
    graphs: list[dict],
    *,
    replay: bool,
) -> tuple[dict, dict]:
    require(report["schema"] == 2, "report schema mismatch")
    configuration = report["configuration"]
    require(stable_hash(configuration) == report["config_sha256"],
            "report configuration hash mismatch")
    require(report["total_graphs"] == report["completed_graphs"] == len(graphs),
            "report graph count mismatch")
    search = configuration["search"]
    require(search["cap"] in CAPS, "unexpected cap")
    require(search == {
        "orders": 4,
        "cap": search["cap"],
        "slices": 24,
        "zero_circle_only": False,
        "include_bulk_order": False,
    }, "unexpected search parameters")
    require(configuration["launch"]["outer_launch_command_user_supplied"] is True,
            "outer launch command not supplied")
    require("taskpolicy -a" in configuration["launch"]["outer_launch_command"],
            "outer launch command omits taskpolicy")
    require(report["session"]["config_sha256"] == report["config_sha256"],
            "report session config mismatch")
    require(report["session"]["workers"] == 11, "benchmark did not use 11 workers")
    require(len(report["results"]) == len(graphs), "result row count mismatch")

    for result, graph in zip(report["results"], graphs, strict=True):
        validate_result_structure(result, graph, configuration)
    statuses = Counter(result["status"] for result in report["results"])
    status_counts = {status: statuses[status] for status in RESULT_STATUSES}
    require(report["status_counts"] == status_counts, "status summary mismatch")
    require(report["certified_killed"] == statuses["KILLED"], "kill count mismatch")
    require(report["unresolved_total"] == len(graphs) - statuses["KILLED"],
            "unresolved count mismatch")
    population_counts = {
        population: sum(
            result["status"] == "KILLED" and result["population"] == population
            for result in report["results"]
        ) for population in ("K7", "K6")
    }
    require(report["killed_by_population"] == population_counts,
            "population kill count mismatch")
    kernel_totals = {
        "KILLED": sum(result["kernel_killed"] for result in report["results"]),
        "SURVIVORS": sum(result["kernel_survivors"] for result in report["results"]),
        "ABORT": sum(result["kernel_aborts"] for result in report["results"]),
    }
    require(report["kernel_status_totals"] == kernel_totals,
            "kernel total mismatch")

    replayed_slices = 0
    replayed_nodes = 0
    killed_witnesses: list[dict] = []
    for result, graph in zip(report["results"], graphs, strict=True):
        if result["status"] != "KILLED":
            continue
        witness = {
            key: result[key] for key in (
                "ordinal", "index", "population", "winning_order",
                "winning_circles", "winning_seed", "winning_placement_order",
                "winning_slice_records",
            )
        }
        killed_witnesses.append(witness)
        if replay:
            for record in result["winning_slice_records"]:
                got = cdriver6.decide6(
                    graph["adjacency"],
                    N,
                    seed=result["winning_seed"],
                    order=result["winning_placement_order"],
                    th0=(record["lo"], record["hi"]),
                    max_nodes=search["cap"],
                )
                expected = (
                    record["status"], record["nodes"], record["unresolved_cells"]
                )
                require(got == expected, f"winning slice replay differs: {result['index']}")
                replayed_slices += 1
                replayed_nodes += got[1]
    checkpoint = validate_checkpoint_directory(report, report_path, graphs)
    return ({
        "cap": search["cap"],
        "report": str(report_path.relative_to(ROOT)),
        "report_sha256": sha256(report_path),
        "config_sha256": report["config_sha256"],
        "checkpoint": checkpoint,
        "status_counts": status_counts,
        "certified_killed": statuses["KILLED"],
        "killed_indices": [
            result["index"] for result in report["results"]
            if result["status"] == "KILLED"
        ],
        "killed_by_population": population_counts,
        "kernel_status_totals": kernel_totals,
        "runtime": report["runtime"],
        "outer_launch_command": configuration["launch"]["outer_launch_command"],
        "winning_witnesses": killed_witnesses,
        "winning_witnesses_sha256": stable_hash(killed_witnesses),
        "witness_replay": {
            "status": "PASS" if replay else "SKIPPED",
            "slices": replayed_slices,
            "nodes": replayed_nodes,
            "exact_status_node_cell_match": replay,
        },
    }, configuration)


def verify_resume(report: dict, graphs: list[dict]) -> dict:
    """Resume a copied checkpoint set while making recomputation fatal."""

    import run_d6_interval_residue as runner

    source = resolve_recorded(report["checkpoint_directory"])
    original_tsv = resolve_recorded(report["decisions"]["path"]).read_bytes()
    with tempfile.TemporaryDirectory() as temporary:
        destination = Path(temporary) / "checkpoints"
        shutil.copytree(source, destination)
        with patch.object(
            runner,
            "analyze_graph_safe",
            side_effect=AssertionError("resume attempted graph recomputation"),
        ) as analyze:
            resumed, has_infra = runner.run_campaign(
                graphs,
                configuration=report["configuration"],
                workers=1,
                checkpoint_every=1,
                progress_every=64,
                retry_infra_errors=False,
                checkpoint_dir=destination,
                report_path=Path(temporary) / "resumed.json",
                decisions_path=Path(temporary) / "resumed.tsv",
            )
        require(analyze.call_count == 0, "resume recomputed a completed graph")
        require(not has_infra, "resume introduced infrastructure error")
        require(resumed["runtime"]["newly_completed_graphs"] == 0,
                "resume completed a graph anew")
        require((Path(temporary) / "resumed.tsv").read_bytes() == original_tsv,
                "resume changed decision TSV bytes")
        return {
            "status": "PASS",
            "newly_completed_graphs": 0,
            "graph_recomputations": 0,
            "decision_TSV_byte_identity": True,
        }


def aggregate(
    report_paths: Sequence[Path],
    *,
    replay: bool,
    resume: bool,
) -> dict:
    reports = [load_json(path) for path in report_paths]
    require([report["configuration"]["search"]["cap"] for report in reports]
            == list(CAPS), "reports are not ordered by expected caps")
    common_selection = reports[0]["configuration"]["selection"]
    common_inputs = reports[0]["configuration"]["inputs"]
    for report in reports[1:]:
        require(report["configuration"]["selection"] == common_selection,
                "reports use different samples")
        require(report["configuration"]["inputs"] == common_inputs,
                "reports use different corpus inputs")
    graphs = reconstruct_sample(reports[0]["configuration"])
    graph_summary = [
        {"ordinal": graph["ordinal"], "index": graph["index"],
         "population": graph["population"]}
        for graph in graphs
    ]

    audited: list[dict] = []
    configurations: list[dict] = []
    provenances: list[dict] = []
    for path, report in zip(report_paths, reports, strict=True):
        audit, configuration = validate_report(report, path, graphs, replay=replay)
        audited.append(audit)
        configurations.append(configuration)
        provenances.append(validate_sources(configuration))
    require(all(provenance == provenances[0] for provenance in provenances),
            "source/git provenance differs across caps")
    controls = configurations[0]["kernel_controls"]
    require(all(configuration["kernel_controls"] == controls
                for configuration in configurations),
            "control records differ across caps")
    control_audit = validate_controls(controls, replay=replay)

    killed_sets = [set(entry["killed_indices"]) for entry in audited]
    require(killed_sets[0] <= killed_sets[1] <= killed_sets[2],
            "KILLED sets are not nested by cap")
    require([len(values) for values in killed_sets] == [3, 7, 10],
            "expected 3/7/10 KILLED counts not reproduced")
    restart_audits = [
        verify_resume(report, graphs) for report in reports
    ] if resume else [{"status": "SKIPPED"} for _ in reports]

    trust = configurations[0]["trust_assumptions"]
    require(all(configuration["trust_assumptions"] == trust
                for configuration in configurations),
            "trust assumptions differ across caps")
    source_hashes = configurations[0]["sources"]
    require(all(configuration["sources"] == source_hashes
                for configuration in configurations),
            "source hashes differ across caps")
    report_artifacts = [
        [entry["report"], entry["report_sha256"]] for entry in audited
    ]
    return {
        "schema": 2,
        "description": (
            "Independent structural, provenance, checkpoint, control, and "
            "exact winning-slice replay audit of the hardened bounded "
            "dimension-six interval benchmark. No new residue search was run."
        ),
        "scope": {
            "population": (
                "Base union of 17,764 K7 enhanced-rank survivors and 1,098 "
                "K6 support survivors; this benchmark predates and does not "
                "apply later strict-H/support/bipartite residue filters."
            ),
            "sample_graphs": 64,
            "sample_seed": 6003,
            "K7": 60,
            "K6": 4,
            "sample_indices_sha256": common_selection["indices_sha256"],
            "sample": graph_summary,
        },
        "verifier": {
            "file": Path(__file__).name,
            "sha256": sha256(Path(__file__).resolve()),
            "runner_focused_tests": {
                "file": "test_d6_interval_residue.py",
                "sha256": sha256(ROOT / "test_d6_interval_residue.py"),
            },
        },
        "provenance": {
            **provenances[0],
            "sources": source_hashes,
            "inputs": {
                "corpus": common_inputs,
                "selection_layers": common_selection["selection_layers"],
            },
            "report_artifact_set_sha256": stable_hash(report_artifacts),
        },
        "common_search": {
            "workers": 11,
            "orders": 4,
            "slices": 24,
            "include_circle_orders": True,
            "include_bulk_order": False,
            "sample_seed": 6003,
        },
        "kernel_controls": controls,
        "control_audit": control_audit,
        "benchmarks": audited,
        "restart_audits": [
            {"cap": cap, **audit} for cap, audit in zip(CAPS, restart_audits)
        ],
        "checks": {
            "configuration_hashes": "PASS",
            "committed_source_provenance": "PASS",
            "corpus_and_selection_reconstruction": "PASS",
            "atomic_checkpoint_and_TSV_indexes": "PASS",
            "controls": control_audit["status"],
            "exact_winning_slice_replay": "PASS" if replay else "SKIPPED",
            "restart_without_recomputation": (
                "PASS" if resume else "SKIPPED"
            ),
            "kill_counts": "PASS: 3, 7, 10",
            "killed_sets_nested_with_cap": "PASS",
            "infrastructure_errors": 0,
        },
        "classification": {
            "KILLED": (
                "Certified for this graph: one recorded placement order was "
                "KILLED on every one of the 24 theta slices."
            ),
            "ABORT": "Node cap reached in at least one call; no claim.",
            "UNRESOLVED": "Survivor boxes or no eligible order; no claim.",
            "INFRA_ERROR": "Infrastructure failure; no claim.",
        },
        "trust_assumptions": trust,
        "production_campaign_launched": False,
    }


def run_focused_unit_tests() -> dict:
    command = [sys.executable, "-m", "unittest", "-q", "test_d6_interval_residue.py"]
    completed = subprocess.run(
        command,
        cwd=ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    require(completed.returncode == 0, "focused interval-runner unit tests failed")
    require("Ran 5 tests" in completed.stderr, "unexpected focused-test count")
    return {
        "status": "PASS",
        "tests": 5,
        "command": " ".join(command),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path, default=ROOT / "d6_interval_residue_benchmark.json"
    )
    parser.add_argument("--no-replay", action="store_true")
    parser.add_argument("--no-resume", action="store_true")
    parser.add_argument("--no-unit-tests", action="store_true")
    args = parser.parse_args()
    report_paths = [ROOT / REPORT_TEMPLATE.format(cap=cap) for cap in CAPS]
    summary = aggregate(
        report_paths, replay=not args.no_replay, resume=not args.no_resume
    )
    focused_tests = (
        {"status": "SKIPPED"}
        if args.no_unit_tests else run_focused_unit_tests()
    )
    summary["focused_unit_tests"] = focused_tests
    summary["checks"]["focused_runner_unit_tests"] = focused_tests["status"]
    rendered = json.dumps(summary, indent=2, sort_keys=True) + "\n"
    args.output.write_text(rendered, encoding="utf-8")
    print(
        "verified hardened interval benchmarks: KILLED 3/7/10; "
        f"wrote {args.output}"
    )


if __name__ == "__main__":
    main()
