#!/usr/bin/env python3
"""Run the audited interval kernel on the exact v5 dimension-six residue.

This is a thin, provenance-pinned selection wrapper.  The search and
checkpoint semantics live in :mod:`run_d6_interval_residue`; this wrapper
only reconstructs the 911 current graphs from the independently verified v5
manifest and records that manifest in the campaign configuration.

As in the underlying runner, only ``KILLED`` is a mathematical conclusion.
``ABORT``, ``UNRESOLVED``, and ``INFRA_ERROR`` reject nothing.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import random
import shlex
import subprocess
import sys
from collections import Counter
from pathlib import Path

import run_d6_interval_residue as engine


ROOT = Path(__file__).resolve().parent
EXPECTED_MANIFEST_SHA256 = (
    "164b2a12813845ef1f6d6ca5e3713ed1d2edac4be58563c1648a39ff8580c0b5"
)
EXPECTED_VERIFICATION_SHA256 = (
    "e871431b8b28fc04f0e58922ff3a6386054d15d37de5c9f9e67601ede47e4b46"
)
EXPECTED_COUNTS = {"K7": 155, "K6_only": 756}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _object_hash(value: object) -> str:
    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


def load_v5_graphs(manifest_path: Path, verification_path: Path) -> tuple[list[dict], dict]:
    manifest_path = manifest_path.resolve()
    verification_path = verification_path.resolve()
    manifest_hash = engine.sha256(manifest_path)
    verification_hash = engine.sha256(verification_path)
    require(manifest_hash == EXPECTED_MANIFEST_SHA256, "unexpected v5 manifest hash")
    require(
        verification_hash == EXPECTED_VERIFICATION_SHA256,
        "unexpected v5 verification hash",
    )
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    verification = json.loads(verification_path.read_text(encoding="utf-8"))
    require(manifest.get("status") == "COMPLETE_EXACT_FILTER_UNION", "v5 status")
    require(manifest.get("class_order") == ["K7", "K6_only"], "v5 class order")
    require(verification.get("status") == "PASS", "v5 verification status")
    require(
        verification.get("manifest", {}).get("sha256") == manifest_hash,
        "verification is not bound to v5 manifest",
    )
    require(verification.get("counts", {}).get("combined") == 911, "verification count")

    graphs: list[dict] = []
    seen: set[int] = set()
    ordered_indices: list[int] = []
    for class_name in manifest["class_order"]:
        block = manifest.get("classes", {}).get(class_name)
        require(isinstance(block, dict), f"missing class {class_name}")
        records = block.get("graphs")
        indices = block.get("indices")
        require(isinstance(records, list) and isinstance(indices, list), "class payload")
        require(len(records) == len(indices) == EXPECTED_COUNTS[class_name], "class count")
        require(indices == [int(record["index"]) for record in records], "class ordering")
        require(_object_hash(indices) == block.get("indices_sha256"), "class indices hash")
        require(_object_hash(records) == block.get("graphs_sha256"), "class graphs hash")
        population = "K7" if class_name == "K7" else "K6"
        for record in records:
            index = int(record["index"])
            require(index not in seen, "duplicate graph index")
            seen.add(index)
            adjacency = tuple(map(int, record["adjacency"]))
            engine.validate_adjacency(adjacency)
            graphs.append(
                {"index": index, "population": population, "adjacency": adjacency}
            )
            ordered_indices.append(index)

    combined = manifest.get("combined", {})
    require(len(graphs) == combined.get("count") == 911, "combined count")
    require(combined.get("cross_class_overlap") == 0, "cross-class overlap")
    require(
        _object_hash(ordered_indices) == combined.get("ordered_indices_sha256"),
        "combined ordering hash",
    )
    require(
        _object_hash(sorted(ordered_indices)) == combined.get("sorted_indices_sha256"),
        "combined sorted hash",
    )
    provenance = {
        "kind": "exact_d6_current_residue_v5",
        "manifest": {"path": str(manifest_path), "sha256": manifest_hash},
        "verification": {
            "path": str(verification_path),
            "sha256": verification_hash,
            "status": "PASS",
        },
        "counts": EXPECTED_COUNTS,
        "combined": len(graphs),
        "ordered_indices_sha256": _object_hash(ordered_indices),
    }
    return graphs, provenance


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest", type=Path, default=ROOT / "d6_current_residue_manifest_v5.json"
    )
    parser.add_argument(
        "--verification",
        type=Path,
        default=ROOT / "d6_current_residue_manifest_v5_verification.json",
    )
    parser.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 1) - 1))
    parser.add_argument("--orders", type=int, default=12)
    parser.add_argument("--cap", type=int, default=20_000)
    parser.add_argument("--slices", type=int, default=1)
    parser.add_argument("--population", choices=("all", "K7", "K6"), default="all")
    parser.add_argument("--sample", type=int)
    parser.add_argument("--sample-seed", type=int, default=600_019_005)
    parser.add_argument("--shards", type=int, default=1)
    parser.add_argument("--shard", type=int, default=0)
    parser.add_argument("--include-circle-orders", action="store_true")
    parser.add_argument("--include-bulk-order", action="store_true")
    parser.add_argument(
        "--checkpoint-dir",
        type=Path,
        default=ROOT / ".runs/d6_interval_v5_checkpoints",
    )
    parser.add_argument(
        "--output", type=Path, default=ROOT / ".runs/d6_interval_v5_report.json"
    )
    parser.add_argument(
        "--decisions", type=Path, default=ROOT / ".runs/d6_interval_v5_decisions.tsv"
    )
    parser.add_argument("--checkpoint-every", type=int, default=1)
    parser.add_argument("--progress-every", type=int, default=25)
    parser.add_argument("--retry-infra-errors", action="store_true")
    parser.add_argument("--selection-only", action="store_true")
    parser.add_argument("--outer-launch-command")
    args = parser.parse_args()
    if (
        args.workers < 1
        or args.orders < 1
        or args.cap < 1
        or args.slices < 1
        or args.checkpoint_every < 1
        or args.progress_every < 1
    ):
        parser.error("worker/search/checkpoint values must be positive")
    if args.shards < 1 or not 0 <= args.shard < args.shards:
        parser.error("require shards >= 1 and 0 <= shard < shards")
    if args.sample is not None and args.sample < 1:
        parser.error("sample must be positive")

    graphs, selection_provenance = load_v5_graphs(args.manifest, args.verification)
    if args.population != "all":
        graphs = [graph for graph in graphs if graph["population"] == args.population]
    graphs = [
        graph for position, graph in enumerate(graphs) if position % args.shards == args.shard
    ]
    if args.sample is not None and args.sample < len(graphs):
        random.Random(args.sample_seed).shuffle(graphs)
        graphs = graphs[: args.sample]
    graphs.sort(key=lambda graph: graph["index"])
    for ordinal, graph in enumerate(graphs):
        graph["ordinal"] = ordinal
    counts = Counter(graph["population"] for graph in graphs)
    selected_summary = {
        "graphs": len(graphs),
        "population_counts": dict(counts),
        "indices_sha256": engine.stable_hash([graph["index"] for graph in graphs]),
        "selection_layer": selection_provenance,
    }
    if args.selection_only:
        print(json.dumps(selected_summary, indent=2, sort_keys=True))
        return

    with engine._KERNEL_LOCK:
        engine.cdriver6._kernel()  # type: ignore[attr-defined]
    kernel_binary = Path(engine.cdriver6.HERE) / engine.cdriver6.LIBNAME
    controls = engine.run_kernel_controls(args.slices)
    try:
        compiler = subprocess.run(
            ["cc", "--version"],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        ).stdout.splitlines()[0]
    except (OSError, subprocess.SubprocessError, IndexError) as error:
        compiler = f"unavailable: {type(error).__name__}: {error}"

    search = {
        "orders": args.orders,
        "cap": args.cap,
        "slices": args.slices,
        "zero_circle_only": not args.include_circle_orders,
        "include_bulk_order": args.include_bulk_order,
    }
    configuration = {
        "schema": 3,
        "sources": {
            "wrapper": engine.sha256(Path(__file__).resolve()),
            "interval_runner": engine.sha256(ROOT / "run_d6_interval_residue.py"),
            "cdriver6.py": engine.sha256(ROOT / "cdriver6.py"),
            "ckernel6.c": engine.sha256(ROOT / "ckernel6.c"),
            "ival.py": engine.sha256(ROOT / "ival.py"),
            engine.cdriver6.LIBNAME: engine.sha256(kernel_binary),
        },
        "selection": selected_summary,
        "search": search,
        "sampling": {
            "sample": args.sample,
            "sample_seed": args.sample_seed,
            "population": args.population,
            "shards": args.shards,
            "shard": args.shard,
        },
        "trust_assumptions": {
            "basic_arithmetic": (
                "IEEE-754 binary64 basic operations and sqrt are correctly rounded; "
                "each interval endpoint is expanded with nextafter."
            ),
            "transcendentals": (
                "macOS libm cos endpoint values are assumed within 8 ulps; the "
                "kernel pads both directions by 8 nextafter steps and includes all "
                "interior extrema."
            ),
            "candidate_nonedges": "unconstrained and allowed to be unit distance",
        },
        "environment": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "cpu_count": os.cpu_count(),
            "compiler": compiler,
            "kernel_compile_command": (
                f"cc -O2 -shared -o {engine.cdriver6.LIBNAME} ckernel6.c -lm"
            ),
        },
        "git": engine.git_provenance(),
        "launch": {
            "python_argv": [sys.executable, *sys.argv],
            "python_command": shlex.join([sys.executable, *sys.argv]),
            "outer_launch_command": args.outer_launch_command,
            "outer_launch_command_user_supplied": args.outer_launch_command is not None,
        },
        "kernel_controls": controls,
    }
    print(
        f"interval v5: {len(graphs)} graphs, {args.workers} native threads, "
        f"up to {args.orders} orders, cap {args.cap}; checkpoint {args.checkpoint_dir}",
        file=sys.stderr,
        flush=True,
    )
    report, has_infra = engine.run_campaign(
        graphs,
        configuration=configuration,
        workers=args.workers,
        checkpoint_every=args.checkpoint_every,
        progress_every=args.progress_every,
        retry_infra_errors=args.retry_infra_errors,
        checkpoint_dir=args.checkpoint_dir.resolve(),
        report_path=args.output.resolve(),
        decisions_path=args.decisions.resolve(),
    )
    status_counts = report["status_counts"]
    print(
        "complete: "
        + "; ".join(f"{name} {status_counts[name]}" for name in engine.RESULT_STATUSES),
        file=sys.stderr,
        flush=True,
    )
    if has_infra:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
