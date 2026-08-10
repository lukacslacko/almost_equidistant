#!/usr/bin/env python3
"""Run the audited interval kernel on the 19-graph v6 K7 residue.

This provenance-pinned wrapper reconstructs all 644 graphs from the
independently verified mixed-certificate v6 manifest, selects exactly its 19
K7 graphs, and delegates search/checkpoint mechanics to
``run_d6_interval_residue``.  Only ``KILLED`` is a mathematical conclusion;
all other result statuses remain unresolved.
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
    "c513bcc40c037af23340dce3cc91d593a101dd1d073f93feab82abd556efadc9"
)
EXPECTED_VERIFICATION_SHA256 = (
    "3d11cbd20cb59b42a8c8fd5252fb3ef3dea649576982ea8017622fad3279a7bc"
)
EXPECTED_COUNTS = {"K7": 19, "K6_only": 625}
EXPECTED_CLASS_HASHES = {
    "K7": "af25b842ca9eaa6e9721f51d0400be2436d930b72cc33c72bc366f0ac5e41290",
    "K6_only": "04d9475217e294efed0c6eac34a7fad88616c90ad7d37a2525776a713e6dd0e1",
}
EXPECTED_COMBINED_ORDER_SHA256 = (
    "907b6fe460846ab613e7d35a94dbfc6eb1eb94322a36f10b4500a8ce4fcc6ccf"
)
EXPECTED_TRUST = {
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
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def object_hash(value: object) -> str:
    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


def load_v6_graphs(
    manifest_path: Path, verification_path: Path
) -> tuple[list[dict], dict]:
    manifest_path = manifest_path.resolve()
    verification_path = verification_path.resolve()
    manifest_hash = engine.sha256(manifest_path)
    verification_hash = engine.sha256(verification_path)
    require(manifest_hash == EXPECTED_MANIFEST_SHA256, "unexpected v6 manifest hash")
    require(
        verification_hash == EXPECTED_VERIFICATION_SHA256,
        "unexpected v6 verification hash",
    )
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    verification = json.loads(verification_path.read_text(encoding="utf-8"))
    require(
        manifest.get("schema") == "d6-current-certified-residue-v6",
        "v6 manifest schema",
    )
    require(
        manifest.get("status") == "COMPLETE_MIXED_CERTIFICATE_UNION",
        "v6 manifest status",
    )
    require(manifest.get("class_order") == ["K7", "K6_only"], "v6 class order")
    require(
        verification.get("schema")
        == "d6-current-certified-residue-v6-verification-v1"
        and verification.get("status") == "PASS",
        "v6 verification status",
    )
    require(
        verification.get("manifest", {}).get("sha256") == manifest_hash,
        "verification is not bound to v6 manifest",
    )
    require(
        verification.get("counts", {}).get("combined") == 644,
        "v6 verification count",
    )
    require(
        verification.get("counts", {}).get("exact_rejections_from_v5") == 264
        and verification.get("counts", {}).get("interval_incremental_rejections")
        == 3,
        "v6 verification trust-tier counts",
    )

    accounting = manifest.get("certificate_accounting", {})
    exact = accounting.get("exact_algebra_and_graph_logic", {})
    interval = accounting.get("interval", {})
    require(
        exact.get("total_rejections") == 264
        and exact.get("floating_point_enters_rejection") is False,
        "v6 exact trust tier",
    )
    require(
        interval.get("incremental_rejections") == 3
        and interval.get("non_KILLED_statuses_used_for_rejection") is False,
        "v6 interval trust tier",
    )
    trust_assumptions = interval.get("trust_assumptions")
    require(
        trust_assumptions == EXPECTED_TRUST,
        "v6 interval trust assumptions",
    )

    graphs: list[dict] = []
    ordered_indices: list[int] = []
    seen: set[int] = set()
    for class_name in manifest["class_order"]:
        block = manifest.get("classes", {}).get(class_name)
        require(isinstance(block, dict), f"missing class {class_name}")
        records = block.get("graphs")
        indices = block.get("indices")
        require(isinstance(records, list) and isinstance(indices, list), "class payload")
        require(
            len(records) == len(indices) == EXPECTED_COUNTS[class_name],
            "class count",
        )
        require(indices == [int(record["index"]) for record in records], "class ordering")
        require(object_hash(indices) == block.get("indices_sha256"), "class indices hash")
        require(
            block.get("indices_sha256") == EXPECTED_CLASS_HASHES[class_name],
            "frozen class indices hash",
        )
        require(object_hash(records) == block.get("graphs_sha256"), "class graphs hash")
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
    require(len(graphs) == combined.get("count") == 644, "combined count")
    require(combined.get("class_counts") == EXPECTED_COUNTS, "combined class counts")
    require(combined.get("cross_class_overlap") == 0, "cross-class overlap")
    require(
        object_hash(ordered_indices) == combined.get("ordered_indices_sha256")
        == EXPECTED_COMBINED_ORDER_SHA256,
        "combined ordering hash",
    )
    provenance = {
        "kind": "mixed_certificate_d6_current_residue_v6",
        "manifest": {"path": str(manifest_path), "sha256": manifest_hash},
        "verification": {
            "path": str(verification_path),
            "sha256": verification_hash,
            "status": "PASS",
        },
        "counts": EXPECTED_COUNTS,
        "combined": 644,
        "ordered_indices_sha256": object_hash(ordered_indices),
        "certificate_accounting": {
            "exact_rejections_from_v5": 264,
            "interval_incremental_rejections": 3,
            "interval_non_KILLED_used": False,
            "interval_trust_assumptions": trust_assumptions,
        },
    }
    return graphs, provenance


def select_k7(
    graphs: list[dict],
    *,
    shards: int,
    shard: int,
    sample: int | None,
    sample_seed: int,
    indices: frozenset[int] | None = None,
) -> list[dict]:
    selected = [graph.copy() for graph in graphs if graph["population"] == "K7"]
    if indices is not None:
        available = {int(graph["index"]) for graph in selected}
        missing = sorted(indices - available)
        require(not missing, f"requested indices absent from v6 K7 class: {missing}")
        selected = [graph for graph in selected if int(graph["index"]) in indices]
    selected = [
        graph
        for position, graph in enumerate(selected)
        if position % shards == shard
    ]
    if sample is not None and sample < len(selected):
        random.Random(sample_seed).shuffle(selected)
        selected = selected[:sample]
    selected.sort(key=lambda graph: graph["index"])
    for ordinal, graph in enumerate(selected):
        graph["ordinal"] = ordinal
    return selected


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        type=Path,
        default=ROOT / "d6_current_residue_manifest_v6.json",
    )
    parser.add_argument(
        "--verification",
        type=Path,
        default=ROOT / "d6_current_residue_manifest_v6_verification.json",
    )
    parser.add_argument("--workers", type=int, default=9)
    parser.add_argument("--orders", type=int, default=4)
    parser.add_argument("--cap", type=int, default=500_000)
    parser.add_argument("--slices", type=int, default=24)
    parser.add_argument("--sample", type=int)
    parser.add_argument("--sample-seed", type=int, default=600_019_006)
    parser.add_argument(
        "--index",
        dest="indices",
        action="append",
        type=int,
        help=(
            "restrict to this v6 K7 graph index; repeat for multiple indices "
            "(selection remains sorted and provenance-bound)"
        ),
    )
    parser.add_argument("--shards", type=int, default=1)
    parser.add_argument("--shard", type=int, default=0)
    parser.add_argument(
        "--checkpoint-dir",
        type=Path,
        default=ROOT / ".runs/d6_interval_v6_k7_cap500000",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / ".runs/d6_interval_v6_k7_cap500000.json",
    )
    parser.add_argument(
        "--decisions",
        type=Path,
        default=ROOT / ".runs/d6_interval_v6_k7_cap500000.tsv",
    )
    parser.add_argument("--checkpoint-every", type=int, default=1)
    parser.add_argument("--progress-every", type=int, default=1)
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

    all_graphs, selection_provenance = load_v6_graphs(
        args.manifest, args.verification
    )
    graphs = select_k7(
        all_graphs,
        shards=args.shards,
        shard=args.shard,
        sample=args.sample,
        sample_seed=args.sample_seed,
        indices=frozenset(args.indices) if args.indices else None,
    )
    selected_summary = {
        "graphs": len(graphs),
        "population_counts": dict(Counter(graph["population"] for graph in graphs)),
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
        "zero_circle_only": False,
        "include_bulk_order": False,
    }
    configuration = {
        "schema": 4,
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
            "population": "K7",
            "shards": args.shards,
            "shard": args.shard,
            "requested_indices": sorted(set(args.indices or [])),
        },
        "trust_assumptions": EXPECTED_TRUST,
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
        f"interval v6 K7: {len(graphs)} graphs, {args.workers} native threads, "
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
        + "; ".join(
            f"{name} {status_counts[name]}" for name in engine.RESULT_STATUSES
        ),
        file=sys.stderr,
        flush=True,
    )
    if has_infra:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
