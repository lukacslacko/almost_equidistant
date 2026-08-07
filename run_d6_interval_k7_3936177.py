#!/usr/bin/env python3
"""Source-bound interval replay of K7 residue graph 3936177 at cap 2,000,000.

The wrapper reconstructs the target from the independently verified v6
manifest, runs the audited interval engine with one worker, four generated
orders and 24 theta slices, and writes atomic restartable checkpoints.  Only
``KILLED`` is a mathematical conclusion; every other status is unresolved.

Production launch is intentionally impossible before this file is committed:
every tracked source/input blob must equal the working file and the tracked
worktree must be clean.  Untracked files are recorded but do not block launch.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import shlex
import subprocess
import sys
from collections import Counter
from pathlib import Path

import run_d6_interval_residue as engine
import run_d6_interval_v6 as v6


ROOT = Path(__file__).resolve().parent
TARGET_INDEX = 3_936_177
TARGET_INDICES_SHA256 = (
    "26e5f40303a16a85deb1f94cadaf3af5826c55d87cd42a01cde5e897a23e0106"
)
ORDERS = 4
CAP = 2_000_000
SLICES = 24
WORKERS = 1
EXPECTED_TRACKED_FILES = {
    "v6_wrapper": (
        "run_d6_interval_v6.py",
        "509ec8055542bb60863756d1ce64a24cccfe57ccc9f56f13925edf3290ac38be",
    ),
    "interval_runner": (
        "run_d6_interval_residue.py",
        "58fbbbc548ec3cc6b9d5b231a6790249b02a75d7bfd82f1afcd7fc77fce59404",
    ),
    "cdriver6.py": (
        "cdriver6.py",
        "8f65948c8bc95242d0fde50a8cbcb798cb0e225b78a8ce626382d00135b95173",
    ),
    "ckernel6.c": (
        "ckernel6.c",
        "383ef7a17c328869c338f64a16f2b90875ff968fe2fc20e9a5c3094c830e0874",
    ),
    "ival.py": (
        "ival.py",
        "524e41e0d0637a5352c59ec998009b968f2c9f5e7f31b7f7927e8bc36c2e6fa0",
    ),
}
EXPECTED_INPUT_FILES = {
    "v6_manifest": (
        "d6_current_residue_manifest_v6.json",
        v6.EXPECTED_MANIFEST_SHA256,
    ),
    "v6_verification": (
        "d6_current_residue_manifest_v6_verification.json",
        v6.EXPECTED_VERIFICATION_SHA256,
    ),
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def select_target(graphs: list[dict]) -> list[dict]:
    """Select exactly the frozen K7 target, independently of list position."""

    selected = [graph.copy() for graph in graphs if graph["index"] == TARGET_INDEX]
    require(len(selected) == 1, "v6 manifest does not contain the target exactly once")
    require(selected[0]["population"] == "K7", "target is not in the K7 class")
    selected[0]["ordinal"] = 0
    indices = [graph["index"] for graph in selected]
    require(engine.stable_hash(indices) == TARGET_INDICES_SHA256, "target index hash")
    return selected


def source_hashes(kernel_binary: Path) -> dict[str, str]:
    hashes = {"focused_wrapper": engine.sha256(Path(__file__).resolve())}
    for key, (filename, expected) in EXPECTED_TRACKED_FILES.items():
        actual = engine.sha256(ROOT / filename)
        require(actual == expected, f"unexpected frozen source hash: {filename}")
        hashes[key] = actual
    hashes[engine.cdriver6.LIBNAME] = engine.sha256(kernel_binary)
    return hashes


def validate_committed_launch(git: dict, sources: dict[str, str]) -> dict:
    """Require a clean tracked tree and exact committed blobs for all inputs."""

    require(git.get("available") is True, "git provenance unavailable")
    require(git.get("branch") == "codex/dimension6", "wrong launch branch")
    lines = git.get("porcelain_lines")
    require(isinstance(lines, list), "missing porcelain lines")
    require(
        all(isinstance(line, str) and line.startswith("?? ") for line in lines),
        "production launch has a modified tracked file",
    )
    porcelain = "\n".join(lines)
    require(
        hashlib.sha256(porcelain.encode("utf-8")).hexdigest()
        == git.get("porcelain_sha256"),
        "porcelain hash mismatch",
    )
    require(git.get("dirty") == bool(lines), "dirty flag mismatch")

    commit = git.get("commit")
    require(isinstance(commit, str) and len(commit) == 40, "invalid launch commit")
    committed = {
        "focused_wrapper": (Path(__file__).name, sources["focused_wrapper"]),
        **EXPECTED_TRACKED_FILES,
        **EXPECTED_INPUT_FILES,
    }
    for key, (filename, expected) in committed.items():
        if key in sources:
            require(sources[key] == expected, f"recorded source hash: {key}")
        require(engine.sha256(ROOT / filename) == expected, f"working file hash: {filename}")
        blob = subprocess.run(
            ["git", "show", f"{commit}:{filename}"],
            cwd=ROOT,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        ).stdout
        require(
            hashlib.sha256(blob).hexdigest() == expected,
            f"commit blob hash: {filename}",
        )
    return {
        "status": "PASS",
        "commit": commit,
        "branch": git["branch"],
        "tracked_files_clean": True,
        "committed_source_and_input_blobs_match": True,
        "porcelain_sha256": git["porcelain_sha256"],
    }


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
    parser.add_argument(
        "--checkpoint-dir",
        type=Path,
        default=ROOT / ".runs/d6_interval_k7_3936177_cap2000000",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / ".runs/d6_interval_k7_3936177_cap2000000.json",
    )
    parser.add_argument(
        "--decisions",
        type=Path,
        default=ROOT / ".runs/d6_interval_k7_3936177_cap2000000.tsv",
    )
    parser.add_argument("--retry-infra-errors", action="store_true")
    parser.add_argument("--selection-only", action="store_true")
    parser.add_argument("--outer-launch-command")
    args = parser.parse_args()

    all_graphs, selection_provenance = v6.load_v6_graphs(
        args.manifest, args.verification
    )
    graphs = select_target(all_graphs)
    selected_summary = {
        "graphs": 1,
        "target_index": TARGET_INDEX,
        "population_counts": dict(
            Counter(graph["population"] for graph in graphs)
        ),
        "indices_sha256": engine.stable_hash(
            [graph["index"] for graph in graphs]
        ),
        "selection_layer": selection_provenance,
    }
    if args.selection_only:
        print(json.dumps(selected_summary, indent=2, sort_keys=True))
        return
    if args.outer_launch_command is None or "caffeinate" not in args.outer_launch_command:
        parser.error("production launch must record a caffeinate outer command")

    with engine._KERNEL_LOCK:
        engine.cdriver6._kernel()  # type: ignore[attr-defined]
    kernel_binary = Path(engine.cdriver6.HERE) / engine.cdriver6.LIBNAME
    sources = source_hashes(kernel_binary)
    git = engine.git_provenance()
    source_boundary = validate_committed_launch(git, sources)
    controls = engine.run_kernel_controls(SLICES)
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

    configuration = {
        "schema": 1,
        "sources": sources,
        "source_boundary": source_boundary,
        "inputs": {
            key: {"path": str((ROOT / filename).resolve()), "sha256": expected}
            for key, (filename, expected) in EXPECTED_INPUT_FILES.items()
        },
        "selection": selected_summary,
        "selector": {
            "kind": "exact_v6_K7_corpus_index",
            "index": TARGET_INDEX,
            "population": "K7",
        },
        "search": {
            "orders": ORDERS,
            "cap": CAP,
            "slices": SLICES,
            "zero_circle_only": False,
            "include_bulk_order": False,
        },
        "trust_assumptions": v6.EXPECTED_TRUST,
        "environment": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "cpu_count": os.cpu_count(),
            "compiler": compiler,
            "kernel_compile_command": (
                f"cc -O2 -shared -o {engine.cdriver6.LIBNAME} ckernel6.c -lm"
            ),
        },
        "git": git,
        "launch": {
            "python_argv": [sys.executable, *sys.argv],
            "python_command": shlex.join([sys.executable, *sys.argv]),
            "outer_launch_command": args.outer_launch_command,
            "outer_launch_command_user_supplied": True,
        },
        "kernel_controls": controls,
    }
    print(
        f"focused interval replay: index {TARGET_INDEX}, {WORKERS} native thread, "
        f"up to {ORDERS} orders, {SLICES} slices, cap {CAP}; "
        f"checkpoint {args.checkpoint_dir}",
        file=sys.stderr,
        flush=True,
    )
    report, has_infra = engine.run_campaign(
        graphs,
        configuration=configuration,
        workers=WORKERS,
        checkpoint_every=1,
        progress_every=1,
        retry_infra_errors=args.retry_infra_errors,
        checkpoint_dir=args.checkpoint_dir.resolve(),
        report_path=args.output.resolve(),
        decisions_path=args.decisions.resolve(),
    )
    print(
        "complete: "
        + "; ".join(
            f"{name} {report['status_counts'][name]}"
            for name in engine.RESULT_STATUSES
        ),
        file=sys.stderr,
        flush=True,
    )
    if has_infra:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
