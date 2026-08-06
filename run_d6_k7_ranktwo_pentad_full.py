#!/usr/bin/env python3
"""Full exact coefficientwise-pentad scan of all 53 K7 no-near covers."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import subprocess
import sys
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from datetime import UTC, datetime
from fractions import Fraction
from pathlib import Path

import d6_k7_rank_reference as prior
import d6_k7_ranktwo_pentad as pentad
from verify_profile_d6 import lower_bound_18_graph


ROOT = Path(__file__).resolve().parent
Q = Fraction
EXPECTED_TARGETS_SHA256 = (
    "d6a1719ec12a5461c019dc01ea9f95549fefd16c71beb42991d6f2416597f313"
)


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


def utc_now() -> str:
    return datetime.now(UTC).isoformat()


def atomic_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    temporary.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    os.replace(temporary, path)


def git_provenance() -> dict:
    def git(*arguments: str) -> str:
        return subprocess.run(
            ["git", *arguments], cwd=ROOT, check=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        ).stdout.strip()

    try:
        status = git("status", "--porcelain=v1", "--untracked-files=all")
        return {
            "available": True,
            "commit": git("rev-parse", "HEAD"),
            "branch": git("branch", "--show-current"),
            "dirty": bool(status),
            "porcelain_sha256": hashlib.sha256(status.encode()).hexdigest(),
            "porcelain_lines": status.splitlines(),
        }
    except (OSError, subprocess.SubprocessError) as error:
        return {"available": False, "error": f"{type(error).__name__}: {error}"}


def initialize_worker() -> None:
    for name in (
        "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
    ):
        os.environ[name] = "1"


def serialized_polynomial(value: pentad.Polynomial) -> list[dict]:
    return [
        {
            "monomial": list(exponent),
            "coefficient": pentad.fraction_json(coefficient),
        }
        for exponent, coefficient in sorted(value.items())
    ]


def equation_summary(index: int, vertices: tuple[int, ...], equation: pentad.Polynomial) -> dict:
    if not equation:
        classification = "ZERO"
    elif all(value > 0 for value in equation.values()):
        classification = "POSITIVE"
    elif all(value < 0 for value in equation.values()):
        classification = "NEGATIVE"
    else:
        classification = "MIXED"
    encoded = serialized_polynomial(equation)
    return {
        "pentad_index": index,
        "vertices": list(vertices),
        "classification": classification,
        "terms": len(equation),
        "degree": pentad.polynomial_degree(equation),
        "polynomial_sha256": stable_hash(encoded),
    }


def make_certificate(
    system: pentad.RankTwoSystem, index: int, sign: int
) -> dict:
    signed = {
        exponent: Q(sign) * value
        for exponent, value in system.pentad_equations[index].items()
    }
    certificate = {
        "schema": 1,
        "kind": "rank_two_schur_pentad_coefficientwise_positive",
        "clique": list(system.clique),
        "remainder": list(system.remainder),
        "basis_neighbour_masks": [list(mask) for mask in system.masks],
        "pair_labels": [list(pair) for pair in system.pair_labels],
        "pentad_index": index,
        "pentad_vertices": list(system.pentad_labels[index]),
        "sign": sign,
        "polynomial_degree": pentad.polynomial_degree(signed),
        "positive_polynomial": serialized_polynomial(signed),
    }
    return certificate


def evaluate_cover(payload: tuple[int, dict]) -> dict:
    ordinal, target = payload
    started = time.perf_counter()
    graph_n = tuple(int(row) for row in target["graph_n"])
    clique = tuple(int(vertex) for vertex in target["fixed_maximum_clique"])
    clique_mask = sum(1 << vertex for vertex in clique)
    system = pentad.rank_two_system(graph_n, clique_mask)
    summaries = [
        equation_summary(index, vertices, equation)
        for index, (vertices, equation) in enumerate(zip(
            system.pentad_labels, system.pentad_equations
        ))
    ]
    one_sign = [
        item for item in summaries
        if item["classification"] in ("POSITIVE", "NEGATIVE")
    ]
    certificate = None
    if one_sign:
        index = int(one_sign[0]["pentad_index"])
        sign = 1 if one_sign[0]["classification"] == "POSITIVE" else -1
        certificate = make_certificate(system, index, sign)
        pentad.verify_coefficientwise_certificate(
            graph_n, clique_mask, certificate
        )
    return {
        "ordinal": ordinal,
        "graph_index": int(target["graph_index"]),
        "seed": target["seed"],
        "zmask": int(target["zmask"]),
        "nvertices": target["nvertices"],
        "graph_n": target["graph_n"],
        "rank_upper": int(target["rank_upper"]),
        "fixed_maximum_clique": list(clique),
        "pentads_tested": len(summaries),
        "equation_summaries": summaries,
        "classification_counts": dict(Counter(
            item["classification"] for item in summaries
        )),
        "one_sign_pentad_indices": [
            int(item["pentad_index"]) for item in one_sign
        ],
        "status": "REJECTED" if certificate is not None else "UNRESOLVED",
        "certificate": certificate,
        "elapsed_seconds": time.perf_counter() - started,
    }


def evaluate_cover_subprocess(payload: tuple[int, dict]) -> dict:
    """Run exactly one cover in a fresh process and parse its atomic record."""

    ordinal, target = payload
    encoded = json.dumps(
        [ordinal, target], sort_keys=True, separators=(",", ":")
    )
    environment = os.environ.copy()
    for name in (
        "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
    ):
        environment[name] = "1"
    result = subprocess.run(
        [sys.executable, str(Path(__file__).resolve()), "--worker-payload", encoded],
        cwd=ROOT,
        env=environment,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if result.returncode:
        raise RuntimeError(
            f"cover ordinal {ordinal}, graph {target['graph_index']} worker "
            f"exit {result.returncode}: {result.stderr.strip()}"
        )
    try:
        record = json.loads(result.stdout)
    except json.JSONDecodeError as error:
        raise RuntimeError(
            f"cover ordinal {ordinal}, graph {target['graph_index']} emitted "
            f"invalid JSON: {result.stdout[:500]!r}; stderr={result.stderr[:500]!r}"
        ) from error
    if int(record.get("ordinal", -1)) != ordinal:
        raise RuntimeError("fresh cover worker returned the wrong ordinal")
    return record


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--targets", type=Path,
        default=Path("d6_k7_ranktwo_pentad_targets.json"),
    )
    parser.add_argument(
        "--union", type=Path,
        default=Path("d6_k7_rankone_pattern_union.json"),
    )
    parser.add_argument(
        "--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1)
    )
    parser.add_argument("--batch-size", type=int, default=16)
    parser.add_argument(
        "--checkpoint", type=Path,
        default=Path(".runs/d6_k7_ranktwo_pentad_full_checkpoint.json"),
    )
    parser.add_argument("--resume", action="store_true")
    parser.add_argument(
        "--benchmark-prefix", type=int, default=0,
        help="evaluate only this many target covers and emit no artifact",
    )
    parser.add_argument(
        "--benchmark-start", type=int, default=0,
        help="zero-based target offset for benchmark-prefix mode",
    )
    parser.add_argument(
        "--output", type=Path,
        default=Path("d6_k7_ranktwo_pentad_full_report.json"),
    )
    parser.add_argument(
        "--worker-payload", help=argparse.SUPPRESS,
    )
    args = parser.parse_args()
    if args.worker_payload is not None:
        ordinal, target = json.loads(args.worker_payload)
        initialize_worker()
        print(json.dumps(
            evaluate_cover((int(ordinal), target)),
            sort_keys=True, separators=(",", ":"),
        ))
        return
    if (
        args.workers <= 0 or args.batch_size <= 0
        or args.benchmark_prefix < 0 or args.benchmark_start < 0
    ):
        parser.error("workers must be positive and benchmark bounds nonnegative")
    if args.benchmark_start and not args.benchmark_prefix:
        parser.error("benchmark-start requires benchmark-prefix")
    started = time.perf_counter()
    started_at = utc_now()

    targets_hash = sha256(args.targets)
    if targets_hash != EXPECTED_TARGETS_SHA256:
        raise ValueError("exact pentad target manifest hash mismatch")
    target_manifest = json.loads(args.targets.read_text(encoding="utf-8"))
    if (
        target_manifest.get("schema") != 1
        or target_manifest.get("kind")
        != "d6_k7_ranktwo_pentad_exact_target_manifest"
        or target_manifest.get("status") != "COMPLETE"
    ):
        raise ValueError("exact pentad target manifest schema/status mismatch")
    targets = target_manifest["targets"]
    graph_indices = list(dict.fromkeys(
        int(item["graph_index"]) for item in targets
    ))
    provenance = {
        **target_manifest["provenance"],
        "target_manifest": {
            "path": str(args.targets), "sha256": targets_hash,
        },
    }
    if len(targets) != 53 or len({item["graph_index"] for item in targets}) != 50:
        raise ValueError("exact no-near target population is not 53 covers/50 graphs")
    if any(
        len(item["graph_n"]) != 12
        or int(item["rank_upper"]) != 7
        or int(item["maximum_clique_size"]) != 5
        or int(item["zmask"]) != 0
        for item in targets
    ):
        raise ValueError("target population does not share n=12/U=7/omega=5/z=0")

    run_targets = (
        targets[
            args.benchmark_start:args.benchmark_start + args.benchmark_prefix
        ] if args.benchmark_prefix else targets
    )
    workers = min(args.workers, len(run_targets))
    # macOS ProcessPoolExecutor repeatedly terminated at pool-lifecycle
    # boundaries despite every cover succeeding alone.  A bounded thread
    # scheduler launches one fresh subprocess per cover: at most `workers`
    # covers are live, every child exits after one atomic JSON result, and an
    # infrastructure failure is attributed to an exact cover key.  Closing
    # and reopening the tiny scheduler after each 16-cover batch avoids a
    # separate macOS process-lifecycle resource accumulation seen at 53 jobs.
    source_hash = sha256(Path(__file__).resolve())
    configuration_hash = stable_hash({
        "targets_sha256": targets_hash,
        "source_sha256": source_hash,
        "workers": workers,
        "batch_size": args.batch_size,
    })
    records: list[dict] = []
    if args.resume and args.checkpoint.exists() and not args.benchmark_prefix:
        checkpoint = json.loads(args.checkpoint.read_text(encoding="utf-8"))
        if checkpoint.get("configuration_sha256") != configuration_hash:
            raise ValueError("pentad checkpoint configuration mismatch")
        records = checkpoint["records"]
        if [int(item["ordinal"]) for item in records] != list(range(len(records))):
            raise ValueError("pentad checkpoint ordinals are not a prefix")
    start_ordinal = len(records)
    enumerated = list(enumerate(run_targets))[start_ordinal:]
    for batch_start in range(0, len(enumerated), args.batch_size):
        batch = enumerated[batch_start:batch_start + args.batch_size]
        with ThreadPoolExecutor(max_workers=min(workers, len(batch))) as pool:
            records.extend(pool.map(evaluate_cover_subprocess, batch))
        if not args.benchmark_prefix:
            atomic_json(args.checkpoint, {
                "schema": 1,
                "kind": "d6_k7_ranktwo_pentad_full_checkpoint",
                "configuration_sha256": configuration_hash,
                "completed": len(records),
                "records": records,
            })
    records.sort(key=lambda item: int(item["ordinal"]))
    if args.benchmark_prefix:
        print(json.dumps({
            "benchmark_prefix": len(run_targets),
            "benchmark_start": args.benchmark_start,
            "workers": workers,
            "status_counts": dict(Counter(item["status"] for item in records)),
            "aggregate_worker_seconds": sum(
                float(item["elapsed_seconds"]) for item in records
            ),
            "serialized_record_bytes_total": sum(len(json.dumps(item)) for item in records),
            "serialized_record_bytes_max": max(
                (len(json.dumps(item)) for item in records), default=0
            ),
            "wall_seconds_including_target_replay": time.perf_counter() - started,
        }, indent=2, sort_keys=True))
        return

    union = json.loads(args.union.read_text(encoding="utf-8"))
    profiles = {
        int(item["index"]): item
        for item in union["cover_structure"]["residue_profiles"]
    }
    cover_counts = Counter(int(item["graph_index"]) for item in records)
    for index, count in cover_counts.items():
        if count != int(profiles[index]["no_saturating_without_near_clique"]):
            raise ValueError("per-graph no-near cover count differs from union profile")
    rejected_by_graph: dict[int, int] = Counter(
        int(item["graph_index"]) for item in records if item["status"] == "REJECTED"
    )
    all_no_near_rejected = sorted(
        index for index, count in cover_counts.items()
        if rejected_by_graph[index] == count
    )
    marginal_graph_rejections = sorted(
        index for index in all_no_near_rejected
        if not profiles[index]["has_saturating_cover"]
        and not profiles[index]["has_tetrad_resistant_near_clique_cover"]
    )

    positive = tuple(int(row) for row in lower_bound_18_graph())
    prior.validate_graph(positive)
    positive_k7 = sum(1 for _ in prior.clique_masks(positive, 7))
    if positive_k7:
        raise AssertionError("known realizable 18-point control unexpectedly has K7")

    status_counts = Counter(item["status"] for item in records)
    classification_counts = Counter()
    for item in records:
        classification_counts.update(item["classification_counts"])
    source_paths = [
        Path(__file__).resolve(),
        ROOT / "d6_k7_ranktwo_pentad.py",
        ROOT / "run_d6_k7_arbitrary_basis_extension_pilot.py",
        ROOT / "d6_k7_positive_polynomial_dual.py",
        ROOT / "d6_k7_rankone_tetrad_pilot.py",
        ROOT / "d6_k7_rank_reference.py",
        ROOT / "verify_profile_d6.py",
    ]
    report = {
        "schema": 1,
        "kind": "d6_k7_ranktwo_pentad_full_coefficientwise_scan",
        "status": "COMPLETE",
        "claim": (
            "Every REJECTED cover has a stored exact rank-two Schur pentad "
            "whose signed Sherman expansion is coefficientwise strictly "
            "positive. UNRESOLVED is not a realizability claim. Marginal graph "
            "rejections use only exact union profiles whose other surviving "
            "cover mechanisms are absent."
        ),
        "started_at": started_at,
        "finished_at": utc_now(),
        "configuration": {
            "workers": workers,
            "logical_cpus": os.cpu_count(),
            "command": " ".join(sys.argv),
            "labeling_policy": (
                "One increasing labeling per five-subset. The independent "
                "checker proves all 120 labelings equal this pentad up to sign."
            ),
            "locator": "none",
            "execution": (
                "bounded fresh subprocess per cover; at most workers live; "
                "scheduler recycled every batch-size covers"
            ),
            "batch_size": args.batch_size,
            "configuration_sha256": configuration_hash,
        },
        "provenance": provenance,
        "source_hashes": {path.name: sha256(path) for path in source_paths},
        "checkpoint": {
            "path": str(args.checkpoint),
            "sha256": sha256(args.checkpoint),
            "completed": len(records),
        },
        "git": git_provenance(),
        "machine": {
            "platform": platform.platform(),
            "python": sys.version,
            "implementation": platform.python_implementation(),
        },
        "positive_18_control": {
            "vertices": len(positive),
            "K7_seeds": positive_k7,
            "status": "PASS_NOT_APPLICABLE_NO_K7",
        },
        "population": {
            "graphs": len(graph_indices),
            "covers": len(targets),
            "cover_keys_sha256": stable_hash([
                [item["graph_index"], item["seed"], item["zmask"]]
                for item in targets
            ]),
            "shared_parameters": {
                "n": 12, "rank_upper": 7, "maximum_clique": 5,
                "remainder": 7, "pentads_per_cover": 21, "zmask": 0,
            },
        },
        "summary": {
            "cover_status_counts": dict(status_counts),
            "pentads_tested": sum(int(item["pentads_tested"]) for item in records),
            "pentad_classification_counts": dict(classification_counts),
            "covers_with_multiple_one_sign_pentads": sum(
                len(item["one_sign_pentad_indices"]) > 1 for item in records
            ),
            "graphs_with_all_no_near_covers_rejected": len(all_no_near_rejected),
            "all_no_near_covers_rejected_indices": all_no_near_rejected,
            "marginal_graph_rejections": len(marginal_graph_rejections),
            "marginal_graph_rejected_indices": marginal_graph_rejections,
            "aggregate_worker_seconds": sum(
                float(item["elapsed_seconds"]) for item in records
            ),
            "wall_seconds": time.perf_counter() - started,
        },
        "covers": records,
        "trust_scope": {
            "exact": (
                "Hash-pinned input selection and union; exact K7 cover replay; "
                "integer/Fraction Sherman expansion and pentad identities."
            ),
            "floating_point": "No floating-point locator is used.",
            "unresolved": "No mathematical conclusion.",
        },
    }
    args.output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "output": str(args.output),
        "sha256": sha256(args.output),
        "summary": report["summary"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
