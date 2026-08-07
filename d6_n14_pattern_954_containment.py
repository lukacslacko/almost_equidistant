#!/usr/bin/env python3
"""Restartable exact containment scan for the pattern-954 obstruction core.

The C kernel searches for a non-induced subgraph monomorphism from the
13-vertex required-unit-edge obstruction into each 19-vertex target.  A HIT
is accepted only after this Python layer independently checks injectivity and
every required edge.  TIMEOUT and infrastructure failures remain unresolved.

The default action only compiles the hash-pinned kernel and runs controls.
Passing ``--run`` is required to scan targets.  Fixed-size atomic chunk files
make a scan restartable; the configuration hash prevents incompatible chunks
from being reused.
"""

from __future__ import annotations

import argparse
import ctypes
import hashlib
import json
import os
import platform
import shlex
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Sequence


ROOT = Path(__file__).resolve().parent
EXPECTED_INPUT_SHA256 = (
    "e302b818747a08e07fbd137211c05886a7e290cd4b561c1e977d4e0830b953a0"
)
EXPECTED_C_SOURCE_SHA256 = (
    "038488956021093cfed039973b1f6ad3e4e2e73c8beb7bc0b94c9619ad8021f6"
)
EXPECTED_TARGET_SHA256 = (
    "9355854172c324f9d93cc4085020974fb9622a010b2e8fe19f5a954d17a7277d"
)
EXPECTED_TARGET_INDICES_SHA256 = (
    "320a68221ef597a1252d0c16cba6b30b4bef2bb6c13b73ce529e860ed431d497"
)
EXPECTED_TARGET_POPULATION = 12_839
KNOWN_HIT_INDEX = 3_862_619
PATTERN_ORDER = 13
TARGET_ORDER = 19
STATUS_NAMES = {
    0: "NO_HIT",
    1: "HIT",
    2: "TIMEOUT",
    3: "INVALID",
}
BASE_COMPILE_FLAGS = (
    "-O3", "-std=c11", "-Wall", "-Wextra", "-pedantic", "-shared"
)
DARWIN_LINK_FLAGS = (
    "-Wl,-no_uuid",
    "-Wl,-install_name,d6_n14_pattern_954_containment.dylib",
)


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


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def atomic_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    with temporary.open("w", encoding="ascii", newline="") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def validate_graph(adjacency: Sequence[int], order: int) -> tuple[int, ...]:
    if len(adjacency) != order:
        raise ValueError(f"expected graph order {order}, got {len(adjacency)}")
    full = (1 << order) - 1
    answer = tuple(adjacency)
    for vertex, row in enumerate(answer):
        if not isinstance(row, int) or isinstance(row, bool):
            raise ValueError("noninteger adjacency row")
        if row & ~full or row & (1 << vertex):
            raise ValueError("out-of-range adjacency bit or loop")
        for other in range(vertex):
            if bool(row & (1 << other)) != bool(
                answer[other] & (1 << vertex)
            ):
                raise ValueError("asymmetric adjacency")
    return answer


def verify_hit_mapping(
    pattern: Sequence[int], target: Sequence[int], mapping: Sequence[int]
) -> None:
    if len(mapping) != len(pattern):
        raise ValueError("HIT mapping has the wrong length")
    if any(not isinstance(vertex, int) or isinstance(vertex, bool)
           for vertex in mapping):
        raise ValueError("HIT mapping contains a noninteger")
    if any(vertex < 0 or vertex >= len(target) for vertex in mapping):
        raise ValueError("HIT mapping contains an out-of-range vertex")
    if len(set(mapping)) != len(mapping):
        raise ValueError("HIT mapping is not injective")
    for first in range(len(pattern)):
        neighbours = pattern[first]
        while neighbours:
            bit = neighbours & -neighbours
            second = bit.bit_length() - 1
            neighbours ^= bit
            if not (target[mapping[first]] & (1 << mapping[second])):
                raise ValueError(
                    f"HIT mapping loses pattern edge ({first},{second})"
                )


class Matcher:
    """Thin ctypes wrapper whose HIT results are independently verified."""

    def __init__(self, library: Path, pattern: Sequence[int]) -> None:
        self.pattern = validate_graph(pattern, PATTERN_ORDER)
        self._pattern_array = (ctypes.c_uint32 * PATTERN_ORDER)(*self.pattern)
        self._library = ctypes.CDLL(str(library))
        self._function = self._library.pattern954_match
        self._function.argtypes = (
            ctypes.POINTER(ctypes.c_uint32),
            ctypes.POINTER(ctypes.c_uint32),
            ctypes.c_uint64,
            ctypes.POINTER(ctypes.c_uint64),
            ctypes.POINTER(ctypes.c_int8),
        )
        self._function.restype = ctypes.c_int

    def call(self, target: Sequence[int], node_limit: int) -> dict:
        if len(target) != TARGET_ORDER:
            raise ValueError("target has the wrong order")
        if node_limit < 0 or node_limit >= 1 << 64:
            raise ValueError("node limit is out of uint64 range")
        target_array = (ctypes.c_uint32 * TARGET_ORDER)(*target)
        nodes = ctypes.c_uint64(0)
        mapping_array = (ctypes.c_int8 * PATTERN_ORDER)(*([-1] * PATTERN_ORDER))
        code = int(self._function(
            self._pattern_array,
            target_array,
            ctypes.c_uint64(node_limit),
            ctypes.byref(nodes),
            mapping_array,
        ))
        status = STATUS_NAMES.get(code, "UNKNOWN")
        result = {"status": status, "nodes": int(nodes.value)}
        if status == "HIT":
            mapping = [int(value) for value in mapping_array]
            verify_hit_mapping(self.pattern, target, mapping)
            result["mapping"] = mapping
        return result


def compiler_identity(cc: str) -> dict:
    executable = shutil.which(cc)
    if executable is None:
        raise RuntimeError(f"compiler {cc!r} was not found")
    completed = subprocess.run(
        [executable, "--version"], check=True, capture_output=True, text=True
    )
    first_line = (completed.stdout or completed.stderr).splitlines()[0]
    return {"requested": cc, "path": executable, "version": first_line}


def compile_kernel(source: Path, output: Path, cc: str) -> dict:
    if sha256(source) != EXPECTED_C_SOURCE_SHA256:
        raise ValueError("containment C source hash mismatch")
    identity = compiler_identity(cc)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + f".tmp.{os.getpid()}")
    flags = list(BASE_COMPILE_FLAGS)
    postprocess_commands: list[list[str]] = []
    if platform.system() == "Darwin":
        # ld embeds the output name and an ad-hoc signature by default.  A
        # stable install name, no UUID, and a fixed re-signing identity make
        # equivalent rebuilds byte-identical and therefore safely resumable.
        flags.extend(DARWIN_LINK_FLAGS)
    command = [identity["path"], *flags, "-o", str(temporary), str(source)]
    try:
        completed = subprocess.run(
            command, check=True, capture_output=True, text=True
        )
        if platform.system() == "Darwin":
            remove_signature = ["codesign", "--remove-signature", str(temporary)]
            sign = [
                "codesign", "--force", "--sign", "-",
                "--identifier", "org.openai.d6.pattern954",
                "--timestamp=none", str(temporary),
            ]
            subprocess.run(remove_signature, check=True, capture_output=True, text=True)
            subprocess.run(sign, check=True, capture_output=True, text=True)
            postprocess_commands.extend((remove_signature, sign))
        os.replace(temporary, output)
    finally:
        if temporary.exists():
            temporary.unlink()
    return {
        "compiler": identity,
        "flags": flags,
        "command": shlex.join(command),
        "postprocess_commands": [shlex.join(value) for value in postprocess_commands],
        "stdout": completed.stdout,
        "stderr": completed.stderr,
        "binary_sha256": sha256(output),
        "binary_bytes": output.stat().st_size,
    }


def load_inputs(input_path: Path) -> tuple[dict, tuple[int, ...], list[dict]]:
    if sha256(input_path) != EXPECTED_INPUT_SHA256:
        raise ValueError("pattern-954 input JSON hash mismatch")
    payload = json.loads(input_path.read_text(encoding="utf-8"))
    core = payload["obstruction_core"]
    pattern = validate_graph(core["adjacency"], PATTERN_ORDER)
    if core.get("nonedges_used_as_distance_constraints") is not False:
        raise ValueError("input does not state required-edge-only semantics")
    target_metadata = payload["current_K7_targets"]
    target_path = ROOT / target_metadata["transient_tsv"]
    if (
        target_metadata.get("population") != EXPECTED_TARGET_POPULATION
        or target_metadata.get("transient_tsv_sha256") != EXPECTED_TARGET_SHA256
        or target_metadata.get("indices_sha256")
        != EXPECTED_TARGET_INDICES_SHA256
    ):
        raise ValueError("target metadata differs from pinned selection")
    if sha256(target_path) != EXPECTED_TARGET_SHA256:
        raise ValueError("target TSV hash mismatch")
    rows: list[dict] = []
    with target_path.open("r", encoding="ascii") as stream:
        for line_number, line in enumerate(stream, 1):
            fields = tuple(map(int, line.split()))
            if len(fields) != TARGET_ORDER + 1:
                raise ValueError(f"malformed target TSV row {line_number}")
            rows.append({
                "index": fields[0],
                "adjacency": validate_graph(fields[1:], TARGET_ORDER),
            })
    indices = [row["index"] for row in rows]
    if len(rows) != EXPECTED_TARGET_POPULATION:
        raise ValueError("target TSV population mismatch")
    if len(set(indices)) != len(indices):
        raise ValueError("target TSV repeats an index")
    if stable_hash(indices) != EXPECTED_TARGET_INDICES_SHA256:
        raise ValueError("target index hash mismatch")
    return payload, pattern, rows


def embedded_pattern_target(pattern: Sequence[int]) -> tuple[int, ...]:
    return validate_graph(tuple(pattern) + (0,) * (TARGET_ORDER - len(pattern)), TARGET_ORDER)


def complete_target() -> tuple[int, ...]:
    full = (1 << TARGET_ORDER) - 1
    return tuple(full ^ (1 << vertex) for vertex in range(TARGET_ORDER))


def run_controls(
    matcher: Matcher, pattern: Sequence[int], rows: Sequence[dict], node_limit: int
) -> dict:
    by_index = {row["index"]: row for row in rows}
    if KNOWN_HIT_INDEX not in by_index:
        raise AssertionError("known containment control is absent from selection")
    controls = {
        "identity_embedding": matcher.call(embedded_pattern_target(pattern), node_limit),
        "complete_target": matcher.call(complete_target(), node_limit),
        "empty_target": matcher.call((0,) * TARGET_ORDER, node_limit),
        "known_corpus_hit": matcher.call(
            by_index[KNOWN_HIT_INDEX]["adjacency"], node_limit
        ),
    }
    invalid = list((0,) * TARGET_ORDER)
    invalid[0] = 1
    controls["invalid_loop"] = matcher.call(invalid, node_limit)
    expected = {
        "identity_embedding": "HIT",
        "complete_target": "HIT",
        "empty_target": "NO_HIT",
        "known_corpus_hit": "HIT",
        "invalid_loop": "INVALID",
    }
    observed = {name: row["status"] for name, row in controls.items()}
    if observed != expected:
        raise AssertionError(f"containment controls failed: {observed!r}")
    return {
        "status": "PASS",
        "known_hit_index": KNOWN_HIT_INDEX,
        "cases": controls,
    }


def evaluate_target(matcher: Matcher, row: dict, node_limit: int) -> dict:
    started = time.perf_counter()
    try:
        result = matcher.call(row["adjacency"], node_limit)
        if result["status"] in ("INVALID", "UNKNOWN"):
            raise RuntimeError(f"kernel returned {result['status']}")
    except Exception as error:  # retained as an auditable infrastructure result
        result = {
            "status": "INFRA_ERROR",
            "nodes": 0,
            "error": f"{type(error).__name__}: {error}",
        }
    result["index"] = row["index"]
    result["wall_seconds"] = time.perf_counter() - started
    return result


def validate_chunk(
    payload: dict,
    config_hash: str,
    expected_rows: Sequence[dict],
    pattern: Sequence[int],
) -> list[dict]:
    if payload.get("schema") != 1 or payload.get("kind") != "pattern954_chunk":
        raise ValueError("checkpoint chunk schema/kind mismatch")
    if payload.get("config_sha256") != config_hash:
        raise ValueError("checkpoint chunk configuration mismatch")
    results = payload.get("results")
    if not isinstance(results, list) or len(results) != len(expected_rows):
        raise ValueError("checkpoint chunk result count mismatch")
    if [result.get("index") for result in results] != [
        row["index"] for row in expected_rows
    ]:
        raise ValueError("checkpoint chunk index coverage/order mismatch")
    allowed = {"HIT", "NO_HIT", "TIMEOUT", "INFRA_ERROR"}
    for result, row in zip(results, expected_rows):
        if result.get("status") not in allowed:
            raise ValueError("checkpoint chunk contains an invalid status")
        if result["status"] == "HIT":
            verify_hit_mapping(pattern, row["adjacency"], result.get("mapping", []))
    return results


def chunks(values: Sequence[dict], size: int) -> Iterable[tuple[int, Sequence[dict]]]:
    for start in range(0, len(values), size):
        yield start, values[start:start + size]


def run_scan(
    args: argparse.Namespace,
    compilation: dict,
    controls: dict,
    pattern: Sequence[int],
    all_rows: Sequence[dict],
) -> dict:
    selected_rows = list(all_rows[:args.limit] if args.limit is not None else all_rows)
    selected_indices = [row["index"] for row in selected_rows]
    config = {
        "schema": 1,
        "kind": "d6_n14_pattern_954_containment_config",
        "required_edge_only": True,
        "input_json_sha256": EXPECTED_INPUT_SHA256,
        "c_source_sha256": EXPECTED_C_SOURCE_SHA256,
        "binary_sha256": compilation["binary_sha256"],
        "target_tsv_sha256": EXPECTED_TARGET_SHA256,
        "full_target_indices_sha256": EXPECTED_TARGET_INDICES_SHA256,
        "full_target_population": EXPECTED_TARGET_POPULATION,
        "selected_indices_sha256": stable_hash(selected_indices),
        "selected_population": len(selected_rows),
        "node_limit": args.node_limit,
        "chunk_size": args.chunk_size,
        "kernel_status_semantics": {
            "HIT": "mapping independently edge-checked; exact rejection",
            "NO_HIT": "exhaustive kernel search ended within node limit",
            "TIMEOUT": "unresolved; never a rejection",
            "INFRA_ERROR": "unresolved; never a rejection",
        },
    }
    config_hash = stable_hash(config)
    checkpoint_dir = args.checkpoint_root / config_hash
    checkpoint_dir.mkdir(parents=True, exist_ok=True)
    campaign_path = checkpoint_dir / "campaign.json"
    campaign = {
        "schema": 1,
        "kind": "d6_n14_pattern_954_containment_campaign",
        "config": config,
        "config_sha256": config_hash,
        "started_utc": utc_now(),
        "pid": os.getpid(),
        "workers": args.workers,
        "command": shlex.join(sys.argv),
        "machine": {
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python": platform.python_version(),
            "logical_cpus": os.cpu_count(),
        },
        "compilation": compilation,
        "controls": controls,
    }
    if campaign_path.exists():
        previous = json.loads(campaign_path.read_text(encoding="utf-8"))
        if previous.get("config_sha256") != config_hash:
            raise ValueError("existing campaign has an incompatible configuration")
    else:
        atomic_json(campaign_path, campaign)

    matcher = Matcher(args.library, pattern)
    aggregate: list[dict] = []
    started = time.perf_counter()
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        for start, rows in chunks(selected_rows, args.chunk_size):
            stop = start + len(rows)
            path = checkpoint_dir / f"chunk_{start:06d}_{stop:06d}.json"
            if path.exists():
                payload = json.loads(path.read_text(encoding="utf-8"))
                aggregate.extend(validate_chunk(payload, config_hash, rows, pattern))
                continue
            results = list(pool.map(
                lambda row: evaluate_target(matcher, row, args.node_limit), rows
            ))
            payload = {
                "schema": 1,
                "kind": "pattern954_chunk",
                "config_sha256": config_hash,
                "start": start,
                "stop": stop,
                "created_utc": utc_now(),
                "results": results,
            }
            atomic_json(path, payload)
            aggregate.extend(validate_chunk(payload, config_hash, rows, pattern))

    counts = {
        status: sum(result["status"] == status for result in aggregate)
        for status in ("HIT", "NO_HIT", "TIMEOUT", "INFRA_ERROR")
    }
    hits = [result for result in aggregate if result["status"] == "HIT"]
    hit_indices = [result["index"] for result in hits]
    unresolved_indices = [
        result["index"] for result in aggregate
        if result["status"] in ("TIMEOUT", "INFRA_ERROR")
    ]
    if sum(counts.values()) != len(selected_rows):
        raise AssertionError("final containment accounting mismatch")
    report = {
        "schema": 1,
        "kind": "d6_n14_pattern_954_containment_report",
        "status": (
            "SAMPLE" if len(selected_rows) != len(all_rows)
            else "COMPLETE" if not unresolved_indices
            else "COMPLETE_WITH_UNRESOLVED"
        ),
        "claim_scope": (
            "Only independently edge-checked HIT mappings are exact graph "
            "rejections. TIMEOUT and INFRA_ERROR are unresolved."
        ),
        "required_edge_only": True,
        "config": config,
        "config_sha256": config_hash,
        "checkpoint_directory": str(checkpoint_dir),
        "campaign_sha256": sha256(campaign_path),
        "processed": len(aggregate),
        "counts": counts,
        "certified_hit_indices": hit_indices,
        "certified_hit_indices_sha256": stable_hash(hit_indices),
        "certified_hits": hits,
        "unresolved_indices": unresolved_indices,
        "unresolved_indices_sha256": stable_hash(unresolved_indices),
        "total_kernel_nodes": sum(result["nodes"] for result in aggregate),
        "sum_target_wall_seconds": sum(
            result["wall_seconds"] for result in aggregate
        ),
        "orchestrator_wall_seconds": time.perf_counter() - started,
        "finished_utc": utc_now(),
        "workers": args.workers,
        "command": shlex.join(sys.argv),
        "compilation": compilation,
        "controls": controls,
    }
    atomic_json(args.report, report)
    return report


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input", type=Path, default=ROOT / "d6_n14_pattern_954_input.json"
    )
    parser.add_argument(
        "--source", type=Path, default=ROOT / "d6_n14_pattern_954_containment.c"
    )
    parser.add_argument(
        "--library",
        type=Path,
        default=ROOT / ".runs/d6_n14_pattern_954_containment.dylib",
    )
    parser.add_argument(
        "--checkpoint-root",
        type=Path,
        default=ROOT / ".runs/d6_n14_pattern_954_containment_checkpoints",
    )
    parser.add_argument(
        "--report",
        type=Path,
        default=ROOT / ".runs/d6_n14_pattern_954_containment_report.json",
    )
    parser.add_argument("--cc", default="cc")
    parser.add_argument("--node-limit", type=int, default=50_000)
    parser.add_argument("--chunk-size", type=int, default=256)
    parser.add_argument(
        "--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1)
    )
    parser.add_argument(
        "--limit", type=int, help="scan only this initial target prefix"
    )
    parser.add_argument(
        "--run", action="store_true", help="perform target scan after controls"
    )
    args = parser.parse_args()
    if args.node_limit <= 0:
        parser.error("--node-limit must be positive")
    if args.chunk_size <= 0 or args.workers <= 0:
        parser.error("--chunk-size and --workers must be positive")
    if args.limit is not None and not (0 <= args.limit <= EXPECTED_TARGET_POPULATION):
        parser.error("--limit is outside the target population")
    return args


def main() -> None:
    args = parse_args()
    _, pattern, rows = load_inputs(args.input)
    compilation = compile_kernel(args.source, args.library, args.cc)
    matcher = Matcher(args.library, pattern)
    controls = run_controls(matcher, pattern, rows, args.node_limit)
    if not args.run:
        print(json.dumps({
            "status": "PREFLIGHT_PASS",
            "targets_scanned": 0,
            "target_population": len(rows),
            "target_indices_sha256": stable_hash(
                [row["index"] for row in rows]
            ),
            "c_source_sha256": sha256(args.source),
            "binary_sha256": compilation["binary_sha256"],
            "controls": controls,
            "next_command": (
                "python3 d6_n14_pattern_954_containment.py --run"
            ),
        }, indent=2, sort_keys=True))
        return
    report = run_scan(args, compilation, controls, pattern, rows)
    print(json.dumps({
        "status": report["status"],
        "processed": report["processed"],
        "counts": report["counts"],
        "report": str(args.report),
        "config_sha256": report["config_sha256"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
