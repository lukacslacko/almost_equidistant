#!/usr/bin/env python3
"""Verify the archived independent full Python K7 decision stream.

The compact gzip is a proof-audit artifact, not a replacement for the
mathematical reference.  This checker pins all sources, validates every row
and aggregate in the Python stream, reproduces the gzip byte-for-byte, and
optionally compares every graph decision with the exhaustive C decision TSV.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import subprocess
from collections import Counter
from contextlib import ExitStack
from pathlib import Path


ROOT = Path(__file__).resolve().parent
PYTHON_FIELDS = (
    "index",
    "decision",
    "first_failing_seed",
    "seeds_checked",
    "covers_checked",
    "clique_failures",
    "degree_failures",
    "mask_failures",
    "basis_failures",
    "support_failures",
    "subspace_K_failures",
    "component_B_failures",
)
EXPECTED_SHA256 = {
    "d6_k7_rank_reference.py": (
        "e94e2fd92ad03a921b26440827b8365cb4740d9a58399379860f7106fec37db0"
    ),
    "run_d6_k7_rank_full_reference.py": (
        "66858a15f8a55a1830ca1cc3d0e2c0dd7e739375faf2cdef6eb107217b6e7c3c"
    ),
    "profile_d6_k7_rank.c": (
        "f41cfd25565d766fd53b8f0e08de1007de91edae1a29047e8f2744b308abfb92"
    ),
    "d6_k7_rank_c_full_report.json": (
        "ac3e73e728ca95624d0bf1da7b62777239289e052d08c62b1b863f4e0f3aa971"
    ),
    "d6_k7_rank_python_full_decisions.tsv.gz": (
        "d331016042c14ba412a29e42b1f2ee06ee101a7d4376067b64e99e249677f10f"
    ),
}
EXPECTED_UNCOMPRESSED_SHA256 = (
    "2ee0d1d490b97b48668e119005a4487bd89c5f7ccc9b5582181c68c57edcd5f7"
)
EXPECTED_C_DECISIONS_SHA256 = (
    "eaa4d5061b44cc86e6ab98ac6315b29a509539e5ac22ab33d1b7a03a34358391"
)
EXPECTED = {
    "graphs": 113_136,
    "rejected": 95_372,
    "survivors": 17_764,
    "first_index": 87,
    "last_index": 3_971_777,
    "seeds_checked": 124_664,
    "covers_checked": 5_078_810,
    "clique_failures": 4_150_349,
    "degree_failures": 93_695,
    "mask_failures": 741_589,
    "basis_failures": 2_055,
    "support_failures": 9,
    "subspace_K_failures": 0,
    "component_B_failures": 61_821,
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def require_equal(label: str, observed: object, expected: object) -> None:
    if observed != expected:
        raise SystemExit(f"{label}: observed {observed!r}, expected {expected!r}")


def audit_sources(archive: Path) -> dict[str, str]:
    observed: dict[str, str] = {}
    for name, expected in EXPECTED_SHA256.items():
        path = archive if name.endswith(".tsv.gz") else ROOT / name
        actual = sha256(path)
        require_equal(f"SHA-256 of {path}", actual, expected)
        observed[name] = actual
    return observed


def c_header(reader: csv.reader) -> dict[str, int]:
    header = next(reader, None)
    if header is None:
        raise SystemExit("empty C decision TSV")
    positions = {name: position for position, name in enumerate(header)}
    needed = {"index", "joint_rejected", "first_joint_seed", "internal_error"}
    missing = needed - positions.keys()
    if missing:
        raise SystemExit(f"C decision TSV lacks columns {sorted(missing)}")
    return positions


def audit_rows(archive: Path, c_decisions: Path | None) -> dict[str, object]:
    totals: Counter[str] = Counter()
    digest = hashlib.sha256()
    previous = -1
    first_index: int | None = None
    mismatches = 0
    with ExitStack() as stack:
        raw = stack.enter_context(gzip.open(archive, "rb"))
        header_raw = raw.readline()
        digest.update(header_raw)
        require_equal(
            "Python decision header",
            tuple(header_raw.decode("ascii").rstrip("\n").split("\t")),
            PYTHON_FIELDS,
        )
        if c_decisions is not None:
            c_stream = stack.enter_context(
                c_decisions.open(newline="", encoding="ascii")
            )
            c_reader = csv.reader(c_stream, delimiter="\t")
            c_positions = c_header(c_reader)
        else:
            c_reader = None
            c_positions = {}

        for position, raw_line in enumerate(raw):
            digest.update(raw_line)
            if not raw_line.endswith(b"\n"):
                raise SystemExit(f"Python decision row {position} lacks newline")
            fields = raw_line.decode("ascii").rstrip("\n").split("\t")
            require_equal(f"field count at Python row {position}", len(fields), 12)
            index = int(fields[0])
            decision = fields[1]
            first_seed = int(fields[2])
            values = [int(value) for value in fields[3:]]
            if index <= previous:
                raise SystemExit(
                    f"non-increasing Python index {index} after {previous}"
                )
            previous = index
            if first_index is None:
                first_index = index
            if decision not in {"REJECTED", "SURVIVOR"}:
                raise SystemExit(f"bad decision {decision!r} at index {index}")
            require_equal(
                f"first-seed consistency at index {index}",
                first_seed != 0,
                decision == "REJECTED",
            )
            totals["graphs"] += 1
            totals[decision.lower()] += 1
            for name, value in zip(PYTHON_FIELDS[3:], values, strict=True):
                totals[name] += value

            if c_reader is not None:
                c_row = next(c_reader, None)
                if c_row is None:
                    raise SystemExit(f"C decision TSV ends before index {index}")
                c_index = int(c_row[c_positions["index"]])
                c_decision = (
                    "REJECTED"
                    if int(c_row[c_positions["joint_rejected"]])
                    else "SURVIVOR"
                )
                c_seed = int(c_row[c_positions["first_joint_seed"]])
                internal_error = int(c_row[c_positions["internal_error"]])
                if (index, decision, first_seed) != (
                    c_index,
                    c_decision,
                    c_seed,
                ):
                    mismatches += 1
                if internal_error:
                    raise SystemExit(f"C internal error at index {c_index}")

        if c_reader is not None and next(c_reader, None) is not None:
            raise SystemExit("C decision TSV has trailing rows")

    require_equal("uncompressed decision SHA-256", digest.hexdigest(),
                  EXPECTED_UNCOMPRESSED_SHA256)
    summary: dict[str, object] = {
        "graphs": totals["graphs"],
        "rejected": totals["rejected"],
        "survivors": totals["survivor"],
        "first_index": first_index,
        "last_index": previous,
        "seeds_checked": totals["seeds_checked"],
        "covers_checked": totals["covers_checked"],
        "clique_failures": totals["clique_failures"],
        "degree_failures": totals["degree_failures"],
        "mask_failures": totals["mask_failures"],
        "basis_failures": totals["basis_failures"],
        "support_failures": totals["support_failures"],
        "subspace_K_failures": totals["subspace_K_failures"],
        "component_B_failures": totals["component_B_failures"],
    }
    for name, expected in EXPECTED.items():
        require_equal(name, summary[name], expected)
    summary["uncompressed_sha256"] = digest.hexdigest()
    summary["C_rows_compared"] = totals["graphs"] if c_decisions else 0
    summary["C_row_mismatches"] = mismatches if c_decisions else None
    if c_decisions is not None:
        require_equal("C row mismatches", mismatches, 0)
    return summary


def audit_reproducible_gzip(archive: Path) -> None:
    payload = gzip.decompress(archive.read_bytes())
    replay = subprocess.run(
        ["gzip", "-n", "-9", "-c"],
        input=payload,
        stdout=subprocess.PIPE,
        check=True,
    ).stdout
    require_equal("deterministic gzip replay", replay, archive.read_bytes())


def audit_c_report() -> None:
    report = json.loads((ROOT / "d6_k7_rank_c_full_report.json").read_text())
    require_equal("C report graphs", report["graphs"], EXPECTED["graphs"])
    require_equal(
        "C report enhanced rejections",
        report["graph_rejections"]["enhanced_joint_existential"],
        EXPECTED["rejected"],
    )
    require_equal("C report survivors", report["survivors"], EXPECTED["survivors"])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--archive",
        type=Path,
        default=ROOT / "d6_k7_rank_python_full_decisions.tsv.gz",
    )
    parser.add_argument(
        "--c-decisions",
        type=Path,
        help="optional exhaustive C TSV for all-row decision/seed comparison",
    )
    args = parser.parse_args()

    source_hashes = audit_sources(args.archive)
    audit_reproducible_gzip(args.archive)
    audit_c_report()
    if args.c_decisions is not None:
        require_equal(
            "C decision TSV SHA-256",
            sha256(args.c_decisions),
            EXPECTED_C_DECISIONS_SHA256,
        )
    summary = audit_rows(args.archive, args.c_decisions)
    print(json.dumps({
        "status": "PASS",
        "source_hashes": source_hashes,
        "gzip_reproduction": "PASS: gzip -n -9 is byte-identical",
        **summary,
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
