#!/usr/bin/env python3
"""Stream the complete current K7 residue into a compact reference input.

This is the full-corpus companion to ``extract_d6_k7_rank_sample.py``.  It
uses the same exact predecessor flags, checks the fixed source hashes and
counts, and streams JSON so that extraction itself has bounded memory.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from extract_d6_k7_rank_sample import (
    EXACT_K7_COLUMNS,
    EXPECTED_CORPUS_SHA256,
    EXPECTED_DECISIONS_SHA256,
    EXPECTED_GRAPHS,
    EXPECTED_K7_RESIDUE,
    EXPECTED_KILL_LOG_SHA256,
)


def load_killed(path: Path) -> tuple[bytearray, str]:
    killed = bytearray(EXPECTED_GRAPHS)
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for line_number, raw in enumerate(stream, 1):
            digest.update(raw)
            fields = raw.split()
            if len(fields) != 2:
                raise ValueError(f"malformed kill-log line {line_number}")
            index = int(fields[0])
            if not 0 <= index < EXPECTED_GRAPHS:
                raise ValueError(f"bad kill-log index {index}")
            killed[index] = 1
    actual = digest.hexdigest()
    if actual != EXPECTED_KILL_LOG_SHA256:
        raise ValueError(f"unexpected kill-log SHA-256 {actual}")
    return killed, actual


def load_residue(
    path: Path, killed: bytearray
) -> tuple[bytearray, str, int]:
    residue = bytearray(EXPECTED_GRAPHS)
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        header_raw = stream.readline()
        digest.update(header_raw)
        header = header_raw.decode("ascii").rstrip("\n").split("\t")
        position = {name: i for i, name in enumerate(header)}
        required = {"index", "omega", *EXACT_K7_COLUMNS}
        if missing := required - position.keys():
            raise ValueError(f"decision header missing {sorted(missing)}")
        count = 0
        for expected_index, raw in enumerate(stream):
            digest.update(raw)
            fields = raw.decode("ascii").rstrip("\n").split("\t")
            index = int(fields[position["index"]])
            if index != expected_index:
                raise ValueError(f"decision index {index}, expected {expected_index}")
            if (
                not killed[index]
                and int(fields[position["omega"]]) == 7
                and all(
                    int(fields[position[column]]) == 0
                    for column in EXACT_K7_COLUMNS
                )
            ):
                residue[index] = 1
                count += 1
    if expected_index + 1 != EXPECTED_GRAPHS:
        raise ValueError(
            f"decisions have {expected_index + 1}, expected {EXPECTED_GRAPHS}"
        )
    actual = digest.hexdigest()
    if actual != EXPECTED_DECISIONS_SHA256:
        raise ValueError(f"unexpected decisions SHA-256 {actual}")
    if count != EXPECTED_K7_RESIDUE:
        raise ValueError(f"K7 residue is {count}, expected {EXPECTED_K7_RESIDUE}")
    return residue, actual, count


def stream_output(
    corpus: Path,
    output: Path,
    residue: bytearray,
    kill_hash: str,
    decisions_hash: str,
) -> tuple[str, int]:
    digest = hashlib.sha256()
    written = 0
    with corpus.open("rb") as source, output.open("w", encoding="utf-8") as sink:
        sink.write('{"schema":1,"description":"complete current exact K7 residue",')
        sink.write('"sources":')
        json.dump(
            {
                "corpus": corpus.name,
                "corpus_sha256": EXPECTED_CORPUS_SHA256,
                "old_kill_log_sha256": kill_hash,
                "decisions_sha256": decisions_hash,
            },
            sink,
            separators=(",", ":"),
            sort_keys=True,
        )
        sink.write(',"graphs":[')
        first = True
        for index, raw in enumerate(source):
            digest.update(raw)
            if residue[index]:
                fields = raw.split()
                if len(fields) != 20 or fields[0] != b"19":
                    raise ValueError(f"malformed corpus line {index}")
                if not first:
                    sink.write(",")
                first = False
                json.dump(
                    {
                        "index": index,
                        "stratum": "K7_current_exact_residue",
                        "adjacency": [int(field) for field in fields[1:]],
                    },
                    sink,
                    separators=(",", ":"),
                )
                written += 1
        sink.write("]}\n")
    if index + 1 != EXPECTED_GRAPHS:
        raise ValueError(f"corpus has {index + 1}, expected {EXPECTED_GRAPHS}")
    actual = digest.hexdigest()
    if actual != EXPECTED_CORPUS_SHA256:
        raise ValueError(f"unexpected corpus SHA-256 {actual}")
    if written != EXPECTED_K7_RESIDUE:
        raise ValueError(f"wrote {written}, expected {EXPECTED_K7_RESIDUE}")
    return actual, written


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--corpus", type=Path, default=Path("aeq_d6_n19.txt"))
    parser.add_argument("--kill-log", type=Path, default=Path("killed_d6_n19.log"))
    parser.add_argument(
        "--decisions",
        type=Path,
        default=Path("/tmp/d6_profile_stream_decisions.tsv"),
    )
    parser.add_argument(
        "--output", type=Path, default=Path("/tmp/d6_k7_rank_full_residue.json")
    )
    args = parser.parse_args()
    killed, kill_hash = load_killed(args.kill_log)
    residue, decisions_hash, count = load_residue(args.decisions, killed)
    corpus_hash, written = stream_output(
        args.corpus, args.output, residue, kill_hash, decisions_hash
    )
    print(
        f"wrote {args.output}: {written} graphs; corpus {corpus_hash}; "
        f"decisions {decisions_hash}; classified {count}"
    )


if __name__ == "__main__":
    main()
