#!/usr/bin/env python3
"""Compare individual C and independent-Python K6 Lorentz decisions.

The Python decisions come from ``d6_k6_lorentz_report.json``.  Every sampled
graph is then run as a one-record corpus through the compiled C profiler, so
agreement of aggregate totals cannot hide decisions exchanged between graphs.

Example after building the production profiler::

    cc -O3 -Wall -Wextra -Werror -pthread -o /tmp/profile_d6 profile_d6.c
    python3 crosscheck_d6_k6_lorentz_c.py /tmp/profile_d6 \
        d6_k6_lorentz_sample.json d6_k6_lorentz_report.json --workers 12
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path


C_JSON_KEY = "K6_odd_component_two_light_ray_rejections"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def c_decision(
    profiler: Path, record: dict, directory: Path
) -> tuple[int, bool, int]:
    """Run one graph and return original index, rejection, and seed witness."""

    index = record["index"]
    corpus = directory / f"graph_{index}.txt"
    corpus.write_text(
        "19 " + " ".join(map(str, record["adjacency"])) + "\n",
        encoding="ascii",
    )
    process = subprocess.run(
        [str(profiler), str(corpus), "-", "1", "1"],
        check=True,
        text=True,
        capture_output=True,
    )
    population = json.loads(process.stdout)["populations"]["all"]
    if C_JSON_KEY not in population:
        raise KeyError(
            f"C profiler JSON has no {C_JSON_KEY!r}; rebuild the Lorentz-CSP branch"
        )
    witnesses = population.get("new_rule_witnesses", {}).get(
        "K6_odd_component_two_light_ray", []
    )
    seed = int(witnesses[0]["seed"]) if witnesses else 0
    return index, bool(population[C_JSON_KEY]), seed


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("profiler", type=Path)
    parser.add_argument("sample", type=Path)
    parser.add_argument("python_report", type=Path)
    parser.add_argument("--workers", type=int, default=12)
    args = parser.parse_args()
    profiler = args.profiler.resolve()
    if not profiler.is_file():
        raise SystemExit(f"profiler does not exist: {profiler}")
    if args.workers < 1:
        parser.error("--workers must be positive")

    with args.sample.open(encoding="utf-8") as stream:
        sample = json.load(stream)
    with args.python_report.open(encoding="utf-8") as stream:
        report = json.load(stream)
    bound_sample_hash = report["inputs"]["sample"]["sha256"]
    observed_sample_hash = sha256(args.sample)
    if bound_sample_hash != observed_sample_hash:
        raise SystemExit(
            "Python report is bound to a different sample: "
            f"{bound_sample_hash} != {observed_sample_hash}"
        )
    expected = {
        item["index"]: bool(item["decision"]["rejected"])
        for item in report["graph_results"]
    }
    sample_indices = {item["index"] for item in sample["graphs"]}
    if set(expected) != sample_indices:
        raise SystemExit("Python report and sample have different graph indices")

    with tempfile.TemporaryDirectory(prefix="d6-k6-lorentz-crosscheck-") as tmp:
        directory = Path(tmp)
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            observed_rows = list(
                pool.map(
                    lambda record: c_decision(profiler, record, directory),
                    sample["graphs"],
                )
            )

    mismatches = [
        {
            "index": index,
            "python_rejected": expected[index],
            "C_rejected": rejected,
            "C_seed_witness_mask": seed,
        }
        for index, rejected, seed in observed_rows
        if rejected != expected[index]
    ]
    if mismatches:
        raise SystemExit(
            "C/Python mismatches:\n" + json.dumps(mismatches[:20], indent=2)
        )
    c_rejected = sum(rejected for _, rejected, _ in observed_rows)
    print(
        json.dumps(
            {
                "result": "PASS",
                "graphs_compared_individually": len(observed_rows),
                "rejected_by_both": c_rejected,
                "survived_by_both": len(observed_rows) - c_rejected,
                "C_JSON_key": C_JSON_KEY,
                "sample_sha256": observed_sample_hash,
                "python_report_sha256": sha256(args.python_report),
                "workers": args.workers,
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
