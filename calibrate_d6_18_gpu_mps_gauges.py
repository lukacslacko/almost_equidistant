#!/usr/bin/env python3
"""Unseeded all-K6-gauge calibration on the 14 standard controls."""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import run_d6_18_gpu_mps_screen as screen


ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT / "d6_18_gpu_mps_all_k6_gauge_calibration.json"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--starts-per-gauge", type=int, default=16)
    parser.add_argument("--lm-top-k", type=int, default=8)
    parser.add_argument("--steps", type=int, default=1000)
    parser.add_argument("--cpu-workers", type=int, default=2)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    producer_args = screen.parser().parse_args(
        [
            "--benchmark", "--device", "mps",
            "--total-restarts", str(args.starts_per_gauge),
            "--restarts-per-wave", str(args.starts_per_gauge),
            "--steps", str(args.steps),
            "--initialization", "seed_sphere",
            "--lm-top-k", str(args.lm_top_k),
            "--cpu-workers", str(args.cpu_workers),
            "--cpu-refine-max-nfev", "1000",
            "--lm-collision-distance", "0.05",
            "--lm-collision-weight", "100",
            "--no-seed-standard-controls",
        ]
    )
    producer_hash = screen.file_sha256(Path(screen.__file__))
    config = screen.configuration(producer_args, producer_hash)
    config_hash = screen.stable_hash(config)
    corpus = json.loads((ROOT / screen.CORPUS_NAME).read_text(encoding="utf-8"))

    variants = []
    variant_metadata = []
    for index, raw in enumerate(corpus["unique_deletions"]):
        if not raw["standard18_compatible"]:
            continue
        rows = screen.validate_rows(raw["adjacency"], raw["class_id"])
        seeds = screen.enumerate_cliques(rows, 6)
        for ordinal, seed in enumerate(seeds):
            variants.append(screen.reorder_record_with_seed(index, raw, seed))
            variant_metadata.append(
                {"index": index, "gauge_ordinal": ordinal, "seed": list(seed)}
            )

    executor = screen.make_executor(config)
    started = time.perf_counter()
    try:
        result = screen.process_batch(variants, config, executor)
    finally:
        if executor is not None:
            executor.shutdown(wait=True)
    elapsed = time.perf_counter() - started

    controls = []
    for index in sorted({item["index"] for item in variant_metadata}):
        selected = [
            (metadata, row)
            for metadata, row in zip(variant_metadata, result["rows"])
            if metadata["index"] == index
        ]
        recovered = [
            (metadata, row) for metadata, row in selected if row["lm_candidate"]
        ]
        best = min(selected, key=lambda item: item[1]["best_distinct_lm_edge_rms"])
        controls.append(
            {
                "index": index,
                "gauge_count": len(selected),
                "optimizer_starts": len(selected) * config["total_restarts"],
                "cpu_lm_attempts": sum(row["cpu_lm_attempts"] for _, row in selected),
                "recovered": bool(recovered),
                "recovered_gauge_ordinals": [
                    metadata["gauge_ordinal"] for metadata, _ in recovered
                ],
                "best_lm_edge_rms": best[1]["best_distinct_lm_edge_rms"],
                "best_lm_minimum_distance": best[1]["best_distinct_lm_minimum_distance"],
                "best_gauge_ordinal": best[0]["gauge_ordinal"],
                "best_gauge_seed": best[0]["seed"],
            }
        )

    report = {
        "schema": "d6-18-gpu-mps-all-k6-gauge-calibration-v1",
        "proof_status": "HEURISTIC_ONLY_NO_REJECTIONS",
        "producer_source_sha256": producer_hash,
        "calibration_source_sha256": screen.file_sha256(Path(__file__)),
        "config": config,
        "config_sha256": config_hash,
        "wall_seconds": elapsed,
        "variant_count": len(variants),
        "controls": controls,
        "recovered_controls": sum(item["recovered"] for item in controls),
        "total_controls": len(controls),
        "wave_profile": result["wave_profiles"],
        "lm_profile": result["lm_profile"],
        "mathematical_rejections": 0,
        "mathematical_realizability_conclusions": 0,
    }
    screen.atomic_json(args.output, report)
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
