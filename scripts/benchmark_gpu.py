#!/usr/bin/env python3
"""Sweep GPU distance parameters for a sketch database.

Examples:
  python3 scripts/benchmark_gpu.py 1m
  python3 scripts/benchmark_gpu.py /data/large_db \
      --batch-sizes 2500,5000,10000 --gpu-counts 1:3 \
      --max-in-flight-tiles 2,4,8 --threads 16 --repeats 3 \
      --output results.csv

Integer lists accept comma-separated values and inclusive ranges, for example
``1,3:5`` expands to ``1,3,4,5``. GPU counts select that many IDs from
``--gpu-ids`` in order. Warm-ups and measurements are both recorded; use the
``phase`` column to distinguish them.
"""

import argparse
import csv
from datetime import datetime, timezone
import itertools
from pathlib import Path
import shlex
import subprocess
import time


def integer_values(value):
    values = []
    try:
        for item in value.split(","):
            parts = [int(part) for part in item.split(":")]
            if len(parts) == 1:
                values.append(parts[0])
                continue
            if len(parts) not in (2, 3):
                raise ValueError
            start, stop = parts[:2]
            step = parts[2] if len(parts) == 3 else (1 if stop >= start else -1)
            if step == 0 or (stop - start) * step < 0:
                raise ValueError
            values.extend(range(start, stop + (1 if step > 0 else -1), step))
    except ValueError as error:
        raise argparse.ArgumentTypeError(
            "expected comma-separated integers or inclusive ranges (start:stop[:step])"
        ) from error
    return values


def positive_integer_values(value):
    values = integer_values(value)
    if any(number < 1 for number in values):
        raise argparse.ArgumentTypeError("values must be positive")
    return values


def main():
    repo = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("database", help="sketch database prefix or path")
    parser.add_argument(
        "--binary", type=Path, default=repo / "target/release/sketchlib",
        help="path to the sketchlib executable",
    )
    parser.add_argument(
        "--batch-sizes", type=positive_integer_values, default="5000",
        help="sketch batch sizes to sweep; accept comma-separated values and inclusive ranges",
    )
    parser.add_argument(
        "--gpu-counts", type=positive_integer_values, default="1",
        help="number of GPUs to use for each run; each count selects that many IDs from --gpu-ids",
    )
    parser.add_argument(
        "--gpu-ids", type=integer_values, default="0,1,2,3",
        help="CUDA device IDs to select from, in order; accepts comma-separated values and inclusive ranges",
    )
    parser.add_argument(
        "--max-in-flight-tiles", type=positive_integer_values, default="4",
        help="maximum in-flight tiles per GPU to sweep; accepts comma-separated values and inclusive ranges",
    )
    parser.add_argument(
        "--threads", type=positive_integer_values, default="16",
        help="CPU thread counts to sweep; accepts comma-separated values and inclusive ranges",
    )
    parser.add_argument(
        "--kmers", type=positive_integer_values, default="31",
        help="k-mer sizes to sweep; accepts comma-separated values and inclusive ranges",
    )
    parser.add_argument(
        "--knn", type=positive_integer_values, default="50",
        help="nearest-neighbour counts to sweep; accepts comma-separated values and inclusive ranges",
    )
    parser.add_argument(
        "--distance", choices=("ani", "jaccard"), default="ani",
        help="distance measure to use",
    )
    parser.add_argument(
        "--warmups", type=int, default=1,
        help="number of unscored warm-up runs per parameter combination (also recorded in the CSV)",
    )
    parser.add_argument(
        "--repeats", type=int, default=3,
        help="number of measured runs per parameter combination",
    )
    parser.add_argument(
        "--output", type=Path, default=Path("gpu_benchmark.csv"),
        help="CSV results file; new runs are appended",
    )
    args = parser.parse_args()

    if args.warmups < 0 or args.repeats < 1:
        parser.error("--warmups must be non-negative and --repeats must be at least 1")
    if any(device < 0 for device in args.gpu_ids):
        parser.error("GPU IDs must be non-negative")
    if max(args.gpu_counts) > len(args.gpu_ids):
        parser.error("--gpu-ids must list at least as many devices as the largest --gpu-counts value")

    database = Path(args.database).expanduser().resolve()
    binary = args.binary.expanduser().resolve()
    output = args.output.expanduser()
    output.parent.mkdir(parents=True, exist_ok=True)

    fields = [
        "timestamp_utc", "database", "phase", "iteration", "gpu_count",
        "gpu_ids", "max_sketch_batch", "max_in_flight_tiles", "threads",
        "kmer", "knn", "distance", "elapsed_seconds", "exit_code", "command",
    ]
    write_header = not output.exists() or output.stat().st_size == 0

    with output.open("a", newline="") as result_file:
        writer = csv.DictWriter(result_file, fieldnames=fields)
        if write_header:
            writer.writeheader()

        combinations = itertools.product(
            args.batch_sizes, args.gpu_counts, args.max_in_flight_tiles,
            args.threads, args.kmers, args.knn,
        )
        for batch, gpu_count, in_flight, threads, kmer, knn in combinations:
            gpu_ids = args.gpu_ids[:gpu_count]
            command = [
                str(binary), "dist", "-k", str(kmer), "-v",
                "--max-sketch-batch", str(batch),
                "--max-in-flight-tiles", str(in_flight),
                "--threads", str(threads), "--knn", str(knn),
                "--gpu", ",".join(map(str, gpu_ids)),
            ]
            if args.distance == "ani":
                command.append("--ani")
            command.append(str(database))

            for phase, count in (("warmup", args.warmups), ("measure", args.repeats)):
                for iteration in range(1, count + 1):
                    print(f"{phase} {iteration}/{count}: {shlex.join(command)}", flush=True)
                    timestamp = datetime.now(timezone.utc).isoformat(timespec="seconds")
                    started = time.perf_counter()
                    result = subprocess.run(command, cwd=repo, stdout=subprocess.DEVNULL)
                    elapsed = time.perf_counter() - started
                    writer.writerow({
                        "timestamp_utc": timestamp,
                        "database": database,
                        "phase": phase,
                        "iteration": iteration,
                        "gpu_count": gpu_count,
                        "gpu_ids": ",".join(map(str, gpu_ids)),
                        "max_sketch_batch": batch,
                        "max_in_flight_tiles": in_flight,
                        "threads": threads,
                        "kmer": kmer,
                        "knn": knn,
                        "distance": args.distance,
                        "elapsed_seconds": f"{elapsed:.6f}",
                        "exit_code": result.returncode,
                        "command": shlex.join(command),
                    })
                    result_file.flush()
                    print(f"exit={result.returncode}, elapsed={elapsed:.3f}s", flush=True)

    print(f"Results appended to {output}")


if __name__ == "__main__":
    main()
