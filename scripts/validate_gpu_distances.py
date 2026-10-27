#!/usr/bin/env python3
"""Sequential CUDA smoke/CPU comparisons, run remotely as a non-root user.

Build first with CUDA_HOME=/usr/local/cuda-13.4 cargo oxide build, then run:
  python3 scripts/validate_gpu_distances.py --dataset 1k_subset --gpu 0

Outputs, fixtures, command logs and report.json stay in the printed temporary
directory. No Cargo test harness, progress/quiet checks or precision-boundary
checks are used. Remote commands require the user's approval before launch.
"""

import argparse
import itertools
import json
import math
import os
from pathlib import Path
import shlex
import subprocess
import tempfile
import time
import traceback


TOLERANCE = 1e-6
KNN = 10


def run(command, repo, log, expected_failure=False):
    print("RUN " + shlex.join(map(str, command)), flush=True)
    start = time.monotonic()
    with log.open("w") as stream:
        result = subprocess.run(
            list(map(str, command)), cwd=repo, stdout=stream,
            stderr=subprocess.STDOUT, timeout=180,
        )
    elapsed = time.monotonic() - start
    if expected_failure:
        contents = log.read_text()
        if result.returncode <= 0 or "No space left on device" not in contents:
            raise RuntimeError(
                f"Expected a writer-error exit, got {result.returncode}; log: {log}\n"
                + contents[-6000:]
            )
        print(f"PASS writer failure ({elapsed:.2f}s); log: {log}", flush=True)
        return {"exit_code": result.returncode, "elapsed_seconds": elapsed}
    if result.returncode:
        raise RuntimeError(
            f"Command exited {result.returncode}; log: {log}\n"
            + log.read_text()[-6000:]
        )
    print(f"OK ({elapsed:.2f}s); log: {log}", flush=True)


def read_output(path, dense=False, names=None, ani=False):
    rows = {}
    pairs = {}
    with path.open() as stream:
        for number, line in enumerate(stream, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 3:
                raise ValueError(f"{path}:{number}: expected three fields")
            sample, neighbour, raw = fields
            value = float(raw)
            if not sample or not neighbour or sample == neighbour:
                raise ValueError(f"{path}:{number}: invalid name or self-pair")
            if not math.isfinite(value) or not -TOLERANCE <= value <= 1 + TOLERANCE:
                raise ValueError(f"{path}:{number}: invalid distance {value}")
            if names is not None and not {sample, neighbour} <= names:
                raise ValueError(f"{path}:{number}: unknown sample")
            if dense:
                pair = tuple(sorted((sample, neighbour)))
                if pair in pairs:
                    raise ValueError(f"{path}:{number}: duplicate unordered pair")
                pairs[pair] = value
            else:
                row = rows.setdefault(sample, {})
                if neighbour in row:
                    raise ValueError(f"{path}:{number}: duplicate neighbour")
                row[neighbour] = value
    if dense:
        expected = set(itertools.combinations(sorted(names), 2))
        if set(pairs) != expected:
            raise ValueError(f"{path}: incomplete dense pair coverage ({len(pairs)} records)")
        return pairs
    if not rows or (names is not None and set(rows) != names):
        raise ValueError(f"{path}: wrong sample rows")
    sample_names = set(rows)
    for sample, row in rows.items():
        if len(row) != KNN or not set(row) <= sample_names:
            raise ValueError(f"{path}: invalid neighbour count or names for {sample}")
        scores = [-value if ani else value for value in row.values()]
        if any(a > b + TOLERANCE for a, b in zip(scores, scores[1:])):
            raise ValueError(f"{path}: incorrectly sorted row {sample}")
    return rows


def compare_dense(reference, candidate):
    error = max(abs(reference[pair] - candidate[pair]) for pair in reference)
    if error > TOLERANCE:
        raise ValueError(f"Dense CPU/GPU maximum absolute error {error} > {TOLERANCE}")
    return {"pairs": len(reference), "maximum_absolute_error": error}


def compare_knn(reference, candidate, oracle, ani):
    error = 0.0
    substitutions = 0
    for sample, cpu_row in reference.items():
        gpu_row = candidate[sample]
        cutoff = min(cpu_row.values()) if ani else max(cpu_row.values())
        for neighbour in set(cpu_row) & set(gpu_row):
            error = max(error, abs(cpu_row[neighbour] - gpu_row[neighbour]))
        for neighbour in set(cpu_row) - set(gpu_row):
            if abs(cpu_row[neighbour] - cutoff) > TOLERANCE:
                raise ValueError(f"{sample}: omitted CPU neighbour outside cutoff tie")
        for neighbour in set(gpu_row) - set(cpu_row):
            actual = oracle[tuple(sorted((sample, neighbour)))]
            if abs(actual - cutoff) > TOLERANCE:
                raise ValueError(f"{sample}: alternative neighbour outside CPU cutoff tie")
            error = max(error, abs(actual - gpu_row[neighbour]))
            substitutions += 1
    if error > TOLERANCE:
        raise ValueError(f"kNN CPU/GPU maximum absolute error {error} > {TOLERANCE}")
    return {
        "samples": len(reference), "neighbours_per_sample": KNN,
        "maximum_absolute_error": error,
        "verified_cutoff_tie_substitutions": substitutions,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path.cwd())
    parser.add_argument("--dataset", default="1k_subset")
    parser.add_argument("--gpu", type=int, default=0)
    args = parser.parse_args()
    if os.geteuid() == 0:
        raise RuntimeError("Run this validation as jlees, never root")
    repo = args.repo.resolve()
    binary = repo / "target/release/sketchlib"
    work = Path(tempfile.mkdtemp(prefix="sketchlib-cuda-async-"))
    print(f"Validation outputs and logs: {work}", flush=True)
    report = {"status": "running", "repo": str(repo), "results": {}}

    def distance_run(label, subset=None, ani=False, sparse=False, mixed=False,
                     gpu=False, batch=None, threads=1, expected_failure=False):
        output = Path("/dev/full") if expected_failure else work / f"{label}.tsv"
        command = [binary, "dist", "-o", output, "-k", "21", "-v",
                   "--threads", str(threads)]
        if sparse:
            command += ["--knn", str(KNN)]
        if subset is not None:
            command += ["--subset", subset]
        if ani:
            command += ["--ani"]
        if mixed:
            command += ["--ref-completeness-file", work / "completeness.tsv"]
        if gpu:
            command += ["--gpu", str(args.gpu)]
        if batch is not None:
            command += ["--max-sketch-batch", str(batch)]
        command += [args.dataset]
        result = run(command, repo, work / f"{label}.log", expected_failure)
        if expected_failure:
            return result
        return output

    try:
        smoke_path = distance_run("smoke_gpu", sparse=True, gpu=True)
        smoke = read_output(smoke_path)
        all_names = sorted(smoke)
        if len(all_names) < 256:
            raise ValueError("Smoke dataset must have at least 256 samples")
        report["smoke_samples"] = len(all_names)
        # Spread the subset through the source collection and reverse its order,
        # exercising logical sample indices independently of metadata order.
        step = max(1, len(all_names) // 256)
        selected = all_names[::step][:256][::-1]
        names = set(selected)
        subset = work / "subset.txt"
        subset.write_text("\n".join(selected) + "\n")
        completeness = {
            name: (0.95, 0.85, 0.50, 1.0)[index % 4]
            for index, name in enumerate(selected)
        }
        (work / "completeness.tsv").write_text(
            "".join(f"{name}\t{completeness[name]}\n" for name in selected)
        )
        uncorrected = {}
        for mixed in (False, True):
            for ani in (False, True):
                tag = ("mixed" if mixed else "none") + ("_ani" if ani else "_jaccard")
                cpu_dense = read_output(
                    distance_run(f"{tag}_cpu_dense", subset, ani=ani, mixed=mixed),
                    dense=True, names=names, ani=ani,
                )
                cpu_knn = read_output(
                    distance_run(f"{tag}_cpu_knn", subset, ani=ani, sparse=True, mixed=mixed),
                    names=names, ani=ani,
                )
                if not mixed:
                    uncorrected[ani] = cpu_dense
                else:
                    changed = [
                        pair for pair in cpu_dense
                        if abs(cpu_dense[pair] - uncorrected[ani][pair]) > TOLERANCE
                    ]
                    if not changed:
                        raise ValueError("Mixed completeness did not change any eligible pair")
                    if any(completeness[a] * completeness[b] < 0.64 for a, b in changed):
                        raise ValueError("Correction changed a pair below the cutoff")
                    report[f"{tag}_changed_pairs"] = len(changed)
                # Four blocks produce ten tiles, exercising window refills and
                # both the smallest CLI pool and concurrent CPU consumers.
                for batch, threads, layout in (
                    (None, 1, "single_t1"),
                    (128, 1, "two_blocks_t1"),
                    (64, 1, "four_blocks_t1"),
                    (64, 4, "four_blocks_t4"),
                ):
                    for sparse in (False, True):
                        label = f"{tag}_gpu_{'knn' if sparse else 'dense'}_{layout}"
                        candidate = read_output(
                            distance_run(label, subset, ani=ani, sparse=sparse, mixed=mixed,
                                         gpu=True, batch=batch, threads=threads),
                            dense=not sparse, names=names, ani=ani,
                        )
                        result = (compare_knn(cpu_knn, candidate, cpu_dense, ani) if sparse
                                  else compare_dense(cpu_dense, candidate))
                        report["results"][label] = result
                        print(f"PASS {label}: {json.dumps(result, sort_keys=True)}", flush=True)
                if mixed and not ani:
                    label = "mixed_jaccard_gpu_dense_singleton"
                    candidate = read_output(
                        distance_run(label, subset, mixed=True, gpu=True, batch=255),
                        dense=True, names=names,
                    )
                    result = compare_dense(cpu_dense, candidate)
                    report["results"][label] = result
                    print(f"PASS {label}: {json.dumps(result, sort_keys=True)}", flush=True)
        report["writer_failures"] = {}
        for threads in (1, 4):
            label = f"writer_failure_t{threads}"
            report["writer_failures"][label] = distance_run(
                label, subset, mixed=True, gpu=True, batch=64,
                threads=threads, expected_failure=True,
            )
        report["status"] = "passed"
        print(f"PASS: smoke and {len(report['results'])} CPU/GPU comparisons", flush=True)
    except Exception as error:
        report["status"] = "failed"
        report["error"] = str(error)
        traceback.print_exc()
        raise
    finally:
        (work / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        print(f"Report: {work / 'report.json'}", flush=True)


if __name__ == "__main__":
    main()
