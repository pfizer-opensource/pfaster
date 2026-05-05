#!/usr/bin/env python3
"""Run MashScreen across a directory of genomes and export results to CSV."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import csv
import os
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from predict import Mash


_WORKER_SCREENER: Mash.MashScreen | None = None
_WORKER_SEROTYPES: tuple[str, ...] = ()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Screen each genome in a directory with MashScreen and write a "
            "serotype-by-genome CSV table."
        )
    )
    parser.add_argument(
        "--genomes-dir",
        default="/workspace/data/VRD_Jonathan/training_genomes",
        help="Directory containing input genome FASTA files.",
    )
    parser.add_argument(
        "--sketch",
        default="sketch/sketch_k70_model73.pkl",
        help="Sketch database path (supports repo-relative and ref-relative paths).",
    )
    parser.add_argument(
        "--output",
        default="training_genomes_mashscreen.csv",
        help="Output CSV file path.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=os.cpu_count() or 1,
        help="Number of parallel worker processes to use.",
    )
    parser.add_argument(
        "--chunksize",
        type=int,
        default=0,
        help=(
            "Task chunk size for multiprocessing map. "
            "Use 0 to choose automatically."
        ),
    )
    return parser.parse_args()


def _init_worker(sketch_path: str) -> None:
    """Initialize one MashScreen instance per worker process."""
    global _WORKER_SCREENER, _WORKER_SEROTYPES
    _WORKER_SCREENER = Mash.MashScreen(sketch_path)
    _WORKER_SEROTYPES = tuple(_WORKER_SCREENER.minhash_db.keys())


def _screen_genome_worker(genome_path: str) -> tuple[str, dict[str, int]]:
    if _WORKER_SCREENER is None:
        raise RuntimeError("Worker not initialized")

    _WORKER_SCREENER.screen(genome_path)
    return (
        Path(genome_path).name,
        {serotype: _WORKER_SCREENER.ref_counts[serotype] for serotype in _WORKER_SEROTYPES},
    )


def resolve_sketch_path(sketch_arg: str, repo_root: Path) -> Path:
    sketch = Path(sketch_arg)
    candidates = [
        sketch,
        repo_root / sketch,
        repo_root / "ref" / sketch,
        repo_root / "ref" / "sketch" / sketch.name,
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()
    raise FileNotFoundError(f"Could not find sketch database: {sketch_arg}")


def collect_genomes(genomes_dir: Path) -> list[Path]:
    if not genomes_dir.exists():
        raise FileNotFoundError(f"Genome directory does not exist: {genomes_dir}")

    valid_suffixes = {".fasta", ".fa", ".fna", ".fas"}
    genomes = [
        path
        for path in genomes_dir.rglob("*")
        if path.is_file() and path.suffix.lower() in valid_suffixes
    ]
    return sorted(genomes)


def main() -> None:
    args = parse_args()
    repo_root = REPO_ROOT

    genomes_dir = Path(args.genomes_dir)
    output_path = Path(args.output)
    if not output_path.is_absolute():
        output_path = (repo_root / output_path).resolve()

    sketch_path = resolve_sketch_path(args.sketch, repo_root)
    genomes = collect_genomes(genomes_dir)
    if not genomes:
        raise RuntimeError(f"No FASTA files found in {genomes_dir}")

    workers = max(1, args.workers)

    # Load once in the parent process to determine output columns.
    screener = Mash.MashScreen(str(sketch_path))
    serotypes = list(screener.minhash_db.keys())

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as out_handle:
        writer = csv.DictWriter(out_handle, fieldnames=["genome"] + serotypes)
        writer.writeheader()

        total = len(genomes)
        if workers == 1:
            for idx, genome_path in enumerate(genomes, start=1):
                screener.screen(str(genome_path))
                row = {"genome": genome_path.name}
                for serotype in serotypes:
                    row[serotype] = screener.ref_counts[serotype]
                writer.writerow(row)

                if idx % 100 == 0 or idx == total:
                    print(f"Processed {idx}/{total}: {genome_path.name}")
        else:
            print(f"Processing {total} genomes with {workers} workers")
            chunksize = args.chunksize
            if chunksize <= 0:
                chunksize = max(1, total // (workers * 4))

            with ProcessPoolExecutor(
                max_workers=workers,
                initializer=_init_worker,
                initargs=(str(sketch_path),),
            ) as pool:
                genome_iter = (str(genome_path) for genome_path in genomes)
                results = pool.map(_screen_genome_worker, genome_iter, chunksize=chunksize)

                for idx, (genome_name, counts) in enumerate(results, start=1):
                    row = {"genome": genome_name}
                    for serotype in serotypes:
                        row[serotype] = counts[serotype]
                    writer.writerow(row)

                    if idx % 100 == 0 or idx == total:
                        print(f"Processed {idx}/{total}: {genome_name}")

    print(f"Wrote MashScreen matrix to {output_path}")


if __name__ == "__main__":
    main()
