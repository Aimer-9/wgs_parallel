#!/usr/bin/env python3
"""Recursively symlink FASTQs from a result tree into one raw_link directory."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path


FASTQ_SUFFIXES = (".fq", ".fastq", ".fq.gz", ".fastq.gz")


def is_fastq(path: Path) -> bool:
    lower = path.name.lower()
    return lower.endswith(FASTQ_SUFFIXES)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Find FASTQ files recursively under a result directory and symlink them into one raw_link directory."
    )
    parser.add_argument(
        "paths", nargs=2, type=Path,
        help="RESULT_DIR RAW_LINK",
    )
    parser.add_argument("--dry-run", action="store_true", help="print links without creating files")
    args = parser.parse_args()
    args.result_dir, args.raw_link = args.paths

    source = args.result_dir.resolve()
    if not source.is_dir():
        parser.error(f"result directory not found: {source}")

    rawdata_roots = sorted(
        path for path in source.rglob("*")
        if path.is_dir() and path.name.lower() == "rawdata"
    )
    search_roots = rawdata_roots or [source]
    files = sorted(
        path for root in search_roots
        for path in root.rglob("*")
        if path.is_file() and is_fastq(path)
    )
    if not files:
        parser.error(f"no FASTQ files found under result directory: {source}")
    raw_link = args.raw_link.resolve()
    linked = 0
    if not args.dry_run:
        raw_link.mkdir(parents=True, exist_ok=True)
    for source_file in files:
        destination = raw_link / source_file.name
        print(f"{destination} -> {source_file}")
        if args.dry_run:
            continue
        if destination.exists() or destination.is_symlink():
            if destination.is_symlink() and destination.resolve() == source_file:
                continue
            raise SystemExit(f"refusing to overwrite existing path: {destination}")
        os.symlink(source_file, destination)
        linked += 1

    print(f"found {len(files)} FASTQ files; linked {linked} files into {raw_link}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except OSError as exc:
        print(f"link error: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
