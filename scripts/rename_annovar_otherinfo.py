#!/usr/bin/env python3
"""Replace ANNOVAR Otherinfo headers with generated and original VCF names."""

from __future__ import annotations

import argparse
import gzip
import os
import re
import tempfile
from pathlib import Path


OTHERINFO = re.compile(r"^Otherinfo(?:\d+)?$")
GENERATED_NAMES = ("ANNOVAR_AF", "ANNOVAR_QUAL", "ANNOVAR_DP")


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def vcf_headers(path: Path) -> list[str]:
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#CHROM\t"):
                return line.rstrip("\r\n").split("\t")
    raise ValueError(f"VCF column header not found: {path}")


def avinput_extra_count(path: Path) -> int:
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.strip() and not line.startswith("#"):
                return max(0, len(line.rstrip("\r\n").split("\t")) - 5)
    raise ValueError(f"No variant records found in avinput: {path}")


def replacement_headers(vcf: Path, avinput: Path) -> list[str]:
    original = vcf_headers(vcf)
    generated_count = avinput_extra_count(avinput) - len(original)
    if generated_count < 0:
        raise ValueError("avinput has fewer extra columns than the VCF header")
    generated = list(GENERATED_NAMES[:generated_count])
    if generated_count > len(GENERATED_NAMES):
        generated.extend(
            f"ANNOVAR_EXTRA{i}"
            for i in range(len(GENERATED_NAMES) + 1, generated_count + 1)
        )
    return generated + original


def rename_headers(multianno: Path, names: list[str]) -> None:
    mode = multianno.stat().st_mode
    with multianno.open(encoding="utf-8") as source:
        header_line = source.readline()
        if not header_line:
            raise ValueError(f"Empty multianno file: {multianno}")
        header = header_line.rstrip("\r\n").split("\t")
        indexes = [index for index, name in enumerate(header) if OTHERINFO.fullmatch(name)]
        if not indexes:
            return
        if len(indexes) != len(names):
            raise ValueError(
                f"Otherinfo/header count mismatch: {len(indexes)} columns, "
                f"{len(names)} replacement names"
            )
        for index, name in zip(indexes, names):
            header[index] = name

        with tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", dir=multianno.parent, delete=False
        ) as destination:
            temporary = Path(destination.name)
            destination.write("\t".join(header) + "\n")
            for line in source:
                destination.write(line)
    os.chmod(temporary, mode)
    os.replace(temporary, multianno)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--multianno", required=True, type=Path)
    parser.add_argument("--vcf", required=True, type=Path)
    parser.add_argument("--avinput", required=True, type=Path)
    args = parser.parse_args()
    rename_headers(args.multianno, replacement_headers(args.vcf, args.avinput))


if __name__ == "__main__":
    main()
