#!/usr/bin/env python3
"""Convert ClinVar avinput records to an ANNOVAR annotation database."""

from __future__ import annotations

import argparse
from pathlib import Path


DEFAULT_FIELDS = ("ALLELEID", "CLNDN", "CLNDISDB", "CLNREVSTAT", "CLNSIG")
OUTPUT_FIELDS = ("CLNALLELEID", "CLNDN", "CLNDISDB", "CLNREVSTAT", "CLNSIG")


def clinvaravinput2annovardb(
    avinput_file: str | Path,
    fields: list[str] | tuple[str, ...] = DEFAULT_FIELDS,
    output_file: str | Path | None = None,
) -> str:
    """Write an ANNOVAR ClinVar database and return its absolute path."""
    source = Path(avinput_file).expanduser().resolve()
    destination = (source.with_suffix(source.suffix + ".txt") if output_file is None
                   else Path(output_file).expanduser().resolve())
    destination.parent.mkdir(parents=True, exist_ok=True)

    wanted = set(fields)
    with source.open(encoding="utf-8") as input_handle, destination.open("w", encoding="utf-8") as output:
        output.write("#Chr\tStart\tEnd\tRef\tAlt\t" + "\t".join(OUTPUT_FIELDS) + "\n")
        for line in input_handle:
            if not line.strip() or line.startswith("#"):
                continue
            columns = line.rstrip("\n").split("\t")
            if len(columns) < 13:
                continue
            values = {field: "." for field in OUTPUT_FIELDS}
            for item in columns[12].split(";"):
                key, separator, value = item.partition("=")
                if not separator or key not in wanted:
                    continue
                output_key = "CLNALLELEID" if key == "ALLELEID" else key
                if output_key in values:
                    values[output_key] = value or "."
            annotation = [values[field].replace(",", r"\x2c") for field in OUTPUT_FIELDS]
            output.write("\t".join(columns[:5] + annotation) + "\n")
    return str(destination)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("avinput", type=Path)
    parser.add_argument("--output", type=Path, help="Output ANNOVAR database path")
    args = parser.parse_args()
    print(clinvaravinput2annovardb(args.avinput, DEFAULT_FIELDS, args.output))


if __name__ == "__main__":
    main()
