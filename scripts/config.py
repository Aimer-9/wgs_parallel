#!/usr/bin/env python3
"""Configuration and metadata normalization for the WGS pipeline.

This command has three subcommands:

``load``
    Validate config.yaml and sample.csv, write a headerless normalized TSV
    manifest, and print shell-quoted ``export`` statements consumed by wgs.sh.
``contigs``
    Read a FASTA index/dictionary and emit the ordered primary chromosomes
    1-22, X, Y, and MT while preserving reference-specific names such as chrM.
``merge-machine``
    Merge machine-detected tool/reference paths into the example YAML template
    without discarding the template's documented parameter defaults.

All user-facing validation failures are raised as ConfigError so the CLI can
print concise messages without Python tracebacks.
"""

from __future__ import annotations

import argparse
import csv
import os
import shlex
import sys
from pathlib import Path
from typing import Any

try:
    import yaml
except ModuleNotFoundError as exc:  # pragma: no cover - exercised by setup validation
    raise SystemExit("Config error: PyYAML is required (python3 -m pip install pyyaml)") from exc


class ConfigError(Exception):
    """Expected configuration or sample-input validation failure."""

    pass


def load_yaml(path: Path) -> dict[str, Any]:
    """Load one YAML mapping and reject unreadable or non-mapping roots."""
    try:
        data = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as exc:
        raise ConfigError(f"cannot read YAML {path}: {exc}") from exc
    if not isinstance(data, dict):
        raise ConfigError(f"{path}: root must be a mapping")
    return data


def mapping(parent: dict[str, Any], key: str) -> dict[str, Any]:
    """Return a nested mapping, treating an omitted/null section as empty."""
    value = parent.get(key, {})
    if value is None:
        return {}
    if not isinstance(value, dict):
        raise ConfigError(f"'{key}' must be a mapping")
    return value


def scalar(section: dict[str, Any], key: str, default: Any = "") -> str:
    """Read a scalar as text while rejecting accidental lists/mappings."""
    value = section.get(key, default)
    if value is None:
        return ""
    if isinstance(value, (dict, list)):
        raise ConfigError(f"'{key}' must be a scalar")
    return str(value)


def integer(section: dict[str, Any], key: str, default: int, minimum: int = 1) -> str:
    """Validate an integer parameter and return its shell-exportable text."""
    value = section.get(key, default)
    if isinstance(value, bool):
        raise ConfigError(f"'{key}' must be an integer")
    try:
        parsed = int(value)
    except (TypeError, ValueError) as exc:
        raise ConfigError(f"'{key}' must be an integer") from exc
    if parsed < minimum:
        raise ConfigError(f"'{key}' must be >= {minimum}")
    return str(parsed)


def decimal(section: dict[str, Any], key: str, default: float) -> str:
    """Validate a CPU fraction in the inclusive range (0, 1]."""
    value = section.get(key, default)
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise ConfigError(f"'{key}' must be numeric") from exc
    if not 0 < parsed <= 1:
        raise ConfigError(f"'{key}' must be > 0 and <= 1")
    return str(parsed)


def yes_no(section: dict[str, Any], key: str, default: str = "no") -> str:
    """Normalize a strict YAML yes/no setting to lowercase text."""
    value = scalar(section, key, default).lower()
    if value not in {"yes", "no"}:
        raise ConfigError(f"'{key}' must be 'yes' or 'no'")
    return value


def extra_args(section: dict[str, Any], key: str = "extra_args", default: str = "", forbidden=()) -> str:
    """Validate shell-like extra arguments and block pipeline-owned options.

    The string is parsed with shlex only for validation. Stage scripts parse it
    again into Bash arrays, preserving quoted argument values without eval.
    """
    value = scalar(section, key, default)
    try:
        tokens = shlex.split(value)
    except ValueError as exc:
        raise ConfigError(f"'{key}' contains invalid shell quoting: {exc}") from exc
    for token in tokens:
        option = token.split("=", 1)[0]
        if option in forbidden:
            raise ConfigError(f"'{key}' cannot override pipeline-managed option '{option}'")
    return value


def comma_items(value: str, field: str) -> list[str]:
    """Split a comma list and reject empty entries that shift ANNOVAR columns."""
    items = [item.strip() for item in value.split(",")]
    if not items or any(not item for item in items):
        raise ConfigError(f"'{field}' must be a comma-separated list without empty entries")
    return items


def suggest_annovar_operations(protocols: list[str]) -> str:
    """Suggest ANNOVAR operation types for common gene/region/filter databases."""
    gene_protocols = {"refgene", "knowngene", "ensgene"}
    operations = []
    for protocol in protocols:
        lowered = protocol.lower()
        if lowered in gene_protocols:
            operations.append("g")
        elif lowered == "cytoband":
            operations.append("r")
        else:
            operations.append("f")
    return ",".join(operations)


def normalize_ref(raw: str) -> str:
    """Normalize hg19/GRCh37 to 19 and hg38/GRCh38 to 38."""
    value = raw.strip().lower().replace("grch", "").replace("hg", "")
    if value in {"19", "37"}:
        return "19"
    if value == "38":
        return "38"
    raise ConfigError(f"invalid ref_version '{raw}' (expected 19/hg19/37 or 38/hg38)")


def load_samples(path: Path) -> list[dict[str, str]]:
    """Validate sample.csv and return normalized absolute-path records.

    Required columns are sample_id, sample_prefix, ref_version, and raw_dir.
    normal_id and normal_bam are optional. A run must use one reference build,
    sample IDs must be unique/shell-safe, and every raw_dir must already exist.
    """
    rows: list[dict[str, str]] = []
    try:
        handle = path.open("r", encoding="utf-8-sig", newline="")
    except OSError as exc:
        raise ConfigError(f"cannot read sample CSV {path}: {exc}") from exc
    with handle:
        reader = csv.DictReader(handle)
        fields = set(reader.fieldnames or [])
        required = {"sample_id", "sample_prefix", "ref_version", "raw_dir"}
        missing = sorted(required - fields)
        if missing:
            raise ConfigError(f"{path}: missing required columns: {', '.join(missing)}")
        for line, row in enumerate(reader, start=2):
            sample_id = (row.get("sample_id") or "").strip()
            sample_prefix = (row.get("sample_prefix") or "").strip()
            ref_raw = (row.get("ref_version") or "").strip()
            raw_dir = (row.get("raw_dir") or "").strip()
            normal_bam = (row.get("normal_bam") or "").strip()
            normal_id = (row.get("normal_id") or "").strip()
            if not sample_id or not sample_prefix or not ref_raw or not raw_dir:
                raise ConfigError(f"{path}:{line}: sample_id, sample_prefix, ref_version, and raw_dir are required")
            for field_name, field_value in (("sample_id", sample_id), ("sample_prefix", sample_prefix)):
                if any(ch in field_value for ch in "\t\r\n/"):
                    raise ConfigError(f"{path}:{line}: unsafe {field_name} '{field_value}'")
            if not os.path.isdir(raw_dir):
                raise ConfigError(f"{path}:{line}: raw_dir not found: {raw_dir}")
            rows.append({
                "sample_id": sample_id,
                "sample_prefix": sample_prefix,
                "ref_version": normalize_ref(ref_raw),
                "raw_dir": os.path.abspath(raw_dir),
                "normal_bam": os.path.abspath(normal_bam) if normal_bam else "",
                "normal_id": normal_id,
            })
    if not rows:
        raise ConfigError(f"{path}: no sample rows found")
    ids = [row["sample_id"] for row in rows]
    if len(ids) != len(set(ids)):
        duplicate = next(item for item in ids if ids.count(item) > 1)
        raise ConfigError(f"{path}: duplicate sample_id: {duplicate}")
    builds = sorted({row["ref_version"] for row in rows})
    if len(builds) != 1:
        raise ConfigError(f"sample.csv must contain one ref_version per run; found: {', '.join(builds)}")
    return rows


def q(value: Any) -> str:
    """Quote one value for safe evaluation by the calling Bash process."""
    return shlex.quote(str(value))


def build_environment(config: dict[str, Any], samples: list[dict[str, str]]) -> dict[str, str]:
    """Flatten YAML tools, references, and parameters into shell variables.

    Defaults are applied here so stage scripts receive a complete contract even
    when optional YAML keys are omitted. Extra-argument fields also declare
    options owned by the pipeline that users must not override.
    """
    work_dir = scalar(config, "work_dir")
    if not work_dir:
        raise ConfigError("missing required field: work_dir")
    tools = mapping(config, "tools")
    refs = mapping(config, "references")
    parameters = mapping(config, "parameters")
    build = samples[0]["ref_version"]
    build_refs = refs.get(build, {})
    if not isinstance(build_refs, dict):
        raise ConfigError(f"references.'{build}' must be a mapping")

    global_p = mapping(parameters, "global")
    raw_qc = mapping(parameters, "fastqc_raw")
    trim = mapping(parameters, "trim")
    trim_qc = mapping(parameters, "fastqc_trimmed")
    bwa = mapping(parameters, "bwa")
    sort = mapping(parameters, "sort")
    index = mapping(parameters, "index")
    mosdepth_p = mapping(parameters, "mosdepth")
    markdup = mapping(parameters, "mark_duplicates")
    bqsr = mapping(parameters, "bqsr")
    hc = mapping(parameters, "haplotype_caller")
    genotype = mapping(parameters, "genotype_gvcfs")
    vqsr = mapping(parameters, "vqsr")
    annovar = mapping(parameters, "annovar")
    mutect = mapping(parameters, "mutect2")

    annovar_protocol = scalar(annovar, "protocol", "refGene,cytoBand,avsnp150,exac03,{clinvar}")
    annovar_operation = scalar(annovar, "operation", "g,r,f,f,f")
    protocol_items = comma_items(annovar_protocol, "parameters.annovar.protocol")
    operation_items = comma_items(annovar_operation, "parameters.annovar.operation")
    malformed_placeholders = [
        item for item in protocol_items
        if ("{" in item or "}" in item) and item != "{clinvar}"
    ]
    if malformed_placeholders:
        raise ConfigError(
            "parameters.annovar.protocol contains malformed placeholder(s): "
            + ", ".join(malformed_placeholders)
            + "; only the exact token {clinvar} may contain braces"
        )
    if protocol_items.count("{clinvar}") > 1:
        raise ConfigError("parameters.annovar.protocol may contain {clinvar} only once")
    if len(protocol_items) != len(operation_items):
        suggestion = suggest_annovar_operations(protocol_items)
        raise ConfigError(
            "parameters.annovar protocol/operation count mismatch: "
            f"{len(protocol_items)} protocols but {len(operation_items)} operations; "
            f"set operation to a {len(protocol_items)}-item list, for example: {suggestion}"
        )
    invalid_operations = sorted(set(operation_items) - {"g", "r", "f"})
    if invalid_operations:
        raise ConfigError(
            "parameters.annovar.operation supports only g, r, or f; invalid: "
            + ", ".join(invalid_operations)
        )

    env = {
        "work_dir": os.path.abspath(work_dir),
        "hsa_version": build,
        "MAX_PROCESSOR_USE_PERCENT": decimal(global_p, "max_processor_use_fraction", 0.8),
        "update_clinvar": yes_no(annovar, "update_clinvar", "no"),
        "fastqc": scalar(tools, "fastqc"),
        "trim_path": scalar(tools, "trimmomatic_dir"),
        "java": scalar(tools, "java", "java"),
        "conda_bin": scalar(tools, "conda", "conda"),
        "bwamem": scalar(tools, "bwa_mem2"),
        "samtools": scalar(tools, "samtools"),
        "gatk": scalar(tools, "gatk"),
        "bcftools": scalar(tools, "bcftools", "bcftools"),
        "parallel_bin": scalar(tools, "parallel", "parallel"),
        "multiqc": scalar(tools, "multiqc", "multiqc"),
        "seqkit": scalar(tools, "seqkit", "seqkit"),
        "mosdepth": scalar(tools, "mosdepth", "mosdepth"),
        "annovar_path": scalar(tools, "annovar_dir"),
        "mosdepth_path": scalar(tools, "mosdepth_dir"),
        "ref_fa": scalar(build_refs, "fasta"),
        "known_indels": scalar(build_refs, "known_indels"),
        "dbsnp": scalar(build_refs, "dbsnp"),
        "Mills": scalar(build_refs, "mills"),
        "hapmap": scalar(build_refs, "hapmap"),
        "file_1000G_omni": scalar(build_refs, "omni_1000g"),
        "file_1000G_phase1": scalar(build_refs, "phase1_1000g"),
        "germline_resource": scalar(build_refs, "germline_resource"),
        "panel_of_normals": scalar(build_refs, "panel_of_normals"),
        "MAX_FASTQC_RAW_JOBS": integer(raw_qc, "max_jobs", 10),
        "FASTQC_RAW_Q30_THREADS": integer(raw_qc, "q30_threads", 20),
        "FASTQC_RAW_EXTRA_ARGS": extra_args(raw_qc),
        "MAX_TRIM_JOBS": integer(trim, "max_jobs", 10),
        "TRIM_LEADING": integer(trim, "leading_quality", 10, 0),
        "TRIM_TRAILING": integer(trim, "trailing_quality", 10, 0),
        "TRIM_MINLEN": integer(trim, "minimum_length", 25, 1),
        "TRIM_SLIDINGWINDOW": scalar(trim, "sliding_window", "4:15"),
        "TRIM_EXTRA_ARGS": extra_args(trim),
        "MAX_FASTQC_TRIM_JOBS": integer(trim_qc, "max_jobs", 10),
        "FASTQC_TRIM_Q30_THREADS": integer(trim_qc, "q30_threads", 20),
        "FASTQC_TRIM_EXTRA_ARGS": extra_args(trim_qc),
        "MAX_BWA_JOBS": integer(bwa, "max_jobs", 3),
        "BWA_EXTRA_ARGS": extra_args(bwa, default="-v 1 -M", forbidden=("-R", "-t")),
        "MAX_SORT_JOBS": integer(sort, "max_jobs", 2),
        "SORT_THREADS_PER_JOB": integer(sort, "threads_per_job", 2),
        "SORT_MEMORY_PER_THREAD": scalar(sort, "memory_per_thread", "1G"),
        "SORT_EXTRA_ARGS": extra_args(sort, forbidden=("-o", "-@", "-m", "-O")),
        "MAX_INDEX_JOBS": integer(index, "max_jobs", 10),
        "MAX_MOSDEPTH_JOBS": integer(mosdepth_p, "max_jobs", 10),
        "MOSDEPTH_THREADS": integer(mosdepth_p, "threads", 4),
        "MOSDEPTH_EXTRA_ARGS": extra_args(mosdepth_p, default="-n", forbidden=("-t",)),
        "MAX_MARKDUP_JOBS": integer(markdup, "max_jobs", 10),
        "MARKDUP_EXTRA_ARGS": extra_args(markdup, forbidden=("-I", "-O", "-M")),
        "MAX_BQSR_JOBS": integer(bqsr, "max_jobs", 10),
        "BQSR_MIN_OUTPUT_BYTES": integer(bqsr, "minimum_output_bytes", 10485760, 1),
        "BQSR_EXTRA_ARGS": extra_args(bqsr, forbidden=("-I", "-O", "-R", "--known-sites")),
        "MAX_HC_JOBS": integer(hc, "max_jobs", 10),
        "HC_NATIVE_THREADS": integer(hc, "native_threads", 1),
        "HC_KEEP_SCATTER": yes_no(hc, "keep_scatter_gvcfs", "no"),
        "HC_EXTRA_ARGS": extra_args(hc, forbidden=("-I", "-O", "-R", "-L", "--intervals", "-XL", "--exclude-intervals")),
        "MAX_GENOTYPE_JOBS": integer(genotype, "max_jobs", 10),
        "GENOTYPE_EXTRA_ARGS": extra_args(genotype, forbidden=("-V", "--variant", "-O", "-R")),
        "MAX_VQSR_JOBS": integer(vqsr, "max_jobs", 10),
        "VQSR_SNP_TRUTH_SENSITIVITY": scalar(vqsr, "snp_truth_sensitivity", "99.5"),
        "VQSR_INDEL_TRUTH_SENSITIVITY": scalar(vqsr, "indel_truth_sensitivity", "99.0"),
        "VQSR_EXTRA_ARGS": extra_args(vqsr, forbidden=("-V", "-O", "-R", "-mode")),
        "MAX_ANNOVAR_JOBS": integer(annovar, "max_jobs", 10),
        "ANNOVAR_PROTOCOL": annovar_protocol,
        "ANNOVAR_OPERATION": annovar_operation,
        "ANNOVAR_EXTRA_ARGS": extra_args(annovar),
        "MAX_MUTECT2_JOBS": integer(mutect, "max_jobs", 10),
        "MUTECT2_EXTRA_ARGS": extra_args(mutect, forbidden=("-I", "-O", "-R", "-L", "--intervals", "-XL", "--exclude-intervals", "-normal", "--normal-sample")),
        "FILTER_MUTECT_EXTRA_ARGS": extra_args(mutect, "filter_extra_args", forbidden=("-V", "-O", "-R")),
    }
    return env


def command_load(args: argparse.Namespace) -> None:
    """Implement ``load``: write the TSV manifest and print Bash exports."""
    config_path = Path(args.config).resolve()
    sample_path = Path(args.samples).resolve()
    config = load_yaml(config_path)
    samples = load_samples(sample_path)
    env = build_environment(config, samples)
    manifest = Path(args.manifest).resolve()
    manifest.parent.mkdir(parents=True, exist_ok=True)
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        for row in samples:
            handle.write("\t".join(row[key] for key in ("sample_id", "sample_prefix", "ref_version", "raw_dir", "normal_id", "normal_bam")) + "\n")
    env.update({
        "config_yaml": str(config_path),
        "sample_csv": str(sample_path),
        "sample_manifest_tsv": str(manifest),
    })
    for key, value in env.items():
        print(f"export {key}={q(value)}")


def read_reference_contigs(fasta: Path) -> list[tuple[str, int]]:
    """Read contig names/lengths from FASTA .fai, falling back to .dict."""
    fai = Path(f"{fasta}.fai")
    dictionary = fasta.with_suffix(".dict")
    contigs: list[tuple[str, int]] = []
    if fai.is_file():
        for line in fai.read_text(encoding="utf-8").splitlines():
            fields = line.split("\t")
            if len(fields) >= 2:
                contigs.append((fields[0], int(fields[1])))
    elif dictionary.is_file():
        for line in dictionary.read_text(encoding="utf-8").splitlines():
            if not line.startswith("@SQ"):
                continue
            values = dict(field.split(":", 1) for field in line.split("\t")[1:] if ":" in field)
            if "SN" in values and "LN" in values:
                contigs.append((values["SN"], int(values["LN"])))
    else:
        raise ConfigError(f"reference index not found: expected {fai} or {dictionary}")
    return contigs


def command_contigs(args: argparse.Namespace) -> None:
    """Emit the ordered 25-contig table using names present in the reference."""
    contigs = read_reference_contigs(Path(args.fasta))
    by_name = {name: length for name, length in contigs}
    logical_names = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]
    for order, logical in enumerate(logical_names, start=1):
        aliases = [logical, f"chr{logical}"]
        if logical == "MT":
            aliases = ["MT", "chrM", "M", "chrMT"]
        matches = [name for name in aliases if name in by_name]
        if len(matches) != 1:
            raise ConfigError(
                f"primary contig {logical}: expected one of {', '.join(aliases)}, found {matches or 'none'}"
            )
        actual = matches[0]
        print(f"{order:02d}\t{logical}\t{actual}\t{by_name[actual]}")


def command_merge_machine(args: argparse.Namespace) -> None:
    """Overlay detected machine paths on the distributable YAML template."""
    output = Path(args.output).resolve()
    base_path = output if output.is_file() else Path(args.template).resolve()
    base = load_yaml(base_path)
    detected = load_yaml(Path(args.detected).resolve())
    detected_tools = detected.get("tools")
    detected_refs = detected.get("references")
    if not isinstance(detected_tools, dict) or not isinstance(detected_refs, dict):
        raise ConfigError("detected config must contain tools and references mappings")
    base["tools"] = detected_tools
    existing_refs = base.get("references", {})
    if not isinstance(existing_refs, dict):
        existing_refs = {}
    existing_refs.update(detected_refs)
    base["references"] = existing_refs
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.tmp")
    temporary.write_text(yaml.safe_dump(base, sort_keys=False), encoding="utf-8")
    os.replace(temporary, output)


def parser() -> argparse.ArgumentParser:
    """Build the subcommand-oriented command-line parser."""
    root = argparse.ArgumentParser()
    sub = root.add_subparsers(dest="command", required=True)
    load = sub.add_parser("load")
    load.add_argument("--config", required=True)
    load.add_argument("--samples", required=True)
    load.add_argument("--manifest", required=True)
    load.set_defaults(func=command_load)
    contigs = sub.add_parser("contigs")
    contigs.add_argument("--fasta", required=True)
    contigs.set_defaults(func=command_contigs)
    merge = sub.add_parser("merge-machine")
    merge.add_argument("--detected", required=True)
    merge.add_argument("--template", required=True)
    merge.add_argument("--output", required=True)
    merge.set_defaults(func=command_merge_machine)
    return root


def main() -> None:
    """Dispatch the CLI and translate ConfigError into a concise exit message."""
    args = parser().parse_args()
    try:
        args.func(args)
    except ConfigError as exc:
        print(f"Config error: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc


if __name__ == "__main__":
    main()
