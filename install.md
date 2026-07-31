# Installation Guide

This guide describes the software required by the WGS pipeline and how to use
`setup.sh`. The setup script installs software that can be downloaded without
registration, creates the shared Conda environment, downloads reference data,
and writes `config.yaml`.

ANNOVAR is not downloaded automatically. Its license requires users to register
and download the package themselves.

## System Requirements

Use a 64-bit Linux system with enough storage for tools, reference bundles, and
WGS output. Reference downloads and intermediate BAM/VCF files can require
hundreds of gigabytes.

Install these commands before running `setup.sh`:

- Bash
- Java 17 or newer
- Python 3
- Miniconda, Anaconda, or another compatible `conda` installation
- GNU Parallel
- Google Cloud CLI with `gsutil`
- Git
- Wget
- `tar`, `unzip`, and bzip2 support
- GCC and Make, required when SAMtools is built

For Ubuntu or Debian, most system packages can be installed with:

```bash
sudo apt update
sudo apt install -y \
  build-essential bzip2 git openjdk-17-jre-headless \
  parallel python3 tar unzip wget
```

Install Miniconda or Anaconda separately and make sure `conda` is available in
the shell. Install the Google Cloud CLI separately and confirm that `gsutil` is
available. Package names and Google Cloud installation steps vary by Linux
distribution.

Check the required commands:

```bash
java -version
python3 --version
conda --version
parallel --version
gsutil version
git --version
```

## Software Managed by setup.sh

The default tools directory is `tools/` under the repository. `setup.sh` can
download or install:

| Software | Installation method |
| --- | --- |
| FastQC 0.12.1 | Downloaded into `tools/FastQC/` |
| Trimmomatic 0.39 | Downloaded into `tools/Trimmomatic-0.39/` |
| BWA-MEM2 2.3 | Downloaded into `tools/bwa-mem2-2.3_x64-linux/` |
| SAMtools 1.22.1 | Downloaded, compiled, and installed under `tools/` |
| bcftools | Installed into the Conda base environment |
| GATK 4.6.2.0 | Downloaded into `tools/gatk-4.6.2.0/` |
| mosdepth plotting script | Downloaded into `tools/mosdepth/scripts/` |

The shared `wgs_parallel` Conda environment is created from
`wgs_parallel.yml`. It contains:

- Python 3.12
- MultiQC
- SeqKit
- mosdepth
- Requests, urllib3, and certifi for ClinVar updates

The mosdepth executable comes from the Conda environment. The separate
`tools/mosdepth/scripts/` directory contains the coverage plotting helper.

## ANNOVAR: Manual Academic Download

ANNOVAR cannot be installed directly by `setup.sh` because its download
requires registration and acceptance of its academic-use terms.

1. Register at <https://annovar.openbioinformatics.org/en/latest/>.
2. Download `annovar.latest.tar.gz` using the instructions supplied by ANNOVAR.
3. Extract the archive into the setup tools directory.

For the default tools directory, the final layout must be:

```text
tools/annovar/annotate_variation.pl
tools/annovar/convert2annovar.pl
tools/annovar/index_annovar.pl
tools/annovar/table_annovar.pl
tools/annovar/humandb/
```

For example, after downloading the archive into the repository directory:

```bash
mkdir -p tools
tar -xzf annovar.latest.tar.gz -C tools
mkdir -p tools/annovar/humandb
```

Do not redistribute the downloaded ANNOVAR archive or package through this
repository. If `--tools-dir /opt/wgs-tools` is used, extract ANNOVAR into
`/opt/wgs-tools/annovar/` instead.

During annotation, the pipeline can download missing ANNOVAR databases such as
`refGene`, `cytoBand`, `avsnp150`, and `exac03` into `humandb/`. When
`parameters.annovar.update_clinvar` is `yes`, the integrated updater downloads,
converts, and indexes ClinVar in the same directory. These operations require
internet access after ANNOVAR itself has been installed manually.

## Full Setup

Run the complete interactive setup from the repository root:

```bash
bash setup.sh
```

The script checks dependencies and prompts before each optional tool or
reference download. Answer `y` for the components required on the machine.
The default full setup performs three tasks:

1. Installs downloadable command-line software under `tools/`.
2. Creates the `wgs_parallel` Conda environment.
3. Downloads hg19 and hg38 reference bundles under `ref/`.

It then merges detected paths with `example/config/config.yaml` and writes
`config.yaml`. Existing documented pipeline parameters and an existing
`work_dir` are preserved when the file is updated.

Prompts require an interactive terminal. In a non-interactive shell, prompted
downloads are skipped rather than accepted automatically.

## Partial Setup Commands

Install only the downloadable command-line tools:

```bash
bash setup.sh download-tools
```

Create only the shared Conda environment:

```bash
bash setup.sh --conda-only
```

Download only hg38 references:

```bash
bash setup.sh download-ref --build 38
```

Download only hg19 references:

```bash
bash setup.sh download-ref --build 19
```

Download both reference builds:

```bash
bash setup.sh download-ref --build all
```

Equivalent legacy selectors are available:

```bash
bash setup.sh --software-only
bash setup.sh --conda-only
bash setup.sh --refs-only --build 38
```

## Custom Installation Directories

Set a different tools or reference root with:

```bash
bash setup.sh \
  --tools-dir /opt/wgs-tools \
  --ref-dir /data/wgs-references \
  --build 38
```

The aliases `--prefix` and `--refs-dir` are also accepted. The selected
directories must be writable. If ANNOVAR is installed manually, place it under
the selected tools root as `<tools-dir>/annovar/`.

## Reference Data

Reference bundles are downloaded with `gsutil`. Depending on the selected
build, the generated configuration includes paths for:

- Reference FASTA
- Known indels
- dbSNP
- Mills indels
- HapMap
- 1000 Genomes Omni
- 1000 Genomes Phase 1
- Mutect2 germline population resource
- Mutect2 panel of normals

Use one reference build for every sample in a pipeline invocation. External
normal BAMs must use the same reference dictionary and chromosome naming.

## Validation

Validate tools without installing or downloading anything:

```bash
bash setup.sh validate-only --validate-tools
```

Validate references for one build:

```bash
bash setup.sh validate-only --validate-refs --build 38
```

Validate both tools and references:

```bash
bash setup.sh validate-only --build all
```

After manually extracting ANNOVAR, rerun tool validation. ANNOVAR should be
reported as `OK` rather than `MANUAL` or `MISSING`.

Also verify the shared environment directly:

```bash
conda run -n wgs_parallel python --version
conda run -n wgs_parallel multiqc --version
conda run -n wgs_parallel seqkit version
conda run -n wgs_parallel mosdepth --version
```

If the environment already exists but `wgs_parallel.yml` has changed, update it
manually because `setup.sh` does not replace an existing environment:

```bash
conda env update -n wgs_parallel -f wgs_parallel.yml --prune
```

## Configure and Run the Pipeline

Review the generated configuration and sample manifest:

```bash
vim config.yaml
vim example/config/sample.csv
```

At minimum, confirm:

- `work_dir` is writable.
- Tool paths point to installed executables.
- The selected reference paths exist.
- `tools.annovar_dir` points to the manually installed ANNOVAR directory.
- Tumor-normal rows use compatible normal BAMs and indexes.

Run quality control, variant processing, or both:

```bash
bash wgs.sh qc  example/config/sample.csv config.yaml
bash wgs.sh run example/config/sample.csv config.yaml
bash wgs.sh all example/config/sample.csv config.yaml
```

## Restart and Troubleshooting

`setup.sh` is restart-aware. Files already installed at the expected paths are
skipped, so it is safe to rerun after an interrupted download or after ANNOVAR
has been installed manually.

Common installation problems:

- `required command not found: parallel`: install GNU Parallel.
- `required command not found: gsutil`: install Google Cloud CLI.
- `ANNOVAR ... MANUAL`: register, download, and extract ANNOVAR manually.
- `conda env 'wgs_parallel' already exists`: use `conda env update` if package
  definitions changed.
- Missing reference validation: rerun `download-ref` for the selected build or
  edit `config.yaml` to point to an existing compatible bundle.
- Non-interactive setup skipped downloads: rerun from an interactive terminal.

Do not start a production run until the relevant `validate-only` checks pass.
