# WGS Parallel Pipeline

A Bash workflow for WGS QC, germline variant calling, optional matched-normal
Mutect2 calling, and ANNOVAR annotation.

## Quick Start

```bash
# Install/detect tools and references and write ./config.yaml.
bash setup.sh

# Edit the generated machine/run configuration and sample manifest.
vim config.yaml
vim example/config/sample.csv

bash wgs.sh qc  example/config/sample.csv config.yaml
bash wgs.sh run example/config/sample.csv config.yaml
# Or run both phases:
bash wgs.sh all example/config/sample.csv config.yaml
```

The CLI defaults to `example/config/sample.csv` and `./config.yaml` when those
arguments are omitted. `run` expects the trimmed FASTQs produced by `qc`.

## Inputs

`sample.csv` has five columns:

```csv
sample_id,sample_prefix,ref_version,raw_dir,normal_id,normal_bam
TUMOR_001,TUMOR_001,38,/data/linked_raw/TUMOR_001,NORMAL_001,/data/normals/NORMAL_001.bam
GERMLINE_001,GERMLINE_001,38,/data/linked_raw/GERMLINE_001,,
```

Sequencing companies often return FASTQs under different nested directory
layouts. Link them into the per-sample `raw_dir` locations by giving the
company's result directory as the source:

```bash
bash scripts/link_fastqc.sh \
  /data/sequencing_company/result_dir \
  /data/project/raw_link
```

The linker finds every FASTQ file recursively and creates flat symlinks directly
under `raw_link/`. It does not read or write `sample.csv`; maintain that
manifest separately. `link_fastq.sh` remains an alias.

- All rows run through QC, alignment, BQSR, and germline calling.
- A nonblank `normal_bam` requests matched-normal Mutect2 for that row.
- `normal_id` is the expected BAM read-group sample name and is checked before calling.
- External normal BAMs must be coordinate sorted, indexed, contain one sample
  name in their read groups, and use the selected reference dictionary.
- Invalid normal BAMs are skipped and reported in
  `output/10_Mutect2/skipped_pairs.tsv`; they do not stop germline work.
- One invocation must use a single reference build.

`config.yaml` is the only YAML configuration. It contains `work_dir`, tool
paths, build-specific reference paths, and a `parameters` section for every
pipeline stage. Each step provides typed concurrency/tuning settings and an
optional `extra_args` string. Arguments use shell quoting but are parsed into
an array and are never evaluated as shell code.

`setup.sh` updates `tools` and `references` atomically while preserving the
existing `work_dir` and `parameters`. The complete example is at
`example/config/config.yaml`.

## Primary-Contig Calling

Before calling, the pipeline resolves the reference's exact names for
chromosomes 1-22, X, Y, and mitochondrial DNA. Both HaplotypeCaller and Mutect2
are restricted to those 25 contigs. Aliases such as `1`/`chr1` and `MT`/`chrM`
are supported.

HaplotypeCaller scatters each sample into 25 interval jobs that read the same
indexed BQSR BAM. It gathers the shards in canonical chromosome order into the
existing `${sample_id}.g.vcf.gz` output. Missing shards resume independently;
successful shards are removed after gathering unless
`parameters.haplotype_caller.keep_scatter_gvcfs` is `yes`.

## Commands And Layout

```bash
bash wgs.sh qc   [sample.csv] [config.yaml]
bash wgs.sh run  [sample.csv] [config.yaml]
bash wgs.sh all  [sample.csv] [config.yaml]
bash wgs.sh stop
```

The supported implementation is `wgs.sh`, `scripts/00_util.sh`, and numbered
step scripts in `scripts/`. Older duplicated entry points are retained only for
reference under `legacy/` and are not maintained.

Important generated files in the work directory include:

```text
sample_list.txt
primary_contigs.tsv
output/logs/<command>_<timestamp>_<pid>.log
output/logs/parallel/<timestamp>_<pid>/<step>.joblog
output/01_fastqc/
output/02_trim/
output/03_fastqc_trim/
output/04_bwa/
output/05_markduplicated/
output/06_bqsr/
output/07_HaplotypeCaller/
output/08_VQSR/
output/09_annovar/
output/10_Mutect2/{pairs,skipped_pairs,failed_pairs}.tsv
output/10_Mutect2/success_samples.txt
output/11_annovar_somatic/
```

Parallel task checkpoints are stored under:

```text
output/logs/checkpoints/<step>/<scope>/running/
output/logs/checkpoints/<step>/<scope>/done/
```

Each worker writes a `running` marker before starting and a `done` marker only
after successful completion. On restart, leftover `running` markers identify
interrupted tasks. The pipeline removes only those tasks' partial outputs and
retries them; completed sample and chromosome outputs are retained.

To limit `samtools sort` memory, configure its concurrency, thread count, and
per-thread memory independently. For example, the following uses approximately
4 GB for sort buffers when two jobs are active, plus process overhead:

```yaml
parameters:
  sort:
    max_jobs: 2
    threads_per_job: 2
    memory_per_thread: 1G
```
