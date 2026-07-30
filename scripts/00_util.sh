#!/bin/bash

# Shared runtime library for wgs.sh and every numbered pipeline stage.
#
# This file is sourced rather than executed. It centralizes:
#   - config.yaml/sample.csv parsing through config.py;
#   - primary-contig discovery (1-22, X, Y, and MT only);
#   - CPU allocation for local, Slurm, SGE, and PBS environments;
#   - recursive FASTQ lookup using sample_id/sample_prefix metadata;
#   - GNU Parallel execution with durable running/done checkpoints;
#   - run logging, single-pipeline locking, and recursive process shutdown.
#
# Several functions read globals exported by config.py, especially work_dir,
# sample_manifest_tsv, ref_fa, parallel_bin, and MAX_PROCESSOR_USE_PERCENT.
# Worker functions passed to parallel_run_sample_list must accept one task
# string. A matching cleanup_<worker_function> must remove partial outputs.

# Load and validate YAML/CSV configuration, then evaluate shell-safe exports
# emitted by config.py. CONFIG_YAML and SAMPLE_CSV can override default paths.
load_pipeline_config() {
  local repo_dir=$1
  local config_yaml="${CONFIG_YAML:-${repo_dir}/config.yaml}"
  local sample_csv="${SAMPLE_CSV:-${repo_dir}/example/config/sample.csv}"
  local manifest_tsv="${repo_dir}/sample_manifest.tsv"

  if [ ! -f "$config_yaml" ]; then
    echo "Config error: config.yaml not found: $config_yaml" >&2
    echo "Run setup.sh or pass a config path to wgs.sh." >&2
    return 1
  fi
  if [ ! -f "$sample_csv" ]; then
    echo "Config error: sample.csv not found: $sample_csv" >&2
    return 1
  fi

  local exports
  exports="$(python3 "${repo_dir}/scripts/config.py" load \
    --config "$config_yaml" --samples "$sample_csv" --manifest "$manifest_tsv")" || return 1

  eval "$exports"
}

# Convert an extra_args string into a Bash array using Python shlex rules.
# This preserves quoted values without using unsafe shell word splitting.
read_extra_args() {
  local value=$1
  local -n target_ref=$2
  local -a parsed=()
  if [ -n "$value" ]; then
    mapfile -t parsed < <(python3 -c 'import shlex,sys; print("\n".join(shlex.split(sys.argv[1])))' "$value") || return 1
  fi
  target_ref=("${parsed[@]}")
}

# Write the canonical 25-contig task table used by germline and somatic callers.
# The generated columns are: order, logical_name, reference_contig, length.
prepare_primary_contigs() {
  local output="${work_dir}/primary_contigs.tsv"
  python3 "${REPO_DIR}/scripts/config.py" contigs --fasta "$ref_fa" > "$output" || return 1
  if [ "$(wc -l < "$output")" -ne 25 ]; then
    echo "Error: expected 25 primary contigs in $output" >&2
    return 1
  fi
  primary_contigs_tsv=$output
  export primary_contigs_tsv
}

# Return the CPU allocation visible to this job. Explicit WGS_CPUS wins, then
# common scheduler variables, with nproc used for an unmanaged local process.
available_cpus() {
  if [ -n "${WGS_CPUS:-}" ]; then
    echo "$WGS_CPUS"
  elif [ -n "${SLURM_CPUS_PER_TASK:-}" ]; then
    echo "$SLURM_CPUS_PER_TASK"
  elif [ -n "${NSLOTS:-}" ]; then
    echo "$NSLOTS"
  elif [ -n "${PBS_NP:-}" ]; then
    echo "$PBS_NP"
  else
    nproc
  fi
}

# Convert a configured CPU fraction (for example 0.8) into an integer budget.
calculate_cpu_budget() {
  local fraction=$1
  awk -v cpus="$(available_cpus)" -v fraction="$fraction" \
    'BEGIN { value=int(cpus*fraction); print value < 1 ? 1 : value }'
}

# Cap GNU Parallel workers by the configured maximum, task count, and CPU
# budget. At least one worker is retained so small jobs can still execute.
limit_parallel_jobs() {
  local fraction=$1
  local requested_jobs=$2
  local task_count=$3
  local cpu_budget
  cpu_budget=$(calculate_cpu_budget "$fraction")
  [ "$requested_jobs" -gt "$task_count" ] && requested_jobs=$task_count
  [ "$requested_jobs" -gt "$cpu_budget" ] && requested_jobs=$cpu_budget
  [ "$requested_jobs" -lt 1 ] && requested_jobs=1
  echo "$requested_jobs"
}

# Calculate uniform threads per worker without exceeding the CPU budget.
calculate_threads() {
  local fraction=$1
  local jobs=$2
  local cpu_budget
  cpu_budget=$(calculate_cpu_budget "$fraction")
  awk -v cpus="$cpu_budget" -v jobs="$jobs" \
    'BEGIN { value=int(cpus/jobs); print value < 1 ? 1 : value }'
}

# Redirect stdout/stderr through tee into output/logs while preserving terminal
# output. Each invocation gets a timestamp/PID-specific main and Parallel log.
init_run_logging() {
  local command_name=$1
  local timestamp
  timestamp=$(date +"%Y%m%d_%H%M%S")
  pipeline_logs_dir="${work_dir}/output/logs"
  parallel_logs_dir="${pipeline_logs_dir}/parallel/${timestamp}_$$"
  mkdir -p "$parallel_logs_dir"
  run_log="${pipeline_logs_dir}/${command_name}_${timestamp}_$$.log"
  export pipeline_logs_dir parallel_logs_dir run_log
  exec > >(tee -a "$run_log") 2>&1
  echo "Run log: $run_log"
}

# Find one R1 or R2 FASTQ recursively below raw_dir. sample_prefix is matched
# literally, while common _R1/_R2 and _1/_2 naming conventions are accepted.
find_read_in_dir() {
  local sample_id=$1
  local raw_dir=$2
  local read_tag=$3
  local sample_prefix=${4:-$sample_id}
  local found=""

  if [ "$read_tag" = "R1" ]; then
    found=$(find "$raw_dir" \( -type f -o -type l \) 2>/dev/null \
      | grep -F "$sample_prefix" \
      | grep -E '(_R1|_1)([_\.-]|\.)' \
      | grep -E '\.(fq|fastq)(\.gz)?$' \
      | sort | head -n 1)
  else
    found=$(find "$raw_dir" \( -type f -o -type l \) 2>/dev/null \
      | grep -F "$sample_prefix" \
      | grep -E '(_R2|_2)([_\.-]|\.)' \
      | grep -E '\.(fq|fastq)(\.gz)?$' \
      | sort | head -n 1)
  fi
  echo "$found"
}

# Resolve sample_prefix and raw_dir from the normalized TSV manifest, then call
# find_read_in_dir. Errors if sample_id has no manifest entry.
find_sample_read() {
  local sample_id=$1
  local manifest_tsv=$2
  local read_tag=$3
  local sample_prefix raw_dir

  IFS=$'\t' read -r sample_prefix raw_dir < <(
    awk -F '\t' -v sample="$sample_id" '$1 == sample { print $2 "\t" $4; exit }' "$manifest_tsv"
  )
  if [ -z "$sample_prefix" ] || [ -z "$raw_dir" ]; then
    echo "Error: sample '${sample_id}' not found in manifest: ${manifest_tsv}" >&2
    return 1
  fi
  find_read_in_dir "$sample_id" "$raw_dir" "$read_tag" "$sample_prefix"
}

# Rebuild sample_list.txt and verify that every manifest row resolves to a pair
# of FASTQs. Reads stay in their configured raw_dir; no staging directory is
# created by the main pipeline.
prepare_sample_inputs() {
  if [ ! -f "${sample_manifest_tsv}" ]; then
    echo "Error: sample manifest not found: ${sample_manifest_tsv}"
    return 1
  fi

  mkdir -p "$work_dir"
  : > "${work_dir}/sample_list.txt"

  while IFS= read -r manifest_row; do
    local normalized_row
    normalized_row=${manifest_row//$'\t'/$'\x1f'}
    IFS=$'\x1f' read -r sample_id sample_prefix sample_ref raw_dir normal_id normal_bam <<< "$normalized_row"
    [ -z "$sample_id" ] && continue
    echo "$sample_id" >> "${work_dir}/sample_list.txt"

    local read_1 read_2
    read_1=$(find_read_in_dir "$sample_id" "$raw_dir" "R1" "$sample_prefix")
    read_2=$(find_read_in_dir "$sample_id" "$raw_dir" "R2" "$sample_prefix")
    if [ -z "$read_1" ] || [ -z "$read_2" ]; then
      echo "Error: cannot find paired reads for sample '${sample_id}' in ${raw_dir}"
      return 1
    fi

  done < "${sample_manifest_tsv}"
}

# Accept GNU Parallel either from PATH or as an executable configured path.
require_gnu_parallel() {
  if ! command -v "${parallel_bin:-parallel}" >/dev/null 2>&1 && [ ! -x "${parallel_bin:-}" ]; then
    echo "Error: GNU parallel is required but not found in PATH"
    echo "Install it and re-run (for example: conda install -c conda-forge parallel)"
    return 1
  fi
}

# Execute a one-task-per-line file with GNU Parallel.
#
# Arguments:
#   1 task file; 2 maximum workers; 3 exported worker function;
#   4 checkpoint/job-log tag; 5 GNU Parallel --halt policy (optional).
#
# Before execution, stale running markers are recovered with the matching
# cleanup_<worker> function. run_worker.sh writes a running marker before the
# tool starts and replaces it with a done marker only after a zero exit status.
parallel_run_sample_list() {
  local sample_list_file=$1
  local max_jobs=$2
  local worker_func=$3
  local job_tag=${4:-$worker_func}
  local halt_policy=${5:-soon,fail=1}
  local v=""

  require_gnu_parallel || return 1

  if [ ! -s "$sample_list_file" ]; then
    echo "Error: sample list is missing or empty: $sample_list_file"
    return 1
  fi
  if [ -z "$max_jobs" ] || [ "$max_jobs" -lt 1 ]; then
    max_jobs=1
  fi

  local parallel_log_dir="${parallel_logs_dir:-${work_dir}/output/logs/parallel}"
  mkdir -p "$parallel_log_dir"

  local checkpoint_scope checkpoint_dir cleanup_func marker task
  checkpoint_scope=$(printf '%s' "$sample_list_file" | sha256sum | awk '{print $1}')
  checkpoint_dir="${work_dir}/output/logs/checkpoints/${job_tag}/${checkpoint_scope}"
  mkdir -p "${checkpoint_dir}/running" "${checkpoint_dir}/done"
  cleanup_func="cleanup_${worker_func}"
  for marker in "${checkpoint_dir}/running/"*; do
    [ -f "$marker" ] || continue
    task=$(sed -n '2p' "$marker")
    echo "Checkpoint recovery: ${job_tag} task '${task}' was interrupted"
    if declare -F "$cleanup_func" >/dev/null; then
      "$cleanup_func" "$task" || return 1
      rm -f "$marker"
      echo "Checkpoint recovery: partial outputs cleared; task will be retried"
    else
      echo "Error: recovery function not found: $cleanup_func" >&2
      return 1
    fi
  done

  # GNU parallel starts fresh shells; export function and current scalar vars.
  export -f "$worker_func"
  export -f read_extra_args
  export -f available_cpus
  export -f calculate_cpu_budget
  export -f limit_parallel_jobs
  export -f calculate_threads
  export -f find_read_in_dir
  export -f find_sample_read
  export checkpoint_dir
  while IFS= read -r v; do
    case "$v" in
      BASHOPTS|BASHPID|BASH_ALIASES|BASH_ARGC|BASH_ARGV|BASH_CMDS|BASH_COMMAND|BASH_LINENO|BASH_SOURCE|BASH_VERSINFO|DIRSTACK|EUID|GROUPS|FUNCNAME|LINENO|PIPESTATUS|PPID|SHELLOPTS|UID)
        continue
        ;;
    esac
    export "$v" 2>/dev/null || true
  done < <(compgen -v)

  local worker_runner="${REPO_DIR}/scripts/run_worker.sh"
  if [ ! -f "$worker_runner" ]; then
    echo "Error: worker runner not found: $worker_runner" >&2
    return 1
  fi

  "${parallel_bin:-parallel}" --will-cite \
    --jobs "$max_jobs" \
    --line-buffer \
    --halt "$halt_policy" \
    --joblog "${parallel_log_dir}/${job_tag}.joblog" \
    bash "$worker_runner" "$worker_func" {} "$checkpoint_dir" :::: "$sample_list_file"
}

# Create the process lock used to prevent two top-level pipelines from sharing
# one repository. A dead PID identifies a stale lock and is recovered safely.
init_pipeline_control_files() {
  local wd=$1
  local control_dir="${wd}/.wgs_parallel_control"
  local pid_file="${control_dir}/pipeline.pid"
  local lock_file="${control_dir}/pipeline.lock"

  mkdir -p "$control_dir"
  if [ -e "$lock_file" ]; then
    local existing_pid=""
    [ -s "$pid_file" ] && existing_pid=$(cat "$pid_file")
    if [ -n "$existing_pid" ] && kill -0 "$existing_pid" 2>/dev/null; then
      echo "Error: pipeline is already running with PID ${existing_pid}"
      return 1
    fi
    echo "Checkpoint recovery: removing stale pipeline lock"
    rm -f "$lock_file" "$pid_file"
  fi

  printf 'pid=%s\nstarted=%s\nwork_dir=%s\n' "$$" "$(date -Is)" "${work_dir:-}" > "$lock_file" || return 1
  echo "$$" > "$pid_file"
}

# Remove the transient top-level PID and lock after normal exit or handled stop.
cleanup_pipeline_control_files() {
  local wd=$1
  local control_dir="${wd}/.wgs_parallel_control"
  rm -f "${control_dir}/pipeline.pid" "${control_dir}/pipeline.lock"
}

# Walk the process tree depth-first so stop can terminate grandchildren created
# by GNU Parallel, Conda, Java, GATK, samtools, and shell pipelines.
collect_child_pids_recursive() {
  local parent_pid=$1
  local child=""
  while IFS= read -r child; do
    [ -z "$child" ] && continue
    child=$(echo "$child" | tr -d '[:space:]')
    [ -z "$child" ] && continue
    collect_child_pids_recursive "$child"
    echo "$child"
  done < <(ps -o pid= --ppid "$parent_pid" 2>/dev/null)
}

# Terminate the recorded pipeline and all descendants, escalating from TERM to
# KILL only for processes that remain alive after the grace period.
stop_pipeline_process_group() {
  local wd=$1
  local control_dir="${wd}/.wgs_parallel_control"
  local pid_file="${control_dir}/pipeline.pid"
  local pid=""
  local child_pids=""

  [ -f "$pid_file" ] && pid=$(cat "$pid_file")

  if [ -z "$pid" ]; then
    echo "No running pipeline metadata found in: $control_dir"
    return 1
  fi

  if ! kill -0 "$pid" 2>/dev/null; then
    echo "Pipeline pid not running: $pid"
    return 0
  fi

  child_pids="$(collect_child_pids_recursive "$pid" | tr '\n' ' ' | xargs)"
  if [ -n "$child_pids" ]; then
    kill -TERM $child_pids 2>/dev/null || true
  fi
  kill -TERM "$pid" 2>/dev/null || true
  sleep 2
  if [ -n "$child_pids" ]; then
    kill -KILL $child_pids 2>/dev/null || true
  fi
  kill -KILL "$pid" 2>/dev/null || true
  echo "Stop signal sent to pipeline pid: $pid"
}
