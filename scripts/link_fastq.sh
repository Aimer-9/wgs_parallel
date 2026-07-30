#!/bin/bash
# Thin compatibility entry point for link_fastq.py.
#
# Usage:
#   bash scripts/link_fastq.sh RESULT_DIR RAW_LINK [--dry-run]
#
# exec replaces this shell with Python so signals and the final exit status are
# delivered directly to the caller without an extra wrapper process.
set -euo pipefail
repo_dir="$(cd "$(dirname "$0")/.." && pwd)"
exec python3 "${repo_dir}/scripts/link_fastq.py" "$@"
