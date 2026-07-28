#!/bin/bash
set -u

worker_func=${1:-}
task=${2:-}
checkpoint_dir=${3:-}
if [ -z "$worker_func" ]; then
  echo "Error: worker function name is empty" >&2
  exit 2
fi
if ! declare -F "$worker_func" >/dev/null; then
  echo "Error: worker function is not available: $worker_func" >&2
  exit 2
fi

if [ -z "$checkpoint_dir" ]; then
  "$worker_func" "$task"
  exit $?
fi

checkpoint_key=$(printf '%s' "$task" | sha256sum | awk '{print $1}')
running_marker="${checkpoint_dir}/running/${checkpoint_key}"
done_marker="${checkpoint_dir}/done/${checkpoint_key}"
mkdir -p "${checkpoint_dir}/running" "${checkpoint_dir}/done"
printf 'started=%s\n%s\n' "$(date -Is)" "$task" > "$running_marker"

if "$worker_func" "$task"; then
  printf 'completed=%s\n%s\n' "$(date -Is)" "$task" > "$done_marker"
  rm -f "$running_marker"
  exit 0
fi
exit 1
