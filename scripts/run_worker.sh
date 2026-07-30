#!/bin/bash
# Execute one exported Bash worker function inside a GNU Parallel child shell.
#
# Arguments:
#   1 worker function name; 2 opaque task string; 3 checkpoint directory.
#
# 00_util.sh exports the worker function and scalar environment before GNU
# Parallel launches this script. The task is hashed only for a filesystem-safe
# marker name; its original text is stored on marker line 2 for recovery. A
# failed or killed worker intentionally leaves its running marker behind. The
# next pipeline invocation calls cleanup_<worker> before retrying that task.
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

# Direct/manual use can omit checkpoint_dir, but normal pipeline work always
# supplies it so incomplete outputs can be recovered deterministically.
if [ -z "$checkpoint_dir" ]; then
  "$worker_func" "$task"
  exit $?
fi

checkpoint_key=$(printf '%s' "$task" | sha256sum | awk '{print $1}')
running_marker="${checkpoint_dir}/running/${checkpoint_key}"
done_marker="${checkpoint_dir}/done/${checkpoint_key}"
mkdir -p "${checkpoint_dir}/running" "${checkpoint_dir}/done"
printf 'started=%s\n%s\n' "$(date -Is)" "$task" > "$running_marker"

# Write done before deleting running so there is never a state with neither
# marker after successful work. A crash in between is recovered conservatively.
if "$worker_func" "$task"; then
  printf 'completed=%s\n%s\n' "$(date -Is)" "$task" > "$done_marker"
  rm -f "$running_marker"
  exit 0
fi
exit 1
