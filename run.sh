#!/bin/bash

export JAVA_HOME=/home/jp2045/software/jdk-23.0.1
export PATH=/home/jp2045/software/jdk-23.0.1/bin:$PATH

mkdir -p data/tmp
export NXF_TEMP="$PWD/data/tmp"

retry=0
max_attempts=5
args=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --retry)
      retry=1
      shift
      ;;
    --retry=*)
      retry=1
      max_attempts="${1#--retry=}"
      shift
      ;;
    *)
      args+=("$1")
      shift
      ;;
  esac
done

run_nextflow() {
  ./nextflow -log /dev/null run main.nf -resume "${args[@]}"
}

if [[ "$retry" -eq 0 ]]; then
  run_nextflow
  exit $?
fi

attempt=1
while true; do
  run_nextflow
  status=$?
  if [[ $status -eq 0 ]]; then
    exit 0
  fi
  # 130 = Ctrl-C, 143 = SIGTERM. Leave an intentional stop alone.
  if [[ $status -eq 130 || $status -eq 143 ]]; then
    exit "$status"
  fi
  if [[ $attempt -ge $max_attempts ]]; then
    echo "Nextflow exited ${status} after ${attempt} attempts; giving up." >&2
    exit "$status"
  fi
  echo "Nextflow exited ${status}; resuming (attempt $((attempt + 1))/${max_attempts}) after 60s." >&2
  attempt=$((attempt + 1))
  sleep 60
done
