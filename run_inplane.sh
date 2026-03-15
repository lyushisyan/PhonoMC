#!/usr/bin/env bash

set -uo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN_PATH="${ROOT_DIR}/build/ntmc"
INPUT_DIR="${ROOT_DIR}/calculation-inplane"

Z_VALUES=(10 20 50 100 200 500 1000 2000 5000)
ROUGH_LABELS=(r0nm r01nm r1nm)

if [[ ! -x "${BIN_PATH}" ]]; then
    echo "Error: binary not found or not executable: ${BIN_PATH}" >&2
    exit 1
fi

threads="${OMP_NUM_THREADS:-36}"
echo "Using OMP_NUM_THREADS=${threads}"

total=0
failed=0

for z in "${Z_VALUES[@]}"; do
    for r in "${ROUGH_LABELS[@]}"; do
        input_file="${INPUT_DIR}/input_inplane_x1000nm_z${z}nm_${r}.toml"
        total=$((total + 1))

        if [[ ! -f "${input_file}" ]]; then
            echo "[MISS] ${input_file}" >&2
            failed=$((failed + 1))
            continue
        fi

        echo "[RUN ] z=${z}nm, rough=${r}, input=$(basename "${input_file}")"
        start_ts="$(date +%s)"
        OMP_NUM_THREADS="${threads}" "${BIN_PATH}" "${input_file}"
        rc=$?
        end_ts="$(date +%s)"
        elapsed=$((end_ts - start_ts))

        if [[ "${rc}" -eq 0 ]]; then
            echo "[ OK ] z=${z}nm, rough=${r}, elapsed=${elapsed}s"
        else
            echo "[FAIL] z=${z}nm, rough=${r}, exit=${rc}, elapsed=${elapsed}s" >&2
            failed=$((failed + 1))
        fi
    done
done

echo "Completed ${total} runs, failed ${failed}."
if [[ "${failed}" -ne 0 ]]; then
    exit 1
fi
