#!/usr/bin/env bash

set -uo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN_PATH="${ROOT_DIR}/build/ntmc"
INPUT_PATH="${ROOT_DIR}/input_cross_100nm.toml"
OUT_ROOT="${ROOT_DIR}/test_threads/cross_100nm_test_threads"
SUMMARY_CSV="${OUT_ROOT}/run_summary.csv"
THREADS=(1 2 4 8 16 32 64)

if [[ ! -x "${BIN_PATH}" ]]; then
    echo "Error: binary not found or not executable: ${BIN_PATH}" >&2
    exit 1
fi

if [[ ! -f "${INPUT_PATH}" ]]; then
    echo "Error: input file not found: ${INPUT_PATH}" >&2
    exit 1
fi

mkdir -p "${OUT_ROOT}"
echo "threads,exit_code,elapsed_seconds,results_folder,log_file" > "${SUMMARY_CSV}"

for n in "${THREADS[@]}"; do
    run_base="${OUT_ROOT}/thread_${n}"
    tmp_input="${OUT_ROOT}/input_cross_100nm_thread_${n}.toml"
    log_file="${OUT_ROOT}/thread_${n}.log"

    if grep -qE '^[[:space:]]*output_folder[[:space:]]*=' "${INPUT_PATH}"; then
        sed -E "s|^[[:space:]]*output_folder[[:space:]]*=.*$|output_folder = \"${run_base}/\"|" \
            "${INPUT_PATH}" > "${tmp_input}"
    else
        cp "${INPUT_PATH}" "${tmp_input}"
        cat >> "${tmp_input}" <<EOF

[io]
output_folder = "${run_base}/"
EOF
    fi

    start_ts="$(date +%s)"
    OMP_NUM_THREADS="${n}" "${BIN_PATH}" "${tmp_input}" > "${log_file}" 2>&1
    exit_code=$?
    end_ts="$(date +%s)"
    elapsed="$((end_ts - start_ts))"

    results_folder="$(grep -m1 '^Results folder:' "${log_file}" | sed 's/^Results folder:[[:space:]]*//')"
    printf "%s,%s,%s,%s,%s\n" "${n}" "${exit_code}" "${elapsed}" "${results_folder}" "${log_file}" >> "${SUMMARY_CSV}"

    if [[ "${exit_code}" -eq 0 ]]; then
        echo "[OK] OMP_NUM_THREADS=${n}, elapsed=${elapsed}s, results=${results_folder}"
    else
        echo "[FAIL] OMP_NUM_THREADS=${n}, exit=${exit_code}, log=${log_file}" >&2
    fi
done

echo "All runs completed. Summary: ${SUMMARY_CSV}"
