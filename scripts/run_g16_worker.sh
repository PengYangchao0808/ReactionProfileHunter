#!/bin/bash
# ==========================================================================
# Gaussian Worker Wrapper (self-cleaning)
# ==========================================================================
# Usage: ./run_g16_worker.sh <input_file> <output_file>
#
# Goals:
# - Use a job-specific scratch directory
# - Always cleanup scratch on success/failure/signals
# - Refuse to start if scratch filesystem is low on space
# - Preserve existing GAUSS_PROFILE_PATH/GAUSS_EXEDIR discovery
# ==========================================================================

set -e

if [ $# -lt 2 ]; then
    echo "Usage: $0 <input_file> <output_file>"
    exit 1
fi

INPUT_FILE="$1"
OUTPUT_FILE="$2"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ------------------------------
# 1) Scratch setup (isolation)
# ------------------------------

# If parent process already provided GAUSS_SCRDIR (e.g. sandbox scratch),
# create a per-job subdir so we can safely clean only our own files.
if [ -n "$GAUSS_SCRDIR" ]; then
    JOB_SCRDIR="$GAUSS_SCRDIR/gau_${USER}_$$"
else
    JOB_SCRDIR="/tmp/gau_${USER}_$$"
fi

mkdir -p "$JOB_SCRDIR"
export GAUSS_SCRDIR="$JOB_SCRDIR"

cleanup() {
    echo ">> [Wrapper] Cleaning up scratch: $GAUSS_SCRDIR" >&2
    rm -rf "$GAUSS_SCRDIR" || true
}

trap cleanup EXIT INT TERM

# Disk guard: require at least 2GB free on scratch filesystem
FREE_KB="$(df -k "$GAUSS_SCRDIR" | awk 'NR==2 {print $4}')"
if [ -n "$FREE_KB" ] && [ "$FREE_KB" -lt 2000000 ]; then
    echo "!! [Wrapper] CRITICAL: Low disk space for scratch ($GAUSS_SCRDIR). Aborting." >&2
    exit 1
fi

# ------------------------------
# 2) Environment setup
# ------------------------------

if [ -z "$GAUSS_PROFILE_PATH" ]; then
    GAUSS_PROFILE_PATHS=(
        "/root/g16/g16.profile"
        "/opt/g16/g16.profile"
        "/usr/local/g16/g16.profile"
        "$SCRIPT_DIR/../../g16/g16.profile"
    )

    for PROFILE_PATH in "${GAUSS_PROFILE_PATHS[@]}"; do
        if [ -f "$PROFILE_PATH" ]; then
            GAUSS_PROFILE_PATH="$PROFILE_PATH"
            break
        fi
    done

    if [ -z "$GAUSS_PROFILE_PATH" ]; then
        GAUSS_ROOT_PATHS=(
            "/root/g16"
            "/opt/g16"
            "/usr/local/g16"
        )

        for ROOT_PATH in "${GAUSS_ROOT_PATHS[@]}"; do
            if [ -d "$ROOT_PATH" ]; then
                GAUSS_EXEDIR="$ROOT_PATH"
                export GAUSS_EXEDIR
                export PATH="$ROOT_PATH:$PATH"
                break
            fi
        done
    fi
fi

if [ -f "$GAUSS_PROFILE_PATH" ]; then
    # shellcheck disable=SC1090
    source "$GAUSS_PROFILE_PATH"
fi

ulimit -s unlimited 2>&1 || true

# ------------------------------
# 3) Run Gaussian
# ------------------------------

if command -v g16 >/dev/null 2>&1; then
    g16 < "$INPUT_FILE" > "$OUTPUT_FILE"
elif [ -n "$GAUSS_EXEDIR" ] && [ -f "$GAUSS_EXEDIR/g16" ]; then
    "$GAUSS_EXEDIR/g16" < "$INPUT_FILE" > "$OUTPUT_FILE"
else
    echo "ERROR: g16 not found in PATH and GAUSS_EXEDIR not set" >&2
    exit 1
fi
