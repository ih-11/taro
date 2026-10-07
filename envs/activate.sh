#!/usr/bin/env bash
# source envs/activate.sh  — activate env, set project paths
#
# Two halves:
#   TARO_CODE   ~/Code/taro                          git repo (ext4)
#   TARO_STORE  /mnt/f/RA/Downstream/Project5_TARO   heavy files (F drive)
#
# /tmp is tmpfs (RAM-backed, 45 GB). Large sorts there would consume memory
# and fail confusingly, so TMPDIR is pinned to an ext4 scratch directory.

export TARO_CODE="$HOME/Code/taro"
export TARO_STORE=/mnt/f/RA/Downstream/Project5_TARO

export TARO_REFS="$TARO_STORE/refs"
export TARO_DATA="$TARO_STORE/data"
export TARO_WORK="$TARO_STORE/work"
export TARO_LOGS="$TARO_STORE/logs"
export TARO_RES="$TARO_CODE/results"

export TMPDIR="$TARO_CODE/.tmp"
export TARO_THREADS=24          # of 28, leaves headroom

mkdir -p "$TMPDIR"

if [ "${CONDA_DEFAULT_ENV:-}" != "taro" ]; then
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate taro
fi

echo "taro env active"
echo "  code    $TARO_CODE"
echo "  store   $TARO_STORE"
echo "  tmpdir  $TMPDIR"
echo "  threads $TARO_THREADS"
