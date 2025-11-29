#!/usr/bin/env bash
# Merge EMBEDDING outputs per pThat bin (one file per job_ptX_Y).
# Reads:  <BASE>/job_ptX_Y/production/*.root
# Writes: <OUT_DIR>/embedding_ptX_Y.root
#
# Default behavior:
#   - Assumes this script lives in:  <REPO>/trees
#   - Looks under:                   <REPO>/data/submit/<YYYY-MM-DD>/
#   - Picks the newest <YYYY-MM-DD> directory
#
# Usage:
#   ./merge_pThatbins.sh [BASE] [OUT_DIR]
#
# Defaults:
#   BASE    = <REPO>/data/submit/<newest_YYYY-MM-DD>
#   OUT_DIR = <REPO>/trees/merged_all
#
# Env:
#   MIN_SIZE=5000    # ignore tiny files (bytes)
#   BATCH=200        # files per shard
#   CLEAN_SHARDS=1   # remove shard_*.root after final merge

set -euo pipefail

log() {
  printf '[%s] %s\n' "$(date '+%F %T')" "$*"
}

########################
# Locate repo + defaults
########################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
SUBMIT_DIR="${REPO_DIR}/data/submit"

latest_run=""
if [[ -d "$SUBMIT_DIR" ]]; then
  latest_run="$(find "$SUBMIT_DIR" -maxdepth 1 -mindepth 1 -type d -printf '%f\n' \
                | sort | tail -n 1 || true)"
fi

if [[ -z "$latest_run" ]]; then
  echo "ERROR: No run directories found under: $SUBMIT_DIR"
  echo "       Expected something like: $SUBMIT_DIR/2025-11-28"
  exit 1
fi

DEFAULT_BASE="${SUBMIT_DIR}/${latest_run}"
DEFAULT_OUT_DIR="${SCRIPT_DIR}/merged_all"

########################
# Arguments
########################

BASE="${1:-$DEFAULT_BASE}"
OUT_DIR="${2:-$DEFAULT_OUT_DIR}"

MIN_SIZE="${MIN_SIZE:-5000}"
BATCH="${BATCH:-200}"
CLEAN_SHARDS="${CLEAN_SHARDS:-1}"

# Fixed list of EMBEDDING bins (no data)
BINS=( job_pt3_5 job_pt5_7 job_pt7_9 job_pt9_11 job_pt11_15 job_pt15_20
       job_pt20_25 job_pt25_30 job_pt30_40 job_pt40_50 job_pt50_-1 )

mkdir -p "$OUT_DIR"

# Detect hadd flags (older ROOT may not support -k)
HADD_FLAGS="-f"
if hadd -h 2>&1 | grep -q -- '-k'; then
  HADD_FLAGS="-fk"
fi

merge_one_bin () {
  local bin="$1"
  local prod="$BASE/$bin/production"
  local tag="${bin#job_}"                 # e.g. pt3_5
  local final="$OUT_DIR/embedding_${tag}.root"

  if [[ ! -d "$prod" ]]; then
    log "  [skip] $bin: no production/ directory at $prod"
    return 0
  fi

  local tmpdir="$OUT_DIR/.${bin}_tmp"
  mkdir -p "$tmpdir"
  local list="$tmpdir/inlist.txt"
  : > "$list"

  # Collect ROOT files above size threshold (sorted)
  local has=0
  while IFS= read -r -d '' f; do
    local sz
    sz=$(stat -c %s "$f" 2>/dev/null || echo 0)
    if [[ "$sz" -gt "$MIN_SIZE" ]]; then
      echo "$f" >> "$list"
      has=1
    fi
  done < <(find "$prod" -type f -name '*.root' -print0 | sort -z)

  if [[ "$has" -eq 0 ]]; then
    log "  [note] $bin: no usable ROOT files"
    rm -rf "$tmpdir"
    return 0
  fi

  local N
  N=$(wc -l < "$list" | tr -d ' ')

  if [[ "$N" -eq 1 ]]; then
    cp -f "$(head -n1 "$list")" "$final"
    log "  [ok]   $bin -> $(basename "$final") (1 file)"
    rm -rf "$tmpdir"
    return 0
  fi

  # Shard to avoid argv-too-long; old-ROOT friendly
  rm -f "$tmpdir"/shard_*.root "$tmpdir"/inlist.batch.* || true
  split -l "$BATCH" -d -a 4 "$list" "$tmpdir/inlist.batch."

  local shard_files=( "$tmpdir"/inlist.batch.* )
  local num_shards="${#shard_files[@]}"

  log "  [bin $bin] merging $N files in $num_shards shard(s)..."
  local idx=0
  for sub in "${shard_files[@]}"; do
    local shard
    shard=$(printf "%s/shard_%04d.root" "$tmpdir" "$idx")
    local files_in_shard
    files_in_shard=$(wc -l < "$sub" | tr -d ' ')
    log "    [shard $((idx+1))/$num_shards] $files_in_shard files"
    hadd $HADD_FLAGS "$shard" $(tr '\n' ' ' < "$sub") >/dev/null
    idx=$((idx+1))
  done

  # Final per-bin merge (trees + histos)
  hadd $HADD_FLAGS "$final" "$tmpdir"/shard_*.root >/dev/null

  # Cleanup
  if [[ "$CLEAN_SHARDS" -eq 1 ]]; then
    rm -rf "$tmpdir"
  else
    rm -f "$tmpdir"/inlist.txt "$tmpdir"/inlist.batch.* || true
  fi

  log "  [ok]   $bin -> $(basename "$final") ($N files)"
}

log "Merging per pThat bin"
log "  BASE       : $BASE"
log "  OUT_DIR    : $OUT_DIR"
log "  Newest run : $latest_run"
echo

for bin in "${BINS[@]}"; do
  merge_one_bin "$bin"
done

log "Done."
