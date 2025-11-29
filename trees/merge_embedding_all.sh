#!/usr/bin/env bash
# Merge embedding outputs (trees + histograms) from a fixed set of pThat bins.
#
# By default it:
#   - Assumes this script lives in:  <REPO>/trees
#   - Looks under:                   <REPO>/data/submit/<YYYY-MM-DD>/
#   - Picks the newest <YYYY-MM-DD> directory
#   - Merges from:                   job_pt*/production under that directory
#
# Usage:
#   ./merge_embedding_all.sh [BASE] [OUT_DIR] [OUT_NAME]
#
# Defaults:
#   BASE     = <REPO>/data/submit/<newest_YYYY-MM-DD>
#   OUT_DIR  = <REPO>/trees/merged_all
#   OUT_NAME = embedding_merged.root
#
# Env:
#   MIN_SIZE=5000   # ignore tiny files (bytes)
#   BATCH=200       # files per shard
#   CLEAN_SHARDS=1  # remove shard_*.root after final merge

set -euo pipefail

log() {
  printf '[%s] %s\n' "$(date '+%F %T')" "$*"
}

########################
# Locate repo + defaults
########################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="${SCRIPT_DIR}/.."
SUBMIT_DIR="${REPO_DIR}/data/submit"

# Find newest dated subdirectory in data/submit/ (by name: YYYY-MM-DD)
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
OUT_NAME="${3:-embedding_merged.root}"

MIN_SIZE="${MIN_SIZE:-5000}"
BATCH="${BATCH:-200}"
CLEAN_SHARDS="${CLEAN_SHARDS:-1}"

# Fixed list of embedding bins
BINS=( job_pt3_5 job_pt5_7 job_pt7_9 job_pt9_11 job_pt11_15 job_pt15_20
       job_pt20_25 job_pt25_30 job_pt30_40 job_pt40_50 job_pt50_-1 )

mkdir -p "$OUT_DIR"

# Detect hadd flags (old ROOT may lack -k)
HADD_FLAGS="-f"
if hadd -h 2>&1 | grep -q -- '-k'; then
  HADD_FLAGS="-fk"
fi

FINAL="$OUT_DIR/$OUT_NAME"
LIST="$OUT_DIR/inlist.txt"
: > "$LIST"

log "Merging embedding files"
log "  BASE        : $BASE"
log "  OUT_DIR     : $OUT_DIR"
log "  OUT_NAME    : $OUT_NAME"
log "  Newest run  : $latest_run"
echo

log "Scanning embedding production files in fixed pThat bins..."
found_any=0
for bin in "${BINS[@]}"; do
  prod="$BASE/$bin/production"
  if [[ ! -d "$prod" ]]; then
    log "  [skip] $bin: no production/ directory at $prod"
    continue
  fi

  # Collect ROOT files above size threshold
  has_files=0
  while IFS= read -r -d '' f; do
    sz=$(stat -c %s "$f" 2>/dev/null || echo 0)
    if [[ "$sz" -gt "$MIN_SIZE" ]]; then
      echo "$f" >> "$LIST"
      has_files=1
      found_any=1
    fi
  done < <(find "$prod" -type f -name '*.root' -print0 | sort -z)

  if [[ "$has_files" -eq 0 ]]; then
    log "  [note] $bin: no usable ROOT files"
  else
    log "  [ok]   $bin"
  fi
done

if [[ "$found_any" -eq 0 ]]; then
  log "No embedding files found in the specified bins under: $BASE"
  exit 1
fi

N=$(wc -l < "$LIST" | tr -d ' ')
log "Total files to merge: $N  (batch size: $BATCH)"

# Single file fast path
if [[ "$N" -eq 1 ]]; then
  src=$(head -n1 "$LIST")
  cp -f "$src" "$FINAL"
  rm -f "$LIST"
  log "Merged 1 file -> $FINAL"
  exit 0
fi

# Shard to avoid argv-too-long; old-ROOT friendly
rm -f "$OUT_DIR"/shard_*.root "$OUT_DIR"/inlist.batch.* || true
split -l "$BATCH" -d -a 4 "$LIST" "$OUT_DIR/inlist.batch."

shard_files=( "$OUT_DIR"/inlist.batch.* )
num_shards="${#shard_files[@]}"

log "Creating $num_shards shard(s)..."
shard_idx=0
for sub in "${shard_files[@]}"; do
  shard=$(printf "%s/shard_%04d.root" "$OUT_DIR" "$shard_idx")
  files_in_shard=$(wc -l < "$sub" | tr -d ' ')
  log "  [shard $((shard_idx+1))/$num_shards] merging $files_in_shard files -> $(basename "$shard")"
  hadd $HADD_FLAGS "$shard" $(tr '\n' ' ' < "$sub") >/dev/null
  shard_idx=$((shard_idx + 1))
done

log "Final merge of $num_shards shard(s) -> $FINAL"
hadd $HADD_FLAGS "$FINAL" "$OUT_DIR"/shard_*.root >/dev/null

# Cleanup
rm -f "$OUT_DIR"/inlist.txt "$OUT_DIR"/inlist.batch.* || true
if [[ "$CLEAN_SHARDS" -eq 1 ]]; then
  rm -f "$OUT_DIR"/shard_*.root || true
fi

log "Done. Merged $N files -> $FINAL"

