#!/usr/bin/env bash
set -euo pipefail

########################
# Configuration
########################

# Base path as seen *inside the container*
BASE="/gpfs/mnt/gpfs01/star/pwg/svomich/JetsTrees"

SIF="${BASE}/analysis/unfolding/roounfold.sif"
MACRO="${BASE}/analysis/unfolding/unfold.cxx"

########################
# Arguments
########################

# 1st arg: input (either basename in trees/merged_all or an absolute path)
if [[ $# -ge 1 ]]; then
  if [[ "$1" = /* ]]; then
    INPUT="$1"
  else
    INPUT="${BASE}/trees/merged_all/$1"
  fi
else
  INPUT="${BASE}/trees/merged_all/embedding_merged.root"
fi

# 2nd arg: output directory (default: analysis/unfolding/out under BASE)
OUT_DIR="${2:-${BASE}/analysis/unfolding/out}"

########################
# Checks
########################

echo "----------------------------------------"
echo "Running unfolding"
echo "SIF        : $SIF"
echo "Macro      : $MACRO"
echo "Input      : $INPUT"
echo "Output dir : $OUT_DIR"
echo "----------------------------------------"

[[ -f "$SIF"   ]] || { echo "ERROR: SIF not found:   $SIF";   exit 1; }
[[ -f "$MACRO" ]] || { echo "ERROR: MACRO not found: $MACRO"; exit 1; }
[[ -f "$INPUT" ]] || { echo "ERROR: Input not found: $INPUT"; exit 1; }

mkdir -p "$OUT_DIR"

########################
# Run inside container
########################

apptainer exec -e -B /gpfs/mnt/gpfs01 \
  "$SIF" \
  root -l -b <<EOF
gSystem->Load("libRooUnfold");
.x ${MACRO}+("${INPUT}","${OUT_DIR}");
.q
EOF

echo "----------------------------------------"
echo "Done."
echo "----------------------------------------"

