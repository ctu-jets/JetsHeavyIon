#!/bin/bash

# Stop on errors
set -e

# Directory where *this* script lives
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Default file names (can be overridden via arguments)
IN_BASENAME="${1:-embedding_merged.root}"
OUT_BASENAME="${2:-hists.root}"

# Paths relative to the script location
INPUT="${SCRIPT_DIR}/../../trees/merged_all/${IN_BASENAME}"
OUTPUT="${SCRIPT_DIR}/${OUT_BASENAME}"
MACRO="${SCRIPT_DIR}/make_hists.C"

echo "----------------------------------------"
echo "Running histogram production"
echo "Script dir : $SCRIPT_DIR"
echo "Input      : $INPUT"
echo "Output     : $OUTPUT"
echo "Macro      : $MACRO"
echo "----------------------------------------"

root -l -b -q "$MACRO(\"$INPUT\",\"$OUTPUT\")"

echo "----------------------------------------"
echo "Done."
echo "----------------------------------------"
