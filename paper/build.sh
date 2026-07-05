#!/usr/bin/env bash
# Build paper.pdf from paper.md using pandoc + tectonic
# Run from the paper/ directory: bash build.sh
#
# Requirements (all in conda env tls_spatial):
#   conda activate tls_spatial
#   conda install -c conda-forge tectonic  # if not already installed

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

OUTPUT="paper.pdf"
LOG="build.log"

echo "Building $OUTPUT ..."

pandoc paper.md \
  --template=template.tex \
  --pdf-engine=tectonic \
  --citeproc \
  --bibliography=references.bib \
  --syntax-highlighting=tango \
  -o "$OUTPUT" \
  2>&1 | tee "$LOG"

if [ -f "$OUTPUT" ]; then
  SIZE=$(du -sh "$OUTPUT" | cut -f1)
  echo "Done: $OUTPUT ($SIZE)"
else
  echo "Build failed. Check $LOG for details."
  exit 1
fi
