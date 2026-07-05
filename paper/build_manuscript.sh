#!/usr/bin/env bash
# Build manuscript.pdf from manuscript.md using pandoc + tectonic
# Run from the paper/ directory: bash build_manuscript.sh
#
# Requirements:
#   conda activate tls_spatial
#   conda install -c conda-forge tectonic  # if not already installed
#
# Before building, copy required figures to paper/figures/:
#   cp ../outputs/stage2_clustering/umap_725.png figures/
#   cp ../outputs/validation_tcga/TCGA_BRCA_TLS_survival.png figures/
#   cp ../outputs/validation_tcga/METABRIC_TLS_survival.png figures/
#   cp ../outputs/validation_geo/shiao_TLS_ICB.png figures/
#   (fig1_overview, fig2_clustering, fig3_coherence, fig4_survival, fig5_shiao
#    are composite figures that need assembly — see manuscript.md figure captions)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

OUTPUT="manuscript.pdf"
LOG="manuscript_build.log"

echo "Building $OUTPUT ..."

pandoc manuscript.md \
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
