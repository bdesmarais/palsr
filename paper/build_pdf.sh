#!/usr/bin/env bash
# Build a JOSS-styled PDF (and .tex) from paper.md using the official JOSS
# (Open Journals / whedon) LaTeX template.
#
# This reproduces, locally, what the JOSS submission system does: it compiles the
# canonical Markdown source `paper.md` through the JOSS LaTeX template. JOSS itself
# ingests `paper.md` --- authors do not submit raw .tex --- but this script lets you
# preview the typeset article and produces a LaTeX version (paper.tex).
#
# Requirements: pandoc (>= 2.11) and a LaTeX engine (xelatex). The JOSS template and
# logo are vendored under paper/.joss/ (downloaded once from openjournals/whedon).
#
# Usage:  cd paper && ./build_pdf.sh
set -euo pipefail
cd "$(dirname "$0")"

# Locate pandoc: prefer one on PATH, else the RStudio/quarto bundle.
if command -v pandoc >/dev/null 2>&1; then
  PANDOC=pandoc
else
  PANDOC="/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64/pandoc"
fi

COMMON_OPTS=(
  paper.md
  --template .joss/latex.template
  --citeproc
  --bibliography paper.bib
  -V journal_name="Journal of Open Source Software"
  -V graphics=true
  -V logo_path=".joss/logo.png"
  -V year=2026
  -V citation_author="Kim et al."
  -V volume="VV" -V issue="II" -V page="NNNN"
  -V submitted="13 June 2026" -V published="13 June 2026"
  -V review_issue_url="https://github.com/openjournals/joss-reviews/issues/XXXX"
  -V repository="https://github.com/bdesmarais/palsr"
  -V archive_doi="https://doi.org/XX.XXXXX/zenodo.XXXXXXX"
  -V formatted_doi="XX.XXXXX/joss.0XXXX"
  -V link-citations=true
)

echo ">> Generating paper.tex"
"$PANDOC" "${COMMON_OPTS[@]}" --pdf-engine=xelatex -o paper.tex

echo ">> Generating paper.pdf"
"$PANDOC" "${COMMON_OPTS[@]}" --pdf-engine=xelatex -o paper.pdf

echo ">> Done: $(ls -1 paper.pdf paper.tex 2>/dev/null)"
