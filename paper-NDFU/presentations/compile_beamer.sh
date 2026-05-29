#!/bin/bash
# Compile Beamer presentation to PDF
# Usage: ./compile_beamer.sh [filename without .tex]

FILENAME=${1:-presentation-cpin2}

echo "Compiling $FILENAME.tex to PDF..."

# Run pdflatex twice to resolve references
pdflatex -interaction=nonstopmode "$FILENAME.tex" && \
pdflatex -interaction=nonstopmode "$FILENAME.tex"

# Clean up auxiliary files
rm -f "$FILENAME.aux" "$FILENAME.log" "$FILENAME.nav" "$FILENAME.out" "$FILENAME.snm" "$FILENAME.toc"

if [ -f "$FILENAME.pdf" ]; then
  echo "✓ Successfully compiled: $FILENAME.pdf"
else
  echo "✗ Compilation failed. Check $FILENAME.log for details."
  exit 1
fi
