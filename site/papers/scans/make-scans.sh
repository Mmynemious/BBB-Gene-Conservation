#!/usr/bin/env bash
# Convert an annotated PDF scan into page images for a paper page.
#
#   ./make-scans.sh <scan.pdf> <paper-slug>
#   ./make-scans.sh ~/Downloads/"Adobe Scan Sep 13, 2026.pdf" kumabe2025
#
# Writes site/papers/scans/<paper-slug>/p-01.jpg, p-02.jpg, ...
# Needs poppler:  brew install poppler

set -euo pipefail

if [ $# -ne 2 ]; then
  echo "usage: $0 <scan.pdf> <paper-slug>" >&2
  exit 1
fi

PDF="$1"
SLUG="$2"
DIR="$(cd "$(dirname "$0")" && pwd)/$SLUG"

[ -f "$PDF" ] || { echo "no such file: $PDF" >&2; exit 1; }

mkdir -p "$DIR"
rm -f "$DIR"/p-*.jpg

# 110 dpi keeps handwriting legible while staying light enough to scroll.
pdftoppm -jpeg -r 110 -jpegopt quality=72 "$PDF" "$DIR/p"

# pdftoppm numbers as p-1.jpg / p-01.jpg depending on page count; normalise to 2 digits.
for f in "$DIR"/p-*.jpg; do
  n=$(basename "$f" .jpg | sed 's/^p-//')
  printf -v padded "%02d" "$((10#$n))"
  [ "$n" = "$padded" ] || mv "$f" "$DIR/p-$padded.jpg"
done

echo "wrote $(ls "$DIR"/p-*.jpg | wc -l | tr -d ' ') pages to $DIR"
du -sh "$DIR"
