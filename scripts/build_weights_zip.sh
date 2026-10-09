#!/usr/bin/env bash
# Assemble the BLAST + ProtiCelli weights zip consumed by run_inference.py.
#
# Layout produced:
#   <out_zip>
#       checkpoint/unet/                ProtiCelli model weights (loaded for inference)
#       checkpoint/unet_ema/            ProtiCelli EMA weights (kept pending use_ema question)
#       vae/                            ProtiCelli VAE weights
#       reference/human_blast_db.zip    from build_human_blast_db.sh
#       reference/human_gene_map.zip    from build_human_gene_map.sh
#
# Usage:
#   scripts/build_weights_zip.sh <proticelli_dir> <blast_db_zip> <gene_map_zip> <out_zip>
#
# <proticelli_dir> is the directory holding the extracted checkpoint/ and vae/
# folders (the proticelli package directory after download_checkpoints()).
# The two reference zips are renamed to the fixed names run_inference.py expects.
#
# Only the folders listed in INCLUDE are packaged. Training state in
# checkpoint/ (optimizer.bin, scheduler.bin, random_states_*.pkl) is left
# out: inference does not read it, and optimizer.bin alone is ~3.5 GB.
# An include list rather than an exclude list means new files in a future
# ProtiCelli release cannot end up in the zip unnoticed.

set -euo pipefail

if [[ $# -ne 4 ]]; then
    sed -n '2,16p' "$0" | sed 's/^# \{0,1\}//'
    exit 1
fi

PROTICELLI_DIR=$1
BLAST_DB_ZIP=$2
GENE_MAP_ZIP=$3
OUT_ZIP=$(cd "$(dirname "$4")" && pwd)/$(basename "$4")

INCLUDE=(checkpoint/unet checkpoint/unet_ema vae)

for d in "${INCLUDE[@]}"; do
    [[ -d "$PROTICELLI_DIR/$d" ]] || { echo "ERROR: $PROTICELLI_DIR/$d not found" >&2; exit 1; }
done
for f in "$BLAST_DB_ZIP" "$GENE_MAP_ZIP"; do
    [[ -f "$f" ]] || { echo "ERROR: $f not found" >&2; exit 1; }
done
[[ ! -e "$OUT_ZIP" ]] || { echo "ERROR: $OUT_ZIP already exists; remove it first" >&2; exit 1; }

STAGE=$(mktemp -d)
trap 'rm -rf "$STAGE"' EXIT
mkdir -p "$STAGE/reference"
cp "$BLAST_DB_ZIP" "$STAGE/reference/human_blast_db.zip"
cp "$GENE_MAP_ZIP" "$STAGE/reference/human_gene_map.zip"

# -X drops macOS extended attributes; -x skips Finder metadata files.
(cd "$PROTICELLI_DIR" && zip -r -X -q "$OUT_ZIP" "${INCLUDE[@]}" -x '*.DS_Store')
(cd "$STAGE" && zip -r -X -q "$OUT_ZIP" reference)

echo "Wrote $OUT_ZIP"
unzip -l "$OUT_ZIP" | awk 'NR>3 && $4 ~ /^(checkpoint\/[^\/]+|vae|reference)\/[^\/]*\/?$/ {printf "  %12s  %s\n", $1, $4}'
shasum -a 256 "$OUT_ZIP"