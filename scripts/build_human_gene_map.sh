#!/usr/bin/env bash
# build_human_gene_map.sh
#
# Downloads accession -> Ensembl cross-reference data for the reviewed
# human UniProt proteome, and compacts it into a simple two-column
# accession -> gene_id(s) table for use as LOCAL_GENE_MAP with
# protein_predictor.py's MAPPING_BACKEND=local. Replaces per-sequence
# mygene.info calls with an in-memory dict lookup.
#
# Requests xref_ensembl_full rather than xref_ensembl: the short form
# omits Ensembl gene ids (ENSG), listing only transcripts and isoforms.
# The parser extracts every ENSG token by regex rather than assuming an
# exact field format.
#
# The UniProt release is taken from the X-UniProt-Release and
# X-UniProt-Release-Date response headers of the download itself.
#
# Usage:
#   ./build_human_gene_map.sh [output_dir]

set -euo pipefail

OUT_DIR="${1:-./human_gene_map}"
mkdir -p "$OUT_DIR"

RAW_TSV_GZ="$OUT_DIR/human_reviewed_ensembl_raw.tsv.gz"
RAW_TSV="$OUT_DIR/human_reviewed_ensembl_raw.tsv"
MAP_TSV="$OUT_DIR/human_gene_map.tsv"
HEADERS="$OUT_DIR/download_headers.txt"

echo "Downloading reviewed human UniProt accession -> symbol -> Ensembl cross-references..."
curl -sS -D "$HEADERS" -o "$RAW_TSV_GZ" \
  "https://rest.uniprot.org/uniprotkb/stream?query=%28organism_id%3A9606%29%20AND%20%28reviewed%3Atrue%29&fields=accession%2Cgene_primary%2Cxref_ensembl_full&format=tsv&compressed=true"

echo "Decompressing..."
gunzip -f "$RAW_TSV_GZ"

ROW_COUNT=$(($(wc -l < "$RAW_TSV") - 1))
echo "Downloaded $ROW_COUNT accession rows to $RAW_TSV"

if [ "$ROW_COUNT" -lt 1000 ]; then
  echo "WARNING: expected roughly 20000 rows (same ballpark as the BLAST db build), got $ROW_COUNT." >&2
  echo "This usually means the query or fields string needs adjusting, check the file contents before proceeding." >&2
fi

echo "Building compact accession -> symbol -> Ensembl gene id map at $MAP_TSV ..."
python3 - "$RAW_TSV" "$MAP_TSV" << 'PYEOF'
import sys
import re
import csv

raw_path, map_path = sys.argv[1], sys.argv[2]

with open(raw_path, newline="") as fin, open(map_path, "w", newline="") as fout:
    reader = csv.reader(fin, delimiter="\t")
    writer = csv.writer(fout, delimiter="\t")
    next(reader, None)  # header, whatever UniProt actually names these columns
    writer.writerow(["accession", "hgnc_symbol", "ensembl_gene_ids"])

    rows_with_gene = 0
    for row in reader:
        if len(row) < 3:
            continue
        accession, symbol, ensembl_blob = row[0], row[1], row[2]
        # Extract every ENSG-prefixed token regardless of surrounding
        # punctuation, robust to formatting details this script's author
        # couldn't verify against a live response.
        gene_ids = sorted(set(re.findall(r"ENSG\d+", ensembl_blob)))
        if gene_ids:
            rows_with_gene += 1
        writer.writerow([accession, symbol.strip(), ";".join(gene_ids)])

    print(f"  {rows_with_gene} of the parsed rows had at least one ENSG id extracted")
PYEOF

PROVENANCE_FILE="$OUT_DIR/PROVENANCE.txt"
{
  echo "Reference data: reviewed human UniProt accession -> symbol -> Ensembl gene id map"
  echo "Generated (UTC): $(date -u +"%Y-%m-%dT%H:%M:%SZ")"
  echo "Source query URL: https://rest.uniprot.org/uniprotkb/stream?query=%28organism_id%3A9606%29%20AND%20%28reviewed%3Atrue%29&fields=accession%2Cgene_primary%2Cxref_ensembl_full&format=tsv&compressed=true"
  echo "Accession rows downloaded: $ROW_COUNT"
  grep -i '^x-uniprot-release' "$HEADERS" | tr -d '\r' || echo "UniProt release: not reported in response headers"
} > "$PROVENANCE_FILE"
echo "Wrote provenance to $PROVENANCE_FILE"

ZIP_FILE="$OUT_DIR/human_gene_map.zip"
(cd "$OUT_DIR" && zip -q "$(basename "$ZIP_FILE")" human_gene_map.tsv PROVENANCE.txt)
echo "Wrote $ZIP_FILE"

echo ""
echo "Done. Package $ZIP_FILE into the weights file with scripts/build_weights_zip.sh."