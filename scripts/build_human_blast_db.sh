#!/usr/bin/env bash
# build_human_blast_db.sh
#
# Downloads the reviewed (Swiss-Prot) human proteome from UniProt's REST
# API and builds a local blastp database from it, for use as
# LOCAL_BLAST_DB with protein_predictor.py's SEARCH_BACKEND=local.
#
# The UniProt release is taken from the X-UniProt-Release and
# X-UniProt-Release-Date response headers of the download itself, so the
# recorded release is the one the data actually came from.
#
# Usage:
#   ./build_human_blast_db.sh [output_dir]
#
# Requires: curl, gunzip, makeblastdb (from ncbi-blast+, e.g.
#   conda install -c bioconda blast   or   brew install blast)

set -euo pipefail

OUT_DIR="${1:-./human_blast_db}"
mkdir -p "$OUT_DIR"

FASTA_GZ="$OUT_DIR/human_reviewed.fasta.gz"
FASTA="$OUT_DIR/human_reviewed.fasta"
DB_PREFIX="$OUT_DIR/human_reviewed_swissprot"
HEADERS="$OUT_DIR/download_headers.txt"

echo "Downloading reviewed human proteome from UniProt..."
curl -sS -D "$HEADERS" -o "$FASTA_GZ" \
  "https://rest.uniprot.org/uniprotkb/stream?query=%28organism_id%3A9606%29%20AND%20%28reviewed%3Atrue%29&format=fasta&compressed=true"

echo "Decompressing..."
gunzip -f "$FASTA_GZ"

SEQ_COUNT=$(grep -c "^>" "$FASTA" || true)
echo "Downloaded $SEQ_COUNT sequences to $FASTA"

if [ "$SEQ_COUNT" -lt 1000 ]; then
  echo "WARNING: expected roughly 20000 reviewed human sequences, got $SEQ_COUNT." >&2
  echo "This usually means the query string needs adjusting, check the file contents before proceeding." >&2
fi

echo "Building blastp database at $DB_PREFIX ..."
makeblastdb -in "$FASTA" -dbtype prot -out "$DB_PREFIX"

PROVENANCE_FILE="$OUT_DIR/PROVENANCE.txt"
{
  echo "Reference database: reviewed human UniProtKB Swiss-Prot BLAST database"
  echo "Generated (UTC): $(date -u +"%Y-%m-%dT%H:%M:%SZ")"
  echo "Source query URL: https://rest.uniprot.org/uniprotkb/stream?query=%28organism_id%3A9606%29%20AND%20%28reviewed%3Atrue%29&format=fasta&compressed=true"
  echo "Sequences downloaded: $SEQ_COUNT"
  grep -i '^x-uniprot-release' "$HEADERS" | tr -d '\r' || echo "UniProt release: not reported in response headers"
  echo ""
  echo "blastp version:"
  blastp -version
  echo ""
  echo "makeblastdb version:"
  makeblastdb -version
} > "$PROVENANCE_FILE"
echo "Wrote provenance to $PROVENANCE_FILE"

ZIP_FILE="$OUT_DIR/human_blast_db.zip"
(cd "$OUT_DIR" && zip -q "$(basename "$ZIP_FILE")" "$(basename "$DB_PREFIX")".* PROVENANCE.txt)
echo "Wrote $ZIP_FILE"

echo ""
echo "Done. Package $ZIP_FILE into the weights file with scripts/build_weights_zip.sh."