#!/usr/bin/env python3
"""Convert a challenge manifest TSV into the limited JSON consumed by the
BLAST + ProtiCelli inference workflow.

Each manifest row becomes one JSON object, keyed by prediction_id: the
identifier for one prediction (one reference image set paired with one
FASTA). It must be unique within a run and is used as the base of the
output filename. Only the fields inference needs
are carried over; provenance, diagnostic, and ground-truth columns
(poi_image, poi_ensembl_gene_id) are dropped, so the same tool is safe to
run against the firewalled test workspace.

Relative paths in the manifest are joined onto --path-prefix, which should
point at where the files actually live: an absolute local directory for
miniWDL, or a gs:// bucket path for Terra. Cromwell localizes to
/cromwell_root on its own, so that prefix should not be written here.

Usage:
    python manifest_to_json.py manifest.tsv -o inputs.json
    python manifest_to_json.py manifest.tsv -o inputs.json \
        --path-prefix gs://my-bucket/cm4ai
"""

import argparse
import csv
import json
import posixpath
import sys
from pathlib import Path

# JSON key -> manifest column. Edit here if either side is renamed.
FIELD_MAP = {
    "prediction_id": "prediction_id",
    "microtubules_image": "microtubules_image",
    "er_image": "er_image",
    "nucleus_image": "nucleus_image",
    "fasta_file": "aa_sequence_fasta_path",
}

# Keys whose values are file paths (and so get the prefix applied).
PATH_KEYS = {"microtubules_image", "er_image", "nucleus_image", "fasta_file"}


def resolve_path(value, prefix):
    """Join a manifest path onto the prefix unless it is already absolute or a URI."""
    if "://" in value or value.startswith("/"):
        return value
    if prefix is None:
        return value
    return posixpath.join(prefix, value)


def convert(manifest_path, prefix):
    with open(manifest_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing_cols = [c for c in FIELD_MAP.values() if c not in reader.fieldnames]
        if missing_cols:
            sys.exit(f"ERROR: manifest is missing required columns: {missing_cols}")

        records = []
        for line_no, row in enumerate(reader, start=2):  # header is line 1
            record = {}
            for key, col in FIELD_MAP.items():
                value = (row[col] or "").strip()
                if not value:
                    sys.exit(f"ERROR: line {line_no}: empty value in column '{col}'")
                record[key] = resolve_path(value, prefix) if key in PATH_KEYS else value
            records.append(record)

    ids = [r["prediction_id"] for r in records]
    dupes = sorted({i for i in ids if ids.count(i) > 1})
    if dupes:
        sys.exit(f"ERROR: duplicate prediction_id values: {dupes}")
    return records


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("manifest", type=Path, help="Manifest TSV")
    ap.add_argument("-o", "--output", type=Path, help="Output JSON (default: stdout)")
    ap.add_argument(
        "--path-prefix",
        help="Prefix for relative paths (default: absolute directory of the manifest)",
    )
    args = ap.parse_args()

    prefix = args.path_prefix or str(args.manifest.resolve().parent)
    records = convert(args.manifest, prefix)

    text = json.dumps(records, indent=2) + "\n"
    if args.output:
        args.output.write_text(text)
        print(f"Wrote {len(records)} record(s) to {args.output}", file=sys.stderr)
    else:
        sys.stdout.write(text)


if __name__ == "__main__":
    main()