#!/usr/bin/env python3
"""Build a WDL inputs file for blast_proticelli.wdl from a folder of input files.

Every {id}.faa in --dir becomes one record when its three reference images
are also present:

    {id}.faa
    {id}_microtubules.png
    {id}_nucleus.png
    {id}_er.png

A .faa without all three images is not a complete record. It is left out,
and its ID and missing files are listed on stderr, so nothing is dropped
silently. Images without a .faa are ignored.

Local file paths are made absolute (miniWDL requires this); gs:// and other
URIs are passed through unchanged.

Usage:
    python scripts/make_wdl_inputs.py \
        --dir test_inputs \
        --weights-file blast_proticelli_weights_v2.zip \
        -o wdl_inputs.json

--uri-prefix gs://BUCKET/path writes each input file as that prefix plus its
file name, for a run on AnVIL once the folder's files are copied there.

--wdl submission_test writes inputs for submission_test.wdl (the AnVIL test
of submission.wdl) instead of blast_proticelli.wdl (the default).
"""

import argparse
import json
import os
import sys
from glob import glob

# Input names for each WDL: (records array, weights file).
INPUT_NAMES = {
    "blast_proticelli": ("records", "weights_file"),
    "submission_test": ("inputs", "weights"),
}


def resolve(path):
    """Absolute path for local files; URIs such as gs:// are left as they are."""
    return path if "://" in path else os.path.abspath(path)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--dir", required=True, help="Directory containing input files")
    ap.add_argument("--weights-file", required=True, help="Weights zip (local path or gs:// URI)")
    ap.add_argument("-o", "--output", default="wdl_inputs.json")
    ap.add_argument("--wdl", choices=sorted(INPUT_NAMES), default="blast_proticelli",
                    help="Which WDL the inputs are for (default: blast_proticelli)")
    ap.add_argument("--uri-prefix", default=None,
                    help="Write input files as PREFIX/<file name> (e.g. gs://bucket/test_cells) instead of local paths")
    args = ap.parse_args()

    records = []
    incomplete = []
    for faa_path in sorted(glob(os.path.join(args.dir, "*.faa"))):
        prediction_id = os.path.basename(faa_path)[: -len(".faa")]
        images = {
            "microtubules_image": faa_path[: -len(".faa")] + "_microtubules.png",
            "er_image": faa_path[: -len(".faa")] + "_er.png",
            "nucleus_image": faa_path[: -len(".faa")] + "_nucleus.png",
        }
        missing = [os.path.basename(p) for p in images.values() if not os.path.exists(p)]
        if missing:
            incomplete.append((prediction_id, missing))
            continue
        def locate(path):
            if args.uri_prefix:
                return args.uri_prefix.rstrip("/") + "/" + os.path.basename(path)
            return resolve(path)

        records.append({
            "prediction_id": prediction_id,
            **{key: locate(path) for key, path in images.items()},
            "fasta_file": locate(faa_path),
        })

    if incomplete:
        print(f"WARNING: {len(incomplete)} .faa file(s) left out, missing reference images:", file=sys.stderr)
        for prediction_id, missing in incomplete:
            print(f"  {prediction_id}: {', '.join(missing)}", file=sys.stderr)
    if not records:
        sys.exit(f"ERROR: no complete records in {args.dir}")

    records_name, weights_name = INPUT_NAMES[args.wdl]
    weights_file = resolve(args.weights_file)
    inputs = {
        f"{args.wdl}.{records_name}": records,
        f"{args.wdl}.{weights_name}": weights_file,
    }

    missing = [p for p in [weights_file] if "://" not in p and not os.path.isfile(p)]
    if missing:
        sys.exit(f"ERROR: weights file not found: {missing[0]}")

    with open(args.output, "w") as f:
        json.dump(inputs, f, indent=2)
        f.write("\n")
    print(f"Wrote {args.output} ({len(records)} record(s))", file=sys.stderr)


if __name__ == "__main__":
    main()
