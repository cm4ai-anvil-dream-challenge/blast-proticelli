#!/usr/bin/env python3
"""Build a WDL inputs file for blast_proticelli.wdl.

Wraps the records produced by manifest_to_json.py together with the
weights file, Docker image, and optional run settings, using the
"<workflow>.<input>" keys that miniWDL and Terra expect.

Local file paths are made absolute (miniWDL requires this); gs:// and other
URIs are passed through unchanged.

Usage:
    python scripts/make_wdl_inputs.py \
        --records inputs.json \
        --weights_file blast_proticelli_weights_v2.zip \
        --docker blast-proticelli:0.2.0 \
        --seed 1 --num_inference_steps 10 \
        -o wdl_inputs.json
"""

import argparse
import json
import os
import sys
from glob import glob

WORKFLOW = "blast_proticelli"


def resolve(path):
    """Absolute path for local files; URIs such as gs:// are left as they are."""
    return path if "://" in path else os.path.abspath(path)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--dir", help="Directory containing input files")
    ap.add_argument("--records", help="JSON from manifest_to_json.py")
    ap.add_argument("--weights-file", required=True, help="Weights zip (local path or gs:// URI)")
    ap.add_argument("-o", "--output", default="wdl_inputs.json")
    args = ap.parse_args()

    if args.dir:
        records = []
        for faa_path in glob(os.path.join(args.dir, "*.faa")):
            prediction_id = os.path.basename(faa_path).replace(".faa", "")

            er_path = faa_path.replace(".faa", "_er.png")
            microtubule_path = faa_path.replace(".faa", "_microtubule.png")
            nucleus_path = faa_path.replace(".faa", "_nucleus.png")
            if os.path.exists(er_path) and os.path.exists(microtubule_path) and os.path.exists(nucleus_path):
                records.append({
                    "prediction_id": prediction_id,
                    "microtubules_image": resolve(microtubule_path),
                    "er_image": resolve(er_path),
                    "nucleus_image": resolve(nucleus_path),
                    "fasta_file": resolve(faa_path),
                })
    else:
        with open(args.records) as f:
            records = json.load(f)
    if isinstance(records, dict):
        records = [records]

    inputs = {
        f"{WORKFLOW}.records": records,
        f"{WORKFLOW}.weights_file": resolve(args.weights_file)
    }

    missing = [p for p in [inputs[f"{WORKFLOW}.weights_file"]] if "://" not in p and not os.path.isfile(p)]
    if missing:
        sys.exit(f"ERROR: weights file not found: {missing[0]}")

    with open(args.output, "w") as f:
        json.dump(inputs, f, indent=2)
        f.write("\n")
    print(f"Wrote {args.output} ({len(records)} record(s))", file=sys.stderr)


if __name__ == "__main__":
    main()
