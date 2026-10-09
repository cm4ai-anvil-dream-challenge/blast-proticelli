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
        --docker blast-proticelli:0.1.0 \
        --seed 1 --num_inference_steps 10 \
        -o wdl_inputs.json
"""

import argparse
import json
import os
import sys

WORKFLOW = "blast_proticelli"


def resolve(path):
    """Absolute path for local files; URIs such as gs:// are left as they are."""
    return path if "://" in path else os.path.abspath(path)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--records", required=True, help="JSON from manifest_to_json.py")
    ap.add_argument("--weights_file", required=True, help="Weights zip (local path or gs:// URI)")
    ap.add_argument("--docker", required=True, help="Docker image, e.g. blast-proticelli:0.1.0")
    ap.add_argument("--cell_line", default=None)
    ap.add_argument("--seed", type=int, default=None)
    ap.add_argument("--num_inference_steps", type=int, default=None, help="Default in the WDL: 50")
    ap.add_argument("-o", "--output", default="wdl_inputs.json")
    args = ap.parse_args()

    with open(args.records) as f:
        records = json.load(f)
    if isinstance(records, dict):
        records = [records]

    inputs = {
        f"{WORKFLOW}.records": records,
        f"{WORKFLOW}.weights_file": resolve(args.weights_file),
        f"{WORKFLOW}.docker": args.docker,
    }
    optional = {"cell_line": args.cell_line, "seed": args.seed, "num_inference_steps": args.num_inference_steps}
    for key, value in optional.items():
        if value is not None:
            inputs[f"{WORKFLOW}.{key}"] = value

    missing = [p for p in [inputs[f"{WORKFLOW}.weights_file"]] if "://" not in p and not os.path.isfile(p)]
    if missing:
        sys.exit(f"ERROR: weights file not found: {missing[0]}")

    with open(args.output, "w") as f:
        json.dump(inputs, f, indent=2)
        f.write("\n")
    print(f"Wrote {args.output} ({len(records)} record(s))", file=sys.stderr)


if __name__ == "__main__":
    main()
