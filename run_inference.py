# run_inference.py
"""
BLAST + ProtiCelli baseline inference.

Reads the limited JSON produced by manifest_to_json.py (one record per
prediction), identifies each record's protein from its FASTA file with
protein_predictor.py (local BLAST + local gene map only, no network), and
runs ProtiCelli on every record whose predicted HGNC symbol is in the
model's vocabulary.

Inputs are assumed to already meet ProtiCelli's requirements (normalized
single-channel images, matching H/W) and protein_predictor's requirements
(single-record FASTA). Malformed inputs fail loudly rather than being
validated and filtered here.

Everything the model needs at runtime ships in one weights zip
(--weights_zip), built with scripts/build_weights_zip.sh:
    blast_proticelli_weights.zip
        checkpoint/unet/                ProtiCelli model weights (loaded for inference)
        checkpoint/unet_ema/            ProtiCelli EMA weights (not currently loaded)
        vae/                            ProtiCelli VAE weights
        reference/human_blast_db.zip    from build_human_blast_db.sh
        reference/human_gene_map.zip    from build_human_gene_map.sh

Outputs, all written flat into --output_dir:
    {prediction_id}.tiff       one per record that could be predicted
    prediction_results.tsv     one row per input record, predicted or not
"""

import argparse
import csv
import zipfile
from importlib import resources
import json
import pickle
import sys
import uuid
from pathlib import Path

from skimage.io import imsave

import pipeline_checks as check
import protein_predictor as pp
from proticelli import Model

# Per-record outcome, recorded in the status column of prediction_results.tsv.
STATUS_PREDICTED = "predicted"
STATUS_PREDICTOR_ERROR = "predictor_error"          # BLAST or mapping raised
STATUS_NO_HIT = "no_blast_hit"                      # BLAST ran, returned nothing
STATUS_NO_SYMBOL = "no_hgnc_symbol"                 # hit found, no HGNC symbol mapped
STATUS_NOT_IN_VOCAB = "not_in_model_vocabulary"     # symbol unknown to ProtiCelli

# The checkpoint subfolder ProtiCelli loads for inference (see load_weights).
INFERENCE_WEIGHTS_SUBDIR = "unet"

# Reference zips inside the weights zip (see module docstring).
BLAST_DB_MEMBER = "reference/human_blast_db.zip"
GENE_MAP_MEMBER = "reference/human_gene_map.zip"

# Columns carried from protein_predictor's result dict into the results TSV.
PREDICTOR_FIELDS = [
    "predicted_uniprot_accession",
    "predicted_ensembl_gene_id",
    "predicted_hgnc_symbol",
    "percent_identity",
    "e_value",
    "alignment_length",
    "gene_candidate_count",
    "mapping_ambiguous",
    "search_database",
    "predictor_pipeline_version",
]


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    parser.add_argument("--inputs", required=True, help="JSON from manifest_to_json.py (list of records, or one record)")
    parser.add_argument("--cell_line", default=None,
                        help="Optional ProtiCelli cell line name applied to every record. Omit to run unconditioned on cell line.")
    parser.add_argument("--weights_zip", required=True, help="Weights zip with checkpoint/, vae/, and reference/ (see module docstring)")
    parser.add_argument("--weights_dir", default="./weights",
                        help="Where the weights are extracted. Must not already hold checkpoint/ or vae/ from another release.")
    parser.add_argument("--antibody_map", default=None,
                        help="Override path to antibody_map.pkl. Default: proticelli/data/ in the installed package.")
    parser.add_argument("--reference_data_dir", default="./reference_data", help="Where the BLAST DB and gene map are extracted")
    parser.add_argument("--num_inference_steps", type=int, default=50)
    parser.add_argument("--batch_size", type=int, default=4)
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("--output_dir", required=True)
    return parser.parse_args()


def load_records(path):
    """Load the input JSON. Accepts a list of records or a single record."""
    data = json.loads(Path(path).read_text())
    records = [data] if isinstance(data, dict) else data
    if not records:
        raise ValueError(f"{path} contains no records")
    ids = [r["prediction_id"] for r in records]
    dupes = sorted({i for i in ids if ids.count(i) > 1})
    if dupes:
        raise ValueError(f"Duplicate prediction_id values would overwrite each other: {dupes}")
    return records


def configure_protein_predictor(blast_db_zip, gene_map_zip, reference_data_dir):
    """
    Force protein_predictor onto its local backends and point it at the
    reference zips from the weights file. Set on the module directly
    (rather than via environment variables) so this run's configuration is
    explicit and cannot silently fall back to the EBI or mygene.info
    network backends.
    """
    pp.SEARCH_BACKEND = "local"
    pp.MAPPING_BACKEND = "local"
    pp.BLAST_DB_ZIP = str(blast_db_zip)
    pp.GENE_MAP_ZIP = str(gene_map_zip)
    pp.REFERENCE_DATA_DIR = reference_data_dir
    # Clear any LOCAL_BLAST_DB / LOCAL_GENE_MAP picked up from the shell
    # environment, so the paths are always derived from the weights zip.
    # Otherwise an exported variable from development work silently wins.
    pp.LOCAL_BLAST_DB = ""
    pp.LOCAL_GENE_MAP = ""
    pp.ensure_reference_data_ready()


def load_weights(weights_zip, weights_dir):
    """
    Extract the weights zip: checkpoint/ and vae/ via ProtiCelli's own
    loader, and the two reference zips directly.

    download_checkpoints() takes URLs, so the local path is converted to a
    file:// URL here, keeping every user-facing argument a plain path. The
    same URL is passed for both checkpoint and vae: the first extraction
    places both folders, so the second is skipped.

    Extraction goes to a run-local weights_dir rather than the proticelli
    package directory. download_checkpoints() skips extraction whenever
    checkpoint/ already exists, so extracting into the package would
    silently reuse whatever weights were already installed and ignore the
    zip this run was given.
    """
    weights_url = Path(weights_zip).resolve().as_uri()
    paths = Model.download_checkpoints(
        dest_dir=weights_dir,
        checkpoint_url=weights_url,
        vae_url=weights_url,
    )
    for key, d in paths.items():
        if not Path(d).is_dir():
            raise FileNotFoundError(f"{key} not found at {d}; check the layout of {weights_zip}")

    # ProtiCelli's Model._load_model() falls back to a randomly initialized
    # model, with no error, when it finds no weights. Predictions would still
    # be written and reported as successful. Require the folder inference
    # actually loads: checkpoint/unet (Model.model calls _load_model() with
    # the default use_ema=False, so checkpoint/unet_ema is not used).
    model_weights = Path(paths["checkpoint_dir"]) / INFERENCE_WEIGHTS_SUBDIR
    if not model_weights.is_dir():
        raise FileNotFoundError(
            f"{model_weights} not found; ProtiCelli would silently use random weights. "
            f"Check the layout of {weights_zip}"
        )

    # The reference zips are pulled out as-is; protein_predictor extracts them.
    with zipfile.ZipFile(weights_zip) as zf:
        missing = [m for m in (BLAST_DB_MEMBER, GENE_MAP_MEMBER) if m not in zf.namelist()]
        if missing:
            raise FileNotFoundError(f"{weights_zip} is missing {missing}; rebuild it with build_weights_zip.sh")
        for member in (BLAST_DB_MEMBER, GENE_MAP_MEMBER):
            zf.extract(member, weights_dir)

    return paths, Path(weights_dir) / BLAST_DB_MEMBER, Path(weights_dir) / GENE_MAP_MEMBER


def find_antibody_map(explicit_path=None):
    """
    Locate ProtiCelli's antibody_map.pkl (HGNC symbol -> embedding index).
    It ships inside the installed proticelli package (proticelli/data/),
    so it is found relative to the package rather than a hardcoded path.
    That resolves the same way on a laptop and inside the Docker image.
    """
    if explicit_path:
        path = Path(explicit_path)
    else:
        path = Path(str(resources.files("proticelli") / "data" / "antibody_map.pkl"))
    if not path.is_file():
        raise FileNotFoundError(f"antibody_map.pkl not found at {path}; pass --antibody_map explicitly")
    return path


def classify(prediction, vocabulary):
    """Map a protein_predictor result to one of the STATUS_* values."""
    if prediction.get("error"):
        return STATUS_PREDICTOR_ERROR
    if not prediction.get("predicted_uniprot_accession"):
        return STATUS_NO_HIT
    symbol = prediction.get("predicted_hgnc_symbol")
    if not symbol:
        return STATUS_NO_SYMBOL
    if symbol not in vocabulary:
        return STATUS_NOT_IN_VOCAB
    return STATUS_PREDICTED


def write_results(rows, path):
    fieldnames = list(rows[0].keys())
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote results: {path}")


def main():
    args = parse_args()

    # Location management belongs to the calling workflow, so outputs are
    # written flat into output_dir. run_id is kept only as a provenance tag.
    run_id = uuid.uuid4().hex[:8]
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    records = load_records(args.inputs)
    print(f"{len(records)} record(s) loaded from {args.inputs}")

    # ---- 1. Weights, reference data, and model vocabulary ----
    checkpoint_paths, blast_db_zip, gene_map_zip = load_weights(args.weights_zip, args.weights_dir)
    configure_protein_predictor(blast_db_zip, gene_map_zip, args.reference_data_dir)
    antibody_map_path = find_antibody_map(args.antibody_map)
    with open(antibody_map_path, "rb") as f:
        vocabulary = pickle.load(f)
    print(f"Model vocabulary: {len(vocabulary)} symbols from {antibody_map_path}")

    # ---- 2. Protein identification, one record at a time ----
    results_rows = []
    to_predict = []  # (record, hgnc_symbol, results_row)
    for i, rec in enumerate(records, 1):
        seq = pp.read_fasta(rec["fasta_file"])
        prediction = pp.predict_for_sequence(seq)
        status = classify(prediction, vocabulary)
        print(f"[{i}/{len(records)}] {rec['prediction_id']}: "
              f"{prediction.get('predicted_hgnc_symbol')} -> {status}")

        row = {
            "prediction_id": rec["prediction_id"],
            "status": status,
            "output_file": "",
            **{k: prediction.get(k) for k in PREDICTOR_FIELDS},
            "predictor_error": prediction.get("error", ""),
            "cell_line": args.cell_line or "",
            "seed": args.seed,
            "num_inference_steps": args.num_inference_steps,
            "run_id": run_id,
        }
        results_rows.append(row)
        if status == STATUS_PREDICTED:
            to_predict.append((rec, prediction["predicted_hgnc_symbol"], row))

    results_path = output_dir / "prediction_results.tsv"

    # Every record failing inside the predictor points at the environment
    # (missing blastp, bad reference zip), not at the data. Fail the run.
    if all(r["status"] == STATUS_PREDICTOR_ERROR for r in results_rows):
        write_results(results_rows, results_path)
        raise RuntimeError("Protein prediction failed for every record; check the BLAST DB and gene map")

    # ---- 3. ProtiCelli, only on records it can handle ----
    if to_predict:
        images = [
            check.assemble_channels([rec["microtubules_image"], rec["nucleus_image"], rec["er_image"]])
            for rec, _, _ in to_predict
        ]
        model = Model(**checkpoint_paths)
        output = model.predict(
            images=images,
            protein_names=[symbol for _, symbol, _ in to_predict],
            cell_line_names=[args.cell_line] * len(to_predict) if args.cell_line else None,
            num_inference_steps=args.num_inference_steps,
            batch_size=args.batch_size,
            seed=args.seed,
        )
        if len(output.images) != len(to_predict):
            raise RuntimeError(f"Output count mismatch: expected {len(to_predict)}, got {len(output.images)}")

        for (rec, _, row), img in zip(to_predict, output.images):
            save_path = output_dir / f"{rec['prediction_id']}.tiff"
            imsave(save_path, img)
            row["output_file"] = save_path.name
            print(f"Saved: {save_path}")
    else:
        print("No records were in the model vocabulary; no images generated", file=sys.stderr)

    write_results(results_rows, results_path)
    counts = {}
    for r in results_rows:
        counts[r["status"]] = counts.get(r["status"], 0) + 1
    print(f"Summary: {counts}")


if __name__ == "__main__":
    main()