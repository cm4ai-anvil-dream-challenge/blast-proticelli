# run_inference.py
"""
BLAST + ProtiCelli baseline inference.

Reads the limited JSON produced by manifest_to_json.py (one record per
prediction), identifies each record's protein from its FASTA file with
protein_predictor.py (local BLAST + local gene map only, no network), and
runs ProtiCelli on every record. The protein used is the best-scoring BLAST
hit that is in ProtiCelli's vocabulary; the best hit overall is reported
alongside it.

Reference images are the challenge's single-channel 8-bit PNGs (512 x 512,
0 to 255). image_scaling.py normalizes them with ProtiCelli's own
ImageNormalizer and maps each prediction back to 0 to 255 (see its
docstring). FASTA files should hold a single record; if there are several, the first is
used and the count is reported. The header never feeds the
prediction; its fields are copied into the results table for our own checks.
Malformed inputs fail loudly rather than being validated and filtered here.

Everything the model needs at runtime ships in one weights zip
(--weights_zip), built with scripts/build_weights_zip.sh:
    blast_proticelli_weights.zip
        checkpoint/unet/                ProtiCelli model weights (loaded for inference)
        checkpoint/unet_ema/            ProtiCelli EMA weights (not loaded, not extracted)
        vae/                            ProtiCelli VAE weights
        reference/human_blast_db.zip    from build_human_blast_db.sh
        reference/human_gene_map.zip    from build_human_gene_map.sh

Outputs, all written flat into --output_dir:
    {prediction_id}.png        one per record that could be predicted (8-bit, 0 to 255)
    prediction_results.tsv     one row per input record, predicted or not
"""

import argparse
import csv
import hashlib
import zipfile
from importlib import resources
import json
import pickle
import sys
import uuid
from pathlib import Path

import image_scaling as isc
import protein_predictor as pp
from proticelli import Model

# Per-record outcome, recorded in the status column of prediction_results.tsv.
STATUS_PREDICTED = "predicted"
STATUS_PREDICTOR_ERROR = "predictor_error"          # BLAST or mapping raised
STATUS_NO_HIT = "no_blast_hit"                      # BLAST ran, returned nothing
STATUS_NOT_IN_VOCAB = "not_in_model_vocabulary"     # hits found, none ProtiCelli knows

# Written into weights_dir after extraction; holds the weights fingerprint.
WEIGHTS_MARKER = ".weights_fingerprint"

# The checkpoint subfolder ProtiCelli loads for inference (see load_weights).
INFERENCE_WEIGHTS_SUBDIR = "unet"

# Reference zips inside the weights zip (see module docstring).
BLAST_DB_MEMBER = "reference/human_blast_db.zip"
GENE_MAP_MEMBER = "reference/human_gene_map.zip"

# Columns carried from protein_predictor's result dict into the results TSV.
PREDICTOR_FIELDS = [
    # The hit used: best-scoring hit in ProtiCelli's vocabulary.
    "predicted_uniprot_accession",
    "predicted_ensembl_gene_id",
    "predicted_hgnc_symbol",
    "model_protein_key",
    "predicted_hit_rank",
    "percent_identity",
    "e_value",
    "alignment_length",
    "bit_score",
    # The best-scoring hit overall, in the vocabulary or not.
    "best_hit_uniprot_accession",
    "best_hit_hgnc_symbol",
    "best_hit_percent_identity",
    "best_hit_e_value",
    "best_hit_bit_score",
    "blast_task",
    "blast_targets_requested",
    "blast_hits_returned",
    "blast_seconds",
    "predictor_cached",
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
                        help="Where the weights are extracted. Reused by later runs given the same weights zip; a folder holding other weights stops the run.")
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


def weights_fingerprint(zf, members):
    """
    Fingerprint of the model weights inside the zip: a hash over each
    member's name, size and CRC-32. These come from the zip's directory, so
    this takes milliseconds even for a multi-GB zip, yet changes whenever any
    weight file's content changes.
    """
    h = hashlib.sha256()
    for info in sorted((zf.getinfo(m) for m in members), key=lambda i: i.filename):
        h.update(f"{info.filename}|{info.file_size}|{info.CRC}\n".encode())
    return h.hexdigest()


def load_weights(weights_zip, weights_dir):
    """
    Extract what inference needs from the weights zip into weights_dir:
    checkpoint/unet/, vae/, and the two reference zips. checkpoint/unet_ema/
    is left in the zip, since inference does not load it.

    Reuse across runs: after extracting, a fingerprint of the weights (see
    weights_fingerprint) is written to weights_dir/WEIGHTS_MARKER. A later
    run given the same zip and weights_dir skips the extraction. If the
    folder holds weights with a different or missing fingerprint, the run
    stops instead of using them, so a run always uses the zip it was given.

    The zip is read directly rather than through ProtiCelli's
    download_checkpoints(). Upstream's version (since 2026-08) copies the
    archive to a temporary file and extracts it once per asset, which for
    a 3.3 GB zip means two full copies and two full extractions. Upstream's
    own validators (checkpoint_ready, vae_ready) are still used to check
    the result, so "loadable" means the same thing here as upstream.
    """
    from proticelli.utils.download import checkpoint_ready, vae_ready

    weights_dir = Path(weights_dir)
    marker = weights_dir / WEIGHTS_MARKER
    model_prefixes = (f"checkpoint/{INFERENCE_WEIGHTS_SUBDIR}/", "vae/")

    with zipfile.ZipFile(weights_zip) as zf:
        names = zf.namelist()
        missing = [m for m in (BLAST_DB_MEMBER, GENE_MAP_MEMBER) if m not in names]
        if missing:
            raise FileNotFoundError(f"{weights_zip} is missing {missing}; rebuild it with build_weights_zip.sh")
        model_members = [n for n in names if n.startswith(model_prefixes)]
        fingerprint = weights_fingerprint(zf, model_members)

        existing = [d for d in ("checkpoint", "vae") if (weights_dir / d).exists()]
        if existing:
            recorded = marker.read_text().strip() if marker.is_file() else None
            if recorded != fingerprint:
                why = "come from a different weights zip" if recorded else "have no record of their source"
                raise FileExistsError(
                    f"{weights_dir} already holds weights that {why}; point --weights_dir at an empty "
                    f"folder, or delete {weights_dir}, so this run cannot use the wrong weights"
                )
            print(f"Reusing weights in {weights_dir} (match {weights_zip})")
        else:
            print(f"Extracting weights to {weights_dir}")
            zf.extractall(weights_dir, members=model_members)
            weights_dir.mkdir(parents=True, exist_ok=True)
            marker.write_text(fingerprint + "\n")

        # The reference zips are small (a few MB); always refresh them.
        zf.extractall(weights_dir, members=[BLAST_DB_MEMBER, GENE_MAP_MEMBER])

    # ProtiCelli's Model._load_model() falls back to a randomly initialized
    # model, with no error, when it finds no weights. Predictions would still
    # be written and reported as successful. Require the folders inference
    # actually loads: checkpoint/unet (Model.model calls _load_model() with
    # the default use_ema=False) and vae/.
    if not checkpoint_ready(weights_dir):
        raise FileNotFoundError(
            f"{weights_dir}/checkpoint/{INFERENCE_WEIGHTS_SUBDIR} is missing or incomplete; ProtiCelli "
            f"would silently use random weights. Check the layout of {weights_zip}"
        )
    if not vae_ready(weights_dir):
        raise FileNotFoundError(f"{weights_dir}/vae is missing or incomplete; check the layout of {weights_zip}")

    paths = {"checkpoint_dir": str(weights_dir / "checkpoint"), "vae_dir": str(weights_dir / "vae")}
    return paths, weights_dir / BLAST_DB_MEMBER, weights_dir / GENE_MAP_MEMBER


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


def classify(prediction):
    """Map a protein_predictor result to one of the STATUS_* values."""
    if prediction.get("error"):
        return STATUS_PREDICTOR_ERROR
    if not prediction.get("best_hit_uniprot_accession"):
        return STATUS_NO_HIT
    if not prediction.get("model_protein_key"):
        return STATUS_NOT_IN_VOCAB
    return STATUS_PREDICTED


# FASTA header fields copied into the results table (reporting only).
HEADER_FIELDS = {
    "fasta_hpa_antibody_id": "hpa_antibody_id",
    "fasta_ensembl_protein_id": "ensembl_protein_id",
    "fasta_ensembl_release": "ensembl",
}


def header_columns(fasta_path):
    header = pp.read_fasta_header(fasta_path)
    fields = pp.parse_header_fields(header)
    return {
        "fasta_record_count": pp.count_fasta_records(fasta_path),
        "fasta_header": header,
        **{col: fields.get(key, "") for col, key in HEADER_FIELDS.items()},
    }


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
        vocabulary = pp.Vocabulary(pickle.load(f).keys())
    print(f"Model vocabulary: {len(vocabulary.keys)} entries from {antibody_map_path}")

    # ---- 2. Protein identification, one record at a time ----
    results_rows = []
    to_predict = []  # (record, model_protein_key, results_row)
    for i, rec in enumerate(records, 1):
        seq = pp.read_fasta(rec["fasta_file"])
        prediction = pp.predict_for_sequence(seq, vocabulary=vocabulary)
        status = classify(prediction)
        print(f"[{i}/{len(records)}] {rec['prediction_id']}: best {prediction.get('best_hit_hgnc_symbol')}, "
              f"used {prediction.get('model_protein_key')} (rank {prediction.get('predicted_hit_rank')}) -> {status}")

        row = {
            "prediction_id": rec["prediction_id"],
            "status": status,
            "output_file": "",
            **{k: prediction.get(k) for k in PREDICTOR_FIELDS},
            "predictor_error": prediction.get("error", ""),
            **header_columns(rec["fasta_file"]),
            # Filled in for predicted records (see step 3).
            "nucleus_max": "",
            "gain_microtubules": "",
            "gain_nucleus": "",
            "gain_er": "",
            "prediction_max": "",
            "output_max": "",
            "cell_line": args.cell_line or "",
            "seed": args.seed,
            "num_inference_steps": args.num_inference_steps,
            "run_id": run_id,
        }
        results_rows.append(row)
        if status == STATUS_PREDICTED:
            to_predict.append((rec, prediction["model_protein_key"], row))

    results_path = output_dir / "prediction_results.tsv"

    # Every record failing inside the predictor points at the environment
    # (missing blastp, bad reference zip), not at the data. Fail the run.
    if all(r["status"] == STATUS_PREDICTOR_ERROR for r in results_rows):
        write_results(results_rows, results_path)
        raise RuntimeError("Protein prediction failed for every record; check the BLAST DB and gene map")

    # ---- 3. ProtiCelli, only on records it can handle ----
    if to_predict:
        # Read the PNGs and normalize with ProtiCelli's ImageNormalizer,
        # keeping each image's nucleus maximum to reconstruct the output.
        images, scales = [], []
        for rec, _, row in to_predict:
            raw = isc.load_reference_channels(rec["microtubules_image"], rec["nucleus_image"], rec["er_image"])
            reference, gains, nucleus_max = isc.normalize_reference(raw)
            images.append(reference)
            scales.append(nucleus_max)
            row["nucleus_max"] = nucleus_max
            row["gain_microtubules"], row["gain_nucleus"], row["gain_er"] = (
                round(float(gains[c]), 4) for c in isc.REFERENCE_CHANNELS
            )

        model = Model(**checkpoint_paths)
        # ProtiCelli picks CUDA, then Apple MPS, then CPU. Log the choice so a
        # silent fall back to the CPU shows in every run log.
        print(f"ProtiCelli device: {model.device}")
        output = model.predict(
            images=images,
            protein_names=[key for _, key, _ in to_predict],
            cell_line_names=[args.cell_line] * len(to_predict) if args.cell_line else None,
            num_inference_steps=args.num_inference_steps,
            batch_size=args.batch_size,
            seed=args.seed,
        )
        if len(output.images) != len(to_predict):
            raise RuntimeError(f"Output count mismatch: expected {len(to_predict)}, got {len(output.images)}")

        for (rec, _, row), pred, nucleus_max in zip(to_predict, output.images, scales):
            out = isc.reconstruct_prediction(pred, nucleus_max)
            save_path = output_dir / f"{rec['prediction_id']}.png"
            isc.write_png(save_path, out)
            row["prediction_max"] = round(float(pred.max()), 4)
            row["output_max"] = int(out.max())
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