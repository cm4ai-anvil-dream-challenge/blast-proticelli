# BLAST-ProtiCelli baseline

A baseline submission for the CM4AI DREAM Challenge on predicting protein subcellular localization. Given reference images of a cell (microtubules, nucleus, ER) and the amino acid sequence of a protein, it generates a predicted fluorescence image of where that protein localizes.


## How to run evaluator


Get the code
```
git clone https://github.com/cm4ai-anvil-dream-challenge/blast-proticelli.git
```

### Pull a prebuilt image
```
docker pull quay.io/ohsu-comp-bio/blast-proticelli-example:v0.2.0
```

If you want to try on a Mac laptop (notice: this will be slow)
```
docker pull --platform linux/amd64 quay.io/ohsu-comp-bio/blast-proticelli-example:v0.2.0
```

### Build the image

The image holds code only (Python 3.11, BLAST+ 2.17.0, ProtiCelli at the pinned commit, and the pipeline scripts in `/app`); weights and reference data arrive at runtime in the weights file. Build for x86_64, which Terra uses, including on Apple Silicon:

```bash
docker buildx build --platform linux/amd64 -t blast-proticelli:0.2.0 --load .
```


### Build an input file
```
./scripts/make_wdl_inputs.py --dir ../test_data --weights-file ../blast_proticelli_weights_v2.zip
```

Each prediction is one record: three single-channel reference images plus one FASTA file holding a single amino acid sequence. If a FASTA file holds several sequences, the first is used with a warning, and `fasta_record_count` in the results shows how many there were. (Six files in the 2026-10 sample set list several Ensembl proteins for one antibody; challenge FASTAs are to hold one each.)

```json
[
  {
    "prediction_id": "cell_0",
    "microtubules_image": "/path/to/cell_0_microtubules.png",
    "er_image": "/path/to/cell_0_er.png",
    "nucleus_image": "/path/to/cell_0_nucleus.png",
    "fasta_file": "/path/to/cell_0.faa"
  }
]
```

`prediction_id` must be unique within a run; it becomes the output filename. The same reference images may appear in several records paired with different FASTA files.

Reference images must be single-channel 8-bit PNGs, 512 x 512, as in the challenge specification; anything else stops the run with an error. Each image set is normalized with ProtiCelli's `ImageNormalizer` (see [Output intensity](#output-intensity)). No resampling happens in this pipeline. The FASTA header never feeds the prediction: the baseline identifies the protein from the sequence alone. Its `key=value` fields (`hpa_antibody_id`, `ensembl_protein_id`, `ensembl`) are copied into the results table for our own checks, and are left blank when absent.



### Install MiniWDL
```
pip install miniwdl
```

### Run miniwdl with the weights
```
miniwdl run ./blast_proticelli.wdl -i wdl_inputs.json
```


### Weights file

Everything the model needs at runtime ships in a single zip:

```
blast_proticelli_weights_v2.zip
    checkpoint/unet/                ProtiCelli model weights (loaded for inference)
    checkpoint/unet_ema/            ProtiCelli EMA weights (not loaded or extracted; see open items)
    vae/                            ProtiCelli VAE weights
    reference/human_blast_db.zip    BLAST database
    reference/human_gene_map.zip    UniProt to Ensembl and HGNC map
```

ProtiCelli's training state (`optimizer.bin`, `scheduler.bin`, `random_states_*.pkl`) is not included. Inference does not read it, and leaving it out removes about 3.5 GB; outputs with and without it differ only by run-to-run GPU noise.

`run_inference.py` extracts only `checkpoint/unet/`, `vae/`, and the reference zips, and checks the model folders with ProtiCelli's own validators (`checkpoint_ready`, `vae_ready`) before predicting. Without it, ProtiCelli silently falls back to a randomly initialized model and still writes images, so a malformed weights file would otherwise produce plausible-looking but meaningless predictions.

The BLAST database and gene map are packaged with the model weights because they define which proteins this baseline can recognize, in the same way ProtiCelli's weights and vocabulary do. Each reference zip contains a `PROVENANCE.txt` recording its UniProt query, build date, and, for zips built after 2026-09-28, the UniProt release it was downloaded from.

The weights file is not stored in git. See [Weights builds](#weights-builds) for released versions.

## Building the weights file

1. Build the reference zips with `training/build_human_blast_db.sh` and `training/build_human_gene_map.sh` (usage in each script's header). Both pull from the current UniProt release, so rebuilding later produces different content; keep built zips rather than relying on rebuilding them.
2. Obtain the ProtiCelli `checkpoint/` and `vae/` folders, for example by running ProtiCelli's `download_checkpoints()`, which fetches them from the Lundberg lab. The build script packages only `checkpoint/unet`, `checkpoint/unet_ema`, and `vae`, and stops if any is missing.
3. Assemble the weights file:

```bash
training/build_weights_zip.sh proticelli \
  human_blast_db/human_blast_db.zip \
  human_gene_map/human_gene_map.zip \
  blast_proticelli_weights_v2.zip
```

The script prints the top-level contents and the SHA-256 checksum. Record every build in the table below.

## Weights builds

| Version | Built | UniProt release (date) | ProtiCelli weights | SHA-256 | Location | Notes |
|---|---|---|---|---|---|---|
| v1 | 2026-09-28 | 2026_03 (2026-09-02), inferred |  Lundberg lab default `download_checkpoints()` URLs | `19fc815a99449e80eb2c703602aa58a5c92103225efb8203dca738b9cd9c0566` | local only (not yet uploaded) | Verified end to end on two test records. Includes ProtiCelli training state files (about 6.5 GB). UniProt release inferred: reference zips built 2026-09-23, when 2026_03 was current; not recorded by the v1 build scripts. |
| v2 | 2026-09-28 | 2026_03 (2026-09-02) | Same as v1 | `2f243c8e0945b5b9830fafae1590c92638415501e7aa6c6db71da3161ed07071` | local only (not yet uploaded) | Training state removed, bringing the file to about 3.6 GB. Same model weights and reference zips as v1; outputs match v1 within MPS run-to-run noise. Verified in Docker. |

## Protein search settings

`protein_predictor.py` runs `blastp` with `-evalue 1000` and `-max_target_seqs 25`, ranks the hits by bit score itself (ties broken by e-value), and walks down the list to the first hit whose gene symbol is in ProtiCelli's vocabulary, matching names the way ProtiCelli does (an exact entry, or a symbol listed inside a combined entry). If none of the 25 qualifies, it searches again with 500, then with every protein in the database (25,000, more than the roughly 20,400 it holds). Each wider search runs only if the previous one came back full. If still nothing is usable, the fallback protein is used. Asking for one target is not the same as asking for the best one: with `-max_target_seqs 1`, BLAST can return the first good hit it finds ([Shah et al. 2019, Bioinformatics 35:1613](https://doi.org/10.1093/bioinformatics/bty833)). Sequences under 30 residues use `-task blastp-short`. The limits can be changed with `BLAST_EVALUE`, `BLAST_MAX_TARGET_SEQS`, `BLAST_MAX_TARGET_SEQS_RETRY` and `BLAST_MAX_TARGET_SEQS_ALL`.

Results are cached by sequence, so a sequence that appears in many records is searched once per run.


## Running inference

```bash
python run_inference.py \
  --inputs inputs.json \
  --weights_zip blast_proticelli_weights_v2.zip \
  --output_dir predictions/test_run
```


### Debugging in Docker


The build fails if the ProtiCelli vocabulary files are missing from the installed package, if `blastp` cannot run, or if the pipeline imports fail.

To test locally, mount the working directory at the same path so absolute paths in `inputs.json` still resolve, and extract weights inside the container:

```bash
docker run --rm --platform linux/amd64 -v "$PWD":"$PWD" -w "$PWD" blast-proticelli:0.2.0 \
  python /app/run_inference.py --inputs inputs.json \
  --weights_zip blast_proticelli_weights_v2.zip \
  --weights_dir /tmp/weights --reference_data_dir /tmp/reference_data \
  --output_dir predictions/docker
```

On a Mac this runs under emulation on the CPU and is slow; `--num_inference_steps 10` is useful for a quick functional test, though fewer steps give visibly rougher images.

### With WDL

`blast_proticelli.wdl` runs all records in one task, so the model loads once; splitting records across tasks is left to the caller. Records are declared as a WDL struct, so every image and FASTA file is a `File` that the workflow engine copies into the task.

| Input | Type | Default |
|---|---|---|
| `records` | `Array[PredictionInput]` (the records described above) | required |
| `weights_file` | `File` | required |
| `docker` | `String` | required |
| `cell_line` | `String?` | none |
| `seed` | `Int?` | none |
| `num_inference_steps` | `Int` | 50 |
| `batch_size` | `Int` | 4 |

Outputs are `predictions` (`Array[File]`, the PNGs) and `prediction_results` (the TSV). Task resources default to Terra's g2-standard-8 (8 CPUs, 32 GiB, 75 GB SSD, one NVIDIA L4) and can be overridden as task inputs. Each run first logs `nvidia-smi` and whether PyTorch can use CUDA, so a silent fall back to the CPU shows at the top of the task log.

To run locally with miniWDL:

```bash
python scripts/manifest_to_json.py manifest.tsv -o inputs.json
python scripts/make_wdl_inputs.py \
  --records inputs.json \
  --weights_file blast_proticelli_weights_v2.zip \
  --docker blast-proticelli:0.2.0 \
  -o wdl_inputs.json
MINIWDL__FILE_IO__OUTPUT_HARDLINKS=true miniwdl run blast_proticelli.wdl -i wdl_inputs.json
```

The environment variable makes miniWDL write outputs as ordinary files rather than symlinks; `_LAST/outputs.json` lists where they are. For Terra, generate the records with `--path-prefix gs://<bucket>/<folder>` and pass the weights file and Docker image by their bucket and registry addresses.

For a step-by-step run of the two-record example, see [docs/WORKED_EXAMPLE.md](docs/WORKED_EXAMPLE.md).

## Outputs

All outputs are written flat into `--output_dir`:

- `{prediction_id}.png`, one per record that could be predicted: 8-bit, 512 x 512, 0 to 255
- `prediction_results.tsv`, one row per input record, predicted or not

### Output intensity

ProtiCelli's `ImageNormalizer` scales every channel of an image against that image's nucleus maximum (the nucleus value at the 99.95th percentile, which also clips all channels), then maps to [-1, 1]. ProtiCelli's decoder returns predictions on the scale before that last step, from 0 to 1.3, where 1.0 means as bright as the nucleus maximum. The output PNG is therefore

    output = prediction x nucleus_max, clipped to 0 to 255 and rounded

which puts the prediction back in the input's units. Values above 1.0 saturate at 255 when the nucleus maximum is 255. ProtiCelli's own 8-bit export (`save_prediction`) uses a fixed display scale instead (prediction / 1.3 x 255), which does not reconstruct input units.

This inversion is exact for proteins at least 30% as bright as the nucleus. Dimmer channels were brightened in training (`ImageNormalizer`'s `dim_gamma`), so a dim protein comes back too bright; see open items. On the 200 challenge sample cells, using each cell's real protein channel as a stand-in prediction, the round trip is within 1 count for 198 cells; the 2 dimmer cells return 63 for a true maximum of 52, and 75 for 71.

### Status values

Every record gets an image. The status says which protein ProtiCelli was asked for:

| Status | Protein used |
|---|---|
| `predicted` | The highest-ranked BLAST hit in ProtiCelli's vocabulary, however weak |
| `fallback_no_vocabulary_hit` | Fallback protein: BLAST found hits, but none ProtiCelli was trained on, even searching every protein in the database |
| `fallback_no_blast_hit` | Fallback protein: BLAST found nothing at all |
| `fallback_predictor_error` | Fallback protein: BLAST or gene mapping raised an error for this record (details in `predictor_error`) |

If every record returns `fallback_predictor_error`, the run stops instead, since that points to a broken environment (for example a missing `blastp` binary or a malformed reference zip) rather than to the data.

### Results columns

Gene symbols come from three places, and the column names say which: `predicted_*` and `best_hit_*` are BLAST hits named by the UniProt gene map in the weights file; `model_protein_key` and `protein_passed_to_model` are entries in ProtiCelli's own vocabulary (`antibody_map.pkl`); `fasta_hgnc_symbol` is the protein the FASTA header names, looked up in the same UniProt gene map.

| Column | Content |
|---|---|
| `prediction_id` | From the input record |
| `status` | See above |
| `output_file` | Image filename |
| `protein_passed_to_model` | What ProtiCelli was asked to predict: the hit's vocabulary entry, or the fallback protein |
| `predicted_uniprot_accession`, `predicted_ensembl_gene_id`, `predicted_hgnc_symbol` | The hit used: the highest-ranked hit in ProtiCelli's vocabulary |
| `model_protein_key` | ProtiCelli vocabulary entry for the hit used, from `antibody_map.pkl` (an entry can list several symbols, for example `PLG,PLGLB1,PLGLB2`) |
| `predicted_hit_rank` | Rank of the hit used among all hits (1 means the best hit overall was usable) |
| `percent_identity`, `e_value`, `alignment_length`, `bit_score` | Statistics for the hit used |
| `best_hit_uniprot_accession`, `best_hit_hgnc_symbol`, `best_hit_percent_identity`, `best_hit_e_value`, `best_hit_bit_score` | The best hit overall, in the vocabulary or not |
| `blast_task` | `blastp`, or `blastp-short` for sequences under 30 residues |
| `blast_targets_requested`, `blast_hits_returned` | 25 unless the wider search (500) was needed |
| `blast_seconds`, `predictor_cached` | BLAST time for this sequence, and whether the result was reused from an earlier record with the same sequence |
| `gene_candidate_count`, `mapping_ambiguous` | Whether the accession mapped to more than one Ensembl gene |
| `search_database` | BLAST database path used |
| `predictor_pipeline_version` | Protein predictor version |
| `predictor_error` | Error message, if any |
| `fasta_record_count` | Sequences in the FASTA file; only the first is used |
| `fasta_header`, `fasta_hpa_antibody_id`, `fasta_ensembl_protein_id`, `fasta_ensembl_release` | From the first FASTA header, reporting only |
| `fasta_hgnc_symbol` | HGNC symbol for `fasta_ensembl_protein_id`, from UniProt's Ensembl cross-references in the gene map (reporting only; compare with `best_hit_hgnc_symbol`) |
| `fasta_hgnc_source` | How `fasta_hgnc_symbol` was found: `uniprot_ensembl_xref`; or why it is blank: `not_in_gene_map` (the header's Ensembl release has an id the current UniProt release no longer lists), `ambiguous_in_gene_map` (several symbols, all listed), `gene_map_has_no_protein_ids` (weights files v1 and v2), `no_ensembl_protein_id_in_header` |
| `nucleus_max` | Nucleus maximum in counts, the output scale (predicted records only) |
| `gain_microtubules`, `gain_nucleus`, `gain_er` | `ImageNormalizer` gain per reference channel |
| `prediction_max` | Brightest decoded prediction value (0 to 1.3) |
| `output_max` | Brightest pixel in the output PNG |
| `cell_line`, `seed`, `num_inference_steps` | Run settings |
| `run_id` | Short random identifier shared by all rows of one run |

