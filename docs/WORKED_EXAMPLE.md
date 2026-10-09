# Worked example: BLAST-ProtiCelli with WDL

This runs the baseline end to end on the two example records in the repository, using miniWDL on a local machine. The same WDL runs on Terra.

## What you need

- Docker
- Python 3 (the helper scripts use only the standard library)
- miniWDL: `pip install miniwdl`
- The weights file, `blast_proticelli_weights_v2.zip` (about 3.3 GB): [blast_proticelli_weights_v2.zip](https://drive.google.com/file/d/1PZkakOI8yQAzUowAQi-qK0TUQ5kh8NiT/view?usp=drive_link)

## 1. Clone the repository

```bash
git clone -b feature/blast-proticelli https://github.com/cm4ai-anvil-dream-challenge/ProtiCelliDream.git
cd ProtiCelliDream
```

## 2. Build the Docker image

```bash
docker buildx build --platform linux/amd64 -t blast-proticelli:0.1.0 --load .
```

The first build takes 10 to 15 minutes, mostly installing PyTorch.

## 3. Get the weights file

Download it into the repository folder, then check it is the expected file:

```bash
shasum -a 256 blast_proticelli_weights_v2.zip
# 2f243c8e0945b5b9830fafae1590c92638415501e7aa6c6db71da3161ed07071
```

## 4. Build the WDL inputs

```bash
python3 scripts/manifest_to_json.py manifest.tsv -o inputs.json

python3 scripts/make_wdl_inputs.py \
  --records inputs.json \
  --weights_file blast_proticelli_weights_v2.zip \
  --docker blast-proticelli:0.1.0 \
  --seed 1 --num_inference_steps 10 \
  -o wdl_inputs.json
```

`--num_inference_steps 10` keeps a local CPU run short. Leave it out for the default of 50, which gives smoother images.

## 5. Run

```bash
MINIWDL__FILE_IO__OUTPUT_HARDLINKS=true miniwdl run blast_proticelli.wdl -i wdl_inputs.json
```

To watch progress from a second terminal:

```bash
tail -f _LAST/call-predict/stdout.txt
```

On a laptop without an NVIDIA GPU the task runs on the CPU (under emulation on Apple Silicon), so expect several minutes. The first two log lines report the GPU: locally they read `nvidia-smi not available` and `CUDA available False`, which is expected.

## 6. Check the results

miniWDL prints the output paths when it finishes, and `_LAST/outputs.json` lists them. You should get two images, `cell_0.tiff` and `cell_1.tiff`, and `prediction_results.tsv`:

| prediction_id | status | predicted_uniprot_accession | predicted_hgnc_symbol | percent_identity | alignment_length |
|---|---|---|---|---|---|
| cell_0 | predicted | P55263 | ADK | 99.4 | 343 |
| cell_1 | predicted | P55263 | ADK | 100.0 | 362 |

The two FASTA files are different ADK isoforms, and BLAST resolves both to the same UniProt entry, so ProtiCelli predicts ADK localization for both cells.

See the README for the full input and output descriptions, the weights file layout, and running on Terra.
