# Preparing a submission: from trained weights to a challenge entry

This guide follows the training notebook (`training/finetune_proticelli.ipynb`). It takes the weights file the notebook produced and turns the model under test into a challenge submission, tested on two cells on AnVIL before it is submitted. BLAST-ProtiCelli stands in for a user-defined model throughout; the steps are the same for any model.

The rules come from the challenge's participant documentation, which is the authority if anything here disagrees:

- [Participant quickstart](https://github.com/cm4ai-anvil-dream-challenge/anvil-dream-participant/blob/main/docs/quickstart.md)
- [Submission contract](https://github.com/cm4ai-anvil-dream-challenge/anvil-dream-participant/blob/main/docs/contract.md)
- [Container images](https://github.com/cm4ai-anvil-dream-challenge/anvil-dream-participant/blob/main/docs/container-images.md)

> **Draft.** Items marked **OPEN** are not settled yet.

## What a submission is

A submission is three things, packed into one archive by the participant command line tool and sent to the challenge's Synapse queue:

| Part | What it holds | Where it lives | This example |
|---|---|---|---|
| WDL | One task named `Predict` that runs the model | GitHub, copied to the workspace bucket | `submission.wdl` |
| Weights file | Everything the model reads at run time: weights, vocabularies, lookup tables | Workspace bucket | `blast_proticelli_weights_finetuned_<mode>_<UTC stamp>.zip` |
| Docker image | Code and software only, no weights | Public registry, named in the WDL by digest | `quay.io/ohsu-comp-bio/blast-proticelli-example@sha256:...` |

The image itself is not uploaded. The WDL's `runtime { docker: ... }` line is the only reference to it, and the challenge pulls it from the registry when the submission runs. Optional supporting files (configuration, vocabularies) can also be included; they arrive in the task's working directory at run time.

Only one step happens off AnVIL: building and publishing the image. The weights stay in the workspace bucket from training through submission.

### Rules that shape the steps

- **No network at run time.** Everything the model needs is in the image (software) or the weights file (data). No downloads from Hugging Face, no `pip install`, no web lookups while predicting. This example forces local BLAST and a local gene map; see `configure_protein_predictor()` in `run_inference.py`.
- **The WDL interface is fixed.** One task named `Predict` and no workflow block; inputs exactly `inputs` and `weights`; one output, `predictions`; `docker` is the only runtime setting. The challenge sets CPU, memory, disk, GPU (one NVIDIA L4) and the time limit.
- **Outputs:** one single channel, 8-bit, 512 x 512 PNG per record, named `<prediction_id>.png`, in a `predictions/` folder.
- **The image** is public, built for `linux/amd64`, includes `tar`, and is named by digest (`@sha256:...`), not by a tag.

## Two WDLs in this repository

| File | Purpose | Submitted? |
|---|---|---|
| `submission.wdl` | Meets the contract. Image fixed by digest; returns predictions only | Yes |
| `submission_test.wdl` | The AnVIL test of `submission.wdl`, generated from it by `scripts/make_submission_test.py`. It adds what the organizers add to their copy: the reference resources (8 CPU, 32 GB, 75 GB SSD, one L4) and a workflow block, since Cromwell, the engine behind AnVIL Workflows, runs workflows. Everything else is copied unchanged | No |
| `blast_proticelli.wdl` | Development: image as an input, results table as an output, resources settable | No |

All three run the same `run_inference.py` command. The model also writes `prediction_results.tsv` (BLAST match and status per record); `submission.wdl` keeps it out of `predictions/`, and the task log carries the same information one line per record.

## 1. Commit and tag the code

Merge all code changes, then tag a release so the image, the WDL and the notebook point at the same code:

```bash
git tag v0.3.0
git push origin v0.3.0
```

## 2. Build and publish the image (off AnVIL)

On any machine with Docker (a laptop, a build server), from the tagged commit. On an Apple Silicon Mac, `--platform linux/amd64` is required; Terra runs x86_64.

```bash
IMAGE=quay.io/ohsu-comp-bio/blast-proticelli-example
docker buildx build --platform linux/amd64 --push -t $IMAGE:v0.3.0 .
docker buildx imagetools inspect $IMAGE:v0.3.0
```

The `Digest:` line of the inspect output is the image's permanent name. Put it in `submission.wdl`:

```wdl
runtime {
  docker: "quay.io/ohsu-comp-bio/blast-proticelli-example@sha256:<digest>"
}
```

Then regenerate the AnVIL test copy, so it carries the same digest, and commit both:

```bash
python scripts/make_submission_test.py
git commit -am "Pin submission image to <digest>"
```

`python scripts/make_submission_test.py --check` reports whether the test copy matches `submission.wdl`. The registry repository must allow anonymous pulls (on Quay, set it to public).

Optional check while Docker is at hand: run the image with networking turned off (see [Developer route](#developer-route)).

## 3. Locate the weights file (on AnVIL)

Notebook step 3b copied the weights file to the workspace bucket under `weights/`. In an AnVIL terminal:

```bash
gsutil ls -l "$WORKSPACE_BUCKET/weights/"
```

Record its SHA-256 for the submission record (see [Record what was submitted](#record-what-was-submitted)).

## 4. Test on AnVIL with two cells

This run tests exactly what will be submitted: the same task, command and image digest.

**Test cells.** Copy two complete input sets (three reference images and a FASTA file each) to the workspace bucket. This example uses the two cells in `test_inputs/`. **OPEN:** whether the tutorial uses these or two KOLF held-out cells.

```bash
gsutil cp test_inputs/* "$WORKSPACE_BUCKET/test_cells/"
```

**Inputs file.** Build it from the local copy, writing bucket paths:

```bash
python scripts/make_wdl_inputs.py --dir test_inputs \
    --wdl submission_test \
    --uri-prefix "$WORKSPACE_BUCKET/test_cells" \
    --weights-file "$WORKSPACE_BUCKET/weights/<weights file name>" \
    -o submission_test_inputs.json
```

**Run.** In the workspace's **Workflows** tab:

1. Add `submission_test.wdl` from this repository through Dockstore. The repository's `.dockstore.yml` registers it, and Dockstore follows the repository, so a pushed change appears there. **OPEN:** screenshots.
2. Choose "Run workflow with inputs defined by file paths".
3. Upload `submission_test_inputs.json`.
4. Launch.

Pass criteria:

- [ ] The job succeeds, and the log shows `CUDA available True`
- [ ] The outputs are `cell_0.png` and `cell_1.png`, and nothing else
- [ ] The log shows one line per cell with the protein used and its status

## 5. Check the outputs against the image specification

`cm4ai-evaluate` in its default `strict` mode refuses any image that does not meet the challenge specification, so a clean run doubles as a format check. Build a manifest pairing each prediction with its ground truth (templates in `cm4ai-evaluate/examples/`), then:

```bash
cm4ai-evaluate manifest.csv -o results/submission_test --metrics basic
```

- [ ] Exit code 0: every row was scored, so every output meets the specification

Scores on two cells mean little; this step checks format, not performance.

## 6. Submit (on AnVIL)

**One-time setup.** The participant tool, `anvil-dream-participant`, needs Python 3.12 or newer; AnVIL's Jupyter environment has 3.10. In a Jupyter terminal, from the repository folder:

```bash
bash scripts/setup_submit_env.sh
conda activate ~/dream-py312
```

The script creates a separate Python 3.12 environment in the home folder (on the persistent disk, so it survives pausing the cloud environment), installs the tool version named in the participant quickstart, and checks that it runs. The notebook's Python is not touched. The environment belongs to the cloud environment, not the workspace, so a cloned workspace does not include it; each user runs the script once. In each new terminal, activate it again with `conda activate ~/dream-py312`.

Then log in to Synapse as described in the [participant quickstart](https://github.com/cm4ai-anvil-dream-challenge/anvil-dream-participant/blob/main/docs/quickstart.md#2-configure-authentication). On AnVIL, no Google login is needed; the tool reads the registered email from Terra.

Copy the submission WDL next to the weights, then submit from the bucket:

```bash
gsutil cp submission.wdl "$WORKSPACE_BUCKET/submission/workflow.wdl"
export DREAM_EVALUATION_ID=9620108
anvil-dream-participant submit-anvil \
  --wdl "$WORKSPACE_BUCKET/submission/workflow.wdl" \
  --weights "$WORKSPACE_BUCKET/weights/<weights file name>" \
  --evaluation-id "$DREAM_EVALUATION_ID" \
  --name blast-proticelli-v0.3.0
```

Before uploading, the tool checks the WDL against the contract and checks that the image can be pulled anonymously for `linux/amd64`. Use a new `--name` for every submission. Save the `submission_id` it prints, and check progress with:

```bash
anvil-dream-participant submission-status --submission-id <submission_id>
```

**OPEN:** confirm the queue ID for the tutorial (the quickstart lists `9620108`), and whether the organizers recommend another Python 3.12 setup on AnVIL.

## Record what was submitted

| Item | Value |
|---|---|
| Code | repository tag |
| Image | `quay.io/...@sha256:...` |
| Weights file | name, bucket path, SHA-256 |
| WDL | `submission.wdl` at the tag |
| Test run | AnVIL job ID, date |
| Submission | `submission_id`, `--name` |

## Developer route

For development on a laptop or cluster, without AnVIL.

**Check the WDL syntax** before any run. This checks WDL syntax only; the contract checks run when `anvil-dream-participant` submits, before anything is uploaded:

```bash
miniwdl check submission.wdl
```

**Sandbox test.** Run the image with networking turned off on the two test cells. This needs a local copy of a weights file; the published baseline weights are enough, since downloads come from code, not weights. Mounting the current folder at the same path keeps the absolute paths in the records file valid:

```bash
python scripts/make_wdl_inputs.py --dir test_inputs --weights-file weights.zip -o wdl_inputs.json
python -c "import json; json.dump(json.load(open('wdl_inputs.json'))['blast_proticelli.records'], open('records.json', 'w'), indent=2)"
docker run --rm --network none --platform linux/amd64 \
    -v "$PWD":"$PWD" -w "$PWD" \
    quay.io/ohsu-comp-bio/blast-proticelli-example:v0.3.0 \
    python /app/run_inference.py --inputs records.json \
        --weights_zip weights.zip --output_dir sandbox_test
```

Without an NVIDIA GPU this runs on the CPU, slowly, but finishes for two cells. With one, add `--gpus all`.

**Development WDL with miniwdl** (Linux with an NVIDIA GPU):

```bash
miniwdl run blast_proticelli.wdl -i wdl_inputs.json
```

Its `docker` input defaults to the release tag and can be pointed at any image under test.
