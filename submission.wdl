version 1.0

# BLAST-ProtiCelli: challenge submission WDL.
#
# Follows the participant submission contract (anvil-dream-participant,
# docs/contract.md): one task named Predict, no workflow block, inputs
# exactly `inputs` and `weights`, output exactly `predictions`, and only
# `docker` in runtime, as a public image pinned by digest. The organizers
# set CPU, memory, disk, GPU and the timeout.
#
# For development (an image input, the results table as an output, resource
# settings), use blast_proticelli.wdl, which runs the same command.

struct PredictionInput {
  String prediction_id
  File microtubules_image
  File nucleus_image
  File er_image
  File fasta_file
  File? mask_image
}

task Predict {
  input {
    Array[PredictionInput] inputs
    File weights
  }

  # Masks, when supplied, are accepted and not used by this model.
  # The results table goes in the working directory, outside predictions/,
  # so it is not returned as a prediction; the log carries one line per
  # record (best BLAST hit, protein used, status).
  command <<<
    set -euo pipefail
    nvidia-smi || echo "nvidia-smi not available"
    python -c "import torch; print('torch', torch.__version__, '| CUDA available', torch.cuda.is_available())" || true

    python /app/run_inference.py \
      --inputs "~{write_json(inputs)}" \
      --weights_zip "~{weights}" \
      --output_dir predictions \
      --results_tsv prediction_results.tsv
  >>>

  output {
    Array[File] predictions = glob("predictions/*.png")
  }

  runtime {
    docker: "quay.io/ohsu-comp-bio/blast-proticelli-example@sha256:3f8b3178ae182fcf852163347574f98e7661aa0f0ce76f1f92993fe2b2e697ff"
  }
}
