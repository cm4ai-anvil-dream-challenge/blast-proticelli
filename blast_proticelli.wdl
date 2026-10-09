version 1.0

# BLAST-ProtiCelli baseline for the CM4AI DREAM Challenge.
#
# One task, one container: identifies each record's protein by BLAST, then
# predicts its localization with ProtiCelli. All records in `records` are
# processed in a single task so the model loads once; sharding across tasks
# is left to the calling workflow.
#
# Everything the model needs at runtime (ProtiCelli weights, BLAST database,
# gene map) arrives in `weights_file`. See README.md for its layout.

struct PredictionInput {
    String prediction_id
    File microtubules_image
    File er_image
    File nucleus_image
    File fasta_file
}

workflow blast_proticelli {
    input {
        Array[PredictionInput] records
        File weights_file
        String docker
        String? cell_line
        Int? seed
        Int num_inference_steps = 50
        Int batch_size = 4
    }

    call predict {
        input:
            records = records,
            weights_file = weights_file,
            docker = docker,
            cell_line = cell_line,
            seed = seed,
            num_inference_steps = num_inference_steps,
            batch_size = batch_size
    }

    output {
        Array[File] predictions = predict.predictions
        File prediction_results = predict.prediction_results
    }

    meta {
        description: "BLAST-ProtiCelli baseline: amino acid sequence to HGNC symbol by BLAST, then ProtiCelli image prediction."
    }
}

task predict {
    input {
        Array[PredictionInput] records
        File weights_file
        String docker
        String? cell_line
        Int? seed
        Int num_inference_steps
        Int batch_size

        # Runtime resources, sized for Terra's g2-standard-8 (1x NVIDIA L4).
        Int cpu = 8
        String memory = "32 GiB"
        Int disk_gb = 75
        Int gpu_count = 1
        String gpu_type = "nvidia-l4"
    }

    # write_json runs after localization, so the paths inside it point at
    # the task's local copies of each file, not the original gs:// URIs.
    File records_json = write_json(records)

    command <<<
        set -euo pipefail

        # Report the GPU and driver, and whether PyTorch can use it. On an
        # unsupported driver PyTorch silently falls back to the CPU; this puts
        # the answer at the top of the log. Neither check can fail the task.
        nvidia-smi || echo "nvidia-smi not available"
        python -c "import torch; print('torch', torch.__version__, '| CUDA build', torch.version.cuda, '| CUDA available', torch.cuda.is_available())" || true

        python /app/run_inference.py \
            --inputs ~{records_json} \
            --weights_zip ~{weights_file} \
            --output_dir out \
            --num_inference_steps ~{num_inference_steps} \
            --batch_size ~{batch_size} \
            ~{"--cell_line " + cell_line} \
            ~{"--seed " + seed}
    >>>

    output {
        Array[File] predictions = glob("out/*.tiff")
        File prediction_results = "out/prediction_results.tsv"
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: memory
        disks: "local-disk " + disk_gb + " SSD"
        gpuCount: gpu_count
        gpuType: gpu_type
    }
}
