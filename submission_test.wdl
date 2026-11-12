version 1.0

# Runs submission.wdl unchanged, for a test on AnVIL before submitting.
#
# The submission contract allows one task and no workflow block, while
# Cromwell (AnVIL's workflow engine) runs workflows. This wrapper adds only
# the workflow; the task, its command and its image digest come straight
# from submission.wdl, so a passing run tests exactly what is submitted.
# This file is not submitted.

import "submission.wdl" as submission

workflow submission_test {
  input {
    Array[PredictionInput] inputs
    File weights
  }

  call submission.Predict {
    input:
      inputs = inputs,
      weights = weights
  }

  output {
    Array[File] predictions = Predict.predictions
  }
}
