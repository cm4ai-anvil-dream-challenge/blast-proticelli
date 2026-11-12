#!/usr/bin/env bash
# One-time setup for submitting: a Python 3.12 environment with the
# challenge's participant command line tool (anvil-dream-participant).
#
# The tool needs Python 3.12 or newer; AnVIL's Jupyter environment has 3.10.
# This creates a separate environment and leaves the notebook's Python alone.
# The environment goes in the home folder, which is on the persistent disk,
# so it survives pausing the cloud environment. Running it again is safe: an
# existing environment is reused.
#
# Usage:
#   bash scripts/setup_submit_env.sh
#   conda activate ~/dream-py312
#
# Settings (environment variables):
#   DREAM_ENV                 where the environment goes (default ~/dream-py312)
#   DREAM_PARTICIPANT_VERSION tool version (default 0.1.2, as in the
#                             participant quickstart)

set -euo pipefail

ENV_PREFIX="${DREAM_ENV:-$HOME/dream-py312}"
VERSION="${DREAM_PARTICIPANT_VERSION:-0.1.2}"
PY="$ENV_PREFIX/bin/python"

if ! command -v conda >/dev/null 2>&1; then
    echo "ERROR: conda not found. On AnVIL, run this in a Jupyter terminal." >&2
    exit 1
fi

if [ -x "$PY" ] && "$PY" -c 'import sys; sys.exit(sys.version_info < (3, 12))'; then
    echo "Reusing $ENV_PREFIX ($("$PY" --version))"
else
    echo "Creating $ENV_PREFIX with Python 3.12"
    conda create -y -q -p "$ENV_PREFIX" python=3.12
fi

echo "Installing anvil-dream-participant==$VERSION"
PIP_USER=false "$PY" -m pip install -q "anvil-dream-participant==$VERSION"

# Check: the tool runs and its submit commands are present.
"$ENV_PREFIX/bin/anvil-dream-participant" submit-anvil --help >/dev/null
echo
echo "Ready: anvil-dream-participant $("$PY" -c 'from importlib.metadata import version; print(version("anvil-dream-participant"))') on $("$PY" --version)"
echo "Activate it in each new terminal with:"
echo "    conda activate $ENV_PREFIX"
