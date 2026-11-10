#!/usr/bin/env bash
# install_blast.sh
#
# Installs NCBI BLAST+ (blastp, makeblastdb) for the example model,
# BLAST-ProtiCelli, outside Docker: the training notebook, a cluster or a
# laptop. The Docker image installs the same version itself (see the
# Dockerfile), so the workflow does not need this script.
#
# BLAST+ is a compiled program, not a Python package, so pip cannot install
# it. On Linux x86_64 this downloads NCBI's prebuilt binaries, pinned to the
# version that built the reference BLAST database, and prints the folder to
# add to PATH. Elsewhere (for example a Mac) install BLAST+ yourself:
#   conda install -c bioconda blast   or   brew install blast
#
# Usage:
#   ./install_blast.sh [tools_dir]
# Prints the bin folder as the last line of output.

set -euo pipefail

BLAST_VERSION="2.17.0"
TOOLS_DIR="${1:-./tools}"
BIN_DIR="${TOOLS_DIR}/ncbi-blast-${BLAST_VERSION}+/bin"

if [[ -x "${BIN_DIR}/blastp" && -x "${BIN_DIR}/makeblastdb" ]]; then
    echo "BLAST+ ${BLAST_VERSION} already in ${BIN_DIR}" >&2
    echo "${BIN_DIR}"
    exit 0
fi

if [[ "$(uname -s)" != "Linux" || "$(uname -m)" != "x86_64" ]]; then
    echo "No prebuilt BLAST+ for $(uname -s) $(uname -m) here. Install it with:" >&2
    echo "  conda install -c bioconda blast   or   brew install blast" >&2
    exit 1
fi

NAME="ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz"
URL="https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${BLAST_VERSION}/${NAME}"
mkdir -p "${TOOLS_DIR}"
echo "Downloading ${URL}" >&2
curl -fsSL --retry 3 "${URL}" -o "${TOOLS_DIR}/${NAME}.part"
mv "${TOOLS_DIR}/${NAME}.part" "${TOOLS_DIR}/${NAME}"
tar -xzf "${TOOLS_DIR}/${NAME}" -C "${TOOLS_DIR}"
"${BIN_DIR}/blastp" -version >&2
echo "${BIN_DIR}"
