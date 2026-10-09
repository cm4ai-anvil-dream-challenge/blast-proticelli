# BLAST-ProtiCelli baseline inference image.
#
# Holds code only. Model weights and reference data arrive at runtime in the
# weights zip (see README.md), so the image does not change when weights do.
#
# Build (--platform is required on an Apple Silicon Mac: Terra runs x86_64.
# Omitting it fails at the blastp check below, since BLAST is x86_64 only):
#   docker buildx build --platform linux/amd64 -t blast-proticelli:0.1.0 --load .
#
# Run: see README.md. The WDL task calls
#   python /app/run_inference.py --inputs ... --weights_zip ... --output_dir ...

FROM python:3.11-slim

# BLAST+ from NCBI's prebuilt binaries, pinned to the version that built the
# v1/v2 reference database (see PROVENANCE.txt in the BLAST db zip), so
# search and database come from the same release.
ARG BLAST_VERSION=2.17.0
RUN apt-get update \
 && apt-get install -y --no-install-recommends curl ca-certificates libgomp1 \
 && rm -rf /var/lib/apt/lists/* \
 && curl -fsSL "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${BLAST_VERSION}/ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz" \
    | tar -xz -C /opt \
 && ln -s "/opt/ncbi-blast-${BLAST_VERSION}+/bin/blastp" /usr/local/bin/blastp \
 && ln -s "/opt/ncbi-blast-${BLAST_VERSION}+/bin/makeblastdb" /usr/local/bin/makeblastdb

# Install the ProtiCelli fork as a regular (non-editable) package, then keep
# only the three pipeline scripts. Removing the source tree means
# "import proticelli" can only resolve to the installed package, which is
# what the build-time check below verifies.
COPY . /tmp/build
RUN pip install --no-cache-dir /tmp/build \
 && mkdir -p /app \
 && cp /tmp/build/run_inference.py /tmp/build/protein_predictor.py /tmp/build/pipeline_checks.py /app/ \
 && rm -rf /tmp/build

# Build-time checks: fail the build, not a later Terra run, if a required
# piece is missing.
#  - antibody_map.pkl and cell_line_map.pkl must ship inside the installed
#    package (a non-editable install drops data files pyproject.toml does
#    not declare).
#  - blastp must run (catches missing shared libraries).
RUN python -c "from importlib import resources; d = resources.files('proticelli') / 'data'; \
missing = [f for f in ('antibody_map.pkl', 'cell_line_map.pkl') if not (d / f).is_file()]; \
assert not missing, f'missing from installed package: {missing}'; print('package data OK')" \
 && blastp -version \
 && python -c "import sys; sys.path.insert(0, '/app'); import run_inference; print('imports OK')"

ENV PYTHONUNBUFFERED=1 \
    MPLBACKEND=Agg
