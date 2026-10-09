#!/usr/bin/env python3
"""
protein_predictor.py

Lightweight, "zero development" baseline protein predictor for the CM4AI
DREAM Challenge stress test.

For each unique aa_sequence in a per-cell manifest:
  1. BLASTp the sequence against UniProtKB via the EBI Job Dispatcher REST
     API (top hit only, no identity/e-value filtering).
  2. Map the resulting UniProt accession to an Ensembl gene ID via
     MyGene.info.

Results are cached on disk keyed by sequence, then joined back onto every
row of the manifest that shares that sequence, and written out as a new
TSV alongside the diagnostic and provenance columns.

Usage:
    python protein_predictor.py manifest.tsv predicted_manifest.tsv

Config (edit below or override with an env var):
    EBI_EMAIL   required by EBI's usage policy for every Job Dispatcher job
"""

import argparse
import glob
import hashlib
import json
import os
import subprocess
import sys
import tempfile
import time
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

import pandas as pd
import requests

# ---- Config --------------------------------------------------------------

EBI_EMAIL = os.environ.get("EBI_EMAIL", "lusardi@ohsu.edu")
EBI_BASE = "https://www.ebi.ac.uk/Tools/services/rest/ncbiblast"
EBI_DATABASE = os.environ.get("EBI_DATABASE", "uniprotkb_swissprot")
ORGANISM_FILTER = os.environ.get("ORGANISM_FILTER", "Homo sapiens")
MYGENE_BASE = "https://mygene.info/v3/query"

# "ebi" (default, unchanged) or "local" (see build_local_db / local_blast_search
# below). Local search needs LOCAL_BLAST_DB pointing at a database prefix
# already built with build_local_db.
SEARCH_BACKEND = os.environ.get("SEARCH_BACKEND", "ebi")
LOCAL_BLAST_DB = os.environ.get("LOCAL_BLAST_DB", "")
BLASTP_BIN = os.environ.get("BLASTP_BIN", "blastp")
MAKEBLASTDB_BIN = os.environ.get("MAKEBLASTDB_BIN", "makeblastdb")

# "mygene" (default, unchanged) or "local" (see build_human_gene_map.sh and
# local_gene_lookup below). Local mapping needs LOCAL_GENE_MAP pointing at
# the compact accession -> gene id file that script produces.
MAPPING_BACKEND = os.environ.get("MAPPING_BACKEND", "mygene")
LOCAL_GENE_MAP = os.environ.get("LOCAL_GENE_MAP", "")
_LOCAL_GENE_MAP_CACHE = None  # lazy-loaded once per process, not once per sequence

# Optional: point these at the zips build_human_blast_db.sh / build_human_gene_map.sh
# produce, instead of already-extracted paths. Purely additive: if these
# aren't set, LOCAL_BLAST_DB / LOCAL_GENE_MAP behave exactly as before,
# already-extracted paths given directly. See ensure_reference_data_ready().
BLAST_DB_ZIP = os.environ.get("BLAST_DB_ZIP", "")
GENE_MAP_ZIP = os.environ.get("GENE_MAP_ZIP", "")
REFERENCE_DATA_DIR = os.environ.get("REFERENCE_DATA_DIR", "./reference_data")


def _file_sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _ensure_extracted(zip_path: str, extract_dir: str, expect_glob: str) -> bool:
    """
    Extract zip_path into extract_dir, unless what's already there was
    extracted from this exact zip (by content hash), recorded in a small
    marker file left alongside the extracted data. This means pointing at
    a different release, even while reusing the same extract_dir, still
    triggers a real extraction rather than silently continuing to serve
    whatever was extracted before, existence alone isn't enough to prove
    it's the same data. Returns True if extraction actually happened,
    False if the existing extracted data already matches this zip.
    """
    Path(extract_dir).mkdir(parents=True, exist_ok=True)
    zip_hash = _file_sha256(zip_path)
    marker_path = os.path.join(extract_dir, f".source_sha256_{expect_glob.replace('*', '')}")

    already_current = (
        glob.glob(os.path.join(extract_dir, expect_glob))
        and os.path.exists(marker_path)
        and Path(marker_path).read_text().strip() == zip_hash
    )
    if already_current:
        return False

    with zipfile.ZipFile(zip_path, "r") as zf:
        zf.extractall(extract_dir)
    Path(marker_path).write_text(zip_hash)
    return True


def ensure_reference_data_ready():
    """
    If BLAST_DB_ZIP / GENE_MAP_ZIP are set, extract them (idempotently) into
    REFERENCE_DATA_DIR and point LOCAL_BLAST_DB / LOCAL_GENE_MAP at the
    extracted files, unless those were already explicitly set to something
    else. No-op for either backend not in "local" mode, and a no-op overall
    if neither zip variable is set, existing already-extracted-path usage is
    completely unaffected.
    """
    global LOCAL_BLAST_DB, LOCAL_GENE_MAP

    if SEARCH_BACKEND == "local" and BLAST_DB_ZIP:
        extracted = _ensure_extracted(BLAST_DB_ZIP, REFERENCE_DATA_DIR, "human_reviewed_swissprot.*")
        if not LOCAL_BLAST_DB:
            LOCAL_BLAST_DB = os.path.join(REFERENCE_DATA_DIR, "human_reviewed_swissprot")
        print(f"Local BLAST db ready at {LOCAL_BLAST_DB} (extracted just now: {extracted})")

    if MAPPING_BACKEND == "local" and GENE_MAP_ZIP:
        extracted = _ensure_extracted(GENE_MAP_ZIP, REFERENCE_DATA_DIR, "human_gene_map.tsv")
        if not LOCAL_GENE_MAP:
            LOCAL_GENE_MAP = os.path.join(REFERENCE_DATA_DIR, "human_gene_map.tsv")
        print(f"Local gene map ready at {LOCAL_GENE_MAP} (extracted just now: {extracted})")

CACHE_DIR = Path("./predictor_cache")
POLL_INTERVAL_SEC = 5
POLL_TIMEOUT_SEC = 300  # ceiling per sequence; blastp vs uniprotkb is usually much faster

PIPELINE_VERSION = "protein_predictor_v0.1_local"

# ---- Small helpers --------------------------------------------------------

def sequence_key(seq: str) -> str:
    """
    Short, filesystem safe internal cache key. Never written to the manifest.

    Includes everything that could change the actual prediction, not just
    the sequence: which search backend produced it (ebi vs local) and its
    config, and which mapping backend resolved the gene id (mygene.info
    vs a local table) and its config. Changing any of these yields a
    different key rather than silently replaying a result computed under
    different settings. Old cache files from a prior config just become
    orphaned, harmless dead files.
    """
    if SEARCH_BACKEND == "local":
        search_tag = f"local|{LOCAL_BLAST_DB}"
    else:
        search_tag = f"ebi|{EBI_DATABASE}|{ORGANISM_FILTER}"

    mapping_tag = f"local|{LOCAL_GENE_MAP}" if MAPPING_BACKEND == "local" else "mygene"

    payload = f"{seq.strip().upper()}|{search_tag}|{mapping_tag}"
    return hashlib.sha256(payload.encode()).hexdigest()[:16]


def cache_path(seq: str) -> Path:
    CACHE_DIR.mkdir(exist_ok=True)
    return CACHE_DIR / f"{sequence_key(seq)}.json"


def load_cached(seq: str):
    p = cache_path(seq)
    if p.exists():
        return json.loads(p.read_text())
    return None


def save_cache(seq: str, result: dict):
    cache_path(seq).write_text(json.dumps(result, indent=2))


# ---- Step 1: EBI BLAST, sequence to UniProt accession --------------------

def submit_blast(seq: str) -> str:
    resp = requests.post(
        f"{EBI_BASE}/run",
        data={
            "email": EBI_EMAIL,
            "program": "blastp",
            "stype": "protein",
            "sequence": seq,
            "database": EBI_DATABASE,
        },
        timeout=30,
    )
    resp.raise_for_status()
    return resp.text.strip()  # job id, plain text


def wait_for_blast(job_id: str) -> str:
    elapsed = 0
    while elapsed < POLL_TIMEOUT_SEC:
        resp = requests.get(f"{EBI_BASE}/status/{job_id}", timeout=30)
        resp.raise_for_status()
        status = resp.text.strip()
        if status == "FINISHED":
            return status
        if status in ("FAILURE", "ERROR", "NOT_FOUND"):
            raise RuntimeError(f"BLAST job {job_id} ended with status {status}")
        time.sleep(POLL_INTERVAL_SEC)
        elapsed += POLL_INTERVAL_SEC
    raise TimeoutError(f"BLAST job {job_id} did not finish within {POLL_TIMEOUT_SEC} seconds")


def fetch_blast_xml(job_id: str) -> str:
    resp = requests.get(f"{EBI_BASE}/result/{job_id}/xml", timeout=30)
    resp.raise_for_status()
    return resp.text


def parse_top_hit(xml_text: str, organism_filter: str = None) -> dict:
    """
    Parse the XML that the EBI Job Dispatcher returns for ncbiblast and
    pull out the best hit, optionally restricted to a given organism.

    This is EBI's own EBIApplicationResult schema, not raw NCBI BLAST XML:
    the root element declares a default namespace, so tags must be matched
    with the "{*}tag" wildcard (any namespace) rather than a bare name, or
    ElementTree silently matches nothing. Accession, organism, and hit
    length live as attributes directly on <hit> (id, ac, description),
    not as separately parsed child text. Percent identity is reported
    directly under <alignment> (no division needed), the e-value field is
    named "expectation", and alignment length is derived from the
    querySeq start/end attributes since there's no dedicated field for it.

    Hits are already returned in descending significance order, so when
    organism_filter is set this walks that list and returns the first hit
    whose description (which carries an "OS=<organism>" tag from the
    original UniProt FASTA header) matches, i.e. the best-scoring hit for
    that organism specifically, not just the single best hit overall.

    NOTE: this is a post hoc filter, not a search-time restriction. EBI's
    hosted ncbiblast REST service has no documented parameter for
    restricting the search itself to a taxon (unlike local BLAST+'s
    -taxids), so this does not reduce search time, it only changes which
    of the already-returned hits gets reported. If no hit in the returned
    list matches, this reports no hit at all rather than falling back to
    a non-matching one.
    """
    root = ET.fromstring(xml_text)
    hits = root.findall(".//{*}hit")

    empty = {
        "uniprot_accession": None,
        "percent_identity": None,
        "e_value": None,
        "alignment_length": None,
    }

    for hit in hits:
        description = hit.get("description", "")
        if organism_filter and organism_filter.lower() not in description.lower():
            continue

        accession = hit.get("ac")

        alignment = hit.find(".//{*}alignment")
        if alignment is None:
            continue

        identity_text = alignment.findtext("{*}identity")
        e_value = alignment.findtext("{*}expectation")

        alignment_length = None
        query_seq_el = alignment.find("{*}querySeq")
        if query_seq_el is not None:
            start, end = query_seq_el.get("start"), query_seq_el.get("end")
            if start is not None and end is not None:
                alignment_length = int(end) - int(start) + 1

        return {
            "uniprot_accession": accession,
            "percent_identity": round(float(identity_text), 1) if identity_text else None,
            "e_value": e_value,
            "alignment_length": alignment_length,
        }

    return empty  # no hit at all, or no hit matched organism_filter


# ---- Local BLAST backend --------------------------------------------------

def build_local_db(fasta_path: str, db_prefix: str):
    """
    Build a local blastp-searchable database from a protein FASTA file
    (see the human-reviewed-proteome download instructions in the README
    / accompanying notes). One-time step, run whenever the source FASTA
    changes, not per query.
    """
    subprocess.run(
        [MAKEBLASTDB_BIN, "-in", fasta_path, "-dbtype", "prot", "-out", db_prefix],
        check=True, capture_output=True, text=True,
    )


def local_blast_search(seq: str, db_prefix: str) -> dict:
    """
    Run blastp locally against a database built with build_local_db, and
    return the top hit in the same shape parse_top_hit returns for the
    EBI backend. No organism filter is applied here: if the database was
    built from a human-only FASTA (as intended), every hit already is
    human by construction, there is nothing to filter after the fact.

    Raises CalledProcessError if blastp itself fails (e.g. bad db path);
    a clean "no hit found" is a successful run with empty stdout, not an
    exception, same distinction the EBI path makes.
    """
    empty = {
        "uniprot_accession": None,
        "percent_identity": None,
        "e_value": None,
        "alignment_length": None,
    }

    with tempfile.NamedTemporaryFile(mode="w", suffix=".faa", delete=False) as f:
        f.write(f">query\n{seq}\n")
        query_path = f.name

    try:
        proc = subprocess.run(
            [
                BLASTP_BIN,
                "-query", query_path,
                "-db", db_prefix,
                "-outfmt", "6",
                "-max_target_seqs", "1",
            ],
            capture_output=True, text=True, timeout=60, check=True,
        )
    finally:
        os.unlink(query_path)

    if not proc.stdout.strip():
        return empty

    # outfmt 6 default columns: qseqid sseqid pident length mismatch gapopen
    # qstart qend sstart send evalue bitscore
    fields = proc.stdout.strip().splitlines()[0].split("\t")
    sseqid, pident, length, evalue = fields[1], fields[2], fields[3], fields[10]
    accession = sseqid.split("|")[1] if "|" in sseqid else sseqid

    return {
        "uniprot_accession": accession,
        "percent_identity": round(float(pident), 1),
        "e_value": evalue,
        "alignment_length": int(length),
    }


# ---- Step 2: MyGene.info, UniProt accession to Ensembl gene ID -----------

def map_to_ensembl_gene(accession: str) -> dict:
    resp = requests.get(
        MYGENE_BASE,
        params={"q": f"uniprot:{accession}", "fields": "symbol,ensembl.gene", "species": "human"},
        timeout=30,
    )
    resp.raise_for_status()
    hits = resp.json().get("hits", [])

    gene_ids = set()
    symbols = set()
    for h in hits:
        ens = h.get("ensembl")
        if isinstance(ens, list):
            for e in ens:
                if "gene" in e:
                    gene_ids.add(e["gene"])
        elif isinstance(ens, dict) and "gene" in ens:
            gene_ids.add(ens["gene"])
        if h.get("symbol"):
            symbols.add(h["symbol"])

    gene_ids = sorted(gene_ids)
    symbols = sorted(symbols)
    return {
        "ensembl_gene_id": gene_ids[0] if gene_ids else None,
        "hgnc_symbol": symbols[0] if symbols else None,
        "gene_candidate_count": len(gene_ids),
        "mapping_ambiguous": len(gene_ids) > 1,
    }


def _load_local_gene_map(path: str) -> dict:
    """
    Load the compact accession -> (symbol, gene id(s)) file
    build_human_gene_map.sh produces (three columns: accession, hgnc
    symbol, semicolon-joined Ensembl gene ids) into memory once per
    process, not once per sequence. Subsequent calls reuse the same
    in-memory dict.
    """
    global _LOCAL_GENE_MAP_CACHE
    if _LOCAL_GENE_MAP_CACHE is not None:
        return _LOCAL_GENE_MAP_CACHE

    mapping = {}
    with open(path) as f:
        next(f, None)  # header
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            accession, symbol, gene_blob = parts[0], parts[1], parts[2]
            mapping[accession] = {
                "hgnc_symbol": symbol or None,
                "gene_ids": [g for g in gene_blob.split(";") if g],
            }

    _LOCAL_GENE_MAP_CACHE = mapping
    return mapping


def local_gene_lookup(accession: str) -> dict:
    """
    Look up an accession's HGNC symbol and Ensembl gene id(s) in a
    pre-built local table instead of calling mygene.info. Pure in-memory
    dict lookup once the table is loaded, no network involved at all.
    """
    mapping = _load_local_gene_map(LOCAL_GENE_MAP)
    entry = mapping.get(accession, {"hgnc_symbol": None, "gene_ids": []})
    gene_ids = entry["gene_ids"]
    return {
        "ensembl_gene_id": gene_ids[0] if gene_ids else None,
        "hgnc_symbol": entry["hgnc_symbol"],
        "gene_candidate_count": len(gene_ids),
        "mapping_ambiguous": len(gene_ids) > 1,
    }


# ---- FASTA parsing ---------------------------------------------------------

def read_fasta(path) -> str:
    """
    Parse a single-record FASTA/.faa file into a raw sequence string.
    Header content is never used for anything, just skipped, so it can be
    stripped down to whatever (or nothing) upstream without affecting this.

    Raises if a ">"-delimited file has zero or more than one record: a
    per-image file with multiple records means something unexpected
    happened upstream, and silently taking the first would hide that
    rather than surface it.
    """
    path = Path(path)
    text = path.read_text()

    if ">" not in text:
        # No header line at all, whole file is the sequence.
        return "".join(text.split()).upper()

    records = [r for r in text.split(">") if r.strip()]

    if len(records) == 0:
        raise ValueError(f"{path} has no FASTA records")
    if len(records) > 1:
        raise ValueError(f"{path} has {len(records)} FASTA records, expected exactly 1")

    lines = records[0].splitlines()
    seq = "".join(line.strip() for line in lines[1:])  # skip header line, whatever it contains
    return seq.upper()


# ---- Orchestration ---------------------------------------------------------

def predict_for_sequence(seq: str, force_rerun: bool = False) -> dict:
    if not force_rerun:
        cached = load_cached(seq)
        if cached is not None:
            return cached

    result = {
        "predicted_uniprot_accession": None,
        "predicted_ensembl_gene_id": None,
        "predicted_hgnc_symbol": None,
        "search_tool": "local_blastp" if SEARCH_BACKEND == "local" else "ebi_blastp",
        "search_database": LOCAL_BLAST_DB if SEARCH_BACKEND == "local" else EBI_DATABASE,
        "percent_identity": None,
        "e_value": None,
        "alignment_length": None,
        "mapping_tool": "local_uniprot_ensembl_map" if MAPPING_BACKEND == "local" else "mygene.info",
        "gene_candidate_count": 0,
        "mapping_ambiguous": False,
        "predictor_pipeline_version": PIPELINE_VERSION,
        "predictor_run_timestamp": pd.Timestamp.now("UTC").isoformat(),
    }

    try:
        if SEARCH_BACKEND == "local":
            if not LOCAL_BLAST_DB:
                raise ValueError("SEARCH_BACKEND=local requires LOCAL_BLAST_DB to be set")
            hit = local_blast_search(seq, LOCAL_BLAST_DB)
            result["search_tool"] = "local_blastp"
            result["search_database"] = LOCAL_BLAST_DB
        else:
            job_id = submit_blast(seq)
            wait_for_blast(job_id)
            xml_text = fetch_blast_xml(job_id)

            # Debug aid: dump the raw XML regardless of what parsing finds, so
            # a hit's actual Hit_def format can be inspected directly instead
            # of guessed at.
            CACHE_DIR.mkdir(exist_ok=True)
            (CACHE_DIR / f"{sequence_key(seq)}_raw.xml").write_text(xml_text)

            hit = parse_top_hit(xml_text, organism_filter=ORGANISM_FILTER)

        result["predicted_uniprot_accession"] = hit["uniprot_accession"]
        result["percent_identity"] = hit["percent_identity"]
        result["e_value"] = hit["e_value"]
        result["alignment_length"] = hit["alignment_length"]

        if hit["uniprot_accession"]:
            if MAPPING_BACKEND == "local":
                if not LOCAL_GENE_MAP:
                    raise ValueError("MAPPING_BACKEND=local requires LOCAL_GENE_MAP to be set")
                mapping = local_gene_lookup(hit["uniprot_accession"])
                result["mapping_tool"] = "local_uniprot_ensembl_map"
            else:
                mapping = map_to_ensembl_gene(hit["uniprot_accession"])
                result["mapping_tool"] = "mygene.info"
            result["predicted_ensembl_gene_id"] = mapping["ensembl_gene_id"]
            result["predicted_hgnc_symbol"] = mapping.get("hgnc_symbol")
            result["gene_candidate_count"] = mapping["gene_candidate_count"]
            result["mapping_ambiguous"] = mapping["mapping_ambiguous"]

    except Exception as exc:
        result["error"] = str(exc)
        print(f"  [warn] prediction failed for sequence hash {sequence_key(seq)}: {exc}", file=sys.stderr)
        return result  # do NOT cache failures, so a rerun retries instead of replaying the same failure

    save_cache(seq, result)
    return result


def run(manifest_path: str, output_path: str, fasta_dir: str = ".", force_rerun: bool = False):
    ensure_reference_data_ready()

    if SEARCH_BACKEND == "local":
        print(f"Search backend: local  (db={LOCAL_BLAST_DB!r})")
        if not LOCAL_BLAST_DB:
            print("  [warn] LOCAL_BLAST_DB is empty, every prediction will fail", file=sys.stderr)
    else:
        print(f"Search backend: ebi  (database={EBI_DATABASE!r}, organism_filter={ORGANISM_FILTER!r})")

    if MAPPING_BACKEND == "local":
        print(f"Mapping backend: local  (map={LOCAL_GENE_MAP!r})")
        if not LOCAL_GENE_MAP:
            print("  [warn] LOCAL_GENE_MAP is empty, every mapping will fail", file=sys.stderr)
    else:
        print("Mapping backend: mygene.info")

    if force_rerun:
        print("Force rerun: ON  (ignoring cache, results will still be re-cached)")

    df = pd.read_csv(manifest_path, sep="\t")

    path_col = "aa_sequence_fasta_path"
    if path_col not in df.columns:
        raise ValueError(f"Manifest is missing a {path_col} column")

    # Resolve each row's FASTA file to its actual sequence content up front.
    # Dedup and caching key off this content, not the file path or filename,
    # so a duplicate sequence under two different protein ID filenames still
    # collapses to one prediction.
    df["_aa_sequence"] = df[path_col].apply(lambda p: read_fasta(Path(fasta_dir) / p))

    unique_seqs = df["_aa_sequence"].dropna().unique()
    print(f"{len(df)} manifest rows, {len(unique_seqs)} unique sequences to predict")

    predictions = {}
    for i, seq in enumerate(unique_seqs, 1):
        print(f"[{i}/{len(unique_seqs)}] predicting sequence hash {sequence_key(seq)} (len={len(seq)})")
        predictions[seq] = predict_for_sequence(seq, force_rerun=force_rerun)

    pred_df = pd.DataFrame.from_dict(predictions, orient="index")
    pred_df.index.name = "_aa_sequence"
    pred_df.reset_index(inplace=True)

    merged = df.merge(pred_df, on="_aa_sequence", how="left")
    merged.drop(columns=["_aa_sequence"], inplace=True)
    merged.to_csv(output_path, sep="\t", index=False)
    print(f"Wrote {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", help="Input manifest TSV, must have an aa_sequence_fasta_path column")
    parser.add_argument(
        "output", nargs="?", default="predicted_manifest.tsv",
        help="Output TSV path (default: predicted_manifest.tsv)"
    )
    parser.add_argument("--fasta-dir", default=".", help="Directory the fasta_path column is relative to")
    parser.add_argument(
        "--force-rerun",
        action="store_true",
        help="Ignore cached results and rerun every sequence fresh (results are still written back to cache)",
    )
    args = parser.parse_args()
    run(args.manifest, args.output, args.fasta_dir, args.force_rerun)