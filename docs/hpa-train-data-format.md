# HPA-train_v1.0: data format notes (2026-10-16)

## Source

- Terra workspace bucket, folder `HPA-train_v1.0/` (requester pays; billed to the Terra Google project, not AnVIL)

### Tar locations

Each antibody's files are in one tar, `<antibody_id>.tar`. In the bucket, tars
are grouped into folders named by the first five characters of the antibody ID:

    HPA-train_v1.0/tars/<first 5 characters of antibody_id>/<antibody_id>.tar

For example, `CAB000011.tar` is at `HPA-train_v1.0/tars/CAB00/CAB000011.tar`.
The manifest's `tar_file` column holds only the tar name, so code must add the
folder. Checked for v1.0: every tar named in the manifest is at this path, and
the bucket holds no other tars.

## Manifest (`manifest.csv.gz`, 1,245,397 rows, one row per cell)

Columns: `image_key, antibody_id, ensembl_protein_id, ensembl_gene_id, gene_name, cell_line, locations, fov_key, tar_file, protein_path, nucleus_path, microtubules_path, er_path, mask_path, sequence_path`

| Item | Format | Example |
|---|---|---|
| `image_key` | `<antibody_id>_<cell_line>_<cell key>` | `CAB000011_A-431_2h4spijca3` |
| Channel images | `<antibody_id>/<image_key>_<channel>.png`, channel = protein, nucleus, microtubules, er, mask | `CAB000011/CAB000011_A-431_2h4spijca3_er.png` |
| Sequence | `<antibody_id>/<antibody_id>.faa`, one per antibody, repeated on every row for that antibody | `CAB000011/CAB000011.faa` |
| Tar | `<antibody_id>.tar`, one per antibody (checked: true for every row) | `CAB000011.tar` |
| Tar location in bucket | `HPA-train_v1.0/tars/<first 5 characters of antibody_id>/<tar_file>`; folders `CAB00` to `CAB08` and `HPA00` to `HPA07`; the manifest holds only the tar name. Verified: every manifest tar exists at this path, and the bucket holds no others | `tars/CAB00/CAB000011.tar` |
| `locations` | semicolon separated; may be empty | |
| `fov_key` | shared by several cells from one field of view | `fdtvdd5t4m` |
| `cell_line` | present on every row (so it can be passed to `fit()`) | `A-431` |

**Rule:** take identifiers and paths from manifest columns; don't parse file names (cell line names contain hyphens, and underscores are not guaranteed absent).

## Implications

- **Training notebook:** build training records from manifest rows; FASTA is per antibody, not per cell. Cell line is now available for `fit()` (ProtiCelli developer said it is optional; unknown lines fall back to the unknown tag).
- **cm4ai-evaluate:** takes a manifest of paths, so likely needs only a documented conversion, not code changes.
