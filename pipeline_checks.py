import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tifffile import imread

def _check_array(
    img,
    report_id,
    expected_hw,
    expected_dtype,
    channel_axis=-1,
    valid_channel_counts=(3, 4),
    expected_range=(-1.0, 1.0),
    range_tolerance=0.05,
):
    """
    Shared validation core: the same channel-count / shape / dtype /
    range checks, run against an already-loaded array regardless of
    whether it came from one combined file (validate_image_file) or was
    assembled from three separate channel files (validate_channel_files).
    report_id is just this image's label in the report (a filename for
    the single-file path, a cell_id for the assembled path); it plays no
    part in any check.
    """
    n_channels = img.shape[channel_axis]
    hw_shape = tuple(s for i, s in enumerate(img.shape) if i != (channel_axis % img.ndim))

    img_min, img_max = float(img.min()), float(img.max())
    issues = []

    if n_channels not in valid_channel_counts:
        issues.append(f"channel count {n_channels} not in {valid_channel_counts}")
    if hw_shape != expected_hw:
        issues.append(f"H/W {hw_shape} != expected {expected_hw}")
    if img.dtype != expected_dtype:
        issues.append(f"dtype {img.dtype} != expected {expected_dtype}")
    if img_min < expected_range[0] - range_tolerance:
        issues.append(f"min {img_min:.3f} below expected {expected_range[0]}")
    if img_max > expected_range[1] + range_tolerance:
        issues.append(f"max {img_max:.3f} above expected {expected_range[1]}")

    passed = len(issues) == 0

    return {
        "filename": report_id, "shape": img.shape, "n_channels": n_channels,
        "dtype": img.dtype, "min": img_min, "max": img_max,
        "passed": passed,
        "issues": "; ".join(issues) if issues else "",
    }


def validate_image_file(
    filename,
    data_dir,
    expected_hw,
    expected_dtype,
    channel_axis=-1,
    valid_channel_counts=(3, 4),
    expected_range=(-1.0, 1.0),
    range_tolerance=0.05,
):
    """Load and validate a single multi-channel reference image file against known-good specs."""
    path = os.path.join(data_dir, filename)
    img = imread(path)
    record = _check_array(
        img, filename, expected_hw, expected_dtype,
        channel_axis=channel_axis, valid_channel_counts=valid_channel_counts,
        expected_range=expected_range, range_tolerance=range_tolerance,
    )
    return img, record

def validate_image_files(
    filenames,
    data_dir,
    expected_hw=None,
    expected_dtype=None,
    csv_path=None,
    only_issues_to_csv=False,
    **kwargs,
):
    """Validate a batch of files. Returns (passing_images, report_df) — failed files excluded from images."""
    images, records = [], []

    if expected_hw is None or expected_dtype is None:
        probe = imread(os.path.join(data_dir, filenames[0]))
        axis = kwargs.get("channel_axis", -1)
        if expected_hw is None:
            expected_hw = tuple(s for i, s in enumerate(probe.shape) if i != (axis % probe.ndim))
        if expected_dtype is None:
            expected_dtype = probe.dtype

    for fname in filenames:
        img, record = validate_image_file(fname, data_dir, expected_hw, expected_dtype, **kwargs)
        if record["passed"]:
            images.append(img)
        records.append(record)  # every file still goes in the report, pass or fail

    report_df = pd.DataFrame(records)
    if csv_path is not None:
        out_df = report_df[~report_df["passed"]] if only_issues_to_csv else report_df
        out_df.to_csv(csv_path, index=False)
        print(f"Wrote {len(out_df)} rows to {csv_path}")

    return images, report_df


# ---- Manifest-driven input: reassemble channels, then reuse the exact ----
# ---- same validation core the single-file path already uses.         ----

# Empirically confirmed order (Theresa inspected the example image stack
# channel by channel), matching the compacted form of ProtiCelli's
# documented 4-channel spec (channels 0, 2, 3) with the ignored protein
# channel (1) dropped. Also the order split_reference_channels.py already
# writes its output in.
CHANNEL_ORDER = ("microtubule", "nucleus", "er")


def assemble_channels(channel_paths, channel_axis=-1):
    """
    Load three single-channel image files and stack them into one
    multi-channel array, in CHANNEL_ORDER (microtubule, nucleus, er).
    channel_paths must be given in that order: [microtubule_path,
    nucleus_path, er_path]. Raises if the three loaded images don't all
    share the same H/W shape, rather than deferring to a less legible
    error from np.stack itself.
    """
    if len(channel_paths) != 3:
        raise ValueError(f"Expected exactly 3 channel paths (microtubule, nucleus, er), got {len(channel_paths)}")

    channels = [imread(p) for p in channel_paths]
    shapes = {c.shape for c in channels}
    if len(shapes) != 1:
        detail = ", ".join(f"{os.path.basename(p)}={c.shape}" for p, c in zip(channel_paths, channels))
        raise ValueError(f"Channel images have mismatched shapes: {detail}")

    return np.stack(channels, axis=channel_axis)


def _id_prefix(path):
    """
    <image_id>_<channel_type>.<ext> -> <image_id>, by dropping the
    filename's extension and its final underscore-separated segment.
    Doesn't hardcode any channel-type word, so it's agnostic to
    "microtubule" vs "microtubules" or other naming drift between what
    a file is actually called and CHANNEL_ORDER's own labels.
    """
    stem = os.path.splitext(os.path.basename(path))[0]
    return stem.rsplit("_", 1)[0]


def check_id_prefix_consistency(channel_paths):
    """
    Check whether the three channel files' inferred <image_id> prefixes
    (see _id_prefix) all agree, catching e.g. a manifest row that
    accidentally paired channels from different cells. Returns
    (consensus_id_or_None, is_consistent). Never raises, a naming
    mismatch is logged as a warning, not treated as a reason to abort.
    """
    prefixes = [_id_prefix(p) for p in channel_paths]
    consistent = len(set(prefixes)) == 1
    if not consistent:
        labeled = dict(zip(CHANNEL_ORDER, prefixes))
        print(f"  [warn] channel filename id mismatch: {labeled}")
    return (prefixes[0] if consistent else None), consistent


def validate_channel_files(
    channel_paths,
    report_id,
    expected_hw,
    expected_dtype,
    channel_axis=-1,
    valid_channel_counts=(3, 4),
    expected_range=(-1.0, 1.0),
    range_tolerance=0.05,
):
    """
    Assemble three separate channel files into one image and run it
    through the same checks _check_array applies to a single combined
    file. report_id labels this row in the report (see
    validate_manifest_channel_images: the manifest's cell_id).
    """
    img = assemble_channels(channel_paths, channel_axis=channel_axis)
    record = _check_array(
        img, report_id, expected_hw, expected_dtype,
        channel_axis=channel_axis, valid_channel_counts=valid_channel_counts,
        expected_range=expected_range, range_tolerance=range_tolerance,
    )
    return img, record


def validate_manifest_channel_images(
    manifest,
    microtubule_col,
    nucleus_col,
    er_col,
    id_col="cell_id",
    expected_hw=None,
    expected_dtype=None,
    csv_path=None,
    only_issues_to_csv=False,
    **kwargs,
):
    """
    Manifest-driven counterpart to validate_image_files: for each row,
    reassemble the three channel files named in microtubule_col /
    nucleus_col / er_col into one multi-channel image, validate it with
    the exact same core the single-file path uses, and return
    (passing_images, report_df) in the same shape validate_image_files
    returns, so callers (run_inference.py) can treat either input style
    the same way downstream.

    manifest is a pandas DataFrame the caller has already loaded; like
    validate_image_files taking a list of filenames rather than reading
    a directory itself, this doesn't read a manifest file, keeping this
    module decoupled from manifest file format (CSV vs TSV, etc).

    Each row's three channel filenames are also checked for a
    consistent <image_id> prefix (check_id_prefix_consistency), recorded
    in the report's id_prefix_consistent column and logged as a warning
    on mismatch. This does not affect passed/issues, which stay governed
    only by the pixel-level checks, per how this was scoped.
    """
    images, records = [], []
    axis = kwargs.get("channel_axis", -1)

    if expected_hw is None or expected_dtype is None:
        first = manifest.iloc[0]
        probe = assemble_channels(
            [first[microtubule_col], first[nucleus_col], first[er_col]], channel_axis=axis
        )
        if expected_hw is None:
            expected_hw = tuple(s for i, s in enumerate(probe.shape) if i != (axis % probe.ndim))
        if expected_dtype is None:
            expected_dtype = probe.dtype

    for _, row in manifest.iterrows():
        channel_paths = [row[microtubule_col], row[nucleus_col], row[er_col]]
        report_id = row[id_col]

        _, id_consistent = check_id_prefix_consistency(channel_paths)

        img, record = validate_channel_files(
            channel_paths, report_id, expected_hw, expected_dtype, **kwargs
        )
        record["id_prefix_consistent"] = id_consistent
        if record["passed"]:
            images.append(img)
        records.append(record)

    report_df = pd.DataFrame(records)
    if csv_path is not None:
        out_df = report_df[~report_df["passed"]] if only_issues_to_csv else report_df
        out_df.to_csv(csv_path, index=False)
        print(f"Wrote {len(out_df)} rows to {csv_path}")

    return images, report_df


def show_predictions_for_cell(
        cell_index, 
        ref_filenames, 
        ref_images_b, 
        images, 
        protein_names, 
        cell_line_names, 
        results,
):
    """Display a reference image alongside all its protein/cell-line predictions.

    Parameters
    ----------
    cell_index : int
        Index into the sorted set of unique reference filenames.
    ref_filenames : list of str
        Filenames corresponding to each entry in `images`/`protein_names`/`cell_line_names`.
    ...

    Returns
    -------
    matplotlib.figure.Figure
        The generated comparison figure.
    """
     
    # Which reference filename are we looking at?
    target_filename = sorted(set(ref_filenames))[cell_index]
    ref_img = ref_images_b[cell_index]

    # Find all predictions matching this reference cell
    matches = [
        (p, c, results.images[i])
        for i, (fn, p, c) in enumerate(zip(ref_filenames, protein_names, cell_line_names))
        if fn == target_filename
    ]

    n = len(matches)
    fig, axes = plt.subplots(1, n + 1, figsize=(3 * (n + 1), 3))

    axes[0].imshow(ref_img)
    axes[0].set_title(f"Reference\n{target_filename}")
    axes[0].axis("off")

    for ax, (protein, cell_line, pred_img) in zip(axes[1:], matches):
        ax.imshow(pred_img)
        ax.set_title(f"{protein}\n{cell_line}")
        ax.axis("off")

    plt.tight_layout()
    plt.show()

    return fig
