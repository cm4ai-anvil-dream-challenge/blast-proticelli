#!/usr/bin/env python3
"""
split_reference_channels.py

Splits a single multi-channel reference TIFF (channel order: microtubule,
nucleus, ER, per the source data) into three separate single-channel TIFF
files, matching the manifest's microtubules_image / nucleus_image /
er_image columns. Output is TIFF, not PNG, since that's ProtiCelli's
current actual input format; this may change to PNG later if formatted
challenge data arrives in that shape instead, at which point only the
output extension/writer needs to change, not the channel-splitting logic.

Assumptions (stated explicitly since this has not been run against a
real ProtiCelli reference image, only a synthetic test fixture):
  - The TIFF holds exactly 3 channels/pages, in the stated order
    (microtubule, nucleus, ER).
  - Channels may be stored as the first axis (shape (3, H, W), typical
    for a multi-page TIFF read as a stack) or the last axis (shape
    (H, W, 3), typical for an RGB-style image). Both are handled; any
    other shape raises rather than guessing which axis is which.
  - Original bit depth (e.g. 16-bit grayscale, common for microscopy) is
    preserved on output, not silently downcast to 8-bit.
  - If the image isn't 512x512, this prints a warning but still writes
    the output at whatever size it actually is, rather than resizing
    unasked.

Usage:
    python split_reference_channels.py input.tiff output_dir

Output: <output_dir>/<input_stem>_microtubules.tiff, <input_stem>_nucleus.tiff,
<input_stem>_er.tiff, where <input_stem> is the input filename without its
extension (e.g. cell_42.tiff -> cell_42_microtubules.tiff), so outputs stay
traceable back to their source file.
"""

import argparse
from pathlib import Path

from skimage import io as skio

CHANNEL_ORDER = ["microtubules", "nucleus", "er"]
EXPECTED_SIZE = 512


def split_channels(image):
    """
    Return the 3 channels as separate 2D arrays, in CHANNEL_ORDER, from
    whichever axis actually holds them. Raises rather than guessing if
    the shape doesn't clearly indicate 3 channels somewhere.
    """
    if image.ndim == 2:
        raise ValueError(
            f"Image is already single-channel (shape {image.shape}), nothing "
            "to split, expected a 3-channel stack."
        )
    if image.ndim != 3:
        raise ValueError(f"Expected a 2D or 3D array, got shape {image.shape}")

    if image.shape[0] == 3:
        return [image[i] for i in range(3)]
    if image.shape[-1] == 3:
        return [image[..., i] for i in range(3)]

    raise ValueError(
        f"Could not find a length-3 channel axis in shape {image.shape}, "
        "expected either (3, H, W) or (H, W, 3)."
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_tiff", help="Path to the multi-channel reference TIFF")
    parser.add_argument("output_dir", help="Directory to write the three split TIFFs into")
    args = parser.parse_args()

    input_path = Path(args.input_tiff)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    prefix = input_path.stem

    image = skio.imread(str(input_path))
    print(f"Loaded {input_path} with shape {image.shape}, dtype {image.dtype}")

    channels = split_channels(image)

    for name, channel in zip(CHANNEL_ORDER, channels):
        if channel.shape != (EXPECTED_SIZE, EXPECTED_SIZE):
            print(
                f"  [warn] {name} channel is {channel.shape}, expected "
                f"({EXPECTED_SIZE}, {EXPECTED_SIZE}), writing it as-is anyway"
            )
        out_path = output_dir / f"{prefix}_{name}.tiff"
        skio.imsave(str(out_path), channel, check_contrast=False)
        print(f"  wrote {out_path}  shape={channel.shape} dtype={channel.dtype}")


if __name__ == "__main__":
    main()