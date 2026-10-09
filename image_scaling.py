"""
PNG input and output for BLAST-ProtiCelli.

Challenge images are single-channel 8-bit PNGs, 512 x 512, values 0 to
255: three reference channels in, one predicted protein image out.
ProtiCelli works in its own normalized units instead, so each prediction
goes through three steps:

  1. load_reference_channels: read the three PNGs into one raw
     [H, W, 3] uint8 stack, in ProtiCelli's order (microtubules,
     nucleus, ER).
  2. normalize_reference: ProtiCelli's own ImageNormalizer (upstream,
     unchanged) maps the raw stack to [-1, 1]. Every channel is scaled
     against the nucleus channel's clipped maximum, which is returned so
     step 3 can undo it.
  3. reconstruct_prediction: maps ProtiCelli's decoded prediction back
     to 0 to 255 with that same nucleus maximum.

Why step 3 is pred x nucleus_max
--------------------------------
ImageNormalizer divides every channel by the nucleus maximum (for any
channel at least r_floor = 30% as bright as the nucleus), so a
normalized value v in [0, 1] corresponds to raw = v x nucleus_max. It
then rescales to [-1, 1]. ProtiCelli's decoder undoes that last rescale
itself ((decoded / 2 + 0.5), clamped to [0, 1.3]), so its output is
already on the v scale: 1.0 means "as bright as the nucleus maximum".
Values above 1.0 (up to 1.3) are brighter than the nucleus maximum and
saturate at 255.

Upstream's own 8-bit export (PredictionResult.save_prediction) instead
writes pred / 1.3 x 255, a fixed display scale, not a reconstruction to
the input's units.

Open (ask the Lundberg lab): dim proteins
-----------------------------------------
For channels dimmer than 30% of the nucleus, training lifted the
intensity (dim_gamma), so plain inversion makes a dim protein too bright
(README example: true max 30 with nucleus 255 returns near 56). The
protein channel is all zeros at inference, so its true ratio is
unknown. This module uses plain inversion; see the README open items.
"""

import numpy as np
from skimage.io import imread, imsave

from proticelli.data import ChannelAssembler, ImageNormalizer

BIT_DEPTH = 8
MAX_VALUE = 255
IMAGE_SHAPE = (512, 512)

# ImageNormalizer channel layout: 0 microtubules, 1 protein, 2 nucleus, 3 ER.
NUCLEUS_CHANNEL = 2
REFERENCE_CHANNELS = [0, 2, 3]  # what the model is conditioned on


def read_png(path):
    """Read one single-channel 8-bit PNG as a 2-D uint8 array."""
    img = imread(path)
    if img.ndim != 2:
        raise ValueError(f"{path}: expected a single-channel image, got shape {img.shape}")
    if img.dtype != np.uint8:
        raise ValueError(f"{path}: expected 8-bit (uint8) pixels, got {img.dtype}")
    if img.shape != IMAGE_SHAPE:
        raise ValueError(f"{path}: expected {IMAGE_SHAPE[0]} x {IMAGE_SHAPE[1]}, got {img.shape}")
    return img


def load_reference_channels(microtubules_path, nucleus_path, er_path):
    """
    Read the three reference PNGs and assemble them with upstream's
    ChannelAssembler. Returns the raw [H, W, 4] uint8 stack with an
    all-zero protein channel (channel 1), the layout ImageNormalizer
    expects.

    Files are read here rather than by ChannelAssembler, because its
    file loader is tifffile, which rejects PNG. ChannelAssembler accepts
    arrays, so it still enforces matching shapes and channel order.
    """
    return ChannelAssembler(has_protein=False).transform({
        "microtubules": read_png(microtubules_path),
        "nucleus": read_png(nucleus_path),
        "er": read_png(er_path),
    })


def normalize_reference(raw_stack, normalizer=None):
    """
    Normalize a raw [H, W, 4] stack with upstream ImageNormalizer.

    Returns
    -------
    reference : float32 [H, W, 3] in [-1, 1]
        Microtubules, nucleus, ER: the model's conditioning input.
    gains : float32 [4]
        Gain applied to each channel (see ImageNormalizer). Kept for the
        results table; the protein entry is meaningless (empty channel).
    nucleus_max : float
        The nucleus channel's clipped maximum in raw counts: the scale
        reconstruct_prediction multiplies by.
    """
    normalizer = normalizer or ImageNormalizer(bit_depth=BIT_DEPTH)
    normalized, gains = normalizer.transform(raw_stack, return_gains=True)
    nucleus_max = clipped_reference_max(raw_stack, normalizer)

    # Consistency check against upstream: the nucleus is its own reference,
    # so its gain is exactly 1 and its normalized maximum is exactly 1
    # (stored as [-1, 1], so 1.0 there too). If a future ProtiCelli changes
    # the normalizer, this fails loudly instead of silently mis-scaling
    # every output.
    if nucleus_max > 0:
        nuc_norm_max = float(normalized[..., NUCLEUS_CHANNEL].max())
        if not (np.isclose(gains[NUCLEUS_CHANNEL], 1.0) and np.isclose(nuc_norm_max, 1.0, atol=1e-5)):
            raise RuntimeError(
                f"ImageNormalizer no longer scales the nucleus to 1 (gain {gains[NUCLEUS_CHANNEL]}, "
                f"normalized max {nuc_norm_max}); the output reconstruction needs revisiting"
            )

    return normalized[..., REFERENCE_CHANNELS], gains, nucleus_max


def clipped_reference_max(raw_stack, normalizer):
    """
    The nucleus maximum ImageNormalizer scales by, in raw counts.

    ImageNormalizer computes this internally but does not return it, so
    its steps 1 and 2 are repeated here: clip every channel at the
    nucleus's `percentile` (capped at the bit-depth maximum), then take the
    nucleus channel's maximum. normalize_reference checks the result
    against the normalizer's own output.
    """
    img = np.asarray(raw_stack, dtype=np.float32)
    max_val = float(normalizer._BIT_DEPTH_MAX[normalizer.bit_depth])
    ref = normalizer.ref_channel
    ref_clip = min(float(np.percentile(img[..., ref], normalizer.percentile)), max_val)
    if ref_clip <= 0:
        ref_clip = min(float(img.max()), max_val)
    if ref_clip <= 0:
        return 0.0
    return min(float(img[..., ref].max()), ref_clip)


def reconstruct_prediction(prediction, nucleus_max):
    """
    Map a decoded ProtiCelli prediction ([H, W], 0 to 1.3) back to raw
    0 to 255 counts: raw = prediction x nucleus_max, clipped and rounded
    to uint8. See the module docstring for why.
    """
    raw = np.asarray(prediction, dtype=np.float32) * float(nucleus_max)
    return np.rint(np.clip(raw, 0, MAX_VALUE)).astype(np.uint8)


def write_png(path, img):
    """Write a 2-D uint8 image as PNG."""
    if img.dtype != np.uint8 or img.ndim != 2:
        raise ValueError(f"expected a 2-D uint8 image, got {img.dtype} {img.shape}")
    imsave(path, img, check_contrast=False)
