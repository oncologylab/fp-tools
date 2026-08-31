"""Multiscale footprint feature helpers.

These helpers provide a lightweight first multiscale scoring backend for
`call-footprints --score multiscale`. They intentionally operate on in-memory
1D arrays so the command can reuse the existing bigWig region processing path.
"""

from __future__ import annotations

import numpy as np

DEFAULT_SCALES = (8, 16, 24, 32, 64, 100, 147)


def parse_scales(values: list[int] | tuple[int, ...] | None) -> tuple[int, ...]:
    """Validate and normalize multiscale window sizes."""

    if values is None:
        return DEFAULT_SCALES
    scales = tuple(sorted({int(value) for value in values}))
    if not scales:
        raise ValueError("At least one multiscale window size is required.")
    invalid = [scale for scale in scales if scale < 3]
    if invalid:
        raise ValueError(f"Multiscale window sizes must be >= 3 bp: {invalid}")
    return scales


def _window_mean(signal: np.ndarray, start: int, end: int) -> float:
    start = max(0, start)
    end = min(len(signal), end)
    if end <= start:
        return np.nan
    return float(np.nanmean(signal[start:end]))


def multiscale_depletion(signal: np.ndarray, scales: list[int] | tuple[int, ...] | None = None) -> dict[int, np.ndarray]:
    """Return central-depletion scores for each scale.

    For each position and scale, the score is the mean of the left and right
    flanking windows minus the mean of the central window. Higher values indicate
    stronger local depletion relative to the flanks.
    """

    arr = np.nan_to_num(np.asarray(signal, dtype=float), nan=0.0)
    scales = parse_scales(scales)
    features: dict[int, np.ndarray] = {}
    for scale in scales:
        half = max(1, scale // 2)
        scores = np.zeros(len(arr), dtype=float)
        for idx in range(len(arr)):
            center = _window_mean(arr, idx - half, idx + half + 1)
            left = _window_mean(arr, idx - scale - half, idx - half)
            right = _window_mean(arr, idx + half + 1, idx + scale + half + 1)
            flank_values = [value for value in (left, right) if not np.isnan(value)]
            flank = float(np.mean(flank_values)) if flank_values else 0.0
            scores[idx] = flank - center
        features[scale] = scores
    return features


def summarize_multiscale(features: dict[int, np.ndarray], method: str = "max") -> np.ndarray:
    """Collapse scale-specific arrays into one summary track."""

    if not features:
        return np.array([], dtype=float)
    matrix = np.vstack([features[scale] for scale in sorted(features)])
    if method == "max":
        return np.nanmax(matrix, axis=0)
    if method == "mean":
        return np.nanmean(matrix, axis=0)
    raise ValueError(f"Unsupported multiscale summary method: {method}")


def symmetric_depletion(
    signal: np.ndarray,
    center_width: int = 33,
    flank_width: int = 32,
    noise_floor: float = 1e-3,
) -> np.ndarray:
    """Return a locally standardized, bilaterally symmetric depletion track.

    Each valid position compares one central window with immediately adjacent
    left and right shoulders. The contrast is divided by the square root of
    local mean absolute signal so high-noise regions do not dominate solely
    because of signal magnitude. Invalid edge positions are zero.
    """

    array = np.nan_to_num(np.asarray(signal, dtype=float), nan=0.0)
    center_width = int(center_width)
    flank_width = int(flank_width)
    noise_floor = float(noise_floor)
    if center_width < 3:
        raise ValueError("center_width must be >= 3")
    if flank_width < 1:
        raise ValueError("flank_width must be >= 1")
    if noise_floor <= 0:
        raise ValueError("noise_floor must be > 0")

    left_half = center_width // 2
    right_half = center_width - left_half
    first = left_half + flank_width
    last = len(array) - right_half - flank_width
    scores = np.zeros(len(array), dtype=float)
    if last < first:
        return scores

    centers = np.arange(first, last + 1, dtype=int)
    prefix = np.concatenate(([0.0], np.cumsum(array, dtype=float)))
    absolute_prefix = np.concatenate(([0.0], np.cumsum(np.abs(array), dtype=float)))

    center_start = centers - left_half
    center_end = centers + right_half
    left_start = center_start - flank_width
    right_end = center_end + flank_width
    center_mean = (prefix[center_end] - prefix[center_start]) / center_width
    left_mean = (prefix[center_start] - prefix[left_start]) / flank_width
    right_mean = (prefix[right_end] - prefix[center_end]) / flank_width
    local_width = center_width + 2 * flank_width
    local_absolute_mean = (absolute_prefix[right_end] - absolute_prefix[left_start]) / local_width
    scores[centers] = ((left_mean + right_mean) / 2.0 - center_mean) / np.sqrt(local_absolute_mean + noise_floor)
    return scores


def hybrid_footprint_score(
    signal: np.ndarray,
    legacy_score: np.ndarray,
    center_width: int = 33,
    flank_width: int = 32,
    weight: float = 0.2,
    noise_floor: float = 1e-3,
) -> np.ndarray:
    """Add a low-weight wide symmetric channel to the existing footprint score."""

    legacy = np.asarray(legacy_score, dtype=float)
    if legacy.shape != np.asarray(signal).shape:
        raise ValueError("signal and legacy_score must have identical shapes")
    if weight < 0:
        raise ValueError("weight must be >= 0")
    return legacy + float(weight) * symmetric_depletion(
        signal,
        center_width=center_width,
        flank_width=flank_width,
        noise_floor=noise_floor,
    )

def trim_multiscale_features(features: dict[int, np.ndarray], flank: int) -> dict[int, np.ndarray]:
    """Trim flank bases from every scale-specific feature array."""

    if flank <= 0:
        return {scale: values.copy() for scale, values in features.items()}
    return {scale: values[flank:-flank] for scale, values in features.items()}


def write_multiscale_npz(
    path: str,
    records: list[tuple[tuple[str, int, int], dict[int, np.ndarray]]],
    scales: list[int] | tuple[int, ...],
    summary_method: str,
) -> None:
    """Write multiscale per-region features to a compressed NumPy sidecar.

    The saved tensor has shape ``n_scales x total_positions``. Region-level
    offsets map each output region to its columns in the concatenated tensor.
    """

    scales = parse_scales(scales)
    chroms: list[str] = []
    starts: list[int] = []
    ends: list[int] = []
    offsets = [0]
    matrices = []

    for region, features in records:
        chrom, start, end = region
        matrix = np.vstack([np.asarray(features[scale], dtype=np.float32) for scale in scales])
        chroms.append(str(chrom))
        starts.append(int(start))
        ends.append(int(end))
        offsets.append(offsets[-1] + matrix.shape[1])
        matrices.append(matrix)

    tensor = np.concatenate(matrices, axis=1) if matrices else np.zeros((len(scales), 0), dtype=np.float32)
    np.savez_compressed(
        path,
        tensor=tensor,
        scales=np.asarray(scales, dtype=np.int32),
        chroms=np.asarray(chroms, dtype=str),
        starts=np.asarray(starts, dtype=np.int64),
        ends=np.asarray(ends, dtype=np.int64),
        offsets=np.asarray(offsets, dtype=np.int64),
        summary_method=np.asarray(summary_method),
        format_version=np.asarray("fp-tools-multiscale-npz-v1"),
    )


def load_multiscale_npz(path: str) -> dict[str, np.ndarray]:
    """Load a multiscale NPZ sidecar into plain NumPy arrays."""

    with np.load(path, allow_pickle=False) as data:
        return {key: data[key].copy() for key in data.files}


def aggregate_multiscale_tensor(data: dict[str, np.ndarray], align: str = "center") -> np.ndarray:
    """Return a scale-by-position average across regions from loaded NPZ data."""

    tensor = data["tensor"]
    offsets = data["offsets"]
    lengths = np.diff(offsets)
    if len(lengths) == 0:
        return np.zeros((tensor.shape[0], 0), dtype=float)
    if align not in {"center", "left"}:
        raise ValueError("align must be 'center' or 'left'")

    max_len = int(np.max(lengths))
    stacks = np.full((len(lengths), tensor.shape[0], max_len), np.nan, dtype=float)
    for idx, (start, end) in enumerate(zip(offsets[:-1], offsets[1:])):
        region_tensor = tensor[:, start:end]
        left = 0 if align == "left" else int((max_len - region_tensor.shape[1]) // 2)
        stacks[idx, :, left : left + region_tensor.shape[1]] = region_tensor
    return np.nanmean(stacks, axis=0)
