"""Multiscale footprint feature helpers.

These helpers provide a lightweight first multiscale scoring backend for
`score-footprints --score multiscale`. They intentionally operate on in-memory
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
