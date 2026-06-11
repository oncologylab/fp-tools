"""Shared signal normalization helpers for BINDetect and PlotAggregate."""

from __future__ import annotations

import numpy as np
from scipy import interpolate
from scipy.optimize import curve_fit

from fp_tools.utils.signals import fast_rolling_math


def sigmoid(x, a, b, L, shift):
    exponent = np.clip(-b * (x - a), -700, 700)
    return L / (1 + np.exp(exponent)) + shift


class ArrayNorm:
    """Multiplicative normalization curve fitted from empirical quantiles."""

    def __init__(self, function, popt, value_min, value_max):
        self.func = function
        self.popt = popt
        self.value_min = value_min
        self.value_max = value_max

    def get_norm_factor(self, arr):
        if self.func == "sigmoid":
            return sigmoid(arr, *self.popt)
        if self.func == "constant":
            return arr * 0 + self.popt
        raise ValueError("Unknown normalization func")

    def normalize(self, arr):
        arr = np.asarray(arr, dtype=float)
        arr_capped = np.where(arr > self.value_max, self.value_max, arr)
        arr_capped = np.where(arr_capped < self.value_min, self.value_min, arr_capped)
        return arr * self.get_norm_factor(arr_capped)


def _positive_values(arr):
    arr = np.asarray(arr, dtype=float)
    arr = arr[np.isfinite(arr)]
    arr = arr[arr > 0]
    return arr if arr.size else np.array([0.0], dtype=float)


def fit_quantile_normalizers(list_of_arrays, names, logger=None, quantiles=None):
    """Fit BINDetect-style multiplicative quantile normalizers."""

    if quantiles is None:
        quantiles = np.linspace(0.05, 0.99, 1000, endpoint=True)
    arrays = [_positive_values(arr) for arr in list_of_arrays]
    array_quantiles = [np.quantile(arr, quantiles) for arr in arrays]
    mean_array_quantiles = np.mean(np.vstack(array_quantiles), axis=0)
    norm_objects = {}
    diagnostics = {
        "quantiles": quantiles,
        "array_quantiles": array_quantiles,
        "mean_array_quantiles": mean_array_quantiles,
    }

    for i, name in enumerate(names):
        xdata = np.asarray(array_quantiles[i], dtype=float)
        denominator = np.where(np.isclose(xdata, 0.0), np.nan, xdata)
        ydata = np.asarray(mean_array_quantiles, dtype=float) / denominator
        finite = np.isfinite(xdata) & np.isfinite(ydata)
        xdata = xdata[finite]
        ydata = ydata[finite]
        if xdata.size < 4 or np.isclose(np.max(xdata), np.min(xdata)):
            factor = 1.0 if not np.any(np.isfinite(ydata)) else float(np.nanmean(ydata))
            if not np.isfinite(factor) or factor <= 0:
                factor = 1.0
            norm_objects[name] = ArrayNorm("constant", factor, value_min=0.0, value_max=max(1.0, float(np.max(arrays[i]))))
            continue

        pad = min(50, max(1, xdata.size // 10))
        ydata_pad = np.pad(ydata, pad, "edge")
        ydata_smooth = fast_rolling_math(ydata_pad, pad, "mean")[pad:-pad]

        p = interpolate.interp1d(xdata, ydata_smooth, kind="linear", fill_value="extrapolate")
        xvals = np.linspace(np.min(xdata), np.max(xdata), 1000)
        y_inter = p(xvals)
        mask = np.isfinite(y_inter)
        xvals, y_inter = xvals[mask], y_inter[mask]
        try:
            a_range = (np.min(xdata), np.max(xdata))
            L_range = np.diff(np.percentile(y_inter, [5, 95]))[0]
            bounds = ((a_range[0], -np.inf, 0, 0), (a_range[1], np.inf, max(L_range, 1e-12), np.inf))
            popt, _ = curve_fit(sigmoid, xvals, y_inter, bounds=bounds)
            func = "sigmoid"
        except Exception as exc:
            popt = float(np.mean(y_inter))
            if logger is not None:
                logger.warning(
                    f"Curve-fitting quantile normalization failed for '{name}'. "
                    f"Falling back to constant factor {popt:.2f}. Error was: {exc}"
                )
            func = "constant"
        norm_objects[name] = ArrayNorm(func, popt, value_min=float(np.min(xdata)), value_max=float(np.max(xdata)))
    return norm_objects, diagnostics


def normalize_arrays(list_of_arrays, names, mode="sample-quantile", logger=None):
    """Return normalized arrays and fitted objects for a supported mode."""

    mode = (mode or "none").replace("_", "-")
    if mode == "none":
        return [np.asarray(arr, dtype=float) for arr in list_of_arrays], {}, {}
    if mode not in {"condition-quantile", "sample-quantile"}:
        raise ValueError(f"Unsupported normalization mode: {mode}")
    norm_objects, diagnostics = fit_quantile_normalizers(list_of_arrays, names, logger=logger)
    normalized = [np.maximum(0.0, norm_objects[name].normalize(arr)) for name, arr in zip(names, list_of_arrays)]
    return normalized, norm_objects, diagnostics
