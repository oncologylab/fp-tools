import unittest

import numpy as np

from fp_tools.utils.normalization import fit_quantile_normalizers, normalize_arrays


class SharedNormalizationTest(unittest.TestCase):
    def test_scaled_arrays_normalize_toward_shared_reference(self):
        base = np.linspace(1.0, 10.0, 100)
        scaled = base * 4.0
        normalized, _, _ = normalize_arrays([base, scaled], ["base", "scaled"], mode="sample-quantile")
        before = abs(np.median(base) - np.median(scaled))
        after = abs(np.median(normalized[0]) - np.median(normalized[1]))
        self.assertLess(after, before)

    def test_none_mode_returns_arrays_unchanged(self):
        arr = np.array([1.0, 2.0, 3.0])
        normalized, objects, diagnostics = normalize_arrays([arr], ["sample"], mode="none")
        np.testing.assert_allclose(normalized[0], arr)
        self.assertEqual(objects, {})
        self.assertEqual(diagnostics, {})

    def test_constant_arrays_fall_back_safely(self):
        objects, _ = fit_quantile_normalizers(
            [np.ones(10), np.ones(10) * 2.0],
            ["a", "b"],
        )
        out = objects["a"].normalize(np.array([1.0, 1.0]))
        self.assertTrue(np.isfinite(out).all())
        self.assertGreater(float(out[0]), 0.0)


if __name__ == "__main__":
    unittest.main()
