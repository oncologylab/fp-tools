import unittest
from unittest import mock

from fp_tools.utils import utilities


class CoreHandlingTest(unittest.TestCase):
    def test_respects_valid_requested_core_count(self):
        logger = mock.Mock()
        with mock.patch.object(utilities.mp, "cpu_count", return_value=16):
            self.assertEqual(utilities.check_cores(4, logger), 4)
        logger.warning.assert_not_called()

    def test_caps_requested_cores_to_available_cores(self):
        logger = mock.Mock()
        with mock.patch.object(utilities.mp, "cpu_count", return_value=8):
            self.assertEqual(utilities.check_cores(32, logger), 8)
        logger.warning.assert_called_once()

    def test_invalid_core_count_uses_all_available_cores(self):
        logger = mock.Mock()
        with mock.patch.object(utilities.mp, "cpu_count", return_value=12):
            self.assertEqual(utilities.check_cores(0, logger), 12)
            self.assertEqual(utilities.check_cores(None, logger), 12)
        self.assertEqual(logger.warning.call_count, 2)


if __name__ == "__main__":
    unittest.main()
