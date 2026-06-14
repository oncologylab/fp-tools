import argparse
import unittest
from unittest import mock

from fp_tools.parsers import add_atacorrect_arguments, add_bindetect_arguments, add_scorebigwig_arguments
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

    def test_omitted_core_count_uses_all_available_cores_without_warning(self):
        logger = mock.Mock()
        with mock.patch.object(utilities.mp, "cpu_count", return_value=12):
            self.assertEqual(utilities.check_cores(None, logger), 12)
        logger.warning.assert_not_called()

    def test_invalid_core_count_uses_all_available_cores_with_warning(self):
        logger = mock.Mock()
        with mock.patch.object(utilities.mp, "cpu_count", return_value=12):
            self.assertEqual(utilities.check_cores(0, logger), 12)
        logger.warning.assert_called_once()

    def test_public_parser_core_defaults_are_auto(self):
        for builder in (add_atacorrect_arguments, add_scorebigwig_arguments, add_bindetect_arguments):
            parser = builder(argparse.ArgumentParser())
            with self.subTest(builder=builder.__name__):
                args = parser.parse_args([])
                self.assertIsNone(args.cores)


if __name__ == "__main__":
    unittest.main()
