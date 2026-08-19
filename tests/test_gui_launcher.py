from __future__ import annotations

import unittest
from pathlib import Path

from fp_tools.cli_gui import _access_urls, _startup_messages


class GuiLauncherTests(unittest.TestCase):
    def test_local_startup_explains_browser_and_ssh_access(self):
        text = "\n".join(_startup_messages("127.0.0.1", 8891, Path("runs")))
        self.assertIn("http://127.0.0.1:8891", text)
        self.assertIn("ssh -N -L 8891:127.0.0.1:8891 USER@SERVER", text)
        self.assertIn("Ctrl+C", text)

    def test_network_mode_warns_about_authentication(self):
        text = "\n".join(_startup_messages("0.0.0.0", 8891, Path("runs")))
        self.assertIn("http://SERVER_IP:8891", text)
        self.assertIn("does not add authentication", text)
        self.assertEqual(
            _access_urls("0.0.0.0", 8891),
            ["http://127.0.0.1:8891", "http://SERVER_IP:8891"],
        )


if __name__ == "__main__":
    unittest.main()
