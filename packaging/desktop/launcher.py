"""PyInstaller entry script for fp-tools GUI release bundles."""

import multiprocessing


if __name__ == "__main__":
    multiprocessing.freeze_support()
    from fp_tools.desktop import main

    raise SystemExit(main())
