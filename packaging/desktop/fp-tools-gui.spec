# -*- mode: python ; coding: utf-8 -*-
"""Cross-platform one-file fp-tools GUI bundle."""

from pathlib import Path

from PyInstaller.utils.hooks import collect_all, copy_metadata


root = Path.cwd().resolve()
datas = []
binaries = []
hiddenimports = []

# fp-tools dispatches commands lazily, and Streamlit discovers components at
# runtime.  Collect both packages explicitly so every GUI workflow is present.
for package in ("fp_tools", "streamlit"):
    package_data, package_binaries, package_hidden = collect_all(package)
    datas += package_data
    binaries += package_binaries
    hiddenimports += package_hidden

for distribution in ("fp-tools-bio", "streamlit"):
    datas += copy_metadata(distribution, recursive=True)

hiddenimports = [
    module
    for module in hiddenimports
    if not module.startswith("streamlit.hello")
]

# Streamlit executes the application from a source path.  Keeping this source
# file as data is intentional even though the module is also frozen.
datas.append((str(root / "src" / "fp_tools" / "gui_app.py"), "fp_tools"))

analysis = Analysis(
    [str(root / "packaging" / "desktop" / "launcher.py")],
    pathex=[str(root / "src")],
    binaries=binaries,
    datas=datas,
    hiddenimports=sorted(set(hiddenimports)),
    hookspath=[],
    runtime_hooks=[],
    excludes=[
        "pytest",
        "sphinx",
        "mkdocs",
        "streamlit.hello",
        "plotly",
        "kaleido",
        "numba",
        "llvmlite",
        "tensorflow",
        "torch",
        "cupy",
        "dask",
    ],
    noarchive=False,
)
pyz = PYZ(analysis.pure)
executable = EXE(
    pyz,
    analysis.scripts,
    analysis.binaries,
    analysis.datas,
    [],
    name="fp-tools-gui",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    console=True,
    disable_windowed_traceback=False,
)
