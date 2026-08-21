# -*- mode: python ; coding: utf-8 -*-
"""Cross-platform one-file fp-tools GUI bundle."""

from importlib.util import find_spec
from pathlib import Path
import sys

from PyInstaller.utils.hooks import collect_all, collect_data_files, copy_metadata


root = Path.cwd().resolve()
icon_suffix = ".ico" if sys.platform == "win32" else ".icns"
icon_path = root / "build" / "desktop-icons" / f"fp-tools{icon_suffix}"
if not icon_path.is_file():
    raise RuntimeError(f"Desktop icon has not been built: {icon_path}")
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

# Cython modules are imported lazily by several command targets. Resolve their
# installed extension files explicitly so one-file bundles cannot silently
# fall back to the source-only package tree.
for module in (
    "fp_tools.utils.sequences",
    "fp_tools.utils.ngs",
    "fp_tools.utils.signals",
):
    module_spec = find_spec(module)
    if module_spec is None or not module_spec.origin:
        raise RuntimeError(f"Desktop build is missing compiled extension {module}")
    binaries.append((module_spec.origin, "fp_tools/utils"))
    hiddenimports.append(module)

# Windows uses bamnostic in place of pysam. bamnostic reads its version from a
# package-data file during import, including for command --help.
if find_spec("bamnostic") is not None:
    datas += collect_data_files("bamnostic")
    hiddenimports.append("bamnostic")

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
datas.append(
    (
        str(root / "docs" / "assets" / "fp_tools_logo_icon_1024.png"),
        "fp_tools/resources",
    )
)

analysis = Analysis(
    [str(root / "packaging" / "desktop" / "launcher.py")],
    pathex=[str(root / "src")],
    binaries=binaries,
    datas=datas,
    hiddenimports=sorted(set(hiddenimports)),
    hookspath=[],
    runtime_hooks=[],
    excludes=[
        "fp_tools.tools.prepare_atac",
        "fp_tools.tools.prepare_atac_legacy",
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
    hide_console="hide-early",
    icon=str(icon_path),
    disable_windowed_traceback=False,
)
