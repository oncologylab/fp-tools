#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
export PATH="$ROOT/.venv/bin:$PATH"

mkdir -p dist
find dist -maxdepth 1 -type f -delete
rm -rf wheelhouse
(cd /tmp && "$ROOT/.venv/bin/python" -m build "$ROOT" --outdir "$ROOT/dist")

AUDITWHEEL=""
if [ -x "$ROOT/.venv/bin/auditwheel" ]; then
  AUDITWHEEL="$ROOT/.venv/bin/auditwheel"
elif command -v auditwheel >/dev/null 2>&1; then
  AUDITWHEEL="$(command -v auditwheel)"
fi

if [ -n "$AUDITWHEEL" ]; then
  if ! command -v patchelf >/dev/null 2>&1; then
    echo
    echo "auditwheel is installed, but patchelf is missing from PATH."
    echo "Install it with '.venv/bin/python -m pip install patchelf' or your system package manager, then rerun this script."
    exit 1
  fi
  mkdir -p wheelhouse
  for wheel in dist/*-linux_x86_64.whl; do
    [ -e "$wheel" ] || continue
    "$AUDITWHEEL" repair "$wheel" -w wheelhouse
    rm -f "$wheel"
  done
  if compgen -G "wheelhouse/*.whl" >/dev/null; then
    cp wheelhouse/*.whl dist/
  fi
else
  echo
  echo "auditwheel is not installed; dist/ may contain a linux_x86_64 wheel that PyPI rejects."
  echo "Install auditwheel and patchelf, then rerun this script before publishing Linux wheels."
fi

echo
echo "Built release artifacts:"
ls -lh dist
