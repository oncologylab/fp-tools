#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

mkdir -p dist
find dist -maxdepth 1 -type f -delete
(cd /tmp && "$ROOT/.venv/bin/python" -m build "$ROOT" --outdir "$ROOT/dist")

echo
echo "Built release artifacts:"
ls -lh dist
