# Contributing to fp-tools

Thank you for helping improve fp-tools. Contributions should preserve the
command-first architecture and reproducible scientific outputs.

## Before opening a change

1. Search existing issues and open a focused issue for substantial behavioral
   changes.
2. Create a branch from `main`.
3. Keep CLI wrappers thin and place scientific implementations under
   `src/fp_tools/tools/` or shared helpers under `src/fp_tools/utils/`.
4. Add tests proportional to the risk of the change.

## Local validation

Run from the repository root:

```bash
.venv/bin/python -m pip check
.venv/bin/python -m unittest discover -s tests -v
.venv/bin/python scripts/smoke_console_scripts.py
.venv/bin/mkdocs build --clean --strict
git diff --check
```

Do not commit public datasets, generated run directories, local agent state, or
manuscript files. See `AGENTS.md` and `RELEASE_CHECKLIST.md` for the complete
repository and release conventions.

## Pull requests

Describe the user-visible behavior, tests performed, and compatibility impact.
Documentation and examples must use the public command names declared in
`pyproject.toml`.
