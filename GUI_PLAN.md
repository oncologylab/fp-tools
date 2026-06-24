# GUI Status

## Current State

The browser GUI is implemented and is now part of the normal `fp-tools` user
surface. It remains an optional wrapper around the command-line tools: direct
CLI execution is still the primary interface, and GUI runs save reusable YAML
configs that can be rerun with `run-workflow`.

The current design and layout are accepted as the baseline GUI style:

- clean sidebar navigation grouped by workflow area
- professional full-width content panels with reduced whitespace
- guided tutorial panel for first-time users
- example YAML loading, upload, path loading, editing, and saving
- runnable YAML preview before launch
- background job launch with run-history inspection
- logs, command records, detected outputs, tables, bigWigs, and HTML report links
- static website preview at `docs/demos/gui/fp-tools-gui-static-demo.html`

## Implemented Command Coverage

The GUI exposes the current command-first workflow:

- `atac-correct`
- `call-footprints`
- `match-motifs`
- `diff-footprints`
- `normalize-bigwig`
- `plot-aggregate`
- `plot-aggregate-batch`
- `motif-discovery`
- `motif-summary`
- `fp-tools-score-variants`
- `pseudobulk-fragments`
- `find-signature-fp`
- `pseudobulk-footprints`
- `run-workflow`

Legacy command aliases remain available at the CLI level for compatibility, but
the GUI and public documentation should use the newer command names.

## Design Rules

- Keep the GUI isolated from scientific implementations in `src/fp_tools/tools/`.
- Keep command-line execution and Python tool modules authoritative.
- Keep YAML optional for CLI users, but shared between GUI and `run-workflow`.
- A GUI-saved config must remain runnable from the command line.
- Do not introduce a shared hosted service, user login, or multi-user scheduler
  into the package-level GUI; run one GUI process per user/session.
- Keep examples simple enough for biologists while preserving full command
  traceability.

## Validation Status

Completed validation includes:

- GUI and config module smoke tests through the unittest suite
- command help checks for packaged entry points
- YAML config expansion and dry-run checks for bundled examples
- local MkDocs build for GUI documentation and static GUI preview
- remote-access command documented as `fp-tools-gui --host 0.0.0.0 --port 8891`

The GUI has also been visually reviewed and aligned with the current Fig. 6/GUI
demo style. The current design should be treated as the style baseline for
future GUI changes.

## Remaining Work

Highest-value future improvements:

1. Add deeper GUI-specific integration tests that click through the main pages
   and verify rendered layout states.
2. Add an optional cancel button for long-running background jobs.
3. Expand output preview support for selected HTML, table, and image outputs
   without turning the run history into a separate analysis system.
4. Continue adding advanced CLI options only when they are useful for common
   workflows; avoid crowding pages with rarely used flags.
5. Add richer rerun-from-history support once the saved YAML and command records
   are stable enough to make reruns predictable.

## Maintenance Notes

- Update this file only when the GUI architecture, command coverage, or accepted
  design baseline changes.
- Website deployment instructions belong in `RELEASE_CHECKLIST.md`; do not use
  `mkdocs gh-deploy` for this repository.
