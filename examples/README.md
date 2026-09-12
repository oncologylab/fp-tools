# Examples

Start with [`gui_configs/`](gui_configs/README.md) for small YAML examples you
can run from the command line or load in the GUI. Run repository examples from
the repository root so their relative input paths resolve correctly. For your
own bulk data, follow the
[bulk workflow guide](https://oncologylab.github.io/fp-tools/get-started/workflows/bulk-atac-seq/)
to prepare sample and comparison tables.

Other example inputs and local output locations:

- `atacorrect/`: ignored manual `atac-correct` outputs.
- `scorebigwig/`: ignored manual `call-footprints` outputs.
- `bindetect/`: ignored historical/manual differential-footprint outputs.
- `gui_configs/`: portable YAML examples for `run-yaml-workflow` and the GUI.
- `nutrient_stress_project/`: portable multi-condition project template.
- `plotaggregate_tfbs_dir/` and `plotaggregate_tfbs_grid/`: small plotting inputs.
- `reports/`: standalone generated PDFs or similar reference files.

When creating new manual validation outputs, place them in the matching subdirectory instead of the repository root.
