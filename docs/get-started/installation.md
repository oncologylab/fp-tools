# Installation

fp-tools requires Python 3.12 or later.

## Command-line package

```bash
pip install fp-tools-bio
```

Confirm that the commands are available:

```bash
atac-correct --help
diff-footprints --help
```

## Optional GUI

```bash
pip install "fp-tools-bio[gui]"
fp-tools-gui
```

The GUI saves YAML configurations that can be run directly with
[`run-workflow`](commands/run-workflow.md).

## Source installation

```bash
git clone https://github.com/oncologylab/fp-tools.git
cd fp-tools
python -m pip install -e .
```

See the [API Reference](../api.md) for complete command options.
