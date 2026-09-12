# [`fp-tools-runtime`](../../api.md#fp-tools-runtime)

Check or install the external programs used by `prepare-atac` and
`discover-motifs`. fp-tools manages these programs separately from your other
software. Linux supports read preparation and motif discovery; macOS and
Windows support motif discovery.

## Example command

```bash
fp-tools-runtime status
```

## Primary inputs

The `status` action takes no input files. To install the programs used for motif
discovery before your first run:

```bash
fp-tools-runtime install meme
```

On Linux, use `fp-tools-runtime install core` for the default `prepare-atac`
workflow, or `fp-tools-runtime install homer` for its `homer-atac` profile.

## Main outputs

The command reports each runtime component, platform, installation state, and
cache location. Commands that need a managed component install it on first use,
so manual installation is optional. If an installation is damaged, run
`fp-tools-runtime repair meme` (or substitute the affected component), then
check `fp-tools-runtime status` again.
