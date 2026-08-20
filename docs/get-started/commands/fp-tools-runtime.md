# [`fp-tools-runtime`](../../api.md#fp-tools-runtime)

Inspect, install, or repair the private external-tool runtime used by raw-read
and de novo motif workflows.

## Example command

```bash
fp-tools-runtime status
```

## Primary inputs

The `status` action takes no input files.

## Main outputs

The command reports each runtime component, platform, installation state, and
cache location. `install core` prepares raw-read tools; MEME Suite and HOMER
components are installed only when requested by their workflows.
