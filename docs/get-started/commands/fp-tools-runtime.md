# [`fp-tools-runtime`](../../api.md#fp-tools-runtime)

Inspect, install, or repair the private external-tool runtime. Linux provides
raw-read and de novo motif components; macOS and Windows provide the optional
de novo motif component.

## Example command

```bash
fp-tools-runtime status
```

## Primary inputs

The `status` action takes no input files.

## Main outputs

The command reports each runtime component, platform, installation state, and
cache location. `install core` and `install homer` are Linux-only raw-read
components. The MEME Suite component is installed only when requested by de
novo motif discovery.
