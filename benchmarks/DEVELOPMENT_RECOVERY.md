# Continue development in another environment

The recovery policy is **GitHub code, public-data metadata and compact research
results—not a Box or raw-data backup**. Public inputs can be downloaded again;
derived arrays and models require regeneration. No background backup is needed.
The temporary Box backup created on 19 October was cancelled and removed.
Earlier Box research reports and manuscript folders were not changed.

## Get the code

Run in a new parent directory, for example `~/projects`:

```bash
git clone https://github.com/oncologylab/fp-tools.git fp-tools
git -C fp-tools worktree add ../fp-tools-research research/footprint-improvement-20260830
git -C fp-tools worktree add ../fp-tools-next-round research/footprint-next-round-20260922
git -C fp-tools worktree add ../fp-tools-calibration research/footprint-calibration-20261019
git -C fp-tools worktree add ../fp-tools-nuclear-bias research/footprint-nuclear-bias-20261019
```

Keep these sibling worktree names: historical helpers refer to them. Research
checkpoints are `0569e2f` (original), `60ff647` (next round), `690978f`
(calibration), and `d071da8` (completed nuclear experiment). Later nuclear
commits add recovery documentation and metadata. Main was `1d8631e` at this
handoff. Do not merge research methods into main.

## Recover public inputs and reproduce derived results

1. Read `AGENTS.md`, `DEV_PLAN.md`, `NUCLEAR_BIAS_RESULTS_20261019.md` and the
   frozen plans first. The nuclear candidate was rejected; recovery does not
   authorize a new search or protected evaluation.
2. Use tracked `benchmarks/manifests/` for accessions, URLs, checksums, reference
   identities and selection rules. Important entries include
   `frozen_parametric_local_inputs_v1.json`, `direct_bias_replicate_depths_v1.json`,
   `naked_dna_gse164997.tsv` and study-specific freezes.
   `benchmarks/scripts/download_manifest.py --help` describes the downloader.
   Do not blindly download every manifest: some identify protected holdouts.
3. See `benchmarks/recovery/20261019/` for sanitized input provenance, environment
   versions, compact all-task results and evidence-archive identities. These
   small files contain no BAMs, FASTQs, per-site arrays, protected labels or
   credentials.
4. Rebuild an isolated Python environment from project metadata and recorded
   scientific versions. Install external tools through the managed runtime.
   The pinned ChromBPNet source applies to its historical reference workflow.
5. Recreate inputs at recorded relative paths and verify download hashes. Follow
   the original preparation/extraction and stage runners under
   `benchmarks/scripts/next_round/`, `calibration_round/` and `nuclear_round/`.
   Read their plans and command help: not every stage runs directly from BAMs.
6. Preserve historical freezes. Map `<PROJECTS>` to the new project parent in
   a separate local runtime registry. Regenerated artifacts must satisfy the
   recorded identity, coordinate and operator contracts before their results
   can be reused. Do not silently waive changed hashes: record discrepancies
   and establish semantic parity separately.

Compact tables preserve the conclusions without large arrays. Public availability
allows redownload; it does not promise a one-command or bit-identical rebuild.
The staged pipeline still needs compute, matching tool versions and original
filtering rules.

## Scientific boundaries and retained results

The primary nuclear model improved held-chromosome naked-DNA calibration but
failed chromatin accuracy: +0.012710 AUROC and +1.130% relative average precision
versus same-detector DWM propensities, paired interval [-0.006426, 0.030941]. All
21 tasks remain represented; HepG2 HMGA1 is unavailable. This is not a native
TOBIAS benchmark or promotion. DWM means dinucleotide weight matrix; AUROC means
area under the receiver operating characteristic curve.

SK-N-SH/GM23338 labels and naked-DNA replicate 3 remain unopened; replicate 2
is not authorized for this round. Retrieval metadata do not authorize reading
those contents. Main defaults and the published manuscript remain unchanged.

Complete evidence ZIPs remain on the original host only. Their names, sizes and
hashes are in `archive_index.json`; GitHub does not host them or claim that large
scientific artifacts are backed up. Personal/nonpublic data, manuscript files
and locally modified datasets cannot be recovered from public accession metadata;
retain their existing authoritative copies separately. This is the requested
metadata-only migration policy.
