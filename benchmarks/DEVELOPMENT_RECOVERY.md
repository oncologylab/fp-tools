# Restore fp-tools development in another environment

Code is on GitHub; scientific data are not stored in Git. The migration backup
destination is `yilab:Yaoxiang/fp-tools/development_backup_2026-10-19/`.
Authenticate rclone to your own Box account in the new environment; credentials
are deliberately not included in the backup.

## Check backup completion first

```bash
rclone cat yilab:Yaoxiang/fp-tools/development_backup_2026-10-19/recovery/status.json
```

Only `COMPLETE_VERIFIED` means the entire allowlisted workspace copy passed its
remote checksum check and the source inventory remained unchanged. `RUNNING`,
`VERIFYING`, a missing status, or `FAILED_INCOMPLETE` does not mean a complete
backup. Do not decommission the original host until completion is confirmed.
The `evidence/` archive is separately useful even while the large transfer runs.

The inventory covers five named project trees. It excludes credentials, virtual
environments, caches, rebuildable build/site outputs and protected holdout
contents. Symlinks are inventoried but not followed; review `symlinks.json` when
restoring. `exclusions.json` is the exact omission ledger. This is not a full
machine backup. No protected contents have been opened to construct the list.

## Restore code and working directories

Run these commands in a new parent directory, such as `~/projects`:

```bash
git clone https://github.com/oncologylab/fp-tools.git fp-tools
git -C fp-tools worktree add ../fp-tools-research research/footprint-improvement-20260830
git -C fp-tools worktree add ../fp-tools-next-round research/footprint-next-round-20260922
git -C fp-tools worktree add ../fp-tools-calibration research/footprint-calibration-20261019
git -C fp-tools worktree add ../fp-tools-nuclear-bias research/footprint-nuclear-bias-20261019
rclone copy yilab:Yaoxiang/fp-tools/development_backup_2026-10-19/recovery ./backup_metadata --checksum
rclone copy yilab:Yaoxiang/fp-tools/development_backup_2026-10-19/parts ./backup_parts --checksum
python - <<'PY'
import hashlib, json
from pathlib import Path
status=json.loads(Path('backup_metadata/status.json').read_text())
assert status['state']=='COMPLETE_VERIFIED', status
for record in json.loads(Path('backup_metadata/parts.json').read_text()):
    p=Path('backup_parts')/record['name'];h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
    assert p.stat().st_size==record['bytes'] and h.hexdigest()==record['sha256'],p
print('All downloaded parts verified')
PY
cat backup_parts/workspace.tar.part-* | tar --no-same-owner -xf -
```

The ordered tar parts preserve relative paths, including data, results,
manuscript material and examples. Git metadata and environments are excluded;
recreate environments using the pinned project/research dependency records.
The manuscript is a private backup only and must not be committed to Git.
Old manifests may contain original absolute paths: map the old project parent
to the new parent explicitly, preserving the original manifests and hashes.
Do not silently rewrite scientific provenance or rerun stages with changed inputs.

For exact historical code, the research checkpoints at backup preparation were
`0569e2f` (original research), `60ff647` (next round), `690978f` (calibration), and
`d071da8` (completed nuclear results); subsequent nuclear commits add recovery
instructions only. Main was `1d8631e`. Verify branch identities before running.

## Recover and verify the self-contained evidence first

```bash
rclone copy yilab:Yaoxiang/fp-tools/development_backup_2026-10-19/evidence ./evidence --checksum
unzip evidence/fp_tools_next_round_evidence_2026-10-19_1473f2d.zip -d evidence/unpacked
cd evidence/unpacked/fp_tools_next_round_evidence
python -B code/verify.py .
python -B code/verify_chromatin.py .
```

Install `code/requirements.txt` in a new isolated Python environment first.
The ZIP SHA-256 is
`8953b12dd95b4d90aadf4e6b525fc6afdab9da02dbc4737f0f930bbd7fe1cc7f`;
size is 374,477,584 bytes. This evidence permits primary reanalysis without BAMs.
It does not contain every historical raw input. The nuclear candidate failed
the chromatin screen; restore its negative decision with its results.

## Resume an interrupted workspace backup on the original host

```bash
python benchmarks/scripts/archive_development_backup.py \
  --base /home/exouser/projects \
  --out benchmarks/results/development_backup_20261019 \
  --remote yilab:Yaoxiang/fp-tools/development_backup_2026-10-19
```

Run from the nuclear research worktree. The runner uses a process lock, rclone
checksum comparison, immutable remote objects, retries, and a final remote
part-set/checksum audit. Two-GiB parts avoid hundreds of thousands of individual
Box upload requests; only a verified temporary transport part is removed locally.
Resuming regenerates the deterministic tar stream and reuses identical remote
parts. It never deletes scientific inputs or remote data. If source
files changed after inventory, it refuses to claim completion; preserve the
old attempt and make a new versioned inventory/destination for the changed
snapshot. Logs and final status are copied to Box recovery metadata.
