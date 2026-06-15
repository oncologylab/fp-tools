# Reproducing the fp-tools manuscript

This repository supports two reproducibility paths:

- **Smoke path:** runs on committed fixtures and checks that the package, paper
  scripts, and LaTeX manuscript build correctly.
- **Full public-data path:** reruns the Buenrostro bulk ATAC, 10x PBMC
  pseudobulk, de novo motif validation, and benchmark scaffolds from downloaded
  public data under `data/public/`.

Large public downloads and generated benchmark outputs are not stored in git.
The committed source of truth is the code, manifests, command scripts, and small
source TSVs beside manuscript figures.

## Environment

Use either Conda/Mamba or Docker.

```bash
mamba env create -f environment.yml
mamba activate fp-tools
```

```bash
docker build -t fp-tools:paper .
docker run --rm -it -v "$PWD":/work -w /work fp-tools:paper
```

## Smoke checks

```bash
make test
make paper-smoke
```

The smoke path regenerates selected figure outputs and compiles
`manuscript/main.pdf`.

## Full public-data reruns

The full path expects public data under `data/public/` and may require substantial
CPU, memory, and disk space.

```bash
bash scripts/run_buenrostro_2x2_atac_replicate_demo.sh
bash scripts/run_buenrostro_denovo_motif_validation.sh
.venv/bin/python benchmarks/scripts/prepare_10x_pbmc_pseudobulk.py --write-example-archive
```

After public-data outputs exist, regenerate manuscript figures with the scripts
under `manuscript/scripts/` and compile the paper:

```bash
make paper-pdf
```

## Benchmark manifests

Top-level TSV files in `benchmarks/manifests/` follow the full benchmark manifest
schema documented in `benchmarks/manifests/README.md`. Compact source manifests
used by helper scripts live under `benchmarks/manifests/compact/` and have their
own explicit schemas validated by `benchmarks/scripts/validate_manifests.py`.

## Engineering benchmarks

Use `benchmarks/scripts/run_engineering_benchmark.py` to record runtime and memory
metadata for future fp-tools, TOBIAS, HINT, or other external-tool comparisons:

```bash
python benchmarks/scripts/run_engineering_benchmark.py \
  --label fp-tools-demo \
  --cores 8 \
  --out benchmarks/results/engineering_runtime.tsv \
  -- diff-footprints --help
```
