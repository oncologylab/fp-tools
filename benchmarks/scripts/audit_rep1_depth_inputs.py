#!/usr/bin/env python3
"""Read-only provenance/count audit; never score labels or replace source files."""
import argparse
import json
from pathlib import Path
import sys

import pysam

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_logged_experiment import file_record, write_receipt

DATA = Path('/home/exouser/projects/fp-tools')
SAMPLE = DATA / 'data/public/processed/naked_dna_gse164997/samples/NakedDNA_rep1'
BAMS = {
    'full_filtered': SAMPLE / 'alignment/NakedDNA_rep1.filtered.bam',
    'first10m': DATA / 'benchmarks/results/footprint_functional_v1/naked_dna_10m_s2026/NakedDNA_rep1.first10m.bam',
}


def audit(out):
    out.mkdir(parents=True, exist_ok=False)
    sources = [Path(__file__).resolve(), Path(__file__).with_name('run_logged_experiment.py'),
               SAMPLE / 'state.json', SAMPLE / 'qc/commands.log',
               SAMPLE / 'qc/flagstat.tsv', SAMPLE / 'qc/NakedDNA_rep1.1.fastp.json']
    for bam in BAMS.values():
        sources.extend([bam, Path(str(bam) + '.bai')])
    records = [file_record(p) for p in sources]
    config = write_receipt(out / 'configuration.freeze.json', {
        'inputs': records, 'pysam_version': pysam.__version__,
        'samtools_version': pysam.__samtools_version__,
        'operation': 'quickcheck, indexed counts and complete flagstat; no scoring',
        'flagstat_arguments': ['-@', '4', '-O', 'json'],
        'new_holdout_labels_read': False, 'final_control_read': False,
        'source_mutations': False,
    }, 'configuration_id')
    summaries = {}
    for name, path in BAMS.items():
        print('Auditing ' + name, flush=True)
        pysam.quickcheck(str(path))
        with pysam.AlignmentFile(str(path), 'rb') as bam:
            assert bam.has_index()
            header = bam.header.to_dict()
            assert header['HD']['SO'] == 'coordinate'
            indexed = {'mapped': bam.mapped, 'unmapped': bam.unmapped,
                       'no_coordinate': bam.nocoordinate}
            references = list(zip(bam.references, bam.lengths))
        flags = json.loads(pysam.flagstat('-@', '4', '-O', 'json', str(path)))
        passed = flags['QC-passed reads']
        failed = flags['QC-failed reads']
        assert indexed['mapped'] == passed['mapped'] + failed['mapped']
        assert passed['total'] == passed['read1'] + passed['read2']
        assert passed['secondary'] == passed['supplementary'] == passed['duplicates'] == 0
        summaries[name] = dict(header=header, indexed=indexed,
                               references=references, flagstat=flags)
    assert summaries['full_filtered']['references'] == summaries['first10m']['references']
    state = json.loads((SAMPLE / 'state.json').read_text())
    full = summaries['full_filtered']['flagstat']['QC-passed reads']
    # Record discrepancies; do not silently repair historical metrics.
    metrics = state['metrics']
    discrepancies = {
        'historical_usable_reads': metrics['usable_reads'],
        'actual_mapped_alignment_records': full['mapped'],
        'historical_fragment_count': metrics['fragment_count'],
        'read1_records': full['read1'], 'read2_records': full['read2'],
        'complete_pair_count_verified': False,
    }
    for entry in records:
        assert file_record(entry['path']) == entry, 'input changed during audit'
    result = write_receipt(out / 'audit.json', {
        'configuration_id': config['configuration_id'], 'summaries': summaries,
        'historical_metric_comparison': discrepancies,
        'same_reference_dictionary': True,
        'depth_only_comparison_qualified': False,
        'reason': 'Full alignment includes read-level blacklist and mitochondrial filtering; first10m does not. Deduplication also precedes depth comparison.',
        'final_control_read': False, 'new_holdout_labels_read': False,
        'promotion_eligible': False,
    }, 'stage_id')
    print(json.dumps(discrepancies, sort_keys=True), flush=True)
    print('Result ' + result['stage_id'], flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', type=Path, required=True)
    audit(parser.parse_args().outdir.resolve())
