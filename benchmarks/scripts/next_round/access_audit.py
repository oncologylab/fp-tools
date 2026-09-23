"""Audit existing receipts/traces, never follow traced paths into data files."""
import argparse
from collections import Counter
import json
from pathlib import Path
import re
from access import DENY,digest
from prepare import dump


def execute(out):
    receipts=[json.loads(line) for line in (out/'logs/access.jsonl').read_text().splitlines()]
    forbidden=DENY+['NakedDNA_rep2','naked_rep2','SRR13450941']
    for row in receipts:
        if any(word.lower() in json.dumps(row).lower() for word in forbidden):
            raise AssertionError('forbidden identity in a recorded input read')
    traces=[];violations=[]
    for path in sorted((out/'logs').glob('*.openat.local.log')):
        counts=Counter();successful=0
        with path.open() as handle:
            for line in handle:
                match=re.search(r'openat\([^,]+, "([^"]+)"',line)
                if not match or re.search(r'= -1\b',line):continue
                name=match.group(1);successful+=1
                if any(w.lower() in name.lower() for w in forbidden):violations.append(dict(trace=path.name,path=name))
                if name.endswith(('.bam','.bai','.fa','.fai','.bw','.npz','.bed','.pickle','.tsv')):counts[Path(name).suffix]+=1
        traces.append(dict(trace=path.name,original_sha256=digest(path),successful_or_resumed_open_records=successful,
                           scientific_suffix_counts=dict(counts),packaged=False,reason='raw paths excluded; metadata summary only'))
    if violations:raise AssertionError('protected-path open detected: '+repr(violations))
    dump(out/'provenance/access_scope_audit.json',dict(receipt_count=len(receipts),
        recorded_input_ids=sorted({r['input_id'] for r in receipts}),trace_summaries=traces,
        prohibited_recorded_reads=0,scope='explicit input allowlists and available syscall traces; not a claim of whole-host monitoring',
        labels='K562/HepG2 chr1-15 permitted historical candidate pool only',
        BAM_scope='all reference contigs streamed for complete-pair inventory/depth selection; only training-locus arrays and labels analyzed',
        controls='examined GSE164997 replicate1 only; replicate2 and final replicate3 contents not used',
        protected_content_packaged=False,raw_payloads_packaged=False,
        limitation='syscall tracing applies to recorded extraction/propensity/reference jobs; other analysis uses only generated allowlisted derived files'))
    print('Access receipts and',len(traces),'traces checked; no recorded prohibited reads',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--outdir',type=Path,required=True);execute(p.parse_args().outdir.resolve())
