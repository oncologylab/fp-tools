"""Assemble permitted derived evidence, preserving original versus package hashes."""
import argparse
import gzip
from hashlib import sha256
import json
import lzma
from pathlib import Path
import re
import shutil
import subprocess
import numpy as np
import pandas as pd
from package import digest


TEXT={'.md','.py','.json','.jsonl','.log','.tsv','.txt','.bed','.toml'}


def sanitize(text):
    for name,label in [('fp-tools-next-round','RESEARCH_ROOT'),('fp-tools-research','ORIGINAL_RESEARCH_ROOT'),('fp-tools','DATA_ROOT')]:
        text=text.replace('/home/exouser/projects/'+name,'<'+label+'>')
    text=text.replace('/home/exouser','<USER_HOME>')
    if re.search(r'-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----|\bgh[pousr]_[A-Za-z0-9]{30,}|\bAKIA[A-Z0-9]{16}',text):
        raise ValueError('possible secret in selected evidence')
    return text


def execute(source,destination):
    if destination.exists():raise ValueError('fresh staging directory required')
    root=destination/'fp_tools_next_round_evidence';root.mkdir(parents=True)
    source_code=Path(__file__).parent;repository=source_code.parents[2]
    selected=set(pd.read_csv(source/'tables/replay_sites.tsv.gz',sep='\t').site_id)|set(pd.read_csv(source/'tables/nested_sites.tsv.gz',sep='\t').site_id)
    records=[];omitted=[]
    def copy(path,relative,logical):
        if relative=='tables/shape_curves.tsv.gz':relative='tables/shape_curves.tsv.xz'
        target=root/relative;target.parent.mkdir(parents=True,exist_ok=True)
        original_hash=digest(path);transformation='unchanged'
        if path.name=='shape_curves.tsv.gz':
            with gzip.open(path,'rb') as original,lzma.open(target,'wb',preset=6) as packed:
                shutil.copyfileobj(original,packed,length=8388608)
            transformation='lossless gzip-to-xz recompression of complete curve source table; no rows/precision omitted'
        elif path.suffix=='.npz':
            with np.load(path,allow_pickle=False) as z:
                data={k:z[k] for k in z.files}
            # Only the pool-aligned chromatin arrays are compacted. Controls,
            # complete archived replay and simulation evidence remain intact.
            if relative.startswith('arrays/') and ('.replay_pool.' in path.name or re.match(r'(HepG2|K562)_rep[123]\.(full|10M|propensities)\.npz$',path.name)):
                ids=data.get('site_id');mask=np.isin(ids,list(selected))
                if ids is None:raise ValueError('array alignment missing')
                count=len(ids)
                for key,value in data.items():
                    if value.ndim and value.shape[0]==count:data[key]=value[mask]
                transformation=f'complete headline union retained: {int(mask.sum())}/{count} pool sites; nonselected profiles optional omitted payload'
                omitted.append(dict(original_logical_path=logical,original_sha256=original_hash,omitted_sites=int((~mask).sum()),reason='not used by any headline score/shape; complete prematch covariate/label table retained'))
            if relative.startswith('arrays/') and re.match(r'(HepG2|K562)_rep[23]\.propensities\.npz$',path.name):
                shared=path.with_name(path.name.split('_rep')[0]+'_rep1.propensities.npz')
                with np.load(shared,allow_pickle=False) as z:
                    selected_reference=np.isin(z['site_id'],list(selected))
                    np.testing.assert_array_equal(data['site_id'],z['site_id'][selected_reference])
                    for key in ['LOG21_plus','LOG21_minus']:
                        np.testing.assert_array_equal(data[key],z[key][selected_reference]);del data[key]
                data['LOG21_reference']=np.array('arrays/'+shared.name)
                transformation+='; identical LOG21 arrays losslessly shared with '+shared.name
            for key,value in data.items():
                if value.dtype.kind in 'US':
                    if any('/home/' in str(v) for v in value.flat):raise ValueError('private path in non-pickle array')
            np.savez_compressed(target,**data)
        elif path.suffix in TEXT:
            before=path.read_text();after=sanitize(before);target.write_text(after)
            if before!=after:transformation='personal absolute paths replaced by logical roots'
        else:shutil.copyfile(path,target)
        records.append(dict(packaged_path=relative,packaged_sha256=digest(target),original_logical_path=logical,
                            original_sha256=original_hash,transformation=transformation))
    for path in sorted(source.iterdir()):
        if path.is_file() and path.suffix in {'.md','.json'}:copy(path,path.name,'RUN_OUTPUT/'+path.name)
    for folder in ['tables','arrays','figures','models','experiments','historical_context','provenance','logs']:
        for path in sorted((source/folder).rglob('*')):
            if not path.is_file():continue
            if path.is_symlink():raise ValueError('no symlink payload')
            if path.name.endswith('.openat.local.log'):continue
            if path.suffix not in TEXT|{'.gz','.npz','.pdf','.png','.pptx'}:continue
            rel=path.relative_to(source).as_posix();copy(path,rel,'RUN_OUTPUT/'+rel)
    for name in ['verify.py','reanalyse.py','recalibrate.py','shrinkage.py','package.py','test_shrinkage.py']:
        copy(source_code/name,'code/'+name,'RESEARCH_SOURCE/benchmarks/scripts/next_round/'+name)
    copy(repository/'benchmarks/scripts/conditional_protection_score.py','code/conditional_protection_score.py','RESEARCH_SOURCE/benchmarks/scripts/conditional_protection_score.py')
    for path in sorted(source_code.glob('*.py')):
        copy(path,'code/source/benchmarks/scripts/next_round/'+path.name,'RESEARCH_SOURCE/benchmarks/scripts/next_round/'+path.name)
    # Exact scientific operators needed to inspect original figures/extraction.
    for name in ['conditional_protection_score.py','fixed_shape_normalization.py','caliper_matching.py','build_strand_functional_profiles.py','run_logged_experiment.py']:
        copy(repository/'benchmarks/scripts'/name,'code/source/benchmarks/scripts/'+name,'RESEARCH_SOURCE/benchmarks/scripts/'+name)
    for name in ['functional_footprints.py','parametric_bias.py','atacorrect_functions.py']:
        copy(repository/'src/fp_tools/tools'/name,'code/source/src/fp_tools/tools/'+name,'RESEARCH_SOURCE/src/fp_tools/tools/'+name)
    (root/'code/requirements.txt').write_text('numpy==2.4.6\nscipy==1.17.1\npandas==2.3.3\nmatplotlib==3.10.9\nscikit-learn==1.9.0\n')
    records.append(dict(packaged_path='code/requirements.txt',packaged_sha256=digest(root/'code/requirements.txt'),original_logical_path='GENERATED_FROM_RECORDED_ENVIRONMENT',original_sha256=None,transformation='portable reanalysis dependency subset'))
    # An explicitly selected historical snapshot, not recursive data discovery.
    old=repository.parent/'fp-tools-research/benchmarks/results/deep_research_review_20260922_v1/fp_tools_research_review'
    for name in ['README_FIRST.md','DEEP_RESEARCH_PROMPT.md','reports/CONDITIONAL_PROTECTION_20260922.md','reports/CONDITIONAL_SCORE_CEILING_20260922.md','reports/FULL_REP1_REGION_CALIBRATION_20260922.md','reports/FULL_REP1_COUNT_SCREEN_20260922.md','reports/STRAND_EXTRACTION_AUDIT_20260922.md','reports/CONTROL_TRACK_SCOPE_20260922.md','reports/REP1_DEPTH_INPUT_AUDIT_20260922.md']:
        if (old/name).is_file():copy(old/name,'historical_context/original_packet/'+name,'ORIGINAL_PACKET/'+name)
    index=dict(schema='fp-tools-next-round-origin-v1',source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repository,text=True).strip(),
        files=records,optional_profile_omissions=omitted,raw_inputs='not distributed; identities/hashes in sanitized input manifests; no verifier requires raw data',
        absolute_root_mapping={'RESEARCH_ROOT':'isolated research checkout','ORIGINAL_RESEARCH_ROOT':'original read-only research checkout','DATA_ROOT':'original read-only scientific input root','USER_HOME':'redacted personal home'},
        reference_semantics='packaged_path references must resolve locally; original_logical_path and original-source input hashes are provenance references, not bundled raw payloads',
        size_target='approximately250MiB if lossless headline evidence allows; evidence completeness overrides target',
        scientific_status='partial development qualification; no promotion')
    (root/'ORIGIN_INDEX.json').write_text(json.dumps(index,indent=2,sort_keys=True)+'\n')
    lines=[]
    for path in sorted(root.rglob('*')):
        if path.is_file():lines.append(digest(path)+'  '+path.relative_to(root).as_posix())
    (root/'SHA256SUMS').write_text('\n'.join(lines)+'\n')
    print(root,'files',len(lines),'bytes',sum(p.stat().st_size for p in root.rglob('*') if p.is_file()),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--destination',type=Path,required=True)
    a=p.parse_args();execute(a.source.resolve(),a.destination.resolve())
