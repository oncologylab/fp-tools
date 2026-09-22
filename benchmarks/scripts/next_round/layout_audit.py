"""Isolate native expected-track dependence on the output-region layout."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import pyBigWig
from prepare import dump
from access import digest


def execute(out):
    sites=pd.read_csv(out/'tables/control_sites.tsv.gz',sep='\t')
    support=pd.read_csv(out/'tables/dwm_support.v3.tsv',sep='\t')
    path=out/'scratch/dwm_original_layout_replay/NakedDNA_rep1_all_sites_expected.bw'
    with np.load(out/'arrays/control_production_DWM.v3.npz',allow_pickle=False) as z:
        old=z['old_expected'];expanded=z['original_expected']
    restored=np.full_like(old,np.nan);rows=[]
    with pyBigWig.open(str(path)) as bw:
        for i,row in enumerate(sites.itertuples()):
            if not support.old_positive_overlap.iloc[i]:continue
            values=np.nan_to_num(bw.values(row.TFBS_chr,int(row.center)-100,int(row.center)+101,numpy=True))
            if row.TFBS_strand=='-':values=values[::-1]
            restored[i]=values
            a=old[i]/old[i].sum();b=values/values.sum();c=expanded[i]/expanded[i].sum()
            rows.append(dict(site_id=row.site_id,chrom=row.TFBS_chr,center=row.center,
                restored_layout_max_difference=float(abs(a-b).max()),expanded_layout_max_difference=float(abs(a-c).max()),
                restored_layout_parity=bool(np.allclose(a,b,atol=1e-7,rtol=1e-5)),
                expanded_layout_parity=bool(np.allclose(a,c,atol=1e-7,rtol=1e-5))))
    frame=pd.DataFrame(rows);frame.to_csv(out/'tables/dwm_layout_isolation.tsv',sep='\t',index=False)
    np.savez_compressed(out/'arrays/control_DWM_restored_layout.npz',site_id=sites.site_id.to_numpy(str),expected=restored)
    dump(out/'provenance/dwm_layout_isolation_result.json',dict(code_sha256=digest(__file__),sites=len(frame),
        restored_mismatches=int((~frame.restored_layout_parity).sum()),expanded_mismatches=int((~frame.expanded_layout_parity).sum()),
        inference='Region layout is the isolated intervention; agreement after restoring the historical layout supports layout dependence, not expanded-layout parity.',
        safety_qualification=False))
    print(frame[['restored_layout_parity','expanded_layout_parity']].sum().to_dict())


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--outdir',type=Path,required=True);execute(p.parse_args().outdir.resolve())
