"""Export permitted short sequence contexts for independent prediction checking."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import pysam
from access import Access,digest,dump


def export(out):
    access=Access(out)
    frame=pd.read_csv(access.path('controls'),sep='\t')
    if len(frame)!=4200 or not frame.TFBS_chr.isin([f'chr{i}' for i in range(1,16)]).all():
        raise ValueError('context export scope')
    path=out/'arrays/control_sequence_contexts.npz'
    if path.exists():raise FileExistsError(path)
    sequences=[]
    with pysam.FastaFile(str(access.path('genome'))) as fasta:
        for row in frame.itertuples():
            sequence=fasta.fetch(row.TFBS_chr,int(row.center)-141,int(row.center)+142).upper()
            if len(sequence)!=283 or set(sequence)-set('ACGT'):
                raise ValueError('unexpected sequence support; stop before new scoring')
            sequences.append(sequence)
    np.savez_compressed(path,site_id=frame.site_id.to_numpy(str),sequence=np.array(sequences),
        chromosome=frame.TFBS_chr.to_numpy(str),motif_strand=frame.TFBS_strand.to_numpy(str))
    dump(out/'context_export.json',dict(code_sha256=digest(__file__),output_sha256=digest(path),
        sites=4200,width=283,scope='short permitted derived sequence contexts, not genome payload',
        sequence_support='all ACGT in full283-base contexts'))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--outdir',type=Path,required=True)
    export(p.parse_args().outdir.resolve())
