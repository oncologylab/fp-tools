"""Explicit independent derived-data replay; outputs go outside the archive."""
import argparse
import json
from pathlib import Path
from verify import verify


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('package',type=Path)
    parser.add_argument('--outdir',type=Path,required=True)
    args=parser.parse_args();root=args.package.resolve();out=args.outdir.resolve()
    if out.is_relative_to(root):raise ValueError('reanalysis must not modify evidence package')
    out.mkdir(parents=True,exist_ok=False)
    result=verify(root,return_tables=True)
    tables=result.pop('tables')
    for name,table in tables.items():table.to_csv(out/(name+'.tsv.gz' if name=='per_site_scores' else name+'.tsv'),sep='\t',index=False)
    result['recomputed_tables']=list(tables)
    (out/'reanalysis.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
