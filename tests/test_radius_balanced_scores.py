from pathlib import Path
import sys,json
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'benchmarks/scripts'))
import evaluate_radius_balanced_scores as runner


def test_secondary_evaluation_retains_empty_tasks_and_recomputes_on_same_sites(tmp_path,monkeypatch):
    sites=[];diag=[];scores=[]
    for i in range(21):
        n=1 if i<20 else 0
        diag.append(dict(cell='HepG2',tf=f'TF{i}',original_pairs=2,retained_pairs=n,status='underpowered'))
        for label in range(2*n):
            row=dict(cell='HepG2',tf=f'TF{i}',TFBS_chr='chr1',TFBS_start=i*100+label*10,
                     TFBS_end=i*100+label*10+5,TFBS_strand='+',chip_label=label)
            sites.append(row)
            for arm in runner.ARMS: scores.append(dict(row,arm=arm,score=float(label)))
    def save(name,rows):
        p=tmp_path/(name+'.tsv');pd.DataFrame(rows).to_csv(p,sep='\t',index=False);return runner.record(p)
    outputs=dict(selected_sites=save('sites',sites),diagnostics=save('diag',diag))
    config=dict(configuration_id='synthetic',runtime={},arms=runner.ARMS,inputs=[],code=[],
        cohort_outputs=outputs,scores=[save('scores',scores)])
    extension=runner.cohorts.parent.prior.wide.old.parent.extension
    monkeypatch.setattr(extension,'runtime',lambda:{})
    monkeypatch.setattr(extension,'read_document',lambda *args:config)
    out=tmp_path/'result';out.mkdir();runner.execute(out)
    table=pd.read_csv(out/'metrics.tsv',sep='\t')
    assert len(table)==147 and table.auroc.notna().sum()==140
    assert table.loc[table.retained_pairs.eq(1),'auroc'].eq(1).all()
    assert table.loc[table.retained_pairs.eq(0),'auroc'].isna().all()
    result=json.loads((out/'comparison_manifest.json').read_text())
    assert not result['shape_reassessed'] and not result['primary_result_replaced']
