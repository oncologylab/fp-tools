"""Focused contracts for counting, tied metrics, and round-trip strict calls."""
import tempfile
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score,average_precision_score
from run import components
from verify import metrics,read_table


class PipelineTests(unittest.TestCase):
    def test_independent_metrics_with_ties(self):
        rng=np.random.default_rng(2026)
        for size in [4,20,101]:
            labels=np.arange(size)%2;score=rng.integers(-3,4,size=size)
            a,b=metrics(labels,score)
            self.assertAlmostEqual(a,roc_auc_score(labels,score),places=14)
            self.assertAlmostEqual(b,average_precision_score(labels,score),places=14)

    def test_strict_float_roundtrip(self):
        threshold=.4738012838460274
        values=np.array([np.nextafter(threshold,-np.inf),threshold,np.nextafter(threshold,np.inf)])
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'values.tsv'
            pd.DataFrame(dict(score=values,threshold=threshold)).to_csv(path,sep='\t',index=False)
            table=read_table(path)
            np.testing.assert_array_equal(table.score.to_numpy().view(np.uint64),values.view(np.uint64))
            self.assertEqual((table.score>table.threshold).tolist(),[False,False,True])

    def test_regional_count_boundary_contract(self):
        p=np.zeros((2,201),int);m=p.copy()
        for offset in [-41,-40,-12,-11,-8,-7,0,7,8,11,12,40,41]:p[0,offset+100]=1
        m[0,100]=2
        arrays={'plus_observed':p,'minus_observed':m}
        props={arm+'_'+strand:np.zeros((2,201)) for arm in ['LOG21','DWM'] for strand in ['plus','minus']}
        frame=components(pd.DataFrame(dict(site_id=['a','b'])),arrays,props)
        self.assertEqual(frame.loc[0,'plus_center'],3)
        self.assertEqual(frame.loc[0,'plus_left'],2)
        self.assertEqual(frame.loc[0,'plus_right'],2)
        self.assertEqual(frame.loc[0,'minus_total'],2)
        self.assertEqual(frame.loc[1,'total'],0)
        self.assertAlmostEqual(frame.loc[0,'LOG21_plus_p'],15/73,places=14)


if __name__=='__main__':unittest.main(verbosity=2)
