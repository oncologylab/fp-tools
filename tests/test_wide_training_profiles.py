from pathlib import Path
import sys
from types import SimpleNamespace
import numpy as np
import pandas as pd
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'benchmarks/scripts'))
import build_wide_training_profiles as wide


def test_background_diagnostic_preserves_zeros_and_exposure_adjustment():
    counts=np.zeros((3,401));bias=np.zeros_like(counts)
    counts[0,100:151]=2;counts[0,250:301]=4
    bias[1,100]=np.nan
    result=wide.background_diagnostic(counts,bias,200,100)
    np.testing.assert_array_equal(result['available'],[True,False,True])
    assert result['left_count'][0]==102 and result['right_count'][0]==204
    assert result['disagreement'][0]==pytest.approx(abs(np.log2(102.5/204.5)))
    assert result['disagreement'][2]==0 and np.isnan(result['disagreement'][1])
    shifted=wide.background_diagnostic(counts,bias+100,200,100)
    np.testing.assert_allclose(result['disagreement'],shifted['disagreement'],equal_nan=True)
    with pytest.raises(ValueError): wide.background_diagnostic(counts,bias,200,500)


def test_wider_context_preserves_genomic_orientation_and_shared_sequence(monkeypatch):
    sequence='ACGT'*1500
    class Fasta:
        references=['chr1'];lengths=[len(sequence)]
        def __enter__(self): return self
        def __exit__(self,*args): pass
        def fetch(self,chrom,start,end): return sequence[start:end]
    monkeypatch.setattr(wide.old.strands,'open_fasta',lambda _:Fasta())
    def bias(model,seq,positions):
        values=np.array(['ACGT'.index(seq[p]) for p in positions],float)
        return values,-values,np.ones(len(values),bool)
    monkeypatch.setattr(wide.old.strands,'strand_log_bias',bias)
    sites=pd.DataFrame(dict(TFBS_chr=['chr1','chr1','absent','chr1'],TFBS_start=[2000,2000,2000,0],
                            TFBS_end=[2010,2010,2010,10],TFBS_strand=['+','-','+','+']))
    model=SimpleNamespace(feature_spec=SimpleNamespace(context_length=21))
    small=wide.infer(sites,model,'unused',200);large=wide.infer(sites,model,'unused',600)
    for before,after in zip(small,large):
        assert after.shape==(4,1201)
        np.testing.assert_array_equal(before[:2],after[:2,400:801])
        np.testing.assert_array_equal(after[0],after[1])
        assert np.isnan(after[2:]).all()
    with pytest.raises(ValueError,match='flank'):
        wide.infer(sites,model,'unused',199)
