from collections import Counter
import numpy as np
import pysam
from run import pairs, oracle
from fp_tools.tools.parametric_bias import BiasFeatureSpec, ConditionalSequenceBiasModel


def record(name,start,flag,mate,cigar='20M'):
    r=pysam.AlignedSegment();r.query_name=name;r.reference_id=0;r.reference_start=start
    r.flag=flag;r.mapping_quality=60;r.cigarstring=cigar;r.next_reference_id=0;r.next_reference_start=mate
    return r


class Handle:
    def __init__(self,reads):self.reads=reads
    def fetch(self,**kwargs):return iter(self.reads)


def test_complete_pair_inventory():
    a=record('a',100,99,200);b=record('a',200,147,100)
    orphan=record('b',100,99,200)
    bad=record('c',100,99,200);dup=record('c',200,147|1024,100)
    stats=Counter();got=list(pairs(Handle([a,b,orphan,bad,dup]),stats))
    assert len(got)==1
    assert stats['complete_fragments']==1
    assert stats['singleton_or_ambiguous_records']==2
    assert stats['excluded_records']==1
    assert stats['mapped_records']==5


def test_oracle_known_coordinates():
    assert oracle(record('x',100,0,200,'5S20M'))==99
    assert oracle(record('x',100,16,200,'20M7S'))==122
    assert oracle(record('x',100,16,200,'10M1000N10M'))==1115


def test_small_model_fit_is_finite():
    rng=np.random.default_rng(2026)
    x=rng.integers(0,4,(4,40,21),dtype=np.uint8)
    y=rng.poisson(.2,(4,40)).astype(float)
    model=ConditionalSequenceBiasModel(BiasFeatureSpec.loglinear21())
    model.fit(x,y,epochs=2,batch_windows=2,l2=.01,seed=2026)
    assert np.isfinite(model.conditional_nll(x,y))
