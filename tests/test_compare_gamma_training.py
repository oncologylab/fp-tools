from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'benchmarks/scripts'))
from compare_gamma_training import align_indices, complete_pairs, load_tracks, KEYS


def test_archive_members_are_loaded_once_before_per_site_scoring():
    class Archive:
        def __init__(self): self.reads={}
        def __getitem__(self,key):
            self.reads[key]=self.reads.get(key,0)+1
            return np.arange(12).reshape(3,4)
    archive=Archive()
    tracks=load_tracks(archive,np.array([2,0]))
    for values in tracks.values():
        np.testing.assert_array_equal(values,[[8,9,10,11],[0,1,2,3]])
    assert len(archive.reads)==4 and set(archive.reads.values())=={1}


def test_review_preserves_adjacent_score_values_and_ties(tmp_path):
    from review_gamma_training import read_table
    values=np.array([.1234567890123456,np.nextafter(.1234567890123456,1),.1234567890123456])
    path=tmp_path/'scores.tsv'
    pd.DataFrame(dict(score=values)).to_csv(path,sep='\t',index=False)
    np.testing.assert_array_equal(read_table(path).score.to_numpy(),values)


def test_alignment_requires_exact_unique_coordinates_and_strands():
    source=pd.DataFrame([['HepG2','CTCF','chr1',10,20,'+'],
                         ['HepG2','CTCF','chr1',30,40,'-']],columns=KEYS)
    np.testing.assert_array_equal(align_indices(source.iloc[::-1],source),[1,0])
    with pytest.raises(ValueError,match='duplicate'):
        align_indices(source,pd.concat([source,source]))
    changed=source.copy();changed.loc[0,'TFBS_strand']='-'
    with pytest.raises(ValueError,match='missing'):
        align_indices(changed,source)


def test_support_removes_both_members_without_crossing_task_boundaries():
    pairs=pd.DataFrame(dict(cell=['HepG2']*4,tf=['CTCF','CTCF','MYC','MYC'],
                            pair_id=[1]*4,chip_label=[1,0,1,0]))
    np.testing.assert_array_equal(complete_pairs(pairs,[True,False,True,True]),[False,False,True,True])
    with pytest.raises(ValueError,match='one positive'):
        complete_pairs(pairs.iloc[:3],[True]*3)
    with pytest.raises(ValueError,match='length'):
        complete_pairs(pairs,[True])
