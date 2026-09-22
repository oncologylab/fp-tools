"""Synthetic pair and coordinate oracles, independent of scientific data."""
import sys
from pathlib import Path
import tempfile
import unittest
import numpy as np
import pandas as pd
import pysam
from fragments import complete_pair,oracle_cut

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'src'))
sys.path.insert(0,str(ROOT/'benchmarks/scripts'))
from fp_tools.tools.parametric_bias import cut_position_from_alignment
from build_strand_functional_profiles import extract_strand_cut_profiles


def read(name,start,flag,cigar='20M',mate=100):
    r=pysam.AlignedSegment();r.query_name=name;r.flag=flag;r.reference_id=0
    r.reference_start=start;r.mapping_quality=60;r.cigarstring=cigar
    r.query_sequence='A'*r.infer_query_length();r.next_reference_id=0;r.next_reference_start=mate
    return r


class FragmentTests(unittest.TestCase):
    def test_coordinates(self):
        for cigar in ['20M','5S20M','20M5S','3H5S20M7S2H','10M1000N10M','10M2I10M','10M3D10M']:
            for flag in [0,16]:
                r=read('x',30,flag,cigar)
                self.assertEqual(oracle_cut(r),cut_position_from_alignment(r,(4,-5)))

    def test_both_template_orientations(self):
        for flags in [(99,147),(163,83)]:
            a,b=read('p',10,flags[0],mate=100),read('p',100,flags[1],mate=10)
            self.assertEqual(complete_pair([a,b])[1],'complete')
            self.assertEqual(complete_pair([a])[1],'singleton')
            self.assertEqual(complete_pair([a,b,b])[1],'ambiguous')
            b.flag|=1024
            self.assertEqual(complete_pair([a,b])[1],'singleton')

    def test_corrected_extractor_independent_assignment(self):
        reads=[read('soft',130,0,'40S20M'),read('long',5,16,'800M'),
               read('other',170,16),read('duplicate',90,1024),read('plain',90,0)]
        sites=pd.DataFrame(dict(TFBS_chr=['chr1']*5,TFBS_start=[90,95,160,790,0],
                              TFBS_end=[110,105,180,810,2],TFBS_strand=['+','-','+','-','+']),index=[4,9,20,33,70])
        with tempfile.TemporaryDirectory() as directory:
            bam=Path(directory)/'fixture.bam'
            with pysam.AlignmentFile(bam,'wb',header={'HD':{'SO':'coordinate'},'SQ':[{'SN':'chr1','LN':1000}]}) as out:
                for r in sorted(reads,key=lambda r:r.reference_start):out.write(r)
            pysam.index(str(bam))
            for duplicates in [False,True]:
                plus,minus,valid=extract_strand_cut_profiles(sites,bam,flank=10,read_shift=(4,-5),keep_duplicates=duplicates)
                expected=np.zeros((5,2,21),np.uint32)
                for i,row in enumerate(sites.itertuples()):
                    center=(row.TFBS_start+row.TFBS_end)//2
                    if center<10 or center+11>1000:continue
                    for r in reads:
                        if r.is_duplicate and not duplicates:continue
                        offset=oracle_cut(r)-center+10
                        if 0<=offset<21:expected[i,int(r.is_reverse),offset]+=1
                np.testing.assert_array_equal(plus,expected[:,0]);np.testing.assert_array_equal(minus,expected[:,1])
                self.assertFalse(valid[-1])


if __name__=='__main__':unittest.main(verbosity=2)
