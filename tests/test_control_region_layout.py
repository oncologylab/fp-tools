from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'benchmarks/scripts'))
from isolate_control_region_layout import historical_command


def test_historical_pair_changes_only_output_scope_and_destinations():
    original=['python','-m','fp_tools.cli','--peaks','old.bed','--regions-out','all.bed',
              '--bias-pkl','fixed.pickle','--norm-off','--outdir','old_output','--prefix','old',
              '--read_shift','4','-5','--cores','4']
    before=original.copy();changed=historical_command(original,Path('/tmp/isolated'))
    assert original==before
    assert [i for i,(a,b) in enumerate(zip(original,changed)) if a!=b]==[6,11,13]
    assert changed[6]=='old.bed' and changed[11]=='/tmp/isolated/tracks' and changed[13]=='NakedDNA_rep2_historical'
