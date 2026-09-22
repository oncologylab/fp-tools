from pathlib import Path
import sys

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'benchmarks/scripts'))
from audit_control_blacklist import blacklist_command, overlap_flags


def test_blacklist_intervention_preserves_source_and_other_options():
    command = ['python', '--peaks', 'old.bed', '--regions-out', 'old.bed',
               '--outdir', 'old', '--prefix', 'old', '--norm-off']
    original = command.copy()
    changed = blacklist_command(command, Path('/tmp/new'), Path('/tmp/blacklist.bed'))
    assert command == original
    assert changed[-2:] == ['--blacklist', '/tmp/blacklist.bed']
    assert changed[:5] == command[:5]
    assert '--norm-off' in changed
    with pytest.raises(ValueError):
        blacklist_command(changed, Path('/tmp/next'), Path('/tmp/other.bed'))


def test_blacklist_half_open_overlap_and_window_flags_keep_denominator():
    sites = pd.DataFrame(dict(TFBS_chr=['chr1'] * 4, TFBS_start=[10, 20, 30, 1000], TFBS_end=[20, 30, 40, 1010]))
    intervals = [('chr1', 20, 30)]
    assert overlap_flags(sites, intervals).tolist() == [False, True, False, False]
    assert overlap_flags(sites, intervals, 100).tolist() == [True, True, True, False]
