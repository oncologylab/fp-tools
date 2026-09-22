from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'benchmarks/scripts'))
from compare_blacklisted_control_layouts import expanded_blacklist_command


def test_expanded_regions_are_preserved_when_blacklist_added():
    command = ['python', '--peaks', 'historical.bed', '--regions-out', 'all.bed',
               '--outdir', 'old', '--prefix', 'old', '--norm-off']
    original = command.copy()
    result = expanded_blacklist_command(command, Path('/tmp/new'), Path('/tmp/blacklist.bed'))
    assert command == original
    assert result[:5] == command[:5]
    assert result[6] == '/tmp/new/tracks'
    assert result[-2:] == ['--blacklist', '/tmp/blacklist.bed']
    with pytest.raises(ValueError):
        expanded_blacklist_command(result, Path('/tmp/next'), Path('/tmp/other.bed'))
