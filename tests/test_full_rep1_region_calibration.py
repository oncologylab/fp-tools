from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'benchmarks/scripts'))
import calibrate_full_rep1_regions as study


def test_freeze_preserves_existing_output(tmp_path):
    marker = tmp_path / 'marker'
    marker.write_text('preserve')
    with pytest.raises(ValueError, match='fresh output'):
        study.freeze(tmp_path)
    assert marker.read_text() == 'preserve'


@pytest.mark.parametrize('changed_policy', [False, True])
def test_verify_never_starts_missing_extraction(tmp_path, monkeypatch, changed_policy):
    config = dict(policy=dict(study.POLICY), protocol=study.calibration.PROTOCOL,
                  runtime={}, inputs=[], code=[])
    if changed_policy:
        config['policy']['model_refitted'] = True
    extension = study.calibration.source.extraction.extension
    monkeypatch.setattr(extension, 'read_document', lambda *args: config)
    monkeypatch.setattr(extension, 'runtime', lambda: {})

    def forbidden(*args):
        pytest.fail('verification must not initiate extraction')

    monkeypatch.setattr(study.calibration.source, 'extract', forbidden)
    with pytest.raises(ValueError, match='frozen study changed' if changed_policy else 'missing profiles'):
        study.execute(tmp_path, verify=True)
