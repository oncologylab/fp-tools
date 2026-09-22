from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'benchmarks/scripts'))
import screen_full_rep1_counts as study


def test_freeze_preserves_existing_directory(tmp_path):
    with pytest.raises(ValueError, match='fresh output'):
        study.freeze(tmp_path)
    assert not list(tmp_path.iterdir())


def test_verify_cannot_create_screen(tmp_path, monkeypatch):
    config = dict(policy=study.POLICY, runtime={}, inputs=[], code=[],
                  source_config={'configuration_id': 'old'}, profile_configuration_id='profiles')
    ext = study.screen.controls.extraction.extension
    monkeypatch.setattr(ext, 'read_document', lambda *args: config)
    monkeypatch.setattr(ext, 'runtime', lambda: {})
    monkeypatch.setattr(study.screen.controls.extraction.base, 'verify_stage', lambda *args: {'outputs': {}})
    monkeypatch.setattr(study, 'record', lambda p: {'path': str(p), 'sha256': 'test'})
    with pytest.raises(ValueError, match='missing screen'):
        study.execute(tmp_path, verify=True)
    assert not (tmp_path / 'screen').exists()
