import json
from pathlib import Path
import pytest
from run import check
from access import digest


def test_amendment_requires_original_identity(tmp_path):
    from protocol import POLICY
    root=Path(__file__).resolve().parents[3]
    relative='benchmarks/scripts/nuclear_round/run.py'
    original='a'*64;current=digest(root/relative)
    (tmp_path/'EXPERIMENT_PLAN.freeze.json').write_text(json.dumps(dict(policy=POLICY,code={relative:original})))
    with pytest.raises(ValueError,match='code changed'):check(tmp_path)
    correction=tmp_path/'IMPLEMENTATION_CORRECTION.freeze.json'
    correction.write_text(json.dumps(dict(code_replacement={relative:dict(original_sha256=original,corrected_sha256=current)})))
    check(tmp_path)
    correction.write_text(json.dumps(dict(code_replacement={relative:dict(original_sha256='b'*64,corrected_sha256=current)})))
    with pytest.raises(ValueError,match='amendment origin'):check(tmp_path)
