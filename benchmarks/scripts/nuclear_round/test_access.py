from pathlib import Path
import pytest
from access import deny


@pytest.mark.parametrize('identity', ['SRR13450941','SRR13450942','ENCSR587TRP','ENCSR485TLP','NakedDNA_rep2','GM23338'])
def test_identity_closed(identity):
    with pytest.raises(PermissionError):
        deny('/tmp/innocent.npz',identity)


def test_symlink_closed(tmp_path):
    target=tmp_path/'NakedDNA_rep3';target.mkdir()
    link=tmp_path/'permitted';link.symlink_to(target,target_is_directory=True)
    with pytest.raises(PermissionError):deny(link/'file.npz','rep1')


def test_permitted():
    deny('/tmp/NakedDNA_rep1/file.npz','GSE164997 SRR13450940')
