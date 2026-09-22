"""Indexed fragment inputs must still work when pysam is unavailable."""

import gzip

import numpy as np
import pandas as pd
import pytest

from fp_tools.tools import find_signature_fp as signature


@pytest.mark.parametrize("create_index", [False, True])
@pytest.mark.parametrize("has_index", [False, True])
def test_fragment_scan_without_pysam_preserves_cut_counts(tmp_path, monkeypatch, create_index, has_index):
    fragments = tmp_path / "fragments.tsv.gz"
    with gzip.open(fragments, "wt") as handle:
        handle.write("# example\nchr22\t90\t110\tcell1\t2\n\nchr22\t500\t510\tcell1\t1\n")
    index = tmp_path / "fragments.tsv.gz.tbi"
    if has_index:
        index.write_bytes(b"index unavailable without pysam")
    monkeypatch.setattr(signature, "pysam", None)
    sites = {"TF": [("chr22", 100)]}
    assert not signature.ensure_tabix_index(fragments, create_index)
    profiles = signature.count_fragment_profiles(
        fragments, pd.DataFrame({"barcode": ["cell1"]}), sites, 20, create_index,
    )
    assert profiles.shape == (1, 1, 40)
    assert profiles.sum() == 4
    assert np.count_nonzero(profiles) == 2
    assert index.exists() == has_index
