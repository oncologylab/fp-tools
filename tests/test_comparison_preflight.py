import copy

import pytest

from fp_tools.tools import bulk_footprinting as bulk
from fp_tools.tools.static_comparison_browser import build_static_browser
from fp_tools.utils.project_layout import read_comparison_table


def test_repeated_pairs_remain_available_without_combined_review(tmp_path, monkeypatch):
    samples = tmp_path / "samples.tsv"
    samples.write_text("sample\tcondition\ns1\tA\ns2\tB\n")
    comparisons = tmp_path / "comparisons.tsv"
    comparisons.write_text("comparison\tcond1\tcond2\nfirst\tA\tB\nsecond\tB\tA\n")
    assert len(read_comparison_table(comparisons)) == 2
    args = bulk.build_parser().parse_args([
        "--sample-table", str(samples), "--comparison-table", str(comparisons),
        "--genome", "hg38", "--outdir", str(tmp_path / "project"),
        "--review-format", "none", "--dry-run",
    ])
    assert bulk.run_bulk_footprinting(args) == 0
    assert not (tmp_path / "project").exists()


def test_unique_pairs_and_self_comparison(tmp_path):
    table = tmp_path / "comparisons.tsv"
    table.write_text("comparison\tcond1\tcond2\none\tA\tB\ntwo\tA\tC\n")
    assert len(read_comparison_table(table, unique_pairs=True)) == 2
    table.write_text("comparison\tcond1\tcond2\nself\tA\tA\n")
    with pytest.raises(ValueError, match=r"self.*line 2.*itself"):
        read_comparison_table(table, unique_pairs=True)


@pytest.mark.parametrize("dry_run", [False, True])
@pytest.mark.parametrize("review", ["auto", "bundle", "standalone"])
@pytest.mark.parametrize("pair", [("A", "B"), ("B", "A")])
def test_bulk_rejects_duplicate_pairs_before_reference(tmp_path, monkeypatch, dry_run, review, pair):
    samples = tmp_path / "samples.tsv"
    samples.write_text("sample\tcondition\ns1\tA\ns2\tB\n")
    comparisons = tmp_path / "comparisons.tsv"
    comparisons.write_text(
        "comparison\tcond1\tcond2\nfirst\tA\tB\n\nsecond\t%s\t%s\n" % pair
    )
    args = bulk.build_parser().parse_args([
        "--sample-table", str(samples), "--comparison-table", str(comparisons),
        "--genome", "hg38", "--outdir", str(tmp_path / "project"),
        "--review-format", review,
    ] + (["--dry-run"] if dry_run else []))
    monkeypatch.setattr(bulk, "resolve_analysis_reference", lambda *a, **k: pytest.fail("reference provisioning reached"))
    with pytest.raises(ValueError, match=r"first.*line 2.*second.*line 4"):
        bulk.run_bulk_footprinting(args)
    assert not (tmp_path / "project").exists()


@pytest.mark.parametrize("existing", [False, True])
@pytest.mark.parametrize("problem", ["duplicate", "reverse", "self", "default", "motif", "plots", "collision"])
def test_invalid_browser_preserves_output(tmp_path, existing, problem):
    payload = {"conditions": ["A", "B"], "points": [{"prefix": "TF1"}]}
    second = copy.deepcopy(payload)
    options = {}
    if problem == "reverse":
        second["conditions"] = ["B", "A"]
    elif problem == "self":
        second["conditions"] = ["C", "C"]
    elif problem == "collision":
        payload["conditions"] = ["A B", "C"]
        second["conditions"] = ["A_B", "C"]
    elif problem in {"default", "motif", "plots"}:
        second["conditions"] = ["A", "C"]
        options = {"default": {"default_comparison": ["X", "Y"]},
                   "motif": {"default_motifs": ["unknown"]},
                   "plots": {"default_aggregate_plots": 13}}[problem]
    output = tmp_path / "review"
    if existing:
        (output / "data").mkdir(parents=True)
        (output / "data" / "sentinel").write_bytes(b"original data")
        (output / "index.html").write_bytes(b"original index")
    before = {p.relative_to(output): p.read_bytes() for p in output.rglob("*") if p.is_file()}
    with pytest.raises(ValueError):
        build_static_browser([payload, second], output, title="Test", **options)
    assert output.exists() == existing
    assert {p.relative_to(output): p.read_bytes() for p in output.rglob("*") if p.is_file()} == before
