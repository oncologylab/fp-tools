"""Contract tests that keep documentation in sync with packaged entry points.

These guard against README/MkDocs/PyPI drift where the documented command
surface diverges from the console scripts declared in ``pyproject.toml``.
"""

import csv
import pathlib
import re
import tomllib
import unittest


ROOT = pathlib.Path(__file__).resolve().parents[1]
ENCODE_DOCS = ROOT / "docs" / "demos" / "data" / "encode"

REMOVED_ALIASES = {
    "ATACorrect",
    "FootprintScores",
    "ScoreBigwig",
    "BINDetect",
    "PlotAggregate",
}


def _load_pyproject():
    return tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))


def _verify_help_commands(readme: str):
    """Return the set of commands that have a ``<cmd> --help`` line in the README."""
    return set(re.findall(r"^([\w.-]+) --help$", readme, flags=re.MULTILINE))


class DocsEntryPointContractTest(unittest.TestCase):
    def setUp(self):
        self.data = _load_pyproject()
        self.project_scripts = self.data["project"]["scripts"]
        self.public_scripts = self.data["tool"]["fp-tools"]["public-console-scripts"]
        self.deprecated_scripts = self.data["tool"]["fp-tools"]["deprecated-console-scripts"]
        self.readme = (ROOT / "README.md").read_text(encoding="utf-8")
        self.markdown_paths = sorted((ROOT / "docs").rglob("*.md"))
        self.site_docs = "\n".join(
            path.read_text(encoding="utf-8") for path in self.markdown_paths
        )
        self.api_reference = (ROOT / "docs" / "api.md").read_text(encoding="utf-8")

    def test_setuptools_and_poetry_scripts_match(self):
        poetry_scripts = self.data["tool"]["poetry"]["scripts"]
        self.assertEqual(
            self.project_scripts,
            poetry_scripts,
            "[project.scripts] and [tool.poetry.scripts] have drifted apart.",
        )

    def test_primary_entry_points_are_documented_in_readme_and_manual(self):
        for command in self.public_scripts:
            self.assertIn(command, self.readme, f"{command} is missing from README.md")
            self.assertIn(
                command, self.site_docs, f"{command} is missing from MkDocs pages"
            )

    def test_help_block_exactly_covers_non_alias_commands(self):
        documented = _verify_help_commands(self.readme)
        self.assertTrue(documented <= set(self.public_scripts))

    def test_tobias_compatible_aliases_are_removed(self):
        for alias in REMOVED_ALIASES:
            self.assertNotIn(alias, self.project_scripts)
        documented = _verify_help_commands(self.readme)
        self.assertFalse(documented & REMOVED_ALIASES)

    def test_gui_is_included_in_standard_install(self):
        extras = self.data["project"].get("optional-dependencies", {})
        self.assertIn(
            "gui", extras, "Expected a [project.optional-dependencies] gui extra."
        )
        self.assertTrue(any("streamlit" in dep for dep in self.data["project"]["dependencies"]))
        self.assertIn("python -m pip install fp-tools-bio", self.readme)
        self.assertNotIn("fp-tools-bio[gui]", self.readme)

    def test_api_reference_is_command_manual(self):
        for command in self.public_scripts:
            self.assertIn(f"## `{command}`", self.api_reference)
            self.assertIn(f"usage: {command}", self.api_reference)
            section = self.api_reference.split(f"## `{command}`", 1)[1]
            section = section.split("\n## `", 1)[0]
            for required in (
                "**Example command**",
                "**Primary inputs**",
                "**Main outputs**",
                "**Complete options**",
            ):
                self.assertIn(required, section)
        for command in self.deprecated_scripts:
            self.assertNotIn(f"## `{command}`", self.api_reference)
        self.assertNotIn("::: fp_tools", self.api_reference)

    def test_each_public_command_has_a_concise_get_started_page(self):
        overview = (ROOT / "docs" / "get-started" / "tool-overview.md").read_text(
            encoding="utf-8"
        )
        command_dir = ROOT / "docs" / "get-started" / "commands"
        self.assertEqual(
            {path.stem for path in command_dir.glob("*.md")},
            set(self.public_scripts),
        )
        for command in self.public_scripts:
            page = command_dir / f"{command}.md"
            content = page.read_text(encoding="utf-8")
            self.assertIn(
                f"# [`{command}`](../../api.md#{command})",
                content,
                f"{command} page title should link to its complete API reference",
            )
            self.assertIn("## Example command", content)
            self.assertIn("## Primary inputs", content)
            self.assertIn("## Main outputs", content)
            self.assertIn(f"../../api.md#{command}", content)
            self.assertIn(f"commands/{command}.md", overview)

            example = content.split("## Example command", 1)[1].split(
                "## Primary inputs", 1
            )[0]
            primary_inputs = content.split("## Primary inputs", 1)[1].split(
                "## Main outputs", 1
            )[0]
            example_flags = set(re.findall(r"--[a-z0-9-]+", example))
            documented_flags = set(re.findall(r"--[a-z0-9-]+", primary_inputs))
            self.assertEqual(
                documented_flags,
                example_flags,
                f"{command} should explain exactly the flags shown in its example",
            )

    def test_get_started_navigation_preserves_four_top_level_tabs(self):
        config = (ROOT / "mkdocs.yml").read_text(encoding="utf-8")
        nav = config.split("nav:\n", 1)[1].split("\nmarkdown_extensions:", 1)[0]
        top_level = re.findall(r"^  - ([^:]+):", nav, flags=re.MULTILINE)
        self.assertEqual(
            top_level,
            ["Get Started", "Output Demo", "GUI Demo", "API Reference"],
        )
        for required in [
            "Home: index.md",
            "Installation: get-started/installation.md",
            "Tool overview: get-started/tool-overview.md",
            "Bulk ATAC-seq: get-started/output-examples/bulk-atac-seq.md",
            "Single-cell ATAC-seq: get-started/output-examples/single-cell-atac-seq.md",
        ]:
            self.assertIn(required, nav)
        self.assertIn("font: false", config)
        self.assertIn("toc_depth: 3", config)
        self.assertIn("javascripts/layout.js", config)
        styles = (ROOT / "docs" / "stylesheets" / "extra.css").read_text(
            encoding="utf-8"
        )
        self.assertIn('"Helvetica Neue", Helvetica, Arial, sans-serif', styles)
        self.assertRegex(styles, r"\.md-typeset h1 \{[^}]*font-size: 1\.25rem;")
        self.assertRegex(styles, r"\.md-typeset h2 \{[^}]*font-size: 1rem;")
        self.assertRegex(styles, r"\.md-typeset h3 \{[^}]*font-size: 0\.9rem;")

        reports = (ROOT / "docs" / "reports.md").read_text(encoding="utf-8")
        gui = (ROOT / "docs" / "gui.md").read_text(encoding="utf-8")
        api = (ROOT / "docs" / "api.md").read_text(encoding="utf-8")
        self.assertIn("  - navigation", reports)
        self.assertIn("  - navigation", gui)
        self.assertIn("  - navigation", api)

    def test_core_command_footer_is_a_closed_sequence(self):
        sequence = [
            "prepare-atac",
            "atac-correct",
            "call-footprints",
            "match-motifs",
            "diff-footprints",
            "normalize-bigwig",
        ]
        command_dir = ROOT / "docs" / "get-started" / "commands"
        for index, command in enumerate(sequence):
            content = (command_dir / f"{command}.md").read_text(encoding="utf-8")
            self.assertIn("core_nav:", content)
            if index:
                self.assertIn(f"title: {sequence[index - 1]}", content)
            if index < len(sequence) - 1:
                self.assertIn(f"title: {sequence[index + 1]}", content)
        for command in set(self.public_scripts) - set(sequence):
            content = (command_dir / f"{command}.md").read_text(encoding="utf-8")
            self.assertNotIn("core_nav:", content)
        footer = (ROOT / "overrides" / "partials" / "footer.html").read_text(
            encoding="utf-8"
        )
        self.assertIn("page.meta.core_nav", footer)
        self.assertNotIn("page.previous_page", footer)
        self.assertNotIn("page.next_page", footer)

    def test_encode_example_tables_are_complete_and_portable(self):
        from fp_tools.tools.prepare_atac import read_preprocess_metadata
        from fp_tools.utils.project_layout import read_comparison_table, read_sample_table

        def rows(name):
            with (ENCODE_DOCS / name).open(encoding="utf-8", newline="") as handle:
                return list(csv.DictReader(handle, delimiter="\t"))

        fastq_rows = rows("encode_hepg2_k562_fastq_urls.tsv")
        local_fastq_rows = rows("local_fastq_template.tsv")
        self.assertEqual(len(fastq_rows), 6)
        self.assertEqual({row["condition"] for row in fastq_rows}, {"HepG2", "K562"})
        for row in fastq_rows:
            for field in ("fastq_1", "fastq_2"):
                self.assertTrue(row[field].startswith("https://www.encodeproject.org/files/"))
            for field in ("fastq_1_md5", "fastq_2_md5"):
                self.assertRegex(row[field], r"^[0-9a-f]{32}$")
        self.assertEqual(
            list(local_fastq_rows[0]),
            ["sample", "condition", "fastq_1", "fastq_2"],
        )

        small_samples = rows("encode_hepg2_k562_bams.tsv")
        local_samples = rows("local_bam_peak_template.tsv")
        full_samples = rows("encode_cancer_7line_bams.tsv")
        small_comparisons = rows("encode_hepg2_k562_comparisons.tsv")
        full_comparisons = rows("encode_cancer_7line_comparisons.tsv")
        self.assertEqual(len(small_samples), 6)
        self.assertEqual(len(local_samples), 4)
        self.assertEqual(len(full_samples), 17)
        self.assertEqual(len(small_comparisons), 1)
        self.assertEqual(len(full_comparisons), 21)
        self.assertEqual(len({row["sample"] for row in full_samples}), 17)
        self.assertEqual(
            {row["condition"] for row in full_samples},
            {"A549", "HCT116", "HepG2", "K562", "MCF-7", "PC-3", "Panc1"},
        )
        conditions = {row["condition"] for row in full_samples}
        for row in full_samples:
            self.assertFalse(pathlib.Path(row["bam"]).is_absolute())
            self.assertFalse(pathlib.Path(row["peaks"]).is_absolute())
            self.assertTrue(row["peaks"].endswith(".bed"))
            self.assertFalse(row["peaks"].endswith(".bed.gz"))
            self.assertRegex(row["bam_md5"], r"^[0-9a-f]{32}$")
            self.assertRegex(row["peak_md5"], r"^[0-9a-f]{32}$")
        for row in full_comparisons:
            self.assertIn(row["cond1"], conditions)
            self.assertIn(row["cond2"], conditions)
            self.assertNotEqual(row["cond1"], row["cond2"])

        self.assertEqual(
            len(read_preprocess_metadata(ENCODE_DOCS / "encode_hepg2_k562_fastq_urls.tsv")),
            6,
        )
        self.assertEqual(
            len(read_preprocess_metadata(ENCODE_DOCS / "local_fastq_template.tsv")),
            2,
        )
        self.assertEqual(
            len(read_sample_table(ENCODE_DOCS / "encode_hepg2_k562_bams.tsv")), 6
        )
        self.assertEqual(
            len(read_sample_table(ENCODE_DOCS / "local_bam_peak_template.tsv")), 4
        )
        self.assertEqual(
            len(read_sample_table(ENCODE_DOCS / "encode_cancer_7line_bams.tsv")), 17
        )
        self.assertEqual(
            len(read_comparison_table(ENCODE_DOCS / "encode_cancer_7line_comparisons.tsv")),
            21,
        )

    def test_encode_docs_have_no_machine_specific_paths(self):
        paths = [
            *ENCODE_DOCS.glob("*"),
            *(ROOT / "docs" / "demos" / "qc" / "encode").glob("*.tsv"),
            ROOT / "docs" / "get-started" / "workflows" / "bulk-atac-seq.md",
            ROOT / "docs" / "api.md",
        ]
        text = "\n".join(
            path.read_text(encoding="utf-8") for path in paths if path.is_file()
        )
        for forbidden in ("/home/", "169.254.169.254", "localhost:"):
            self.assertNotIn(forbidden, text)

    def test_bulk_workflow_uses_minimal_beginner_examples(self):
        guide = (
            ROOT / "docs" / "get-started" / "workflows" / "bulk-atac-seq.md"
        ).read_text(encoding="utf-8")
        self.assertIn("sample\tcondition\tfastq_1\tfastq_2", guide)
        self.assertIn("sample\tcondition\tbam\tpeaks", guide)
        self.assertIn("Only the four columns shown above are needed", guide)
        for unnecessary_detail in (
            "Plan storage and runtime before downloading",
            "Why the wrapper output is not byte-for-byte identical",
            "Representative ENCODE QC files",
            "not a reproduction of the pair-specific demo",
        ):
            self.assertNotIn(unnecessary_detail, guide)

    def test_command_guides_use_precise_signal_names_and_file_patterns(self):
        command_dir = ROOT / "docs" / "get-started" / "commands"
        required_patterns = {
            "prepare-atac": "{sample}.rp10m.bw",
            "atac-correct": "{sample}_corrected.bw",
            "call-footprints": "{sample}_footprints.bw",
            "match-motifs": "motif_matches_results.txt",
            "diff-footprints": "{prefix}_results.txt",
            "normalize-bigwig": "{sample}_corrected_q95_scaled.bw",
        }
        for command, pattern in required_patterns.items():
            content = (command_dir / f"{command}.md").read_text(encoding="utf-8")
            self.assertIn(pattern, content)
        prepare = (command_dir / "prepare-atac.md").read_text(encoding="utf-8")
        self.assertNotIn("RP10M coverage", prepare)
        self.assertIn("alignment coverage bigWig", prepare)
        atac = (command_dir / "atac-correct.md").read_text(encoding="utf-8")
        for signal in (
            "Observed base-resolution cut-site signal",
            "Tn5 sequence-bias score",
            "Expected cut-site signal",
            "Bias-corrected cut-site signal",
        ):
            self.assertIn(signal, atac)

    def test_public_site_docs_use_current_wording(self):
        public_docs = "\n".join(
            (ROOT / path).read_text(encoding="utf-8")
            for path in [
                "README.md",
                "docs/index.md",
                "docs/reports.md",
                "docs/api.md",
            ]
        )
        for stale in [
            "BINDetect",
            "BINDetect-style",
            "Compatibility Aliases",
            "Compatibility Aliases",
            "Open static report",
            "GUI Preview",
            "The PyPI package is",
            "The Python import name is",
            "ATACCorrect",
        ]:
            self.assertNotIn(stale, public_docs)

    def test_api_reference_uses_general_purpose_names(self):
        api = self.api_reference.lower()
        for project_specific_term in ("nutrient", "legacy"):
            self.assertNotIn(project_specific_term, api)
        self.assertIn("--profile {modern,homer-atac}", api)
        self.assertIn("--footprint-kernel {fast,reference}", api)

    def test_public_site_uses_branded_assets_and_responsive_demo_contracts(self):
        index = (ROOT / "docs" / "index.md").read_text(encoding="utf-8")
        reports = (ROOT / "docs" / "reports.md").read_text(encoding="utf-8")
        gui = (ROOT / "docs" / "gui.md").read_text(encoding="utf-8")
        bulk_example = (
            ROOT
            / "docs"
            / "get-started"
            / "output-examples"
            / "bulk-atac-seq.md"
        ).read_text(encoding="utf-8")
        single_cell_example = (
            ROOT
            / "docs"
            / "get-started"
            / "output-examples"
            / "single-cell-atac-seq.md"
        ).read_text(encoding="utf-8")
        gui_demo = (
            ROOT / "docs" / "demos" / "gui" / "fp-tools-gui-static-demo.html"
        ).read_text(encoding="utf-8")
        report_demo = (
            ROOT
            / "docs"
            / "demos"
            / "reports"
            / "diff_footprints_K562_HepG2.html"
        ).read_text(encoding="utf-8")
        self.assertIn("fp_tools_logo_horizontal.svg", self.readme)
        self.assertIn("# ATAC-seq footprinting and regulatory motif analysis", index)
        self.assertIn("fp_tools_logo_horizontal.svg", index)
        self.assertLess(
            index.index("fp_tools_logo_horizontal.svg"),
            index.index("# ATAC-seq footprinting and regulatory motif analysis"),
        )
        self.assertNotIn("fp-hero", index)
        self.assertNotIn("Where To Go Next", index)
        self.assertIn(
            'src="../ENCODE-Cancer-Cell-lines-Footprinting/"', reports
        )
        self.assertIn('class="fp-live-demo', reports)
        self.assertNotIn("interface_diff_footprints_html.png", reports)
        self.assertIn(
            'src="../demos/gui/fp-tools-gui-static-demo.html"', gui
        )
        self.assertIn('class="fp-live-demo', gui)
        self.assertNotIn("interface_gui_home.png", gui)
        self.assertNotIn("fp-demo-callout", reports + gui)
        self.assertIn(
            'src="../../../ENCODE-Cancer-Cell-lines-Footprinting/"', bulk_example
        )
        self.assertIn("pbmc5k_signature_heatmap.png", single_cell_example)
        self.assertIn("pbmc5k_eight_marker_umaps.png", single_cell_example)
        for asset in [
            "pbmc5k_signature_heatmap.png",
            "pbmc5k_eight_marker_umaps.png",
        ]:
            self.assertTrue((ROOT / "docs" / "assets" / asset).exists())
        self.assertFalse(
            (
                ROOT
                / "docs"
                / "assets"
                / "pbmc5k_single_cell_footprinting_summary.svg"
            ).exists()
        )
        for marker in [
            "STAT6",
            "FOSB",
            "CEBPA",
            "IRF8",
            "RELA",
            "ZNF683",
            "NR4A1",
            "SMAD3",
        ]:
            self.assertIn(marker, single_cell_example)
        self.assertFalse(
            (ROOT / "docs" / "assets" / "interface_diff_footprints_html.png").exists()
        )
        self.assertFalse(
            (ROOT / "docs" / "assets" / "interface_gui_home.png").exists()
        )
        self.assertIn('class="menu-button"', gui_demo)
        self.assertIn('class="documentation-link" href="../../" target="_top"', gui_demo)
        self.assertIn("aria-current", gui_demo)
        self.assertIn('rel="icon"', gui_demo)
        self.assertIn('rel="icon"', report_demo)

        encode_browser = (
            ROOT / "docs" / "ENCODE-Cancer-Cell-lines-Footprinting" / "index.html"
        ).read_text(encoding="utf-8")
        self.assertIn('class="documentation-link"', encode_browser)
        self.assertIn('href="../"', encode_browser)
        self.assertIn('target="_top"', encode_browser)

    def test_landing_pages_are_concise_and_use_current_workflow_names(self):
        index = (ROOT / "docs" / "index.md").read_text(encoding="utf-8")
        self.assertNotIn("Reproducible regulatory genomics", index)
        self.assertNotIn("Where To Go Next", index)
        self.assertNotIn("fp-hero", index)
        self.assertLess(len(index.split()), 190)
        self.assertIn("ATAC-seq footprinting and regulatory motif analysis", index)
        self.assertIn("Output demo with ENCODE cancer cell lines", self.readme)
        self.assertIn(
            "https://oncologylab.github.io/fp-tools/demos/gui/"
            "fp-tools-gui-static-demo.html",
            self.readme,
        )
        self.assertIn("Single-cell ATAC-seq", self.readme)
        self.assertIn("De novo motifs", self.readme)

    def test_prepare_atac_profiles_are_explained_in_plain_language(self):
        public_atac_docs = "\n".join(
            (ROOT / path).read_text(encoding="utf-8")
            for path in ["README.md", "docs/index.md", "docs/api.md"]
        )
        for required in [
            "fastp",
            "MACS3",
            "Trim Galore",
            "Picard",
            "HOMER",
            "--profile {modern,homer-atac}",
        ]:
            self.assertIn(required, public_atac_docs)
        for rejected in [
            "CUT&Tag branches",
            "`XS:i:` filtering",
            "legacy-compatible",
            "Bowtie2/Picard",
        ]:
            self.assertNotIn(rejected, public_atac_docs)

    def test_local_markdown_links_exist(self):
        markdown_files = [
            ROOT / "README.md",
            ROOT / "AGENTS.md",
            ROOT / "DEV_PLAN.md",
            ROOT / "RELEASE_CHECKLIST.md",
            *sorted((ROOT / "docs").rglob("*.md")),
            *sorted((ROOT / "benchmarks").glob("*.md")),
        ]
        pattern = re.compile(r"\[[^\]]+\]\((?!https?://|mailto:|#)([^)#]+)")
        missing = []
        for source in markdown_files:
            for target in pattern.findall(source.read_text(encoding="utf-8")):
                path = (source.parent / target).resolve()
                if not path.exists():
                    missing.append(f"{source.relative_to(ROOT)} -> {target}")
        self.assertEqual(missing, [])


if __name__ == "__main__":
    unittest.main()
