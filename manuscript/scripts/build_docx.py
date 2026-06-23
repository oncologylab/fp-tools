#!/usr/bin/env python
"""Build a Word DOCX from the manuscript LaTeX source.

The PDF manuscript uses PDF/SVG-oriented figure assets and MDPI-specific LaTeX
metadata commands. For Word, render PDF figures to PNG first and feed Pandoc a
small, standard-LaTeX wrapper so figures are embedded reliably.
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
MANUSCRIPT_DIR = REPO_ROOT / "manuscript"
DEFAULT_TEX = MANUSCRIPT_DIR / "main.tex"
DEFAULT_OUTPUT = MANUSCRIPT_DIR / "main.docx"
BUILD_DIR = MANUSCRIPT_DIR / ".docx_build"
MEDIA_DIR = BUILD_DIR / "figures"
VANCOUVER_CSL = Path("/usr/share/citation-style-language/styles/vancouver.csl")


FIGURE_WIDTHS = {
    "Fig1": "3.0in",
    "Fig2": "6.0in",
    "Fig3": "4.7in",
    "Fig4": "6.0in",
    "Fig5": "3.2in",
    "Fig6": "6.0in",
}

REF_REPLACEMENTS = {
    "fig:workflow_overview": "1",
    "fig:denovo_validation": "2",
    "fig:differential_report": "3",
    "fig:single_cell_footprinting": "4",
    "fig:engineering": "5",
    "fig:gui_interface": "6",
    "tab:features": "1",
    "tab:supp_analysis_parameters": "S1",
    "tab:supp_software_versions": "S2",
}

BACK_MATTER = {
    "authorcontributions": "Author Contributions",
    "funding": "Funding",
    "institutionalreview": "Institutional Review Board Statement",
    "informedconsent": "Informed Consent Statement",
    "dataavailability": "Data Availability Statement",
    "acknowledgments": "Acknowledgments",
    "conflictsofinterest": "Conflicts of Interest",
}

CAPTION_PREFIXES = {
    r"\caption{Overview": r"\caption{Figure 1. Overview",
    r"\caption{De novo motif": r"\caption{Figure 2. De novo motif",
    r"\caption{Interactive differential": r"\caption{Figure 3. Interactive differential",
    r"\caption{Single-cell footprinting": r"\caption{Figure 4. Single-cell footprinting",
    r"\caption{Single-machine runtime": r"\caption{Figure 5. Single-machine runtime",
    r"\caption{Browser GUI": r"\caption{Figure 6. Browser GUI",
    r"\caption{Qualitative scope": r"\caption{Table 1. Qualitative scope",
    r"\caption{Primary analysis": r"\caption{Table S1. Primary analysis",
    r"\caption{Software and resource": r"\caption{Table S2. Software and resource",
}


def run(cmd: list[str], cwd: Path | None = None) -> None:
    subprocess.run(cmd, cwd=cwd, check=True)


def brace_content(text: str, command: str) -> str:
    marker = f"\\{command}"
    start = text.find(marker)
    if start < 0:
        return ""
    brace = text.find("{", start)
    if brace < 0:
        return ""
    depth = 0
    for idx in range(brace, len(text)):
        char = text[idx]
        if char == "{":
            depth += 1
        elif char == "}":
            depth -= 1
            if depth == 0:
                return text[brace + 1 : idx].strip()
    return ""


def replace_command_block(text: str, command: str, replacement_name: str) -> str:
    marker = f"\\{command}"
    pieces: list[str] = []
    pos = 0
    while True:
        start = text.find(marker, pos)
        if start < 0:
            pieces.append(text[pos:])
            break
        brace = text.find("{", start)
        if brace < 0:
            pieces.append(text[pos:])
            break
        depth = 0
        end = None
        for idx in range(brace, len(text)):
            if text[idx] == "{":
                depth += 1
            elif text[idx] == "}":
                depth -= 1
                if depth == 0:
                    end = idx
                    break
        if end is None:
            pieces.append(text[pos:])
            break
        content = text[brace + 1 : end].strip()
        pieces.append(text[pos:start])
        pieces.append(f"\\noindent\\textbf{{{replacement_name}:}} {content}\n\n")
        pos = end + 1
    return "".join(pieces)


def render_docx_figures() -> None:
    MEDIA_DIR.mkdir(parents=True, exist_ok=True)
    for fig in ("Fig1", "Fig2", "Fig4", "Fig5"):
        source = MANUSCRIPT_DIR / "figures" / f"{fig}.pdf"
        target_prefix = MEDIA_DIR / fig
        run(["pdftocairo", "-png", "-singlefile", "-r", "300", str(source), str(target_prefix)])
    for fig in ("Fig3", "Fig6"):
        shutil.copy2(MANUSCRIPT_DIR / "figures" / f"{fig}.png", MEDIA_DIR / f"{fig}.png")


def simplify_table_lines(text: str) -> str:
    lines: list[str] = []
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith(r"\begin{tabularx}"):
            if "*{6}" in stripped:
                line = r"\begin{tabular}{lllllll}"
            elif "0.28" in stripped:
                line = r"\begin{tabular}{lll}"
            else:
                line = r"\begin{tabular}{ll}"
        elif stripped == r"\end{tabularx}":
            line = r"\end{tabular}"
        elif stripped in {r"\begingroup", r"\endgroup", r"\scriptsize", r"\footnotesize"}:
            continue
        elif stripped.startswith((r"\setlength", r"\renewcommand", r"\newcommand", r"\vspace", r"\let")):
            continue
        elif stripped in {r"\toprule", r"\midrule", r"\bottomrule"}:
            line = r"\hline"
        lines.append(line)
    return "\n".join(lines)


def prepare_docx_tex(source: Path) -> Path:
    text = source.read_text()
    title = brace_content(text, "Title")
    author = brace_content(text, "AuthorNames") or brace_content(text, "Author")
    address = brace_content(text, "address").replace("%", "").strip()
    address = re.sub(r"\$\^\{[^{}]*\}\$", "", address)
    address = address.replace(r"\quad", "").strip()
    address = re.sub(r"\s+", " ", address)
    abstract = brace_content(text, "abstract")
    keywords = brace_content(text, "keyword")

    body_start = text.index(r"\begin{document}") + len(r"\begin{document}")
    body_end = text.index(r"\bibliography{references}")
    body = text[body_start:body_end]

    for label, value in REF_REPLACEMENTS.items():
        body = body.replace(rf"\ref{{{label}}}", value)
    for fig, width in FIGURE_WIDTHS.items():
        body = re.sub(
            rf"\\includegraphics\[[^\]]*\]\{{figures/{fig}\.(?:pdf|png)\}}",
            rf"\\includegraphics[width={width}]{{.docx_build/figures/{fig}.png}}",
            body,
        )
    for command, replacement_name in BACK_MATTER.items():
        body = replace_command_block(body, command, replacement_name)
    for old, new in CAPTION_PREFIXES.items():
        body = body.replace(old, new)

    body = body.replace(r"\featureyes", r"\checkmark")
    body = body.replace(r"\featurepartial", r"$\circ$")
    body = re.sub(r"\\shortstack\{([^{}]*)\}", lambda m: m.group(1).replace(r"\\", " "), body)
    body = re.sub(r"\\label\{[^{}]*\}", "", body)
    body = body.replace(r"\appendix", "")
    body = simplify_table_lines(body)

    docx_source = f"""\\documentclass{{article}}
\\usepackage{{graphicx}}
\\usepackage{{booktabs}}
\\usepackage{{amsmath}}
\\usepackage{{amssymb}}
\\title{{{title}}}
\\author{{{author}\\\\{address}}}
\\date{{}}

\\begin{{document}}
\\maketitle
\\begin{{abstract}}
{abstract}
\\end{{abstract}}

\\noindent\\textbf{{Keywords:}} {keywords}

{body}

\\bibliography{{references}}
\\end{{document}}
"""
    BUILD_DIR.mkdir(parents=True, exist_ok=True)
    output = BUILD_DIR / "main_docx.tex"
    output.write_text(docx_source)
    return output


def build_docx(source: Path, output: Path, keep_build_dir: bool = False) -> None:
    shutil.rmtree(BUILD_DIR, ignore_errors=True)
    render_docx_figures()
    docx_tex = prepare_docx_tex(source)
    csl_args = ["--csl", str(VANCOUVER_CSL)] if VANCOUVER_CSL.exists() else []
    try:
        run(
            [
                "pandoc",
                "--from=latex",
                "--to=docx",
                "--citeproc",
                "--bibliography",
                "references.bib",
                *csl_args,
                "--resource-path",
                str(MANUSCRIPT_DIR),
                "--metadata",
                "link-citations=true",
                str(docx_tex.relative_to(MANUSCRIPT_DIR)),
                "-o",
                str(output.relative_to(MANUSCRIPT_DIR)),
            ],
            cwd=MANUSCRIPT_DIR,
        )
    finally:
        if not keep_build_dir:
            shutil.rmtree(BUILD_DIR, ignore_errors=True)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=DEFAULT_TEX)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--keep-build-dir", action="store_true")
    args = parser.parse_args(argv)
    build_docx(args.source, args.output, keep_build_dir=args.keep_build_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
