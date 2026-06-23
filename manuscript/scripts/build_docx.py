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
import tempfile
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path


W_NS = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
ET.register_namespace("w", W_NS)

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

DOCX_STYLE_SPECS = {
    "Normal": {"size": 20},
    "BodyText": {"size": 20},
    "FirstParagraph": {"size": 20},
    "Compact": {"size": 18},
    "Abstract": {"size": 20},
    "Bibliography": {"size": 18},
    "Title": {"size": 28, "bold": True},
    "Subtitle": {"size": 22, "bold": True},
    "Author": {"size": 20},
    "Date": {"size": 20},
    "AbstractTitle": {"size": 20, "bold": True},
    "Heading1": {"size": 24, "bold": True},
    "Heading2": {"size": 22, "bold": True},
    "Heading3": {"size": 20, "bold": True},
    "Heading4": {"size": 20, "bold": True},
    "Heading5": {"size": 20, "bold": True},
    "Heading6": {"size": 20, "bold": True},
    "Heading7": {"size": 20, "bold": True},
    "Heading8": {"size": 20, "bold": True},
    "Heading9": {"size": 20, "bold": True},
    "TOCHeading": {"size": 24, "bold": True},
    "Caption": {"size": 18, "italic": True},
    "ImageCaption": {"size": 18, "italic": True},
    "TableCaption": {"size": 18, "italic": True},
    "Table": {"size": 18},
}

DOCX_STYLE_SPACING = {
    "Normal": {"after": "120", "line": "240"},
    "BodyText": {"after": "120", "line": "240"},
    "FirstParagraph": {"after": "120", "line": "240"},
    "Compact": {"after": "60", "line": "220"},
    "Abstract": {"after": "120", "line": "240"},
    "Bibliography": {"after": "60", "line": "220"},
    "Title": {"after": "160", "line": "240"},
    "Author": {"after": "120", "line": "240"},
    "AbstractTitle": {"before": "120", "after": "60", "line": "240"},
    "Heading1": {"before": "260", "after": "100", "line": "240"},
    "Heading2": {"before": "180", "after": "80", "line": "240"},
    "Heading3": {"before": "140", "after": "60", "line": "240"},
    "TOCHeading": {"before": "260", "after": "100", "line": "240"},
    "Caption": {"before": "80", "after": "140", "line": "220"},
    "ImageCaption": {"before": "80", "after": "140", "line": "220"},
    "TableCaption": {"before": "80", "after": "100", "line": "220"},
}


def w_tag(name: str) -> str:
    return f"{{{W_NS}}}{name}"


def w_attr(name: str) -> str:
    return f"{{{W_NS}}}{name}"


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


def get_or_add(parent: ET.Element, tag: str) -> ET.Element:
    child = parent.find(w_tag(tag))
    if child is None:
        child = ET.SubElement(parent, w_tag(tag))
    return child


def set_or_remove_bool(rpr: ET.Element, tag: str, enabled: bool | None) -> None:
    child = rpr.find(w_tag(tag))
    if enabled is None:
        return
    if enabled:
        if child is None:
            ET.SubElement(rpr, w_tag(tag))
    elif child is not None:
        rpr.remove(child)


def set_arial_font(rpr: ET.Element) -> None:
    fonts = get_or_add(rpr, "rFonts")
    for attr in ("ascii", "hAnsi", "eastAsia", "cs"):
        fonts.set(w_attr(attr), "Arial")
    for attr in ("asciiTheme", "hAnsiTheme", "eastAsiaTheme", "cstheme"):
        fonts.attrib.pop(w_attr(attr), None)


def apply_run_style(rpr: ET.Element, size: int, *, bold: bool | None = None, italic: bool | None = None) -> None:
    set_arial_font(rpr)

    sz = get_or_add(rpr, "sz")
    sz.set(w_attr("val"), str(size))
    sz_cs = get_or_add(rpr, "szCs")
    sz_cs.set(w_attr("val"), str(size))

    color = get_or_add(rpr, "color")
    color.set(w_attr("val"), "000000")

    set_or_remove_bool(rpr, "b", bold)
    set_or_remove_bool(rpr, "bCs", bold)
    set_or_remove_bool(rpr, "i", italic)
    set_or_remove_bool(rpr, "iCs", italic)


def apply_paragraph_spacing(ppr: ET.Element, spacing_spec: dict[str, str]) -> None:
    spacing = get_or_add(ppr, "spacing")
    for attr, value in spacing_spec.items():
        spacing.set(w_attr(attr), value)
    if "line" in spacing_spec:
        spacing.set(w_attr("lineRule"), "auto")


def normalize_styles_xml(xml_bytes: bytes) -> bytes:
    root = ET.fromstring(xml_bytes)

    doc_defaults = get_or_add(root, "docDefaults")
    rpr_default = get_or_add(get_or_add(doc_defaults, "rPrDefault"), "rPr")
    apply_run_style(rpr_default, 20)

    for style in root.findall(w_tag("style")):
        existing_rpr = style.find(w_tag("rPr"))
        if existing_rpr is not None:
            set_arial_font(existing_rpr)

        style_id = style.get(w_attr("styleId"))
        if style_id is None:
            continue
        spec = DOCX_STYLE_SPECS.get(style_id)
        spacing_spec = DOCX_STYLE_SPACING.get(style_id)
        if spec is None and spacing_spec is None:
            continue
        if spec is not None:
            rpr = get_or_add(style, "rPr")
            apply_run_style(
                rpr,
                spec["size"],
                bold=spec.get("bold", False),
                italic=spec.get("italic", False),
            )
        if spacing_spec is not None:
            ppr = get_or_add(style, "pPr")
            apply_paragraph_spacing(ppr, spacing_spec)

    return ET.tostring(root, encoding="utf-8", xml_declaration=True)


def normalize_docx_styles(docx_path: Path) -> None:
    with tempfile.NamedTemporaryFile(delete=False, suffix=".docx") as tmp:
        tmp_path = Path(tmp.name)
    try:
        with zipfile.ZipFile(docx_path, "r") as zin, zipfile.ZipFile(tmp_path, "w", zipfile.ZIP_DEFLATED) as zout:
            for item in zin.infolist():
                data = zin.read(item.filename)
                if item.filename == "word/styles.xml":
                    data = normalize_styles_xml(data)
                zout.writestr(item, data)
        shutil.move(tmp_path, docx_path)
    finally:
        tmp_path.unlink(missing_ok=True)


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
    normalize_docx_styles(output)


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
