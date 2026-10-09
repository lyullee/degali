"""Build a Word manuscript draft from the canonical IJHE Markdown source.

The builder is intentionally author-neutral.  It preserves unresolved author
placeholders until the corresponding author supplies them and can be run with
``--strict`` to refuse a final export while any placeholder remains.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

from docx import Document
from docx.enum.section import WD_SECTION
from docx.enum.style import WD_STYLE_TYPE
from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT, WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, RGBColor


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE = ROOT / "MANUSCRIPT_LH2_VALIDATION.md"

# Repository paths are useful in the reproducibility Markdown, but they are
# distracting and non-portable when the same source is rendered as a journal
# manuscript.  Keep the links in the source package while replacing their
# visible labels in DOCX/PDF output with reviewer-facing descriptions.
JOURNAL_PATH_LABELS = {
    "docs/lh2-paper-baseline-2026-09-20.md": "the frozen baseline record",
    "docs/transient-dense-gas-3d.md": "the transient-extension scope record",
    "outputs/transient-3d-error-study-2026-10-08/error_study.md": "the transient error-study record",
    "docs/ijhe-preslhy-trial10-conditional-evidence-2026-10-09.md": "the conditional Trial 10 evidence record",
    "docs/ijhe-ffi-site-geometry-boundary-2026-10-09.md": "the FFI site-geometry boundary record",
    "docs/ijhe-elvhys-conditional-evidence-2026-10-09.md": "the ELVHYS boundary record",
    "docs/ijhe-smedis-transfer-boundary-2026-10-09.md": "the SMEDIS transfer-boundary record",
    "docs/open-channel-h2-boundary.md": "the open-channel H2 boundary record",
    "docs/ffi-test6-decomposition.md": "the Test 6 decomposition record",
    "docs/applied-energy-evidence-2026-10-03.md": "the paired comparison evidence record",
    "docs/time-aligned-uncertainty-results-2026-10-03.md": "the time-aligned uncertainty record",
    "docs/decision-oriented-scenarios-2026-10-03.md": "the decision-scenario record",
    "tools/audit_lh2_paper_baseline.py": "the baseline audit script",
    "tools/audit_applied_energy_comparison.py": "the paired comparison audit script",
    "tools/audit_time_aligned_uncertainty.py": "the time-aligned audit script",
    "tools/audit_decision_scenarios.py": "the decision-scenario audit script",
    "tools/audit_open_channel_h2_dataset.py": "the open-channel H2 audit script",
}
# The combined figure PDF is an evidence bundle rather than the manuscript's
# figure-number order.  Keep an explicit, review-file mapping so each caption
# is paired with the correct panel when the six main figures are embedded.
FIGURE_IMAGE_DIR = ROOT / "outputs" / "ijhe-main-figure-images-2026-10-09-clean"


def _set_font(style, name: str, size: float, bold: bool = False, italic: bool = False) -> None:
    style.font.name = name
    style.font.size = Pt(size)
    style.font.bold = bold
    style.font.italic = italic
    style.font.color.rgb = RGBColor(0, 0, 0)
    style._element.rPr.rFonts.set(qn("w:ascii"), name)
    style._element.rPr.rFonts.set(qn("w:hAnsi"), name)


def _set_cell_shading(cell, fill: str) -> None:
    tc_pr = cell._tc.get_or_add_tcPr()
    shd = tc_pr.find(qn("w:shd"))
    if shd is None:
        shd = OxmlElement("w:shd")
        tc_pr.append(shd)
    shd.set(qn("w:fill"), fill)


def _set_cell_borders(cell, color: str = "D9D9D9") -> None:
    tc_pr = cell._tc.get_or_add_tcPr()
    borders = tc_pr.first_child_found_in("w:tcBorders")
    if borders is None:
        borders = OxmlElement("w:tcBorders")
        tc_pr.append(borders)
    for edge in ("top", "left", "bottom", "right", "insideH", "insideV"):
        tag = "w:" + edge
        element = borders.find(qn(tag))
        if element is None:
            element = OxmlElement(tag)
            borders.append(element)
        element.set(qn("w:val"), "single")
        element.set(qn("w:sz"), "4")
        element.set(qn("w:space"), "0")
        element.set(qn("w:color"), color)


def _set_cell_margins(cell, top: int = 90, start: int = 100, bottom: int = 90, end: int = 100) -> None:
    tc = cell._tc
    tc_pr = tc.get_or_add_tcPr()
    tc_mar = tc_pr.first_child_found_in("w:tcMar")
    if tc_mar is None:
        tc_mar = OxmlElement("w:tcMar")
        tc_pr.append(tc_mar)
    for name, value in (("top", top), ("start", start), ("bottom", bottom), ("end", end)):
        node = tc_mar.find(qn("w:" + name))
        if node is None:
            node = OxmlElement("w:" + name)
            tc_mar.append(node)
        node.set(qn("w:w"), str(value))
        node.set(qn("w:type"), "dxa")


def _set_row_cant_split(row) -> None:
    """Keep a table row together when Word paginates the document."""

    tr_pr = row._tr.get_or_add_trPr()
    if tr_pr.find(qn("w:cantSplit")) is None:
        tr_pr.append(OxmlElement("w:cantSplit"))


def _set_row_header(row) -> None:
    """Repeat a table header row when a table continues on a new page."""

    tr_pr = row._tr.get_or_add_trPr()
    if tr_pr.find(qn("w:tblHeader")) is None:
        tr_pr.append(OxmlElement("w:tblHeader"))


def _add_page_number(paragraph) -> None:
    run = paragraph.add_run()
    fld_char1 = OxmlElement("w:fldChar")
    fld_char1.set(qn("w:fldCharType"), "begin")
    instr = OxmlElement("w:instrText")
    instr.set(qn("xml:space"), "preserve")
    instr.text = " PAGE "
    fld_char2 = OxmlElement("w:fldChar")
    fld_char2.set(qn("w:fldCharType"), "end")
    run._r.append(fld_char1)
    run._r.append(instr)
    run._r.append(fld_char2)


def _remove_paragraph_borders(paragraph) -> None:
    p_pr = paragraph._p.get_or_add_pPr()
    borders = p_pr.find(qn("w:pBdr"))
    if borders is not None:
        p_pr.remove(borders)


def _remove_style_borders(style) -> None:
    p_pr = style._element.get_or_add_pPr()
    borders = p_pr.find(qn("w:pBdr"))
    if borders is not None:
        p_pr.remove(borders)


def _plain_inline(text: str) -> str:
    text = re.sub(r"Artifact:\s*`[^`]+`", "Figure file supplied separately", text)
    text = re.sub(r"!\[([^]]*)\]\([^)]*\)", r"\1", text)
    text = re.sub(r"\[([^]]+)\]\([^)]*\)", r"\1", text)
    for path, label in JOURNAL_PATH_LABELS.items():
        text = text.replace(path, label)
    text = text.replace("**", "").replace("__", "")
    text = text.replace("`", "")
    text = text.replace("*", "")
    return text.strip()


def _figure_number(text: str) -> int | None:
    match = re.match(r"\*\*Figure\s+(\d+)\.", text.strip())
    return int(match.group(1)) if match else None


def _add_inline_figure(doc: Document, raw: str) -> bool:
    """Insert a locally rendered figure before its caption.

    IJHE's Your Paper Your Way guide allows a single Word/PDF file for review,
    but asks that figures appear near the relevant text.  The canonical SVGs
    remain separate upload artifacts; these cropped PNGs are only the review
    copy embedded in the author-neutral DOCX.
    """

    number = _figure_number(raw)
    if number is None:
        return False
    image = FIGURE_IMAGE_DIR / f"figure-{number}.png"
    caption_text = re.sub(r"\s*Artifact:\s*`[^`]+`\.?", "", raw).strip()
    if image.is_file():
        image_paragraph = doc.add_paragraph()
        image_paragraph.alignment = WD_ALIGN_PARAGRAPH.CENTER
        image_paragraph.paragraph_format.keep_with_next = True
        run = image_paragraph.add_run()
        run.add_picture(str(image), width=Inches(6.35))
    caption = doc.add_paragraph(style="Caption")
    caption.add_run(_plain_inline(caption_text))
    return True


def _parse_table(lines: list[str], index: int) -> tuple[list[list[str]], int]:
    rows: list[list[str]] = []
    while index < len(lines) and lines[index].strip().startswith("|"):
        raw = lines[index].strip().strip("|")
        cells = [_plain_inline(part.strip()) for part in raw.split("|")]
        if not all(re.fullmatch(r":?-+:?", cell) for cell in cells):
            rows.append(cells)
        index += 1
    return rows, index


def _add_table(doc: Document, rows: list[list[str]]) -> None:
    if not rows:
        return
    width = max(len(row) for row in rows)
    table = doc.add_table(rows=1, cols=width)
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    table.style = "Table Grid"
    header = table.rows[0].cells
    _set_row_cant_split(table.rows[0])
    _set_row_header(table.rows[0])
    for col in range(width):
        header[col].text = rows[0][col] if col < len(rows[0]) else ""
        header[col].vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER
        _set_cell_shading(header[col], "1F4E78")
        _set_cell_borders(header[col])
        _set_cell_margins(header[col])
        for run in header[col].paragraphs[0].runs:
            run.font.bold = True
            run.font.color.rgb = RGBColor(255, 255, 255)
            run.font.size = Pt(8.5)
        for paragraph in header[col].paragraphs:
            paragraph.paragraph_format.keep_with_next = True
    for row_index, row in enumerate(rows[1:]):
        cells = table.add_row().cells
        _set_row_cant_split(table.rows[-1])
        for col in range(width):
            cells[col].text = row[col] if col < len(row) else ""
            cells[col].vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER
            if row_index % 2:
                _set_cell_shading(cells[col], "F2F6F9")
            _set_cell_borders(cells[col])
            _set_cell_margins(cells[col])
            for paragraph in cells[col].paragraphs:
                paragraph.paragraph_format.space_after = Pt(1)
                for run in paragraph.runs:
                    run.font.size = Pt(8.5)
    doc.add_paragraph().paragraph_format.space_after = Pt(2)


def build(
    source: Path,
    output: Path,
    author_list: str | None = None,
    corresponding: str | None = None,
    strict: bool = False,
    supplementary: bool = False,
) -> dict[str, object]:
    text = source.read_text(encoding="utf-8")
    placeholders = "[AUTHOR TO CONFIRM" in text or "[CORRESPONDING AUTHOR" in text
    declaration_placeholders = "To be completed by the authors" in text
    unresolved_authors = (not supplementary) and (placeholders or not author_list or not corresponding)
    if strict and (unresolved_authors or declaration_placeholders):
        if declaration_placeholders and not unresolved_authors:
            raise ValueError(
                "final export requires completed author declarations; "
                "omit --strict for a draft"
            )
        raise ValueError(
            "final export requires both author_list and corresponding; "
            "omit --strict for a draft"
        )

    doc = Document()
    section = doc.sections[0]
    section.top_margin = Inches(0.8)
    section.bottom_margin = Inches(0.8)
    section.left_margin = Inches(0.85)
    section.right_margin = Inches(0.85)

    normal = doc.styles["Normal"]
    _set_font(normal, "Arial", 10.5)
    normal.paragraph_format.space_after = Pt(6)
    normal.paragraph_format.line_spacing = 1.08
    for name, size in (("Title", 16), ("Heading 1", 13), ("Heading 2", 11.5), ("Heading 3", 10.5)):
        _set_font(doc.styles[name], "Arial", size, bold=True)
        doc.styles[name].paragraph_format.space_before = Pt(10 if name != "Title" else 0)
        doc.styles[name].paragraph_format.space_after = Pt(5)
        _remove_style_borders(doc.styles[name])
    if "Caption" not in doc.styles:
        caption_style = doc.styles.add_style("Caption", WD_STYLE_TYPE.PARAGRAPH)
    else:
        caption_style = doc.styles["Caption"]
    _set_font(caption_style, "Arial", 9, italic=True)
    caption_style.paragraph_format.space_before = Pt(4)
    caption_style.paragraph_format.space_after = Pt(5)

    footer = section.footer.paragraphs[0]
    footer.alignment = WD_ALIGN_PARAGRAPH.CENTER
    footer_run = footer.add_run("International Journal of Hydrogen Energy — ")
    footer_run.font.name = "Arial"
    footer_run.font.size = Pt(8)
    _add_page_number(footer)

    lines = text.splitlines()
    index = 0
    paragraph_buffer: list[str] = []
    seen_title = False

    def flush() -> None:
        if not paragraph_buffer:
            return
        raw = " ".join(part.strip() for part in paragraph_buffer).strip()
        paragraph_buffer.clear()
        if not raw:
            return
        if raw.startswith("**Figure "):
            if _add_inline_figure(doc, raw):
                return
            p = doc.add_paragraph(style="Caption")
        elif raw.startswith("**Table "):
            p = doc.add_paragraph(style="Caption")
            p.paragraph_format.keep_with_next = raw.startswith("**Table ")
        else:
            p = doc.add_paragraph()
        p.add_run(_plain_inline(raw))

    while index < len(lines):
        line = lines[index]
        if line.startswith("# ") and not seen_title:
            flush()
            p = doc.add_paragraph(style="Title")
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER
            _remove_paragraph_borders(p)
            p.add_run(_plain_inline(line[2:]))
            if not supplementary:
                author = author_list or "[AUTHOR TO CONFIRM: complete author list, order, affiliations, postal addresses and ORCIDs]"
                contact = corresponding or "[CORRESPONDING AUTHOR TO CONFIRM: name, institutional email, postal address and ORCID]"
                meta = doc.add_paragraph()
                meta.alignment = WD_ALIGN_PARAGRAPH.CENTER
                meta.add_run(author + "\n" + contact).italic = True
            seen_title = True
            index += 1
            continue
        heading = re.match(r"^(#{2,4})\s+(.*)$", line)
        if heading:
            flush()
            level = min(len(heading.group(1)) - 1, 3)
            p = doc.add_paragraph(style=f"Heading {level}")
            if "5.9 Decision-boundary" in heading.group(2):
                p.paragraph_format.keep_with_next = True
            p.add_run(_plain_inline(heading.group(2)))
            index += 1
            continue
        if line.strip().startswith("|"):
            flush()
            rows, index = _parse_table(lines, index)
            _add_table(doc, rows)
            continue
        if line.strip().startswith("```"):
            flush()
            index += 1
            code: list[str] = []
            while index < len(lines) and not lines[index].strip().startswith("```"):
                code.append(lines[index])
                index += 1
            if index < len(lines):
                index += 1
            p = doc.add_paragraph()
            p.paragraph_format.left_indent = Inches(0.25)
            run = p.add_run("\n".join(code))
            run.font.name = "Courier New"
            run.font.size = Pt(8.5)
            continue
        if not line.strip():
            flush()
        else:
            paragraph_buffer.append(line)
        index += 1
    flush()

    output.parent.mkdir(parents=True, exist_ok=True)
    doc.save(output)
    return {"output": str(output), "placeholder_state": "present" if unresolved_authors else "resolved_or_not_required"}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--author-list")
    parser.add_argument("--corresponding")
    parser.add_argument("--strict", action="store_true")
    parser.add_argument(
        "--supplementary",
        action="store_true",
        help="build a supplementary-information file without author metadata placeholders",
    )
    args = parser.parse_args()
    try:
        record = build(
            args.source,
            args.output,
            args.author_list,
            args.corresponding,
            args.strict,
            args.supplementary,
        )
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    print(record)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
