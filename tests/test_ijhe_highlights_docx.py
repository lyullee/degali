from docx import Document

from tools.build_ijhe_highlights_docx import build, read_bullets


def test_ijhe_highlights_docx_is_standalone_and_portal_sized(tmp_path):
    output = tmp_path / "Highlights.docx"
    result = build(output)
    assert result["count"] == 5
    paragraphs = [p.text for p in Document(output).paragraphs if p.text.strip()]
    assert paragraphs == read_bullets()
    assert all(len(item) <= 85 for item in paragraphs)
