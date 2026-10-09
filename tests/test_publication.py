from tools import check_publication


def test_publication_snapshot_has_no_broken_local_markdown_links():
    assert check_publication.markdown_link_errors(check_publication.publication_files()) == []


def test_local_markdown_link_must_be_in_publication_snapshot(tmp_path, monkeypatch):
    document = tmp_path / "note.md"
    local_only = tmp_path / "local-only.json"
    document.write_text("[local result](local-only.json)\n", encoding="utf-8")
    local_only.write_text("{}\n", encoding="utf-8")
    monkeypatch.setattr(check_publication, "ROOT", tmp_path)

    assert check_publication.markdown_link_errors([document]) == [
        "broken local Markdown link in note.md: local-only.json"
    ]
    assert check_publication.markdown_link_errors([document, local_only]) == []
