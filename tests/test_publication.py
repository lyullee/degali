from tools.check_publication import markdown_link_errors, publication_files


def test_publication_snapshot_has_no_broken_local_markdown_links():
    assert markdown_link_errors(publication_files()) == []
