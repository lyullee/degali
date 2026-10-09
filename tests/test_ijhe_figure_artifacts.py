from tools.audit_ijhe_figure_artifacts import audit


def test_graphical_abstract_meets_elsevier_minimum_dimensions():
    payload = audit()
    record = next(item for item in payload["figures"] if item["path"] == "IJHE_GRAPHICAL_ABSTRACT.svg")
    assert record["status"] == "valid"
    assert record["pixel_width"] >= 1328
    assert record["pixel_height"] >= 531
