"""Tests for Gene Ontology annotation transformations."""

import pandas as pd

from arabidopsis_stress_recovery_atlas.annotation.go import (
    get_go_term_memberships,
    get_goslim_memberships,
)


def make_go_annotations() -> pd.DataFrame:
    """Create small GO annotation data with duplicated evidence records."""

    return pd.DataFrame(
        {
            "gene_id": [
                "AT1G00010",
                "AT1G00010",
                "AT1G00020",
            ],
            "go_id": [
                "GO:0000001",
                "GO:0000001",
                "GO:0000002",
            ],
            "go_term": [
                "term A",
                "term A",
                "term B",
            ],
            "aspect": [
                "P",
                "P",
                "F",
            ],
            "goslim_category": [
                "response to stress",
                "response to stress",
                "enzyme activity",
            ],
        }
    )


def test_go_memberships_remove_duplicate_evidence_rows():
    """Repeated evidence should not duplicate a gene-term membership."""

    memberships = get_go_term_memberships(
        make_go_annotations()
    )

    assert len(memberships) == 2


def test_go_memberships_can_filter_biological_process():
    """Aspect P should retain Biological Process memberships."""

    memberships = get_go_term_memberships(
        make_go_annotations(),
        aspect="P",
    )

    assert len(memberships) == 1
    assert memberships.iloc[0]["gene_id"] == "AT1G00010"


def test_goslim_memberships_are_unique():
    """Repeated annotations should collapse to unique gene-category pairs."""

    memberships = get_goslim_memberships(
        make_go_annotations(),
        aspect="P",
    )

    assert len(memberships) == 1
    assert (
        memberships.iloc[0]["goslim_category"]
        == "response to stress"
    )