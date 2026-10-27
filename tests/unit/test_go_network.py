"""Tests for GO-term overlap network construction."""

import pandas as pd
import pytest

from arabidopsis_stress_recovery_atlas.analysis.go_network import (
    calculate_bp_term_overlap,
    filter_strong_bp_overlap,
)


def make_membership() -> pd.DataFrame:
    """Create small gene-term memberships with known overlap."""

    return pd.DataFrame(
        {
            "gene_id": [
                "g1",
                "g2",
                "g3",
                "g2",
                "g3",
                "g4",
                "g5",
            ],
            "go_id": [
                "GO:A",
                "GO:A",
                "GO:A",
                "GO:B",
                "GO:B",
                "GO:B",
                "GO:C",
            ],
            "go_term": [
                "term A",
                "term A",
                "term A",
                "term B",
                "term B",
                "term B",
                "term C",
            ],
        }
    )


def test_pairwise_term_overlap_is_calculated():
    """Three GO terms should produce three unique term pairs."""

    overlap = calculate_bp_term_overlap(
        make_membership()
    )

    assert len(overlap) == 3


def test_jaccard_and_overlap_coefficient_are_correct():
    """Overlap measures should use the expected set formulas."""

    overlap = calculate_bp_term_overlap(
        make_membership()
    )

    pair = overlap[
        (
            overlap["go_id_a"]
            == "GO:A"
        )
        & (
            overlap["go_id_b"]
            == "GO:B"
        )
    ].iloc[0]

    assert pair["intersection"] == 2
    assert pair["union"] == 4

    assert pair["jaccard"] == pytest.approx(
        2 / 4
    )

    assert (
        pair["overlap_coefficient"]
        == pytest.approx(
            2 / 3
        )
    )


def test_strong_overlap_filter_uses_threshold():
    """Only pairs meeting the overlap threshold should remain."""

    overlap = calculate_bp_term_overlap(
        make_membership()
    )

    strong = filter_strong_bp_overlap(
        overlap,
        threshold=0.50,
    )

    assert len(strong) == 1

    assert {
        strong.iloc[0]["go_id_a"],
        strong.iloc[0]["go_id_b"],
    } == {
        "GO:A",
        "GO:B",
    }