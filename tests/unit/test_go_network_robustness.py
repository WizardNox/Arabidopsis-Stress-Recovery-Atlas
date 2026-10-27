"""Tests for GO-network threshold robustness."""

import pandas as pd
import pytest

from arabidopsis_stress_recovery_atlas.analysis.go_network_robustness import (
    compare_largest_components,
    summarize_network_thresholds,
)


def make_enrichment() -> pd.DataFrame:
    """Create four artificial significant GO terms."""

    return pd.DataFrame(
        {
            "go_id": [
                "GO:A",
                "GO:B",
                "GO:C",
                "GO:D",
            ],
            "go_term": [
                "term A",
                "term B",
                "term C",
                "term D",
            ],
        }
    )


def make_overlap() -> pd.DataFrame:
    """Create artificial GO-term relationships with known strengths."""

    return pd.DataFrame(
        {
            "go_id_a": [
                "GO:A",
                "GO:B",
                "GO:C",
            ],
            "go_term_a": [
                "term A",
                "term B",
                "term C",
            ],
            "go_id_b": [
                "GO:B",
                "GO:C",
                "GO:D",
            ],
            "go_term_b": [
                "term B",
                "term C",
                "term D",
            ],
            "genes_a": [
                10,
                10,
                10,
            ],
            "genes_b": [
                10,
                10,
                10,
            ],
            "intersection": [
                8,
                7,
                5,
            ],
            "union": [
                12,
                13,
                15,
            ],
            "jaccard": [
                8 / 12,
                7 / 13,
                5 / 15,
            ],
            "overlap_coefficient": [
                0.8,
                0.7,
                0.5,
            ],
        }
    )


def test_threshold_summary_reflects_stricter_network():
    """A stricter threshold should remove weaker edges."""

    summary, largest_components = (
        summarize_network_thresholds(
            make_enrichment(),
            make_overlap(),
            thresholds=[
                0.40,
                0.60,
            ],
        )
    )

    loose = summary[
        summary["overlap_threshold"]
        == 0.40
    ].iloc[0]

    strict = summary[
        summary["overlap_threshold"]
        == 0.60
    ].iloc[0]

    assert loose["n_nodes"] == 4
    assert loose["n_edges"] == 3
    assert loose[
        "largest_component_size"
    ] == 4

    assert strict["n_nodes"] == 4
    assert strict["n_edges"] == 2
    assert strict[
        "largest_component_size"
    ] == 3

    assert strict[
        "n_isolated_terms"
    ] == 1

    assert largest_components[
        0.60
    ] == {
        "GO:A",
        "GO:B",
        "GO:C",
    }


def test_largest_component_similarity_is_calculated():
    """Largest-component comparisons should use Jaccard similarity."""

    largest_components = {
        0.40: {
            "GO:A",
            "GO:B",
            "GO:C",
            "GO:D",
        },
        0.60: {
            "GO:A",
            "GO:B",
            "GO:C",
        },
    }

    comparison = (
        compare_largest_components(
            largest_components
        )
    )

    row = comparison.iloc[0]

    assert row["shared_terms"] == 3
    assert row["union_terms"] == 4

    assert row["jaccard"] == pytest.approx(
        3 / 4
    )