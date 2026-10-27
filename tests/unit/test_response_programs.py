"""Tests for persistent-expression pattern summaries."""

import pandas as pd
import pytest

from arabidopsis_stress_recovery_atlas.analysis.response_programs import (
    ResponseProgramError,
    summarize_persistence_pattern,
)


def make_test_classified_data() -> pd.DataFrame:
    """Create classified gene-condition records with known patterns."""

    return pd.DataFrame(
        {
            "gene_id": [
                "gene_a",
                "gene_a",
                "gene_a",
                "gene_a",
                "gene_b",
                "gene_b",
                "gene_c",
                "gene_c",
                "gene_d",
            ],
            "condition": [
                "cold",
                "drought",
                "heat",
                "high_light",
                "cold",
                "drought",
                "cold",
                "drought",
                "cold",
            ],
            "recovery_shift": [
                1.0,
                0.8,
                1.2,
                0.9,
                -1.0,
                -0.8,
                1.0,
                -1.0,
                1.0,
            ],
            "recovery_outcome": [
                "persistent_same_direction",
                "persistent_same_direction",
                "persistent_same_direction",
                "persistent_same_direction",
                "persistent_same_direction",
                "persistent_same_direction",
                "persistent_same_direction",
                "persistent_same_direction",
                "persistent_same_direction",
            ],
        }
    )


def test_persistence_pattern_counts_breadth():
    """High-breadth and multi-stress genes should be counted correctly."""

    summary = summarize_persistence_pattern(
        make_test_classified_data(),
        high_breadth_minimum=4,
        multi_stress_minimum=2,
    )

    assert summary["high_breadth_genes"] == 1
    assert summary["multi_stress_persistent_genes"] == 3


def test_directional_consistency_is_calculated():
    """Multi-stress genes should be evaluated across persistent conditions."""

    summary = summarize_persistence_pattern(
        make_test_classified_data()
    )

    assert (
        summary["directionally_consistent_genes"]
        == 2
    )

    assert (
        summary["consistently_above_control"]
        == 1
    )

    assert (
        summary["consistently_below_control"]
        == 1
    )

    assert (
        summary["directional_consistency_percent"]
        == pytest.approx(
            2 / 3 * 100
        )
    )


def test_single_stress_genes_are_not_used_for_directional_consistency():
    """Genes persistent in only one stress should not enter the comparison."""

    summary = summarize_persistence_pattern(
        make_test_classified_data()
    )

    assert (
        summary["multi_stress_persistent_genes"]
        == 3
    )


def test_missing_columns_raise_error():
    """Incomplete classified tables should be rejected."""

    data = pd.DataFrame(
        {
            "gene_id": ["gene_a"],
        }
    )

    with pytest.raises(ResponseProgramError):
        summarize_persistence_pattern(data)