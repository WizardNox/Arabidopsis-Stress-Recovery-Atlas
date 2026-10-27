"""Tests for persistence-breadth and robustness summaries."""

import pandas as pd
import pytest

from arabidopsis_stress_recovery_atlas.analysis.response_programs import (
    calculate_persistence_rank_correlation,
    get_persistent_gene_condition_pairs,
    summarize_breadth_distribution,
    summarize_condition_persistence_fractions,
)


def test_exact_breadth_distribution_is_counted():
    """Genes should be counted at their exact persistence breadth."""

    classified = pd.DataFrame(
        {
            "gene_id": [
                "gene_a",
                "gene_a",
                "gene_a",
                "gene_b",
                "gene_b",
                "gene_c",
            ],
            "condition": [
                "cold",
                "drought",
                "heat",
                "cold",
                "drought",
                "cold",
            ],
            "recovery_outcome": [
                "persistent_same_direction",
                "persistent_same_direction",
                "persistent_same_direction",
                "persistent_same_direction",
                "persistent_same_direction",
                "persistent_same_direction",
            ],
        }
    )

    distribution = summarize_breadth_distribution(
        classified,
        maximum_conditions=3,
    )

    assert distribution == {
        1: 1,
        2: 1,
        3: 1,
    }


def test_condition_persistence_fraction_is_calculated():
    """Persistence fraction should use endpoint-eligible genes as denominator."""

    classified = pd.DataFrame(
        {
            "condition": [
                "cold",
                "cold",
                "cold",
                "drought",
                "drought",
            ],
            "endpoint_eligible": [
                True,
                True,
                False,
                True,
                True,
            ],
            "recovery_outcome": [
                "persistent_same_direction",
                "recovered",
                "not_endpoint_eligible",
                "persistent_same_direction",
                "persistent_same_direction",
            ],
        }
    )

    summary = (
        summarize_condition_persistence_fractions(
            classified
        )
        .set_index("condition")
    )

    assert (
        summary.loc[
            "cold",
            "persistence_fraction_percent",
        ]
        == pytest.approx(50.0)
    )

    assert (
        summary.loc[
            "drought",
            "persistence_fraction_percent",
        ]
        == pytest.approx(100.0)
    )


def test_rank_correlation_uses_condition_ranks():
    """Persistence-fraction ranking should be compared across thresholds."""

    summary_a = pd.DataFrame(
        {
            "condition": [
                "cold",
                "drought",
                "heat",
            ],
            "persistence_fraction_percent": [
                10.0,
                30.0,
                20.0,
            ],
        }
    )

    summary_b = pd.DataFrame(
        {
            "condition": [
                "heat",
                "cold",
                "drought",
            ],
            "persistence_fraction_percent": [
                25.0,
                5.0,
                40.0,
            ],
        }
    )

    correlation = (
        calculate_persistence_rank_correlation(
            summary_a,
            summary_b,
        )
    )

    assert correlation == pytest.approx(1.0)


def test_persistent_gene_condition_pairs_are_returned():
    """Only persistent same-direction records should enter the set."""

    classified = pd.DataFrame(
        {
            "gene_id": [
                "gene_a",
                "gene_a",
                "gene_b",
            ],
            "condition": [
                "cold",
                "drought",
                "cold",
            ],
            "recovery_outcome": [
                "persistent_same_direction",
                "recovered",
                "persistent_same_direction",
            ],
        }
    )

    pairs = get_persistent_gene_condition_pairs(
        classified
    )

    assert pairs == {
        ("gene_a", "cold"),
        ("gene_b", "cold"),
    }