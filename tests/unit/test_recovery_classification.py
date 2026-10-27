"""Tests for endpoint/recovery shifts and recovery classification."""

import pandas as pd

from arabidopsis_stress_recovery_atlas.analysis.exploratory import (
    calculate_condition_shifts,
    calculate_persistence_breadth,
    classify_recovery_outcomes,
    summarize_persistent_counts,
)


def make_test_expression() -> pd.DataFrame:
    """Create a small expression table with known recovery outcomes."""

    return pd.DataFrame(
        {
            "gene_id": [
                "gene_persistent",
                "gene_recovered",
                "gene_reversal",
                "gene_ineligible",
            ],
            "CONT 1": [1.0, 4.0, 5.0, 2.0],
            "CONT 2": [1.0, 4.0, 5.0, 2.0],
            "CONT 3": [1.0, 4.0, 5.0, 2.0],
            "endpoint_DROUGHT 1": [3.0, 2.0, 3.0, 2.5],
            "endpoint_DROUGHT 2": [3.0, 2.0, 3.0, 2.5],
            "endpoint_DROUGHT 3": [3.0, 2.0, 3.0, 2.5],
            "recovery_DROUGHT 1": [2.0, 3.5, 6.0, 2.2],
            "recovery_DROUGHT 2": [2.0, 3.5, 6.0, 2.2],
            "recovery_DROUGHT 3": [2.0, 3.5, 6.0, 2.2],
        }
    )


def make_test_manifest() -> pd.DataFrame:
    """Create matching sample metadata."""

    return pd.DataFrame(
        {
            "sample_id": [
                "CONT 1",
                "CONT 2",
                "CONT 3",
                "endpoint_DROUGHT 1",
                "endpoint_DROUGHT 2",
                "endpoint_DROUGHT 3",
                "recovery_DROUGHT 1",
                "recovery_DROUGHT 2",
                "recovery_DROUGHT 3",
            ],
            "state": [
                "control",
                "control",
                "control",
                "endpoint",
                "endpoint",
                "endpoint",
                "recovery",
                "recovery",
                "recovery",
            ],
            "condition": [
                "control",
                "control",
                "control",
                "drought",
                "drought",
                "drought",
                "drought",
                "drought",
                "drought",
            ],
        }
    )


def test_condition_shifts_are_calculated_from_control():
    """Endpoint and recovery shifts should use the control mean."""

    shifts = calculate_condition_shifts(
        make_test_expression(),
        make_test_manifest(),
    )

    persistent = shifts[
        shifts["gene_id"] == "gene_persistent"
    ].iloc[0]

    assert persistent["control_mean"] == 1.0
    assert persistent["endpoint_mean"] == 3.0
    assert persistent["recovery_mean"] == 2.0
    assert persistent["endpoint_shift"] == 2.0
    assert persistent["recovery_shift"] == 1.0


def test_recovery_outcomes_are_classified_correctly():
    """The response-floor rules should distinguish recovery outcomes."""

    shifts = calculate_condition_shifts(
        make_test_expression(),
        make_test_manifest(),
    )

    classified = classify_recovery_outcomes(
        shifts,
        response_floor=0.75,
    )

    outcomes = classified.set_index(
        "gene_id"
    )["recovery_outcome"].to_dict()

    assert (
        outcomes["gene_persistent"]
        == "persistent_same_direction"
    )

    assert (
        outcomes["gene_recovered"]
        == "recovered"
    )

    assert (
        outcomes["gene_reversal"]
        == "substantial_reversal"
    )

    assert (
        outcomes["gene_ineligible"]
        == "not_endpoint_eligible"
    )


def test_persistent_counts_are_summarized_by_condition():
    """Persistent genes should be counted within each stress condition."""

    shifts = calculate_condition_shifts(
        make_test_expression(),
        make_test_manifest(),
    )

    classified = classify_recovery_outcomes(
        shifts,
        response_floor=0.75,
    )

    counts = summarize_persistent_counts(
        classified
    )

    assert counts == {
        "drought": 1
    }


def test_persistence_breadth_counts_conditions_per_gene():
    """Breadth should count persistent conditions for each gene."""

    classified = pd.DataFrame(
        {
            "gene_id": [
                "gene_a",
                "gene_a",
                "gene_b",
                "gene_b",
                "gene_c",
            ],
            "condition": [
                "cold",
                "drought",
                "cold",
                "drought",
                "cold",
            ],
            "recovery_outcome": [
                "persistent_same_direction",
                "persistent_same_direction",
                "persistent_same_direction",
                "recovered",
                "recovered",
            ],
        }
    )

    breadth = calculate_persistence_breadth(
        classified
    ).set_index("gene_id")

    assert (
        breadth.loc[
            "gene_a",
            "persistence_breadth",
        ]
        == 2
    )

    assert (
        breadth.loc[
            "gene_b",
            "persistence_breadth",
        ]
        == 1
    )

    assert (
        breadth.loc[
            "gene_c",
            "persistence_breadth",
        ]
        == 0
    )