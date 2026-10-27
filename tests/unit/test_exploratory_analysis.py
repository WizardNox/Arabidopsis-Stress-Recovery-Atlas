"""Tests for exploratory scientific calculations."""

import pandas as pd
import pytest

from arabidopsis_stress_recovery_atlas.analysis.exploratory import (
    ExploratoryAnalysisError,
    calculate_control_pairwise_differences,
    calculate_response_floor,
    calculate_response_floors,
)


def make_test_expression() -> pd.DataFrame:
    """Create a small canonical expression table."""

    return pd.DataFrame(
        {
            "gene_id": [
                "AT1G00010",
                "AT1G00020",
            ],
            "CONT 1": [0.0, 2.0],
            "CONT 2": [1.0, 2.0],
            "CONT 3": [2.0, 2.0],
        }
    )


def make_test_manifest() -> pd.DataFrame:
    """Create a small sample manifest."""

    return pd.DataFrame(
        {
            "sample_id": [
                "CONT 1",
                "CONT 2",
                "CONT 3",
            ],
            "state": [
                "control",
                "control",
                "control",
            ],
        }
    )


def test_control_pairwise_differences_are_correct():
    """All unique control pairs should be compared for every gene."""

    expression = make_test_expression()
    manifest = make_test_manifest()

    differences = (
        calculate_control_pairwise_differences(
            expression,
            manifest,
        )
    )

    assert len(differences) == 6

    assert sorted(
        differences["absolute_difference"].tolist()
    ) == [
        0.0,
        0.0,
        0.0,
        1.0,
        1.0,
        2.0,
    ]


def test_response_floor_uses_requested_quantile():
    """The response floor should use the supplied quantile."""

    expression = make_test_expression()
    manifest = make_test_manifest()

    differences = (
        calculate_control_pairwise_differences(
            expression,
            manifest,
        )
    )

    response_floor = calculate_response_floor(
        differences,
        0.50,
    )

    assert response_floor == pytest.approx(0.5)


def test_multiple_response_floors_are_labeled():
    """Multiple quantiles should receive percentile labels."""

    expression = make_test_expression()
    manifest = make_test_manifest()

    differences = (
        calculate_control_pairwise_differences(
            expression,
            manifest,
        )
    )

    floors = calculate_response_floors(
        differences,
        [0.50, 0.75],
    )

    assert set(floors) == {
        "p50",
        "p75",
    }


def test_invalid_quantile_raises_error():
    """Quantiles outside the valid range should be rejected."""

    differences = pd.DataFrame(
        {
            "absolute_difference": [
                0.1,
                0.2,
            ]
        }
    )

    with pytest.raises(
        ExploratoryAnalysisError
    ):
        calculate_response_floor(
            differences,
            1.5,
        )