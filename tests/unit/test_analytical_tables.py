"""Tests for canonical analytical-table construction."""

import pandas as pd
import pytest

from arabidopsis_stress_recovery_atlas.transformation.analytical_tables import (
    AnalyticalTableError,
    build_gene_condition_response,
    build_gene_persistence_summary,
)


def make_classified() -> pd.DataFrame:
    """Create artificial classified gene-condition results."""

    return pd.DataFrame(
        {
            "gene_id": [
                "g1",
                "g1",
                "g2",
                "g2",
            ],
            "condition": [
                "cold",
                "drought",
                "cold",
                "drought",
            ],
            "control_mean": [
                5.0,
                5.0,
                4.0,
                4.0,
            ],
            "endpoint_mean": [
                7.0,
                8.0,
                4.1,
                2.0,
            ],
            "recovery_mean": [
                6.0,
                7.0,
                4.0,
                3.0,
            ],
            "endpoint_shift": [
                2.0,
                3.0,
                0.1,
                -2.0,
            ],
            "recovery_shift": [
                1.0,
                2.0,
                0.0,
                -1.0,
            ],
            "endpoint_eligible": [
                True,
                True,
                False,
                True,
            ],
            "recovery_above_floor": [
                True,
                True,
                False,
                True,
            ],
            "same_direction": [
                True,
                True,
                True,
                True,
            ],
            "recovery_outcome": [
                "persistent_same_direction",
                "persistent_same_direction",
                "not_endpoint_eligible",
                "persistent_same_direction",
            ],
        }
    )


def make_breadth() -> pd.DataFrame:
    """Create persistence breadth for every artificial gene."""

    return pd.DataFrame(
        {
            "gene_id": [
                "g1",
                "g2",
                "g3",
            ],
            "persistence_breadth": [
                2,
                1,
                0,
            ],
        }
    )


def test_gene_condition_response_preserves_all_classified_rows():
    """Canonical response table should retain one row per gene-condition pair."""

    result = (
        build_gene_condition_response(
            make_classified(),
            make_breadth(),
            response_floor=0.342263,
        )
    )

    assert len(
        result
    ) == 4

    assert not result.duplicated(
        subset=[
            "gene_id",
            "condition",
        ]
    ).any()

    assert (
        result[
            "response_floor"
        ]
        == 0.342263
    ).all()

    g1_drought = result[
        (
            result[
                "gene_id"
            ]
            == "g1"
        )
        & (
            result[
                "condition"
            ]
            == "drought"
        )
    ].iloc[0]

    assert g1_drought[
        "recovery_outcome"
    ] == (
        "persistent_same_direction"
    )

    assert g1_drought[
        "persistence_breadth"
    ] == 2


def test_gene_condition_response_rejects_duplicate_keys():
    """Duplicate gene-condition rows should fail validation."""

    classified = (
        make_classified()
    )

    classified = pd.concat(
        [
            classified,
            classified.iloc[
                [0]
            ],
        ],
        ignore_index=True,
    )

    with pytest.raises(
        AnalyticalTableError
    ):
        build_gene_condition_response(
            classified,
            make_breadth(),
            response_floor=0.342263,
        )


def test_gene_persistence_summary_includes_zero_breadth_genes():
    """Genes without persistent responses must remain explicitly present."""

    result = (
        build_gene_persistence_summary(
            make_classified(),
            make_breadth(),
        )
    )

    assert len(
        result
    ) == 3

    g1 = result[
        result[
            "gene_id"
        ]
        == "g1"
    ].iloc[0]

    assert g1[
        "n_persistent_conditions"
    ] == 2

    assert g1[
        "n_persistent_up"
    ] == 2

    assert g1[
        "n_persistent_down"
    ] == 0

    g2 = result[
        result[
            "gene_id"
        ]
        == "g2"
    ].iloc[0]

    assert g2[
        "n_persistent_conditions"
    ] == 1

    assert g2[
        "n_persistent_down"
    ] == 1

    g3 = result[
        result[
            "gene_id"
        ]
        == "g3"
    ].iloc[0]

    assert g3[
        "persistence_breadth"
    ] == 0

    assert g3[
        "n_persistent_conditions"
    ] == 0

    assert g3[
        "n_persistent_up"
    ] == 0

    assert g3[
        "n_persistent_down"
    ] == 0


def test_gene_persistence_summary_detects_inconsistent_breadth():
    """Breadth should agree with persistent-same-direction counts."""

    breadth = (
        make_breadth()
    )

    breadth.loc[
        breadth[
            "gene_id"
        ]
        == "g1",
        "persistence_breadth",
    ] = 5

    with pytest.raises(
        AnalyticalTableError
    ):
        build_gene_persistence_summary(
            make_classified(),
            breadth,
        )