"""Regression tests for canonical analytical tables."""

from pathlib import Path

import pandas as pd
import pytest
from pandas.testing import (
    assert_frame_equal,
)

from arabidopsis_stress_recovery_atlas.ingestion.load_expression_data import (
    load_expression_data,
)
from arabidopsis_stress_recovery_atlas.transformation.canonical import (
    build_primary_analytical_tables,
    write_primary_analytical_tables,
)


@pytest.fixture(scope="module")
def canonical_tables():
    """Build canonical tables once for this integration-test module."""

    expression = (
        load_expression_data()
    )

    return (
        build_primary_analytical_tables(
            expression,
            response_floor_quantile=0.95,
        )
    )


def test_primary_analytical_tables_match_stage1(
    canonical_tables,
):
    """Canonical tables should reproduce the established Stage 1 results."""

    sample = canonical_tables[
        "sample"
    ]

    response = canonical_tables[
        "gene_condition_response"
    ]

    gene_summary = canonical_tables[
        "gene_persistence_summary"
    ]

    response_floor = canonical_tables[
        "response_floor"
    ]

    # Experimental design.
    assert len(
        sample
    ) == 57

    assert (
        sample[
            "state"
        ]
        == "control"
    ).sum() == 3

    assert (
        sample[
            "state"
        ]
        == "endpoint"
    ).sum() == 27

    assert (
        sample[
            "state"
        ]
        == "recovery"
    ).sum() == 27

    assert sample.loc[
        sample[
            "state"
        ]
        != "control",
        "condition",
    ].nunique() == 9

    # Primary response floor.
    assert response_floor == pytest.approx(
        0.34226272105,
        abs=1e-12,
    )

    # Complete gene-condition table.
    assert len(
        response
    ) == 235062

    assert response[
        "gene_id"
    ].nunique() == 26118

    assert response[
        "condition"
    ].nunique() == 9

    assert not response.duplicated(
        subset=[
            "gene_id",
            "condition",
        ]
    ).any()

    assert (
        response[
            "response_floor"
        ]
        == response_floor
    ).all()

    # Complete gene-level table.
    assert len(
        gene_summary
    ) == 26118

    assert gene_summary[
        "gene_id"
    ].nunique() == 26118

    assert not gene_summary[
        "gene_id"
    ].duplicated().any()

    # Compare exact Stage 1 persistent
    # gene-condition membership.
    stage1_persistent = pd.read_csv(
        Path(
            "data/interim/exploratory_analysis"
        )
        / "persistent_gene_table.csv"
    )

    stage2_persistent_pairs = {
        (
            row.gene_id,
            row.condition,
        )
        for row in response.loc[
            response[
                "recovery_outcome"
            ]
            == "persistent_same_direction",
            [
                "gene_id",
                "condition",
            ],
        ].itertuples(
            index=False
        )
    }

    stage1_persistent_pairs = {
        (
            row.gene_id,
            row.condition,
        )
        for row in stage1_persistent[
            [
                "gene_id",
                "condition",
            ]
        ].itertuples(
            index=False
        )
    }

    assert (
        stage2_persistent_pairs
        == stage1_persistent_pairs
    )

    # Compare the exact Stage 1 high-breadth
    # persistence summary.
    stage1_high_breadth = pd.read_csv(
        Path(
            "data/interim/exploratory_analysis"
        )
        / "high_breadth_persistent_genes.csv"
    )

    stage2_high_breadth = (
        gene_summary[
            gene_summary[
                "persistence_breadth"
            ]
            >= 4
        ][
            [
                "gene_id",
                "n_persistent_conditions",
                "n_persistent_up",
                "n_persistent_down",
            ]
        ]
        .sort_values(
            "gene_id"
        )
        .reset_index(
            drop=True
        )
    )

    stage1_high_breadth = (
        stage1_high_breadth[
            [
                "gene_id",
                "n_persistent_conditions",
                "n_persistent_up",
                "n_persistent_down",
            ]
        ]
        .sort_values(
            "gene_id"
        )
        .reset_index(
            drop=True
        )
    )

    assert len(
        stage2_high_breadth
    ) == 206

    assert_frame_equal(
        stage2_high_breadth,
        stage1_high_breadth,
        check_dtype=False,
    )


def test_primary_analytical_tables_can_be_written(
    canonical_tables,
    tmp_path,
):
    """Canonical tables should write cleanly to the processed layer."""

    output_directory = (
        tmp_path
        / "processed"
    )

    paths = (
        write_primary_analytical_tables(
            canonical_tables,
            output_directory,
        )
    )

    assert set(
        paths
    ) == {
        "sample",
        "gene_condition_response",
        "gene_persistence_summary",
    }

    for path in paths.values():
        assert path.is_file()

    response_written = pd.read_csv(
        paths[
            "gene_condition_response"
        ]
    )

    gene_written = pd.read_csv(
        paths[
            "gene_persistence_summary"
        ]
    )

    assert len(
        response_written
    ) == 235062

    assert len(
        gene_written
    ) == 26118