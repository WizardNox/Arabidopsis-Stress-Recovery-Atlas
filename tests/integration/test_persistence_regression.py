"""Regression tests for the primary p95 persistence analysis."""

import pytest

from arabidopsis_stress_recovery_atlas.analysis.exploratory import (
    calculate_condition_shifts,
    calculate_control_pairwise_differences,
    calculate_persistence_breadth,
    calculate_response_floor,
    classify_recovery_outcomes,
    summarize_persistent_counts,
)
from arabidopsis_stress_recovery_atlas.config import (
    load_project_config,
    load_scientific_baseline,
)
from arabidopsis_stress_recovery_atlas.ingestion.load_expression_data import (
    load_expression_data,
)
from arabidopsis_stress_recovery_atlas.transformation.expression import (
    get_expression_columns,
    get_informative_expression,
    standardize_gene_id_column,
)
from arabidopsis_stress_recovery_atlas.transformation.metadata import (
    build_sample_manifest,
)


def test_primary_p95_persistence_matches_stage1_results():
    """Stage 2 should reproduce primary Stage 1 persistence results."""

    project_config = load_project_config()
    baseline = load_scientific_baseline()

    expression = load_expression_data()

    standardized = standardize_gene_id_column(
        expression
    )

    informative = get_informative_expression(
        standardized
    )

    sample_columns = get_expression_columns(
        standardized
    )

    manifest = build_sample_manifest(
        sample_columns
    )

    control_differences = (
        calculate_control_pairwise_differences(
            informative,
            manifest,
        )
    )

    primary_quantile = project_config[
        "analysis"
    ]["response_floor"]["primary_quantile"]

    response_floor = calculate_response_floor(
        control_differences,
        primary_quantile,
    )

    assert response_floor == pytest.approx(
        baseline["primary_response_floor"][
            "rlog_units"
        ],
        abs=1e-6,
    )

    shifts = calculate_condition_shifts(
        informative,
        manifest,
    )

    classified = classify_recovery_outcomes(
        shifts,
        response_floor,
    )

    persistent_counts = summarize_persistent_counts(
        classified
    )

    assert (
        persistent_counts
        == baseline[
            "primary_p95_persistence_by_condition"
        ]
    )

    breadth = calculate_persistence_breadth(
        classified
    )

    cumulative_expectations = {
        1: baseline["persistent_gene_counts"][
            "at_least_1_stress"
        ],
        2: baseline["persistent_gene_counts"][
            "at_least_2_stresses"
        ],
        3: baseline["persistent_gene_counts"][
            "at_least_3_stresses"
        ],
        4: baseline["persistent_gene_counts"][
            "at_least_4_stresses"
        ],
        5: baseline["persistent_gene_counts"][
            "at_least_5_stresses"
        ],
    }

    for minimum_conditions, expected_count in (
        cumulative_expectations.items()
    ):
        observed_count = int(
            (
                breadth["persistence_breadth"]
                >= minimum_conditions
            ).sum()
        )

        assert observed_count == expected_count

    high_breadth_minimum = project_config[
        "analysis"
    ]["persistence"][
        "high_breadth_min_conditions"
    ]

    high_breadth_count = int(
        (
            breadth["persistence_breadth"]
            >= high_breadth_minimum
        ).sum()
    )

    assert (
        high_breadth_count
        == baseline["persistent_gene_counts"][
            "at_least_4_stresses"
        ]
    )