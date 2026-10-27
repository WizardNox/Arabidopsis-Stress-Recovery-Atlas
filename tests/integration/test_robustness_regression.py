"""Regression tests for response-floor robustness."""

import pytest

from arabidopsis_stress_recovery_atlas.analysis.exploratory import (
    calculate_condition_shifts,
    calculate_control_pairwise_differences,
    calculate_response_floor,
    classify_recovery_outcomes,
)
from arabidopsis_stress_recovery_atlas.analysis.response_programs import (
    summarize_persistence_pattern,
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


def test_persistence_patterns_are_robust_across_response_floors():
    """Stage 2 should reproduce Stage 1 robustness results."""

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

    shifts = calculate_condition_shifts(
        informative,
        manifest,
    )

    quantiles = project_config[
        "robustness"
    ]["response_floor_quantiles"]

    high_breadth_minimum = project_config[
        "analysis"
    ]["persistence"][
        "high_breadth_min_conditions"
    ]

    for quantile in quantiles:
        percentile = (
            f"p{int(round(quantile * 100))}"
        )

        response_floor = (
            calculate_response_floor(
                control_differences,
                quantile,
            )
        )

        classified = classify_recovery_outcomes(
            shifts,
            response_floor,
        )

        summary = summarize_persistence_pattern(
            classified,
            high_breadth_minimum=(
                high_breadth_minimum
            ),
            multi_stress_minimum=2,
        )

        assert (
            summary["high_breadth_genes"]
            == baseline["robustness"][
                "high_breadth_genes"
            ][percentile]
        )

        assert (
            summary[
                "multi_stress_persistent_genes"
            ]
            == baseline["robustness"][
                "multi_stress_persistent_genes"
            ][percentile]
        )

        assert (
            summary[
                "directional_consistency_percent"
            ]
            == pytest.approx(
                baseline["robustness"][
                    "directional_consistency_percent"
                ][percentile],
                abs=0.01,
            )
        )

        assert (
            summary[
                "consistently_above_percent"
            ]
            == pytest.approx(
                baseline["robustness"][
                    "consistently_above_control_percent"
                ][percentile],
                abs=0.01,
            )
        )

        if percentile == "p95":
            assert (
                summary[
                    "directionally_consistent_genes"
                ]
                == baseline[
                    "directional_consistency"
                ][
                    "directionally_consistent_genes"
                ]
            )

            assert (
                 summary[
                    "consistently_above_control"
                ]
                == baseline[
                    "directional_consistency"
                ][
                    "consistently_above_control"
                ]
            )

            assert (
                summary[
                    "consistently_below_control"
                ]
                == baseline[
                    "directional_consistency"
                ][
                    "consistently_below_control"
                ]
            )    