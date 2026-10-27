"""Regression tests for persistence breadth and threshold robustness."""

import pytest

from arabidopsis_stress_recovery_atlas.analysis.exploratory import (
    calculate_condition_shifts,
    calculate_control_pairwise_differences,
    calculate_response_floor,
    classify_recovery_outcomes,
)
from arabidopsis_stress_recovery_atlas.analysis.response_programs import (
    calculate_persistence_rank_correlation,
    get_persistent_gene_condition_pairs,
    summarize_breadth_distribution,
    summarize_condition_persistence_fractions,
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


def test_persistence_breadth_and_threshold_stability_match_stage1():
    """Stage 2 should reproduce Stage 1 breadth and robustness results."""

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

    classified_by_percentile = {}
    fractions_by_percentile = {}
    pairs_by_percentile = {}

    for quantile in quantiles:
        percentile = (
            f"p{int(round(quantile * 100))}"
        )

        floor = calculate_response_floor(
            control_differences,
            quantile,
        )

        classified = classify_recovery_outcomes(
            shifts,
            floor,
        )

        classified_by_percentile[
            percentile
        ] = classified

        fractions_by_percentile[
            percentile
        ] = summarize_condition_persistence_fractions(
            classified
        )

        pairs_by_percentile[
            percentile
        ] = get_persistent_gene_condition_pairs(
            classified
        )

    p95_distribution = summarize_breadth_distribution(
        classified_by_percentile["p95"],
        maximum_conditions=baseline[
            "experimental_design"
        ]["stress_conditions"],
    )

    expected_distribution = {
        int(breadth): count
        for breadth, count in baseline[
            "persistence_breadth_distribution"
        ].items()
    }

    assert (
        p95_distribution
        == expected_distribution
    )

    assert (
        pairs_by_percentile["p99"]
        <= pairs_by_percentile["p95"]
    )

    assert (
        pairs_by_percentile["p95"]
        <= pairs_by_percentile["p90"]
    )

    rank_expectations = baseline[
        "robustness"
    ]["persistence_rank_spearman"]

    p90_p95 = (
        calculate_persistence_rank_correlation(
            fractions_by_percentile["p90"],
            fractions_by_percentile["p95"],
        )
    )

    p95_p99 = (
        calculate_persistence_rank_correlation(
            fractions_by_percentile["p95"],
            fractions_by_percentile["p99"],
        )
    )

    p90_p99 = (
        calculate_persistence_rank_correlation(
            fractions_by_percentile["p90"],
            fractions_by_percentile["p99"],
        )
    )

    assert p90_p95 == pytest.approx(
        rank_expectations["p90_vs_p95"],
        abs=0.001,
    )

    assert p95_p99 == pytest.approx(
        rank_expectations["p95_vs_p99"],
        abs=0.001,
    )

    assert p90_p99 == pytest.approx(
        rank_expectations["p90_vs_p99"],
        abs=0.001,
    )