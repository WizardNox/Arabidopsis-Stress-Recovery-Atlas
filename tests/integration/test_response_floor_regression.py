"""Regression tests for the control-derived response floor."""

from math import comb

import pytest

from arabidopsis_stress_recovery_atlas.analysis.exploratory import (
    calculate_control_pairwise_differences,
    calculate_response_floors,
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


def test_control_response_floors_match_stage1_results():
    """Stage 2 should reproduce the Stage 1 control-derived floors."""

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

    differences = (
        calculate_control_pairwise_differences(
            informative,
            manifest,
        )
    )

    control_samples = baseline[
        "experimental_design"
    ]["control_samples"]

    expected_difference_count = (
        baseline["dataset"]["informative_genes"]
        * comb(control_samples, 2)
    )

    assert (
        len(differences)
        == expected_difference_count
    )

    quantiles = project_config[
        "robustness"
    ]["response_floor_quantiles"]

    floors = calculate_response_floors(
        differences,
        quantiles,
    )

    expected_floors = baseline[
        "robustness"
    ]["response_floors"]

    for percentile, expected_value in (
        expected_floors.items()
    ):
        assert floors[percentile] == pytest.approx(
            expected_value,
            abs=1e-6,
        )

    assert floors["p95"] == pytest.approx(
        baseline["primary_response_floor"][
            "rlog_units"
        ],
        abs=1e-6,
    )