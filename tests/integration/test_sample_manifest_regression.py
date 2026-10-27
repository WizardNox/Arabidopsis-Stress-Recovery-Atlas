"""Regression tests between Stage 1 and Stage 2 sample metadata."""

import pandas as pd
from pandas.testing import assert_frame_equal

from arabidopsis_stress_recovery_atlas.config import (
    PROJECT_ROOT,
    load_scientific_baseline,
)
from arabidopsis_stress_recovery_atlas.ingestion.load_expression_data import (
    load_expression_data,
)
from arabidopsis_stress_recovery_atlas.transformation.expression import (
    get_expression_columns,
    standardize_gene_id_column,
)
from arabidopsis_stress_recovery_atlas.transformation.metadata import (
    build_sample_manifest,
)


COMPARISON_COLUMNS = [
    "state",
    "condition",
    "replicate",
    "heat",
    "drought",
    "high_light",
    "cold",
    "high_humidity",
]


def test_reconstructed_manifest_matches_stage1_biological_design():
    """Stage 2 should reproduce the biological sample design used in Stage 1."""

    baseline = load_scientific_baseline()

    stage1_manifest_path = (
        PROJECT_ROOT
        / "data"
        / "interim"
        / "sample_manifest.csv"
    )

    stage1 = pd.read_csv(stage1_manifest_path)

    expression = load_expression_data()
    standardized = standardize_gene_id_column(
        expression
    )

    sample_columns = get_expression_columns(
        standardized
    )

    stage2 = build_sample_manifest(
        sample_columns
    )

    expected_samples = baseline["dataset"][
        "expression_samples"
    ]

    assert len(stage1) == expected_samples
    assert len(stage2) == expected_samples

    stage1_keys = stage1[
        ["state", "condition", "replicate"]
    ]

    stage2_keys = stage2[
        ["state", "condition", "replicate"]
    ]

    assert not stage1_keys.duplicated().any()
    assert not stage2_keys.duplicated().any()

    stage1_comparison = (
        stage1[COMPARISON_COLUMNS]
        .sort_values(
            ["state", "condition", "replicate"]
        )
        .reset_index(drop=True)
    )

    stage2_comparison = (
        stage2[COMPARISON_COLUMNS]
        .sort_values(
            ["state", "condition", "replicate"]
        )
        .reset_index(drop=True)
    )

    assert_frame_equal(
        stage2_comparison,
        stage1_comparison,
        check_dtype=False,
    )