"""Tests for expression-matrix transformations."""

import pandas as pd

from arabidopsis_stress_recovery_atlas.transformation.expression import (
    get_all_zero_expression,
    get_expression_columns,
    get_informative_expression,
    standardize_gene_id_column,
)


def make_test_expression_data() -> pd.DataFrame:
    """Create a small expression table for transformation tests."""

    return pd.DataFrame(
        {
            "Unnamed: 0": [
                "AT1G00010",
                "AT1G00020",
                "AT1G00030",
            ],
            "sample_1": [0.0, 1.5, 0.0],
            "sample_2": [0.0, 0.0, 2.0],
            "sample_3": [0.0, 3.0, 0.0],
        }
    )


def test_gene_id_column_is_standardized():
    """The source gene column should be renamed to canonical gene_id."""

    data = make_test_expression_data()

    standardized = standardize_gene_id_column(data)

    assert "gene_id" in standardized.columns
    assert "Unnamed: 0" not in standardized.columns


def test_standardization_does_not_modify_input():
    """Standardization should leave the original DataFrame unchanged."""

    data = make_test_expression_data()

    standardize_gene_id_column(data)

    assert "Unnamed: 0" in data.columns
    assert "gene_id" not in data.columns


def test_expression_columns_exclude_gene_id():
    """Only sample-expression columns should be returned."""

    data = make_test_expression_data()

    standardized = standardize_gene_id_column(data)

    expression_columns = get_expression_columns(
        standardized
    )

    assert expression_columns == [
        "sample_1",
        "sample_2",
        "sample_3",
    ]


def test_informative_expression_removes_all_zero_genes():
    """Genes that are zero in every sample should be excluded."""

    data = make_test_expression_data()

    informative = get_informative_expression(data)

    assert list(informative["gene_id"]) == [
        "AT1G00020",
        "AT1G00030",
    ]


def test_all_zero_expression_contains_only_zero_genes():
    """The all-zero subset should contain only invariant zero rows."""

    data = make_test_expression_data()

    all_zero = get_all_zero_expression(data)

    assert list(all_zero["gene_id"]) == [
        "AT1G00010"
    ]