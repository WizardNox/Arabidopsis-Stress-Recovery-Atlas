"""Validate the processed expression matrix."""

from typing import Any

import numpy as np
import pandas as pd
from pandas.api.types import is_numeric_dtype

from arabidopsis_stress_recovery_atlas.config import (
    load_scientific_baseline,
)


class ExpressionValidationError(ValueError):
    """Raised when the expression matrix violates project assumptions."""


GENE_ID_CANDIDATES = (
    "gene_id",
    "Unnamed: 0",
)


def identify_gene_id_column(data: pd.DataFrame) -> str:
    """Identify the gene identifier column in the source matrix."""

    matches = [
        column
        for column in GENE_ID_CANDIDATES
        if column in data.columns
    ]

    if not matches:
        raise ExpressionValidationError(
            "No recognized gene identifier column was found."
        )

    if len(matches) > 1:
        raise ExpressionValidationError(
            "Multiple possible gene identifier columns were found: "
            + ", ".join(matches)
        )

    return matches[0]


def validate_expression_matrix(
    data: pd.DataFrame,
) -> dict[str, Any]:
    """Validate the expression matrix against project expectations."""

    baseline = load_scientific_baseline()

    expected_rows = baseline["dataset"]["source_gene_rows"]
    expected_samples = baseline["dataset"]["expression_samples"]
    expected_informative = baseline["dataset"]["informative_genes"]
    expected_all_zero = baseline["dataset"]["all_zero_genes"]

    gene_id_column = identify_gene_id_column(data)

    expression_columns = [
        column
        for column in data.columns
        if column != gene_id_column
    ]

    if len(data) != expected_rows:
        raise ExpressionValidationError(
            "Unexpected number of gene rows. "
            f"Expected {expected_rows:,}, found {len(data):,}."
        )

    if len(expression_columns) != expected_samples:
        raise ExpressionValidationError(
            "Unexpected number of expression samples. "
            f"Expected {expected_samples}, "
            f"found {len(expression_columns)}."
        )

    if data[gene_id_column].isna().any():
        raise ExpressionValidationError(
            "Gene identifier column contains missing values."
        )

    if data[gene_id_column].duplicated().any():
        duplicate_count = int(
            data[gene_id_column].duplicated().sum()
        )

        raise ExpressionValidationError(
            "Gene identifier column contains "
            f"{duplicate_count:,} duplicate identifiers."
        )

    non_numeric_columns = [
        column
        for column in expression_columns
        if not is_numeric_dtype(data[column])
    ]

    if non_numeric_columns:
        raise ExpressionValidationError(
            "Non-numeric expression columns found: "
            + ", ".join(non_numeric_columns)
        )

    expression = data[expression_columns]

    missing_values = int(expression.isna().sum().sum())

    if missing_values:
        raise ExpressionValidationError(
            "Expression matrix contains "
            f"{missing_values:,} missing values."
        )

    expression_array = expression.to_numpy(dtype=float)

    if not np.isfinite(expression_array).all():
        raise ExpressionValidationError(
            "Expression matrix contains infinite values."
        )

    all_zero_mask = (expression == 0).all(axis=1)

    all_zero_genes = int(all_zero_mask.sum())
    informative_genes = int((~all_zero_mask).sum())

    if all_zero_genes != expected_all_zero:
        raise ExpressionValidationError(
            "Unexpected number of all-zero genes. "
            f"Expected {expected_all_zero:,}, "
            f"found {all_zero_genes:,}."
        )

    if informative_genes != expected_informative:
        raise ExpressionValidationError(
            "Unexpected number of informative genes. "
            f"Expected {expected_informative:,}, "
            f"found {informative_genes:,}."
        )

    return {
        "gene_id_column": gene_id_column,
        "gene_rows": len(data),
        "expression_samples": len(expression_columns),
        "informative_genes": informative_genes,
        "all_zero_genes": all_zero_genes,
        "missing_expression_values": missing_values,
    }