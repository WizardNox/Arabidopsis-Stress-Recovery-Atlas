"""Transform the expression matrix into canonical project representations."""

import pandas as pd

from arabidopsis_stress_recovery_atlas.validation.expression import (
    identify_gene_id_column,
)


CANONICAL_GENE_ID_COLUMN = "gene_id"


def standardize_gene_id_column(
    data: pd.DataFrame,
) -> pd.DataFrame:
    """Return a copy with the gene identifier column named 'gene_id'."""

    gene_id_column = identify_gene_id_column(data)

    standardized = data.copy()

    if gene_id_column != CANONICAL_GENE_ID_COLUMN:
        standardized = standardized.rename(
            columns={
                gene_id_column: CANONICAL_GENE_ID_COLUMN
            }
        )

    return standardized


def get_expression_columns(
    data: pd.DataFrame,
) -> list[str]:
    """Return expression-sample columns from a canonical expression table."""

    if CANONICAL_GENE_ID_COLUMN not in data.columns:
        raise ValueError(
            "Canonical gene identifier column 'gene_id' was not found."
        )

    return [
        column
        for column in data.columns
        if column != CANONICAL_GENE_ID_COLUMN
    ]


def get_informative_expression(
    data: pd.DataFrame,
) -> pd.DataFrame:
    """Return genes with at least one non-zero expression value."""

    standardized = standardize_gene_id_column(data)

    expression_columns = get_expression_columns(standardized)

    informative_mask = (
        standardized[expression_columns] != 0
    ).any(axis=1)

    informative = standardized.loc[
        informative_mask
    ].copy()

    return informative


def get_all_zero_expression(
    data: pd.DataFrame,
) -> pd.DataFrame:
    """Return genes whose expression is zero in every sample."""

    standardized = standardize_gene_id_column(data)

    expression_columns = get_expression_columns(standardized)

    all_zero_mask = (
        standardized[expression_columns] == 0
    ).all(axis=1)

    all_zero = standardized.loc[
        all_zero_mask
    ].copy()

    return all_zero