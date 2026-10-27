"""Build canonical reference and interpretation tables."""

import numpy as np
import pandas as pd


class ReferenceTableError(ValueError):
    """Raised when a canonical reference table cannot be constructed."""


def _require_columns(
    table: pd.DataFrame,
    required_columns: set[str],
    table_name: str,
) -> None:
    """Require a table to contain specific columns."""

    missing = (
        required_columns
        - set(table.columns)
    )

    if missing:
        raise ReferenceTableError(
            f"{table_name} is missing required columns: "
            + ", ".join(
                sorted(missing)
            )
        )


def build_gene_dimension(
    informative_expression: pd.DataFrame,
    araport_loci: pd.DataFrame,
) -> pd.DataFrame:
    """Build one canonical row for every informative gene."""

    _require_columns(
        informative_expression,
        {"gene_id"},
        "informative_expression",
    )

    _require_columns(
        araport_loci,
        {"gene_id"},
        "araport_loci",
    )

    genes = (
        informative_expression[
            ["gene_id"]
        ]
        .copy()
    )

    genes["gene_id"] = (
        genes["gene_id"]
        .astype(str)
        .str.strip()
    )

    annotation = (
        araport_loci.copy()
    )

    annotation["gene_id"] = (
        annotation["gene_id"]
        .astype(str)
        .str.strip()
    )

    if genes[
        "gene_id"
    ].duplicated().any():
        raise ReferenceTableError(
            "Informative expression contains duplicate gene IDs."
        )

    if annotation[
        "gene_id"
    ].duplicated().any():
        raise ReferenceTableError(
            "Araport locus summary contains duplicate gene IDs."
        )

    annotated_ids = set(
        annotation[
            "gene_id"
        ]
    )

    result = genes.merge(
        annotation,
        on="gene_id",
        how="left",
        validate="one_to_one",
    )

    result.insert(
        1,
        "has_araport_annotation",
        result[
            "gene_id"
        ].isin(
            annotated_ids
        ),
    )

    return (
        result
        .sort_values(
            "gene_id"
        )
        .reset_index(
            drop=True
        )
    )


def build_go_term_enrichment_table(
    go_slim_enrichment: pd.DataFrame,
    bp_enrichment: pd.DataFrame,
    alpha: float = 0.05,
) -> pd.DataFrame:
    """Standardize GO-Slim and detailed-BP enrichment into one table.

    The two analyses retain their different multiple-testing
    corrections rather than being treated as statistically identical.
    """

    go_slim_required = {
        "go_aspect",
        "go_slim_term",
        "foreground_with_term",
        "foreground_total",
        "background_with_term",
        "background_total",
        "foreground_fraction",
        "background_fraction",
        "enrichment_ratio",
        "odds_ratio",
        "p_value",
        "fdr_global",
        "fdr_within_aspect",
    }

    bp_required = {
        "go_id",
        "go_term",
        "foreground_with_term",
        "foreground_total",
        "background_with_term",
        "background_total",
        "foreground_fraction",
        "background_fraction",
        "enrichment_ratio",
        "odds_ratio",
        "p_value",
        "fdr_bh",
    }

    _require_columns(
        go_slim_enrichment,
        go_slim_required,
        "go_slim_enrichment",
    )

    _require_columns(
        bp_enrichment,
        bp_required,
        "bp_enrichment",
    )

    common_numeric_columns = [
        "foreground_with_term",
        "foreground_total",
        "background_with_term",
        "background_total",
        "foreground_fraction",
        "background_fraction",
        "enrichment_ratio",
        "odds_ratio",
        "p_value",
    ]

    go_slim = pd.DataFrame(
        {
            "analysis_type": "go_slim",
            "go_aspect": (
                go_slim_enrichment[
                    "go_aspect"
                ]
            ),
            "go_id": pd.NA,
            "term": (
                go_slim_enrichment[
                    "go_slim_term"
                ]
            ),
        }
    )

    for column in common_numeric_columns:
        go_slim[
            column
        ] = go_slim_enrichment[
            column
        ].values

    go_slim[
        "fdr_primary"
    ] = (
        go_slim_enrichment[
            "fdr_global"
        ].values
    )

    go_slim[
        "fdr_secondary"
    ] = (
        go_slim_enrichment[
            "fdr_within_aspect"
        ].values
    )

    go_slim[
        "is_significant"
    ] = (
        go_slim[
            "fdr_primary"
        ]
        < alpha
    )

    bp = pd.DataFrame(
        {
            "analysis_type": (
                "detailed_bp"
            ),
            "go_aspect": "P",
            "go_id": (
                bp_enrichment[
                    "go_id"
                ]
            ),
            "term": (
                bp_enrichment[
                    "go_term"
                ]
            ),
        }
    )

    for column in common_numeric_columns:
        bp[
            column
        ] = bp_enrichment[
            column
        ].values

    bp[
        "fdr_primary"
    ] = (
        bp_enrichment[
            "fdr_bh"
        ].values
    )

    bp[
        "fdr_secondary"
    ] = np.nan

    bp[
        "is_significant"
    ] = (
        bp[
            "fdr_primary"
        ]
        < alpha
    )

    result = pd.concat(
        [
            go_slim,
            bp,
        ],
        ignore_index=True,
    )

    column_order = [
        "analysis_type",
        "go_aspect",
        "go_id",
        "term",
        "foreground_with_term",
        "foreground_total",
        "background_with_term",
        "background_total",
        "foreground_fraction",
        "background_fraction",
        "enrichment_ratio",
        "odds_ratio",
        "p_value",
        "fdr_primary",
        "fdr_secondary",
        "is_significant",
    ]

    return (
        result[
            column_order
        ]
        .sort_values(
            [
                "analysis_type",
                "fdr_primary",
                "term",
            ],
            na_position="last",
        )
        .reset_index(
            drop=True
        )
    )


def build_published_de_concordance_table(
    direction_validation: pd.DataFrame,
) -> pd.DataFrame:
    """Build canonical shared gene-condition published-DE comparison."""

    required_columns = {
        "gene_id",
        "condition",
        "direction",
        "endpoint_shift",
        "recovery_shift",
        "persistence_breadth",
        "published_log2fc",
        "published_adj_p",
        "published_direction",
        "published_log2fc_treatment_vs_control",
        "published_direction_treatment_vs_control",
        "direction_agrees",
    }

    _require_columns(
        direction_validation,
        required_columns,
        "direction_validation",
    )

    result = (
        direction_validation[
            [
                "gene_id",
                "condition",
                "direction",
                "endpoint_shift",
                "recovery_shift",
                "persistence_breadth",
                "published_log2fc",
                "published_adj_p",
                "published_direction",
                "published_log2fc_treatment_vs_control",
                "published_direction_treatment_vs_control",
                "direction_agrees",
            ]
        ]
        .copy()
        .rename(
            columns={
                "direction": (
                    "project_direction"
                ),
                "published_log2fc": (
                    "published_log2fc_original"
                ),
                "published_direction": (
                    "published_direction_original"
                ),
            }
        )
    )

    result[
        "gene_id"
    ] = (
        result[
            "gene_id"
        ]
        .astype(str)
        .str.strip()
    )

    if result.duplicated(
        subset=[
            "gene_id",
            "condition",
        ]
    ).any():
        raise ReferenceTableError(
            "Published-DE concordance contains duplicate "
            "gene-condition rows."
        )

    return (
        result
        .sort_values(
            [
                "condition",
                "gene_id",
            ]
        )
        .reset_index(
            drop=True
        )
    )