"""Build canonical analytical tables from tested scientific outputs."""

from collections.abc import Sequence

import pandas as pd


class AnalyticalTableError(ValueError):
    """Raised when a canonical analytical table cannot be constructed."""


def _require_columns(
    table: pd.DataFrame,
    required_columns: set[str],
    table_name: str,
) -> None:
    """Require a table to contain a defined set of columns."""

    missing_columns = (
        required_columns
        - set(table.columns)
    )

    if missing_columns:
        raise AnalyticalTableError(
            f"{table_name} is missing required columns: "
            + ", ".join(
                sorted(
                    missing_columns
                )
            )
        )


def _require_unique_key(
    table: pd.DataFrame,
    key_columns: Sequence[str],
    table_name: str,
) -> None:
    """Require one row per declared table key."""

    duplicate_mask = table.duplicated(
        subset=list(
            key_columns
        ),
        keep=False,
    )

    if duplicate_mask.any():
        raise AnalyticalTableError(
            f"{table_name} contains duplicate rows for key: "
            + ", ".join(
                key_columns
            )
        )


def _clean_gene_ids(
    table: pd.DataFrame,
) -> pd.DataFrame:
    """Return a copy with standardized gene-ID strings."""

    result = table.copy()

    result[
        "gene_id"
    ] = (
        result[
            "gene_id"
        ]
        .astype(str)
        .str.strip()
    )

    return result


def build_gene_condition_response(
    classified: pd.DataFrame,
    persistence_breadth: pd.DataFrame,
    response_floor: float,
) -> pd.DataFrame:
    """Build one canonical row per gene and stress condition.

    The scientific values and recovery classifications must already
    have been calculated by the tested analysis layer.

    This function only validates, joins, orders, and standardizes them.
    """

    classified_required = {
        "gene_id",
        "condition",
        "control_mean",
        "endpoint_mean",
        "recovery_mean",
        "endpoint_shift",
        "recovery_shift",
        "endpoint_eligible",
        "recovery_above_floor",
        "same_direction",
        "recovery_outcome",
    }

    breadth_required = {
        "gene_id",
        "persistence_breadth",
    }

    _require_columns(
        classified,
        classified_required,
        "classified",
    )

    _require_columns(
        persistence_breadth,
        breadth_required,
        "persistence_breadth",
    )

    classified_clean = (
        _clean_gene_ids(
            classified
        )
    )

    breadth_clean = (
        _clean_gene_ids(
            persistence_breadth
        )
    )

    _require_unique_key(
        classified_clean,
        [
            "gene_id",
            "condition",
        ],
        "classified",
    )

    _require_unique_key(
        breadth_clean,
        [
            "gene_id",
        ],
        "persistence_breadth",
    )

    result = (
        classified_clean
        .merge(
            breadth_clean[
                [
                    "gene_id",
                    "persistence_breadth",
                ]
            ],
            on="gene_id",
            how="left",
            validate="many_to_one",
        )
    )

    if result[
        "persistence_breadth"
    ].isna().any():
        raise AnalyticalTableError(
            "Some gene-condition rows are missing persistence breadth."
        )

    result[
        "response_floor"
    ] = float(
        response_floor
    )

    column_order = [
        "gene_id",
        "condition",
        "control_mean",
        "endpoint_mean",
        "recovery_mean",
        "endpoint_shift",
        "recovery_shift",
        "response_floor",
        "endpoint_eligible",
        "recovery_above_floor",
        "same_direction",
        "recovery_outcome",
        "persistence_breadth",
    ]

    result = result[
        column_order
    ]

    return (
        result
        .sort_values(
            [
                "gene_id",
                "condition",
            ]
        )
        .reset_index(
            drop=True
        )
    )


def build_gene_persistence_summary(
    classified: pd.DataFrame,
    persistence_breadth: pd.DataFrame,
) -> pd.DataFrame:
    """Build one persistence-summary row for every informative gene.

    Genes with no persistent recovery outcome remain in the table
    with zero persistent conditions.

    Direction counts are descriptive summaries of already-classified
    persistent-same-direction rows; they do not perform a new
    biological classification.
    """

    classified_required = {
        "gene_id",
        "condition",
        "recovery_shift",
        "recovery_outcome",
    }

    breadth_required = {
        "gene_id",
        "persistence_breadth",
    }

    _require_columns(
        classified,
        classified_required,
        "classified",
    )

    _require_columns(
        persistence_breadth,
        breadth_required,
        "persistence_breadth",
    )

    classified_clean = (
        _clean_gene_ids(
            classified
        )
    )

    breadth_clean = (
        _clean_gene_ids(
            persistence_breadth
        )
    )

    _require_unique_key(
        classified_clean,
        [
            "gene_id",
            "condition",
        ],
        "classified",
    )

    _require_unique_key(
        breadth_clean,
        [
            "gene_id",
        ],
        "persistence_breadth",
    )

    persistent = (
        classified_clean[
            classified_clean[
                "recovery_outcome"
            ]
            == "persistent_same_direction"
        ]
        .copy()
    )

    if persistent.empty:
        persistent_counts = pd.DataFrame(
            columns=[
                "gene_id",
                "n_persistent_conditions",
                "n_persistent_up",
                "n_persistent_down",
            ]
        )

    else:
        persistent[
            "persistent_up"
        ] = (
            persistent[
                "recovery_shift"
            ]
            > 0
        )

        persistent[
            "persistent_down"
        ] = (
            persistent[
                "recovery_shift"
            ]
            < 0
        )

        persistent_counts = (
            persistent
            .groupby(
                "gene_id",
                as_index=False,
            )
            .agg(
                n_persistent_conditions=(
                    "condition",
                    "nunique",
                ),
                n_persistent_up=(
                    "persistent_up",
                    "sum",
                ),
                n_persistent_down=(
                    "persistent_down",
                    "sum",
                ),
            )
        )

    result = (
        breadth_clean
        .merge(
            persistent_counts,
            on="gene_id",
            how="left",
            validate="one_to_one",
        )
    )

    count_columns = [
        "n_persistent_conditions",
        "n_persistent_up",
        "n_persistent_down",
    ]

    result[
        count_columns
    ] = (
        result[
            count_columns
        ]
        .fillna(0)
        .astype(int)
    )

    result[
        "persistence_breadth"
    ] = (
        result[
            "persistence_breadth"
        ]
        .astype(int)
    )

    # Persistence breadth and the number of persistent
    # stress conditions are two representations of the
    # same upstream scientific result.
    inconsistent = (
        result[
            "persistence_breadth"
        ]
        != result[
            "n_persistent_conditions"
        ]
    )

    if inconsistent.any():
        raise AnalyticalTableError(
            "Persistence breadth is inconsistent with "
            "persistent-same-direction condition counts."
        )

    column_order = [
        "gene_id",
        "persistence_breadth",
        "n_persistent_conditions",
        "n_persistent_up",
        "n_persistent_down",
    ]

    return (
        result[
            column_order
        ]
        .sort_values(
            [
                "persistence_breadth",
                "gene_id",
            ],
            ascending=[
                False,
                True,
            ],
        )
        .reset_index(
            drop=True
        )
    )