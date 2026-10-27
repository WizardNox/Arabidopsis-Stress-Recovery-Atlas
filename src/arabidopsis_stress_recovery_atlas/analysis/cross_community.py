"""Analyze genes represented across multiple significant BP communities."""

from collections.abc import Iterable

import pandas as pd


class CrossCommunityError(ValueError):
    """Raised when cross-community gene analysis cannot be performed."""


RECOVERY_CATEGORIES = (
    "no_qualifying_endpoint_response",
    "persistent_same_direction",
    "recovered",
    "substantial_reversal",
)


def build_community_gene_sets(
    significant_membership: pd.DataFrame,
    communities: list[set[str]],
) -> dict[int, set[str]]:
    """Collect foreground genes represented in each BP community."""

    required_columns = {
        "gene_id",
        "go_id",
    }

    missing_columns = (
        required_columns
        - set(significant_membership.columns)
    )

    if missing_columns:
        raise CrossCommunityError(
            "Significant BP membership is missing required columns: "
            + ", ".join(
                sorted(
                    missing_columns
                )
            )
        )

    community_gene_sets = {}

    for community_id, community in enumerate(
        communities,
        start=1,
    ):
        genes = set(
            significant_membership.loc[
                significant_membership[
                    "go_id"
                ].isin(
                    community
                ),
                "gene_id",
            ]
        )

        community_gene_sets[
            community_id
        ] = genes

    return community_gene_sets


def build_gene_community_membership(
    community_gene_sets: dict[
        int,
        set[str],
    ],
) -> pd.DataFrame:
    """Create one row for every gene-community membership."""

    records = []

    for community_id, genes in (
        community_gene_sets.items()
    ):
        for gene_id in genes:
            records.append(
                {
                    "gene_id": gene_id,
                    "community_id": (
                        community_id
                    ),
                }
            )

    if not records:
        return pd.DataFrame(
            columns=[
                "gene_id",
                "community_id",
            ]
        )

    return (
        pd.DataFrame(records)
        .drop_duplicates()
        .sort_values(
            [
                "gene_id",
                "community_id",
            ]
        )
        .reset_index(drop=True)
    )


def summarize_gene_community_breadth(
    gene_community_membership: pd.DataFrame,
) -> pd.DataFrame:
    """Count the number of BP communities represented by each gene."""

    required_columns = {
        "gene_id",
        "community_id",
    }

    missing_columns = (
        required_columns
        - set(
            gene_community_membership.columns
        )
    )

    if missing_columns:
        raise CrossCommunityError(
            "Gene-community membership is missing required columns: "
            + ", ".join(
                sorted(
                    missing_columns
                )
            )
        )

    breadth = (
        gene_community_membership
        .groupby(
            "gene_id"
        )["community_id"]
        .nunique()
        .rename(
            "n_communities"
        )
        .reset_index()
    )

    return breadth


def build_gene_community_profile(
    gene_community_membership: pd.DataFrame,
) -> pd.DataFrame:
    """Summarize community count and membership IDs for each gene."""

    required_columns = {
        "gene_id",
        "community_id",
    }

    missing_columns = (
        required_columns
        - set(
            gene_community_membership.columns
        )
    )

    if missing_columns:
        raise CrossCommunityError(
            "Gene-community membership is missing required columns: "
            + ", ".join(
                sorted(
                    missing_columns
                )
            )
        )

    profile = (
        gene_community_membership
        .groupby(
            "gene_id"
        )
        .agg(
            n_communities=(
                "community_id",
                "nunique",
            ),
            communities=(
                "community_id",
                lambda values: ",".join(
                    map(
                        str,
                        sorted(
                            set(
                                values
                            )
                        ),
                    )
                ),
            ),
        )
        .reset_index()
    )

    return profile


def select_cross_community_genes(
    gene_community_profile: pd.DataFrame,
    minimum_communities: int = 3,
) -> pd.DataFrame:
    """Select genes represented in at least the requested communities."""

    if minimum_communities < 1:
        raise CrossCommunityError(
            "Minimum community count must be at least 1."
        )

    required_columns = {
        "gene_id",
        "n_communities",
        "communities",
    }

    missing_columns = (
        required_columns
        - set(
            gene_community_profile.columns
        )
    )

    if missing_columns:
        raise CrossCommunityError(
            "Gene-community profile is missing required columns: "
            + ", ".join(
                sorted(
                    missing_columns
                )
            )
        )

    return (
        gene_community_profile[
            gene_community_profile[
                "n_communities"
            ]
            >= minimum_communities
        ]
        .copy()
        .sort_values(
            [
                "n_communities",
                "gene_id",
            ],
            ascending=[
                False,
                True,
            ],
        )
        .reset_index(drop=True)
    )


def build_cross_community_recovery_matrix(
    recovery_categories: pd.DataFrame,
    selected_gene_ids: Iterable[str],
    stress_conditions: Iterable[str],
) -> pd.DataFrame:
    """Extract recovery-state profiles for selected cross-community genes."""

    matrix = recovery_categories.copy()

    if "gene_id" in matrix.columns:
        matrix = matrix.set_index(
            "gene_id"
        )

    matrix.index = (
        matrix.index.astype(str)
    )

    selected_gene_ids = [
        str(gene_id)
        for gene_id
        in selected_gene_ids
    ]

    stress_conditions = list(
        stress_conditions
    )

    missing_genes = (
        set(selected_gene_ids)
        - set(matrix.index)
    )

    if missing_genes:
        raise CrossCommunityError(
            "Recovery matrix is missing selected genes: "
            + ", ".join(
                sorted(
                    missing_genes
                )
            )
        )

    missing_conditions = (
        set(stress_conditions)
        - set(matrix.columns)
    )

    if missing_conditions:
        raise CrossCommunityError(
            "Recovery matrix is missing stress conditions: "
            + ", ".join(
                sorted(
                    missing_conditions
                )
            )
        )

    selected = (
        matrix.loc[
            selected_gene_ids,
            stress_conditions,
        ]
        .copy()
        .fillna(
            "no_qualifying_endpoint_response"
        )
    )

    selected.index.name = "gene_id"

    return selected


def summarize_cross_community_recovery(
    recovery_matrix: pd.DataFrame,
) -> pd.DataFrame:
    """Summarize recovery states across selected genes for each stress."""

    if recovery_matrix.empty:
        raise CrossCommunityError(
            "Recovery matrix is empty."
        )

    total_genes = len(
        recovery_matrix
    )

    records = []

    for condition in recovery_matrix.columns:
        counts = (
            recovery_matrix[
                condition
            ]
            .value_counts()
        )

        record = {
            "condition": condition,
        }

        for category in RECOVERY_CATEGORIES:
            record[category] = int(
                counts.get(
                    category,
                    0,
                )
            )

        record[
            "endpoint_responders"
        ] = (
            record[
                "persistent_same_direction"
            ]
            + record[
                "recovered"
            ]
            + record[
                "substantial_reversal"
            ]
        )

        record[
            "endpoint_engagement_percent"
        ] = (
            record[
                "endpoint_responders"
            ]
            / total_genes
            * 100
        )

        responders = record[
            "endpoint_responders"
        ]

        if responders > 0:
            record[
                "persistence_among_responders_percent"
            ] = (
                record[
                    "persistent_same_direction"
                ]
                / responders
                * 100
            )

            record[
                "recovery_among_responders_percent"
            ] = (
                record[
                    "recovered"
                ]
                / responders
                * 100
            )

            record[
                "reversal_among_responders_percent"
            ] = (
                record[
                    "substantial_reversal"
                ]
                / responders
                * 100
            )

        else:
            record[
                "persistence_among_responders_percent"
            ] = float(
                "nan"
            )

            record[
                "recovery_among_responders_percent"
            ] = float(
                "nan"
            )

            record[
                "reversal_among_responders_percent"
            ] = float(
                "nan"
            )

        records.append(
            record
        )

    summary = (
        pd.DataFrame(records)
        .set_index(
            "condition"
        )
    )

    # Stage 1 produced an unnamed index after transposing
    # the recovery-count table. That exact output is preserved
    summary.index.name = None

    return summary