"""Gene Ontology enrichment analysis for persistent expression genes."""

from collections.abc import Iterable

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

from arabidopsis_stress_recovery_atlas.annotation.go import (
    get_go_term_memberships,
    get_goslim_memberships,
)


class GOEnrichmentError(ValueError):
    """Raised when GO enrichment cannot be calculated."""


ASPECT_NAMES = {
    "P": "Biological Process",
    "F": "Molecular Function",
    "C": "Cellular Component",
}


def _as_gene_set(
    genes: Iterable[str],
) -> set[str]:
    """Convert gene identifiers to a clean set of strings."""

    return {
        str(gene)
        for gene in genes
    }


def _validate_gene_universe(
    background_ids: set[str],
    foreground_ids: set[str],
) -> None:
    """Check that the enrichment foreground belongs to the background."""

    if not foreground_ids:
        raise GOEnrichmentError(
            "Foreground gene set is empty."
        )

    if not background_ids:
        raise GOEnrichmentError(
            "Background gene set is empty."
        )

    if not foreground_ids.issubset(
        background_ids
    ):
        raise GOEnrichmentError(
            "Foreground genes must be a subset of the background."
        )


def calculate_go_slim_enrichment(
    go_annotations: pd.DataFrame,
    background_genes: Iterable[str],
    foreground_genes: Iterable[str],
    alpha: float = 0.05,
    aspects: tuple[str, ...] = (
        "P",
        "F",
        "C",
    ),
)-> tuple[pd.DataFrame, pd.DataFrame]:
    """Calculate aspect-specific GO-Slim over-representation."""

    background_ids = _as_gene_set(
        background_genes
    )

    foreground_ids = _as_gene_set(
        foreground_genes
    )

    _validate_gene_universe(
        background_ids,
        foreground_ids,
    )

    go_slim_pairs = get_goslim_memberships(
        go_annotations
    ).rename(
        columns={
            "aspect": "go_aspect",
            "goslim_category": "go_slim_term",
        }
    )

    for column in (
        "gene_id",
        "go_aspect",
        "go_slim_term",
    ):
        go_slim_pairs[column] = (
            go_slim_pairs[column]
            .astype(str)
            .str.strip()
        )

    enrichment_records = []

    for aspect in aspects:
        aspect_pairs = go_slim_pairs[
            go_slim_pairs["go_aspect"]
            == aspect
        ]

        aspect_gene_ids = set(
            aspect_pairs["gene_id"]
        )

        background_aspect = (
            background_ids
            & aspect_gene_ids
        )

        foreground_aspect = (
            foreground_ids
            & aspect_gene_ids
        )

        if not foreground_aspect:
            raise GOEnrichmentError(
                f"No foreground genes have GO-Slim annotations "
                f"for aspect '{aspect}'."
            )

        terms = sorted(
            aspect_pairs[
                "go_slim_term"
            ].unique()
        )

        for term in terms:
            term_gene_ids = set(
                aspect_pairs.loc[
                    aspect_pairs[
                        "go_slim_term"
                    ]
                    == term,
                    "gene_id",
                ]
            )

            foreground_with_term = len(
                foreground_aspect
                & term_gene_ids
            )

            foreground_without_term = (
                len(foreground_aspect)
                - foreground_with_term
            )

            background_with_term = len(
                background_aspect
                & term_gene_ids
            )

            other_background_with_term = (
                background_with_term
                - foreground_with_term
            )

            other_background_without_term = (
                len(background_aspect)
                - len(foreground_aspect)
                - other_background_with_term
            )

            odds_ratio, p_value = fisher_exact(
                [
                    [
                        foreground_with_term,
                        foreground_without_term,
                    ],
                    [
                        other_background_with_term,
                        other_background_without_term,
                    ],
                ],
                alternative="greater",
            )

            foreground_fraction = (
                foreground_with_term
                / len(foreground_aspect)
            )

            background_fraction = (
                background_with_term
                / len(background_aspect)
            )

            enrichment_ratio = (
                foreground_fraction
                / background_fraction
                if background_fraction > 0
                else np.nan
            )

            enrichment_records.append(
                {
                    "go_aspect": aspect,
                    "go_slim_term": term,
                    "foreground_with_term": (
                        foreground_with_term
                    ),
                    "foreground_total": (
                        len(foreground_aspect)
                    ),
                    "background_with_term": (
                        background_with_term
                    ),
                    "background_total": (
                        len(background_aspect)
                    ),
                    "foreground_fraction": (
                        foreground_fraction
                    ),
                    "background_fraction": (
                        background_fraction
                    ),
                    "enrichment_ratio": (
                        enrichment_ratio
                    ),
                    "odds_ratio": odds_ratio,
                    "p_value": p_value,
                }
            )

    enrichment = pd.DataFrame(
        enrichment_records
    )

    enrichment["fdr_global"] = (
        multipletests(
            enrichment["p_value"],
            method="fdr_bh",
        )[1]
    )

    enrichment["fdr_within_aspect"] = np.nan

    for aspect in aspects:
        mask = (
            enrichment["go_aspect"]
            == aspect
        )

        enrichment.loc[
            mask,
            "fdr_within_aspect",
        ] = multipletests(
            enrichment.loc[
                mask,
                "p_value",
            ],
            method="fdr_bh",
        )[1]

    significant = enrichment[
        enrichment["fdr_global"] < alpha
    ].copy()

    significant["aspect_name"] = (
        significant["go_aspect"]
        .map(ASPECT_NAMES)
    )

    significant["foreground_percent"] = (
        significant["foreground_fraction"]
        * 100
    )

    significant["background_percent"] = (
        significant["background_fraction"]
        * 100
    )

    significant = (
        significant
        .sort_values(
            [
                "fdr_global",
                "enrichment_ratio",
            ],
            ascending=[
                True,
                False,
            ],
        )
        .reset_index(drop=True)
    )

    return enrichment, significant


def calculate_bp_enrichment(
    go_annotations: pd.DataFrame,
    background_genes: Iterable[str],
    foreground_genes: Iterable[str],
    minimum_term_size: int = 10,
    root_go_id: str = "GO:0008150",
    alpha: float = 0.05,
) -> tuple[
    pd.DataFrame,
    pd.DataFrame,
    set[str],
    set[str],
]:
    """Calculate detailed Biological Process over-representation."""

    background_ids = _as_gene_set(
        background_genes
    )

    foreground_ids = _as_gene_set(
        foreground_genes
    )

    _validate_gene_universe(
        background_ids,
        foreground_ids,
    )

    bp_pairs = get_go_term_memberships(
        go_annotations,
        aspect="P",
    )

    bp_gene_ids = set(
        bp_pairs["gene_id"]
    )

    bp_background_ids = (
        background_ids
        & bp_gene_ids
    )

    bp_foreground_ids = (
        foreground_ids
        & bp_gene_ids
    )

    if not bp_foreground_ids:
        raise GOEnrichmentError(
            "No foreground genes have Biological Process annotations."
        )

    bp_term_sizes = (
        bp_pairs[
            bp_pairs["gene_id"].isin(
                bp_background_ids
            )
        ]
        .groupby(
            [
                "go_id",
                "go_term",
            ]
        )["gene_id"]
        .nunique()
    )

    eligible_bp_terms = (
        bp_term_sizes[
            bp_term_sizes
            >= minimum_term_size
        ]
        .reset_index(
            name="background_gene_count"
        )
    )

    eligible_bp_terms = (
        eligible_bp_terms[
            eligible_bp_terms["go_id"]
            != root_go_id
        ]
        .copy()
    )

    bp_term_gene_sets = (
        bp_pairs[
            bp_pairs["gene_id"].isin(
                bp_background_ids
            )
        ]
        .groupby(
            [
                "go_id",
                "go_term",
            ]
        )["gene_id"]
        .apply(set)
    )

    enrichment_records = []

    for row in eligible_bp_terms.itertuples(
        index=False
    ):
        go_id = row.go_id
        go_term = row.go_term

        term_genes = bp_term_gene_sets.loc[
            (go_id, go_term)
        ]

        foreground_with_term = len(
            foreground_ids
            & term_genes
        )

        foreground_without_term = (
            len(bp_foreground_ids)
            - foreground_with_term
        )

        background_with_term = len(
            bp_background_ids
            & term_genes
        )

        other_background_with_term = (
            background_with_term
            - foreground_with_term
        )

        other_background_without_term = (
            len(bp_background_ids)
            - len(bp_foreground_ids)
            - other_background_with_term
        )

        odds_ratio, p_value = fisher_exact(
            [
                [
                    foreground_with_term,
                    foreground_without_term,
                ],
                [
                    other_background_with_term,
                    other_background_without_term,
                ],
            ],
            alternative="greater",
        )

        foreground_fraction = (
            foreground_with_term
            / len(bp_foreground_ids)
        )

        background_fraction = (
            background_with_term
            / len(bp_background_ids)
        )

        enrichment_ratio = (
            foreground_fraction
            / background_fraction
        )

        enrichment_records.append(
            {
                "go_id": go_id,
                "go_term": go_term,
                "foreground_with_term": (
                    foreground_with_term
                ),
                "foreground_total": (
                    len(bp_foreground_ids)
                ),
                "background_with_term": (
                    background_with_term
                ),
                "background_total": (
                    len(bp_background_ids)
                ),
                "foreground_fraction": (
                    foreground_fraction
                ),
                "background_fraction": (
                    background_fraction
                ),
                "enrichment_ratio": (
                    enrichment_ratio
                ),
                "odds_ratio": odds_ratio,
                "p_value": p_value,
            }
        )

    enrichment = pd.DataFrame(
        enrichment_records
    )

    enrichment["fdr_bh"] = (
        multipletests(
            enrichment["p_value"],
            method="fdr_bh",
        )[1]
    )

    significant = (
        enrichment[
            enrichment["fdr_bh"]
            < alpha
        ]
        .copy()
        .sort_values(
            [
                "fdr_bh",
                "enrichment_ratio",
            ],
            ascending=[
                True,
                False,
            ],
        )
        .reset_index(drop=True)
    )

    return (
        enrichment,
        significant,
        bp_background_ids,
        bp_foreground_ids,
    )


def build_significant_bp_membership(
    go_annotations: pd.DataFrame,
    significant_enrichment: pd.DataFrame,
    bp_foreground_ids: set[str],
) -> pd.DataFrame:
    """Return foreground gene memberships for significant BP terms."""

    bp_pairs = get_go_term_memberships(
        go_annotations,
        aspect="P",
    )

    significant_bp_ids = set(
        significant_enrichment["go_id"]
    )

    membership = (
        bp_pairs[
            bp_pairs["go_id"].isin(
                significant_bp_ids
            )
            & bp_pairs["gene_id"].isin(
                bp_foreground_ids
            )
        ][
            [
                "gene_id",
                "go_id",
                "go_term",
            ]
        ]
        .drop_duplicates()
        .reset_index(drop=True)
    )

    return membership