"""Summarize persistent expression patterns across stress conditions."""

import numpy as np
import pandas as pd

from arabidopsis_stress_recovery_atlas.analysis.exploratory import (
    calculate_persistence_breadth,
)


class ResponseProgramError(ValueError):
    """Raised when persistent-expression summaries cannot be calculated."""


def summarize_persistence_pattern(
    classified: pd.DataFrame,
    high_breadth_minimum: int = 4,
    multi_stress_minimum: int = 2,
) -> dict[str, float | int]:
    """Summarize persistence breadth and directional consistency."""

    required_columns = {
        "gene_id",
        "condition",
        "recovery_shift",
        "recovery_outcome",
    }

    missing_columns = (
        required_columns - set(classified.columns)
    )

    if missing_columns:
        raise ResponseProgramError(
            "Classified table is missing required columns: "
            + ", ".join(sorted(missing_columns))
        )

    breadth = calculate_persistence_breadth(
        classified
    )

    high_breadth_genes = int(
        (
            breadth["persistence_breadth"]
            >= high_breadth_minimum
        ).sum()
    )

    multi_stress_gene_ids = set(
        breadth.loc[
            breadth["persistence_breadth"]
            >= multi_stress_minimum,
            "gene_id",
        ]
    )

    persistent = classified[
        classified["recovery_outcome"]
        == "persistent_same_direction"
    ].copy()

    multi_stress_persistent = persistent[
        persistent["gene_id"].isin(
            multi_stress_gene_ids
        )
    ].copy()

    multi_stress_gene_count = len(
        multi_stress_gene_ids
    )

    if multi_stress_gene_count == 0:
        raise ResponseProgramError(
            "No multi-stress persistent genes were found."
        )

    multi_stress_persistent["direction_sign"] = (
        np.sign(
            multi_stress_persistent[
                "recovery_shift"
            ]
        )
    )

    sign_summary = (
        multi_stress_persistent
        .groupby("gene_id")["direction_sign"]
        .agg(
            unique_signs="nunique",
            direction="first",
        )
    )

    consistent = sign_summary[
        sign_summary["unique_signs"] == 1
    ]

    directionally_consistent_genes = len(
        consistent
    )

    consistently_above_control = int(
        (consistent["direction"] > 0).sum()
    )

    consistently_below_control = int(
        (consistent["direction"] < 0).sum()
    )

    directional_consistency_percent = (
        directionally_consistent_genes
        / multi_stress_gene_count
        * 100
    )

    consistently_above_percent = (
        consistently_above_control
        / directionally_consistent_genes
        * 100
    )

    return {
        "high_breadth_genes": high_breadth_genes,
        "multi_stress_persistent_genes": (
            multi_stress_gene_count
        ),
        "directionally_consistent_genes": (
            directionally_consistent_genes
        ),
        "consistently_above_control": (
            consistently_above_control
        ),
        "consistently_below_control": (
            consistently_below_control
        ),
        "directional_consistency_percent": (
            directional_consistency_percent
        ),
        "consistently_above_percent": (
            consistently_above_percent
        ),
    }

def summarize_breadth_distribution(
    classified: pd.DataFrame,
    maximum_conditions: int,
) -> dict[int, int]:
    """Count genes at each exact persistence breadth."""

    breadth = calculate_persistence_breadth(
        classified
    )

    persistent_breadth = breadth[
        breadth["persistence_breadth"] > 0
    ]

    counts = (
        persistent_breadth["persistence_breadth"]
        .value_counts()
        .reindex(
            range(1, maximum_conditions + 1),
            fill_value=0,
        )
        .sort_index()
    )

    return {
        int(breadth_value): int(count)
        for breadth_value, count in counts.items()
    }


def summarize_condition_persistence_fractions(
    classified: pd.DataFrame,
) -> pd.DataFrame:
    """Calculate the persistent fraction of endpoint-eligible genes."""

    required_columns = {
        "condition",
        "endpoint_eligible",
        "recovery_outcome",
    }

    missing_columns = (
        required_columns - set(classified.columns)
    )

    if missing_columns:
        raise ResponseProgramError(
            "Classified table is missing required columns: "
            + ", ".join(sorted(missing_columns))
        )

    endpoint_eligible = (
        classified[
            classified["endpoint_eligible"]
        ]
        .groupby("condition")
        .size()
    )

    persistent = (
        classified[
            classified["recovery_outcome"]
            == "persistent_same_direction"
        ]
        .groupby("condition")
        .size()
    )

    summary = pd.DataFrame(
        {
            "endpoint_eligible_genes": (
                endpoint_eligible
            ),
            "persistent_same_direction_genes": (
                persistent
            ),
        }
    ).fillna(0)

    summary[
        "persistent_same_direction_genes"
    ] = summary[
        "persistent_same_direction_genes"
    ].astype(int)

    summary["persistence_fraction_percent"] = (
        summary[
            "persistent_same_direction_genes"
        ]
        / summary["endpoint_eligible_genes"]
        * 100
    )

    return summary.reset_index()


def calculate_persistence_rank_correlation(
    summary_a: pd.DataFrame,
    summary_b: pd.DataFrame,
) -> float:
    """Calculate Spearman rank correlation of condition persistence fractions."""

    required_columns = {
        "condition",
        "persistence_fraction_percent",
    }

    for summary in (summary_a, summary_b):
        missing_columns = (
            required_columns - set(summary.columns)
        )

        if missing_columns:
            raise ResponseProgramError(
                "Persistence-fraction table is missing required columns: "
                + ", ".join(sorted(missing_columns))
            )

    merged = summary_a[
        [
            "condition",
            "persistence_fraction_percent",
        ]
    ].merge(
        summary_b[
            [
                "condition",
                "persistence_fraction_percent",
            ]
        ],
        on="condition",
        suffixes=("_a", "_b"),
        validate="one_to_one",
    )

    ranks_a = merged[
        "persistence_fraction_percent_a"
    ].rank(method="average")

    ranks_b = merged[
        "persistence_fraction_percent_b"
    ].rank(method="average")

    correlation = ranks_a.corr(ranks_b)

    return float(correlation)


def get_persistent_gene_condition_pairs(
    classified: pd.DataFrame,
) -> set[tuple[str, str]]:
    """Return persistent gene-condition pairs."""

    persistent = classified[
        classified["recovery_outcome"]
        == "persistent_same_direction"
    ]

    return set(
        zip(
            persistent["gene_id"],
            persistent["condition"],
        )
    )