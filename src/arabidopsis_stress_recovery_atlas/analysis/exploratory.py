"""Exploratory scientific calculations for stress-recovery analysis."""

from itertools import combinations

import numpy as np
import pandas as pd


class ExploratoryAnalysisError(ValueError):
    """Raised when an exploratory calculation cannot be performed."""


def calculate_control_pairwise_differences(
    expression: pd.DataFrame,
    manifest: pd.DataFrame,
) -> pd.DataFrame:
    """Calculate absolute expression differences between control replicates."""

    required_expression_column = "gene_id"

    if required_expression_column not in expression.columns:
        raise ExploratoryAnalysisError(
            "Expression table must contain canonical column 'gene_id'."
        )

    required_manifest_columns = {
        "sample_id",
        "state",
    }

    missing_manifest_columns = (
        required_manifest_columns - set(manifest.columns)
    )

    if missing_manifest_columns:
        raise ExploratoryAnalysisError(
            "Sample manifest is missing required columns: "
            + ", ".join(sorted(missing_manifest_columns))
        )

    control_samples = manifest.loc[
        manifest["state"] == "control",
        "sample_id",
    ].tolist()

    if len(control_samples) < 2:
        raise ExploratoryAnalysisError(
            "At least two control samples are required."
        )

    missing_control_samples = [
        sample
        for sample in control_samples
        if sample not in expression.columns
    ]

    if missing_control_samples:
        raise ExploratoryAnalysisError(
            "Control samples are missing from the expression table: "
            + ", ".join(missing_control_samples)
        )

    pairwise_tables = []

    for sample_a, sample_b in combinations(
        control_samples,
        2,
    ):
        absolute_difference = (
            expression[sample_a]
            - expression[sample_b]
        ).abs()

        pairwise_table = pd.DataFrame(
            {
                "gene_id": expression["gene_id"].to_numpy(),
                "sample_a": sample_a,
                "sample_b": sample_b,
                "absolute_difference": (
                    absolute_difference.to_numpy()
                ),
            }
        )

        pairwise_tables.append(pairwise_table)

    return pd.concat(
        pairwise_tables,
        ignore_index=True,
    )


def calculate_response_floor(
    pairwise_differences: pd.DataFrame,
    quantile: float,
) -> float:
    """Calculate a response floor from control-control differences."""

    if not 0 <= quantile <= 1:
        raise ExploratoryAnalysisError(
            "Quantile must be between 0 and 1."
        )

    if "absolute_difference" not in pairwise_differences.columns:
        raise ExploratoryAnalysisError(
            "Pairwise-difference table must contain "
            "'absolute_difference'."
        )

    response_floor = pairwise_differences[
        "absolute_difference"
    ].quantile(quantile)

    return float(response_floor)


def calculate_response_floors(
    pairwise_differences: pd.DataFrame,
    quantiles: list[float],
) -> dict[str, float]:
    """Calculate multiple percentile-based response floors."""

    floors = {}

    for quantile in quantiles:
        percentile_label = (
            f"p{int(round(quantile * 100))}"
        )

        floors[percentile_label] = (
            calculate_response_floor(
                pairwise_differences,
                quantile,
            )
        )

    return floors



def calculate_condition_shifts(
    expression: pd.DataFrame,
    manifest: pd.DataFrame,
) -> pd.DataFrame:
    """Calculate endpoint and recovery shifts from the control mean."""

    if "gene_id" not in expression.columns:
        raise ExploratoryAnalysisError(
            "Expression table must contain canonical column 'gene_id'."
        )

    required_manifest_columns = {
        "sample_id",
        "state",
        "condition",
    }

    missing_manifest_columns = (
        required_manifest_columns - set(manifest.columns)
    )

    if missing_manifest_columns:
        raise ExploratoryAnalysisError(
            "Sample manifest is missing required columns: "
            + ", ".join(sorted(missing_manifest_columns))
        )

    control_samples = manifest.loc[
        manifest["state"] == "control",
        "sample_id",
    ].tolist()

    if not control_samples:
        raise ExploratoryAnalysisError(
            "No control samples were found."
        )

    missing_samples = [
        sample
        for sample in manifest["sample_id"]
        if sample not in expression.columns
    ]

    if missing_samples:
        raise ExploratoryAnalysisError(
            "Manifest samples are missing from the expression table: "
            + ", ".join(missing_samples)
        )

    control_mean = expression[
        control_samples
    ].mean(axis=1)

    stress_manifest = manifest[
        manifest["state"] != "control"
    ]

    conditions = (
        stress_manifest["condition"]
        .drop_duplicates()
        .tolist()
    )

    condition_tables = []

    for condition in conditions:
        endpoint_samples = manifest.loc[
            (manifest["state"] == "endpoint")
            & (manifest["condition"] == condition),
            "sample_id",
        ].tolist()

        recovery_samples = manifest.loc[
            (manifest["state"] == "recovery")
            & (manifest["condition"] == condition),
            "sample_id",
        ].tolist()

        if not endpoint_samples:
            raise ExploratoryAnalysisError(
                f"No endpoint samples found for condition '{condition}'."
            )

        if not recovery_samples:
            raise ExploratoryAnalysisError(
                f"No recovery samples found for condition '{condition}'."
            )

        endpoint_mean = expression[
            endpoint_samples
        ].mean(axis=1)

        recovery_mean = expression[
            recovery_samples
        ].mean(axis=1)

        condition_table = pd.DataFrame(
            {
                "gene_id": expression["gene_id"].to_numpy(),
                "condition": condition,
                "control_mean": control_mean.to_numpy(),
                "endpoint_mean": endpoint_mean.to_numpy(),
                "recovery_mean": recovery_mean.to_numpy(),
                "endpoint_shift": (
                    endpoint_mean - control_mean
                ).to_numpy(),
                "recovery_shift": (
                    recovery_mean - control_mean
                ).to_numpy(),
            }
        )

        condition_tables.append(
            condition_table
        )

    return pd.concat(
        condition_tables,
        ignore_index=True,
    )


def classify_recovery_outcomes(
    shifts: pd.DataFrame,
    response_floor: float,
) -> pd.DataFrame:
    """Classify recovery outcomes using the exploratory response floor."""

    if response_floor < 0:
        raise ExploratoryAnalysisError(
            "Response floor cannot be negative."
        )

    required_columns = {
        "gene_id",
        "condition",
        "endpoint_shift",
        "recovery_shift",
    }

    missing_columns = (
        required_columns - set(shifts.columns)
    )

    if missing_columns:
        raise ExploratoryAnalysisError(
            "Shift table is missing required columns: "
            + ", ".join(sorted(missing_columns))
        )

    classified = shifts.copy()

    classified["endpoint_eligible"] = (
        classified["endpoint_shift"].abs()
        >= response_floor
    )

    classified["recovery_above_floor"] = (
        classified["recovery_shift"].abs()
        >= response_floor
    )

    classified["same_direction"] = (
        np.sign(classified["endpoint_shift"])
        == np.sign(classified["recovery_shift"])
    )

    classified["recovery_outcome"] = (
        "not_endpoint_eligible"
    )

    recovered_mask = (
        classified["endpoint_eligible"]
        & ~classified["recovery_above_floor"]
    )

    persistent_mask = (
        classified["endpoint_eligible"]
        & classified["recovery_above_floor"]
        & classified["same_direction"]
    )

    reversal_mask = (
        classified["endpoint_eligible"]
        & classified["recovery_above_floor"]
        & ~classified["same_direction"]
    )

    classified.loc[
        recovered_mask,
        "recovery_outcome",
    ] = "recovered"

    classified.loc[
        persistent_mask,
        "recovery_outcome",
    ] = "persistent_same_direction"

    classified.loc[
        reversal_mask,
        "recovery_outcome",
    ] = "substantial_reversal"

    return classified


def summarize_persistent_counts(
    classified: pd.DataFrame,
) -> dict[str, int]:
    """Count persistent same-direction genes within each condition."""

    conditions = (
        classified["condition"]
        .drop_duplicates()
        .tolist()
    )

    persistent = classified[
        classified["recovery_outcome"]
        == "persistent_same_direction"
    ]

    counts = (
        persistent.groupby("condition")["gene_id"]
        .nunique()
        .reindex(
            conditions,
            fill_value=0,
        )
    )

    return {
        condition: int(count)
        for condition, count in counts.items()
    }


def calculate_persistence_breadth(
    classified: pd.DataFrame,
) -> pd.DataFrame:
    """Count the number of stress conditions showing persistence per gene."""

    all_genes = (
        classified["gene_id"]
        .drop_duplicates()
    )

    persistent = classified[
        classified["recovery_outcome"]
        == "persistent_same_direction"
    ]

    breadth = (
        persistent.groupby("gene_id")["condition"]
        .nunique()
        .reindex(
            all_genes,
            fill_value=0,
        )
        .astype(int)
    )

    return pd.DataFrame(
        {
            "gene_id": breadth.index,
            "persistence_breadth": breadth.to_numpy(),
        }
    )