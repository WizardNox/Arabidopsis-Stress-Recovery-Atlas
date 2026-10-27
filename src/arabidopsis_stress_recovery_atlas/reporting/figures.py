"""Generate reproducible figures from database-derived reporting tables."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def plot_persistent_genes_by_condition(
    condition_summary: pd.DataFrame,
    output_path: Path,
) -> Path:
    """Plot persistent same-direction gene counts by stress condition."""

    output_path = Path(
        output_path
    )

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    data = (
        condition_summary
        .sort_values(
            "persistent_genes",
            ascending=False,
        )
        .reset_index(
            drop=True
        )
    )

    figure, axis = plt.subplots(
        figsize=(
            10,
            6,
        )
    )

    axis.bar(
        data[
            "condition"
        ],
        data[
            "persistent_genes"
        ],
    )

    axis.set_title(
        "Persistent same-direction genes by stress condition"
    )

    axis.set_xlabel(
        "Stress condition"
    )

    axis.set_ylabel(
        "Persistent genes"
    )

    axis.tick_params(
        axis="x",
        labelrotation=45,
    )

    for label in (
        axis.get_xticklabels()
    ):
        label.set_horizontalalignment(
            "right"
        )

    figure.tight_layout()

    figure.savefig(
        output_path,
        dpi=200,
        bbox_inches="tight",
    )

    plt.close(
        figure
    )

    return output_path


def plot_recovery_outcomes_by_condition(
    condition_summary: pd.DataFrame,
    output_path: Path,
) -> Path:
    """Plot recovery outcomes among endpoint-eligible genes."""

    output_path = Path(
        output_path
    )

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    data = (
        condition_summary
        .sort_values(
            "condition"
        )
        .reset_index(
            drop=True
        )
    )

    x = np.arange(
        len(
            data
        )
    )

    recovered = data[
        "recovered_genes"
    ].to_numpy()

    persistent = data[
        "persistent_genes"
    ].to_numpy()

    reversed_genes = data[
        "substantial_reversal_genes"
    ].to_numpy()

    figure, axis = plt.subplots(
        figsize=(
            11,
            6,
        )
    )

    axis.bar(
        x,
        recovered,
        label="Recovered",
    )

    axis.bar(
        x,
        persistent,
        bottom=recovered,
        label="Persistent same direction",
    )

    axis.bar(
        x,
        reversed_genes,
        bottom=(
            recovered
            + persistent
        ),
        label="Substantial reversal",
    )

    axis.set_title(
        "Recovery outcomes among endpoint-eligible genes"
    )

    axis.set_xlabel(
        "Stress condition"
    )

    axis.set_ylabel(
        "Genes"
    )

    axis.set_xticks(
        x
    )

    axis.set_xticklabels(
        data[
            "condition"
        ],
        rotation=45,
        ha="right",
    )

    axis.legend()

    figure.tight_layout()

    figure.savefig(
        output_path,
        dpi=200,
        bbox_inches="tight",
    )

    plt.close(
        figure
    )

    return output_path