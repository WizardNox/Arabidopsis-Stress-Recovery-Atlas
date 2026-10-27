"""Build reproducible reporting outputs from the analytical database."""

from pathlib import Path
from typing import Any

import pandas as pd

from arabidopsis_stress_recovery_atlas.database.connection import (
    connect_database,
)
from arabidopsis_stress_recovery_atlas.reporting.figures import (
    plot_persistent_genes_by_condition,
    plot_recovery_outcomes_by_condition,
)
from arabidopsis_stress_recovery_atlas.reporting.results import (
    build_reporting_tables,
    get_project_summary,
)


def _format_value(
    value,
) -> str:
    """Format a value for a compact Markdown table."""

    if pd.isna(
        value
    ):
        return ""

    if isinstance(
        value,
        float,
    ):
        return (
            f"{value:.6g}"
        )

    return str(
        value
    ).replace(
        "|",
        "\\|",
    )


def _markdown_table(
    dataframe: pd.DataFrame,
    columns: list[str],
    limit: int | None = None,
) -> str:
    """Render a small DataFrame as Markdown without extra dependencies."""

    data = dataframe[
        columns
    ]

    if limit is not None:
        data = data.head(
            limit
        )

    header = (
        "| "
        + " | ".join(
            columns
        )
        + " |"
    )

    separator = (
        "| "
        + " | ".join(
            "---"
            for _ in columns
        )
        + " |"
    )

    rows = []

    for row in data.itertuples(
        index=False,
        name=None,
    ):
        rows.append(
            "| "
            + " | ".join(
                _format_value(
                    value
                )
                for value in row
            )
            + " |"
        )

    return "\n".join(
        [
            header,
            separator,
            *rows,
        ]
    )


def render_summary_report(
    summary: dict[str, Any],
    tables: dict[str, pd.DataFrame],
) -> str:
    """Render the database-derived scientific summary as Markdown."""

    condition_table = (
        _markdown_table(
            tables[
                "condition_recovery_summary"
            ],
            [
                "condition",
                "persistent_genes",
                "recovered_genes",
                "substantial_reversal_genes",
            ],
        )
    )

    high_breadth_table = (
        _markdown_table(
            tables[
                "high_breadth_genes"
            ],
            [
                "gene_id",
                "persistence_breadth",
                "n_persistent_up",
                "n_persistent_down",
                "short_description",
            ],
            limit=12,
        )
    )

    published_table = (
        _markdown_table(
            tables[
                "published_de_summary"
            ],
            [
                "condition",
                "shared_gene_records",
                "direction_agreements",
                "direction_agreement_percent",
            ],
        )
    )

    return f"""# Arabidopsis Stress Recovery Atlas — Generated Summary

This report is generated from the canonical SQLite analytical database.
It does not rerun or redefine the scientific classification rules.

## Dataset and analytical layer

- Samples: {summary["sample_count"]}
- Informative genes: {summary["gene_count"]}
- Stress conditions: {summary["stress_condition_count"]}
- Gene-condition records: {summary["gene_condition_rows"]}
- Primary control-derived response floor: {summary["response_floor"]:.6f}

The response floor is an exploratory magnitude threshold derived from
control-to-control variation. It is not a differential-expression
significance threshold.

## Persistent recovery signal

Genes classified as `persistent_same_direction` retained a recovery
RNA-abundance shift beyond the primary response floor in the same
direction as their endpoint shift.

High-breadth genes persistent across at least four stress conditions:
**{summary["high_breadth_genes"]}**

### Recovery outcomes by condition

{condition_table}

![Persistent genes by condition](figures/persistent_genes_by_condition.png)

![Recovery outcomes among endpoint-eligible genes](figures/recovery_outcomes_by_condition.png)

### Highest persistence breadth

{high_breadth_table}

The high-breadth definition is an operational project threshold and
should not be interpreted as a distinct biological gene class.

## Functional enrichment

Significant GO-Slim terms: **{summary["significant_go_slim_terms"]}**

Significant detailed Biological Process terms:
**{summary["significant_detailed_bp_terms"]}**

The complete significant-term table is available in
`tables/significant_go_terms.csv`.

GO enrichment describes over-represented annotations among the selected
gene set. It does not establish causal pathways or prove activation of
the biological processes named by the terms.

## Published recovery-DE comparison

Shared project/published gene-condition records:
**{summary["published_shared_records"]}**

Direction agreements after harmonizing the published fold-change
orientation:
**{summary["published_direction_agreements"]} /
{summary["published_shared_records"]}**

{published_table}

This is a same-study methodological concordance check using available
published recovery-DE results. It is not independent biological
validation.

## Interpretation boundary

The project describes persistent residual RNA-abundance patterns after
stress removal. Persistent rlog-normalized expression should not be
interpreted by itself as evidence of continued transcription, protein
activity, causal stress memory, or a longitudinal response in the same
individual plants.
"""


def build_reporting_outputs(
    database_path: Path,
    output_directory: Path,
) -> dict[str, Any]:
    """Generate all reporting tables, figures, and the summary report."""

    database_path = Path(
        database_path
    )

    output_directory = Path(
        output_directory
    )

    tables_directory = (
        output_directory
        / "tables"
    )

    figures_directory = (
        output_directory
        / "figures"
    )

    output_directory.mkdir(
        parents=True,
        exist_ok=True,
    )

    tables_directory.mkdir(
        parents=True,
        exist_ok=True,
    )

    figures_directory.mkdir(
        parents=True,
        exist_ok=True,
    )

    connection = (
        connect_database(
            database_path,
            read_only=True,
        )
    )

    try:
        summary = (
            get_project_summary(
                connection
            )
        )

        tables = (
            build_reporting_tables(
                connection
            )
        )

    finally:
        connection.close()

    table_paths = {}

    for table_name, dataframe in (
        tables.items()
    ):
        path = (
            tables_directory
            / f"{table_name}.csv"
        )

        dataframe.to_csv(
            path,
            index=False,
        )

        table_paths[
            table_name
        ] = path

    figure_paths = {
        "persistent_genes_by_condition": (
            plot_persistent_genes_by_condition(
                tables[
                    "condition_recovery_summary"
                ],
                figures_directory
                / "persistent_genes_by_condition.png",
            )
        ),
        "recovery_outcomes_by_condition": (
            plot_recovery_outcomes_by_condition(
                tables[
                    "condition_recovery_summary"
                ],
                figures_directory
                / "recovery_outcomes_by_condition.png",
            )
        ),
    }

    summary_path = (
        output_directory
        / "summary.md"
    )

    summary_text = (
        render_summary_report(
            summary,
            tables,
        )
    )

    summary_path.write_text(
        summary_text,
        encoding="utf-8",
    )

    return {
        "summary": summary,
        "tables": tables,
        "table_paths": table_paths,
        "figure_paths": figure_paths,
        "summary_path": summary_path,
    }