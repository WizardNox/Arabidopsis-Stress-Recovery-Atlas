"""Compare exploratory persistence with published recovery-DE results."""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

from arabidopsis_stress_recovery_atlas.config import (
    PROJECT_ROOT,
    load_sources_config,
)


class CrossStudyError(ValueError):
    """Raised when the published-DE comparison cannot be performed."""


PUBLISHED_DE_MANIFEST_RECORDS = [
    (
        "HL ep",
        "endpoint",
        "high_light",
    ),
    (
        "Humidity ep",
        "endpoint",
        "high_humidity",
    ),
    (
        "Drought ep",
        "endpoint",
        "drought",
    ),
    (
        "Heat ep",
        "endpoint",
        "heat",
    ),
    (
        "Heat drought ep",
        "endpoint",
        "heat_drought",
    ),
    (
        "Heat HL ep",
        "endpoint",
        "heat_high_light",
    ),
    (
        "Heat drought hl ep",
        "endpoint",
        "heat_drought_high_light",
    ),
    (
        "Cold ep",
        "endpoint",
        "cold",
    ),
    (
        "Cold HL ep",
        "endpoint",
        "cold_high_light",
    ),
    (
        "HL rec",
        "recovery",
        "high_light",
    ),
    (
        "Humidity rec",
        "recovery",
        "high_humidity",
    ),
    (
        "Drought rec",
        "recovery",
        "drought",
    ),
    (
        "Heat drought rec",
        "recovery",
        "heat_drought",
    ),
    (
        "Cold rec",
        "recovery",
        "cold",
    ),
]


def build_published_de_manifest() -> pd.DataFrame:
    """Return the Stage 1 mapping from workbook sheets to project conditions."""

    return pd.DataFrame(
        PUBLISHED_DE_MANIFEST_RECORDS,
        columns=[
            "sheet_name",
            "state",
            "condition",
        ],
    )


def get_published_de_workbook_path() -> Path:
    """Return the registered mmc16 published-DE workbook path."""

    sources = load_sources_config()

    registered_sources = sources.get(
        "sources",
        {},
    )

    for source in registered_sources.values():
        local_path = source.get(
            "local_path"
        )

        if (
            local_path
            and Path(local_path).name
            == "mmc16.xlsx"
        ):
            return (
                PROJECT_ROOT
                / local_path
            )

    raise CrossStudyError(
        "Published-DE source mmc16.xlsx is not registered "
        "in config/sources.yaml."
    )


def load_published_de_tables(
    workbook_path: Path | None = None,
    manifest: pd.DataFrame | None = None,
) -> dict[str, pd.DataFrame]:
    """Load the endpoint and available recovery DEG sheets."""

    if workbook_path is None:
        workbook_path = (
            get_published_de_workbook_path()
        )

    if manifest is None:
        manifest = (
            build_published_de_manifest()
        )

    if not workbook_path.is_file():
        raise CrossStudyError(
            "Published-DE workbook not found: "
            f"{workbook_path}"
        )

    tables = {}

    for row in manifest.itertuples(
        index=False
    ):
        try:
            table = pd.read_excel(
                workbook_path,
                sheet_name=row.sheet_name,
            )
        except (
            ValueError,
            OSError,
        ) as exc:
            raise CrossStudyError(
                "Could not load published-DE sheet "
                f"'{row.sheet_name}'."
            ) from exc

        # One Stage 1 worksheet used different
        # column names from all the other sheets.
        if row.sheet_name == "Heat HL ep":
            table = table.rename(
                columns={
                    "FC": "log2(FC)",
                    "pval": "P",
                    "FDR": "adj-P",
                }
            )

        required_columns = {
            "gene_ID",
            "log2(FC)",
            "adj-P",
        }

        missing_columns = (
            required_columns
            - set(table.columns)
        )

        if missing_columns:
            raise CrossStudyError(
                f"Published-DE sheet '{row.sheet_name}' "
                "is missing required columns: "
                + ", ".join(
                    sorted(
                        missing_columns
                    )
                )
            )

        tables[
            row.sheet_name
        ] = table.copy()

    return tables


def summarize_published_de_workbook(
    tables: dict[str, pd.DataFrame],
    manifest: pd.DataFrame,
) -> pd.DataFrame:
    """Summarize and validate all registered published DEG sheets."""

    records = []

    for row in manifest.itertuples(
        index=False
    ):
        table = tables[
            row.sheet_name
        ]

        records.append(
            {
                "sheet_name": (
                    row.sheet_name
                ),
                "state": row.state,
                "condition": (
                    row.condition
                ),
                "n_genes": len(
                    table
                ),
                "n_unique_genes": (
                    table[
                        "gene_ID"
                    ].nunique()
                ),
                "n_duplicate_genes": (
                    table[
                        "gene_ID"
                    ]
                    .duplicated()
                    .sum()
                ),
                "n_missing_values": (
                    table
                    .isna()
                    .sum()
                    .sum()
                ),
                "min_abs_log2fc": (
                    table[
                        "log2(FC)"
                    ]
                    .abs()
                    .min()
                ),
                "max_adj_p": (
                    table[
                        "adj-P"
                    ].max()
                ),
                "all_abs_log2fc_ge_1": (
                    (
                        table[
                            "log2(FC)"
                        ].abs()
                        >= 1
                    ).all()
                ),
                "all_adj_p_le_0_05": (
                    (
                        table[
                            "adj-P"
                        ]
                        <= 0.05
                    ).all()
                ),
            }
        )

    return pd.DataFrame(
        records
    )


def build_published_condition_counts(
    workbook_summary: pd.DataFrame,
) -> pd.DataFrame:
    """Summarize published endpoint and recovery DEG counts."""

    counts = (
        workbook_summary[
            [
                "state",
                "condition",
                "n_genes",
            ]
        ]
        .pivot(
            index="condition",
            columns="state",
            values="n_genes",
        )
        .reset_index()
    )

    counts.columns.name = None

    counts[
        "recovery_sheet_present"
    ] = counts[
        "recovery"
    ].notna()

    counts[
        "published_recovery_degs"
    ] = counts[
        "recovery"
    ].astype(
        "Int64"
    )

    counts[
        "endpoint"
    ] = counts[
        "endpoint"
    ].astype(
        int
    )

    counts[
        "recovery_count_status"
    ] = np.where(
        counts[
            "recovery_sheet_present"
        ],
        "exact count available",
        "exact count unavailable; paper reports <10",
    )

    return counts


def build_persistent_sets_from_recovery_categories(
    recovery_categories: pd.DataFrame,
) -> dict[str, set[str]]:
    """Build condition-specific persistent-same-direction gene sets."""

    table = recovery_categories.copy()

    if "gene_id" not in table.columns:
        table = (
            table
            .reset_index()
            .rename(
                columns={
                    table.index.name
                    or "index": "gene_id"
                }
            )
        )

    if "gene_id" not in table.columns:
        raise CrossStudyError(
            "Recovery-category table must contain 'gene_id'."
        )

    conditions = [
        column
        for column in table.columns
        if column != "gene_id"
    ]

    persistent_sets = {}

    for condition in conditions:
        persistent_sets[
            condition
        ] = set(
            table.loc[
                table[
                    condition
                ]
                == "persistent_same_direction",
                "gene_id",
            ]
            .astype(str)
            .str.strip()
        )

    return persistent_sets


def build_published_recovery_sets(
    tables: dict[str, pd.DataFrame],
    manifest: pd.DataFrame,
    conditions: list[str],
) -> dict[
    str,
    set[str] | None,
]:
    """Build published recovery-DEG sets where an exact sheet exists."""

    recovery_sets = {}

    for condition in conditions:
        matching_rows = manifest[
            (
                manifest[
                    "state"
                ]
                == "recovery"
            )
            & (
                manifest[
                    "condition"
                ]
                == condition
            )
        ]

        if len(
            matching_rows
        ) == 1:
            sheet = (
                matching_rows
                .iloc[0][
                    "sheet_name"
                ]
            )

            recovery_sets[
                condition
            ] = set(
                tables[
                    sheet
                ][
                    "gene_ID"
                ]
                .astype(str)
                .str.strip()
            )

        else:
            # Missing recovery worksheets mean the
            # exact count is unavailable, not zero.
            recovery_sets[
                condition
            ] = None

    return recovery_sets


def compare_persistent_to_published(
    persistent_sets: dict[
        str,
        set[str],
    ],
    published_recovery_sets: dict[
        str,
        set[str] | None,
    ],
) -> pd.DataFrame:
    """Compare exploratory persistent sets with published recovery DEGs."""

    records = []

    for condition in persistent_sets:
        our_genes = (
            persistent_sets[
                condition
            ]
        )

        published_genes = (
            published_recovery_sets.get(
                condition
            )
        )

        if published_genes is None:
            continue

        shared_genes = (
            our_genes
            & published_genes
        )

        records.append(
            {
                "condition": (
                    condition
                ),
                "our_persistent_genes": (
                    len(
                        our_genes
                    )
                ),
                "published_recovery_degs": (
                    len(
                        published_genes
                    )
                ),
                "shared_genes": (
                    len(
                        shared_genes
                    )
                ),
                "fraction_of_our_persistent_formally_de": (
                    len(
                        shared_genes
                    )
                    / len(
                        our_genes
                    )
                    if our_genes
                    else np.nan
                ),
                "fraction_of_published_degs_captured": (
                    len(
                        shared_genes
                    )
                    / len(
                        published_genes
                    )
                    if published_genes
                    else np.nan
                ),
            }
        )

    return pd.DataFrame(
        records
    )


def calculate_published_overlap_enrichment(
    persistent_sets: dict[
        str,
        set[str],
    ],
    published_recovery_sets: dict[
        str,
        set[str] | None,
    ],
    informative_gene_ids: set[str],
) -> pd.DataFrame:
    """Test enrichment of overlap with published recovery DEGs."""

    universe_size = len(
        informative_gene_ids
    )

    if universe_size == 0:
        raise CrossStudyError(
            "Informative-gene universe is empty."
        )

    records = []

    for condition in persistent_sets:
        our_genes = (
            persistent_sets[
                condition
            ]
        )

        published_genes = (
            published_recovery_sets.get(
                condition
            )
        )

        if published_genes is None:
            continue

        shared_genes = (
            our_genes
            & published_genes
        )

        a = len(
            shared_genes
        )

        b = len(
            our_genes
            - published_genes
        )

        c = len(
            published_genes
            - our_genes
        )

        d = (
            universe_size
            - a
            - b
            - c
        )

        if d < 0:
            raise CrossStudyError(
                "Published and exploratory gene sets are "
                "incompatible with the supplied universe."
            )

        odds_ratio, p_value = (
            fisher_exact(
                [
                    [
                        a,
                        b,
                    ],
                    [
                        c,
                        d,
                    ],
                ],
                alternative="greater",
            )
        )

        expected_overlap = (
            len(
                our_genes
            )
            * len(
                published_genes
            )
            / universe_size
        )

        records.append(
            {
                "condition": (
                    condition
                ),
                "our_persistent_genes": (
                    len(
                        our_genes
                    )
                ),
                "published_recovery_degs": (
                    len(
                        published_genes
                    )
                ),
                "shared_genes": a,
                "expected_overlap_by_chance": (
                    expected_overlap
                ),
                "overlap_enrichment": (
                    a
                    / expected_overlap
                    if expected_overlap
                    > 0
                    else np.nan
                ),
                "odds_ratio": (
                    odds_ratio
                ),
                "p_value": (
                    p_value
                ),
            }
        )

    results = pd.DataFrame(
        records
    )

    if not results.empty:
        results[
            "fdr_bh"
        ] = multipletests(
            results[
                "p_value"
            ],
            method="fdr_bh",
        )[1]

    return results


def build_published_recovery_long(
    tables: dict[str, pd.DataFrame],
    manifest: pd.DataFrame,
) -> pd.DataFrame:
    """Create a long recovery-DE table while preserving original contrasts."""

    records = []

    recovery_manifest = (
        manifest[
            manifest[
                "state"
            ]
            == "recovery"
        ]
    )

    for row in (
        recovery_manifest
        .itertuples(
            index=False
        )
    ):
        table = tables[
            row.sheet_name
        ].copy()

        table[
            "condition"
        ] = row.condition

        table[
            "published_direction"
        ] = np.where(
            table[
                "log2(FC)"
            ]
            > 0,
            "up",
            "down",
        )

        table = table[
            [
                "gene_ID",
                "condition",
                "log2(FC)",
                "adj-P",
                "published_direction",
            ]
        ].rename(
            columns={
                "gene_ID": (
                    "gene_id"
                ),
                "log2(FC)": (
                    "published_log2fc"
                ),
                "adj-P": (
                    "published_adj_p"
                ),
            }
        )

        records.append(
            table
        )

    result = pd.concat(
        records,
        ignore_index=True,
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

    # The published workbook used the opposite
    # contrast orientation from this project.
    #
    # Preserve the original published value above,
    # and create a separate reoriented copy.
    result[
        "published_log2fc_treatment_vs_control"
    ] = -result[
        "published_log2fc"
    ]

    result[
        "published_direction_treatment_vs_control"
    ] = np.where(
        result[
            "published_log2fc_treatment_vs_control"
        ]
        > 0,
        "up",
        "down",
    )

    return result


def validate_published_directions(
    persistent_gene_table: pd.DataFrame,
    published_recovery_long: pd.DataFrame,
) -> pd.DataFrame:
    """Compare persistent-direction labels after contrast harmonization."""

    required_persistent = {
        "gene_id",
        "condition",
        "direction",
    }

    missing = (
        required_persistent
        - set(
            persistent_gene_table.columns
        )
    )

    if missing:
        raise CrossStudyError(
            "Persistent-gene table is missing required columns: "
            + ", ".join(
                sorted(
                    missing
                )
            )
        )

    persistent = (
        persistent_gene_table
        .copy()
    )

    persistent[
        "gene_id"
    ] = (
        persistent[
            "gene_id"
        ]
        .astype(str)
        .str.strip()
    )

    validation = (
        persistent
        .merge(
            published_recovery_long,
            on=[
                "gene_id",
                "condition",
            ],
            how="inner",
            validate="one_to_one",
        )
    )

    validation[
        "direction_agrees"
    ] = (
        validation[
            "direction"
        ]
        == validation[
            "published_direction_treatment_vs_control"
        ]
    )

    return validation


def summarize_direction_concordance(
    direction_validation: pd.DataFrame,
) -> pd.DataFrame:
    """Summarize harmonized direction agreement by condition."""

    summary = (
        direction_validation
        .groupby(
            "condition"
        )
        .agg(
            shared_genes=(
                "gene_id",
                "size",
            ),
            direction_agreements=(
                "direction_agrees",
                "sum",
            ),
        )
        .reset_index()
    )

    summary[
        "direction_concordance"
    ] = (
        summary[
            "direction_agreements"
        ]
        / summary[
            "shared_genes"
        ]
    )

    return summary