"""Load canonical processed tables into SQLite."""

import sqlite3
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from arabidopsis_stress_recovery_atlas.config import (
    PROJECT_ROOT,
)


TABLE_FILES = {
    "sample": "sample.csv",
    "gene": "gene.csv",
    "gene_condition_response": (
        "gene_condition_response.csv"
    ),
    "gene_persistence_summary": (
        "gene_persistence_summary.csv"
    ),
    "go_term_enrichment": (
        "go_term_enrichment.csv"
    ),
    "published_de_concordance": (
        "published_de_concordance.csv"
    ),
}


BOOLEAN_COLUMNS = {
    "sample": {
        "cold",
        "heat",
        "drought",
        "high_light",
        "high_humidity",
    },
    "gene": {
        "has_araport_annotation",
    },
    "gene_condition_response": {
        "endpoint_eligible",
        "recovery_above_floor",
        "same_direction",
    },
    "go_term_enrichment": {
        "is_significant",
    },
    "published_de_concordance": {
        "direction_agrees",
    },
}


DELETE_ORDER = [
    "published_de_concordance",
    "go_term_enrichment",
    "gene_persistence_summary",
    "gene_condition_response",
    "gene",
    "sample",
]


def _sqlite_value(
    value: Any,
) -> Any:
    """Convert pandas/numpy scalar values into SQLite-compatible values."""

    if value is None:
        return None

    try:
        if pd.isna(
            value
        ):
            return None
    except (
        TypeError,
        ValueError,
    ):
        pass

    if isinstance(
        value,
        np.generic,
    ):
        value = value.item()

    if isinstance(
        value,
        bool,
    ):
        return int(
            value
        )

    return value


def _normalize_boolean_columns(
    table_name: str,
    dataframe: pd.DataFrame,
) -> pd.DataFrame:
    """Normalize canonical boolean columns before SQLite loading."""

    result = dataframe.copy()

    for column in BOOLEAN_COLUMNS.get(
        table_name,
        set(),
    ):
        if column not in result.columns:
            continue

        result[
            column
        ] = (
            result[
                column
            ]
            .map(
                {
                    True: 1,
                    False: 0,
                    1: 1,
                    0: 0,
                    "True": 1,
                    "False": 0,
                    "true": 1,
                    "false": 0,
                }
            )
        )

        if result[
            column
        ].isna().any():
            raise ValueError(
                f"Could not normalize boolean column "
                f"{table_name}.{column}."
            )

        result[
            column
        ] = result[
            column
        ].astype(
            int
        )

    return result


def _database_columns(
    connection: sqlite3.Connection,
    table_name: str,
) -> set[str]:
    """Return column names defined for a SQLite table."""

    rows = connection.execute(
        f"PRAGMA table_info({table_name})"
    ).fetchall()

    return {
        row["name"]
        for row in rows
    }


def load_dataframe(
    connection: sqlite3.Connection,
    table_name: str,
    dataframe: pd.DataFrame,
) -> int:
    """Insert one canonical DataFrame into an existing table."""

    dataframe = (
        _normalize_boolean_columns(
            table_name,
            dataframe,
        )
    )

    database_columns = (
        _database_columns(
            connection,
            table_name,
        )
    )

    unexpected_columns = (
        set(
            dataframe.columns
        )
        - database_columns
    )

    if unexpected_columns:
        raise ValueError(
            f"{table_name} contains columns not present "
            "in the database schema: "
            + ", ".join(
                sorted(
                    unexpected_columns
                )
            )
        )

    columns = list(
        dataframe.columns
    )

    column_sql = ", ".join(
        columns
    )

    placeholders = ", ".join(
        "?"
        for _ in columns
    )

    sql = (
        f"INSERT INTO {table_name} "
        f"({column_sql}) "
        f"VALUES ({placeholders})"
    )

    rows = [
        tuple(
            _sqlite_value(
                value
            )
            for value in row
        )
        for row in dataframe.itertuples(
            index=False,
            name=None,
        )
    ]

    connection.executemany(
        sql,
        rows,
    )

    return len(
        rows
    )


def load_processed_tables(
    connection: sqlite3.Connection,
    processed_directory: Path | None = None,
    clear_existing: bool = True,
) -> dict[str, int]:
    """Load all canonical processed tables into SQLite."""

    if processed_directory is None:
        processed_directory = (
            PROJECT_ROOT
            / "data"
            / "processed"
        )

    processed_directory = Path(
        processed_directory
    )

    counts = {}

    with connection:

        if clear_existing:
            for table_name in (
                DELETE_ORDER
            ):
                connection.execute(
                    f"DELETE FROM {table_name}"
                )

        for table_name, filename in (
            TABLE_FILES.items()
        ):
            path = (
                processed_directory
                / filename
            )

            if not path.is_file():
                raise FileNotFoundError(
                    f"Canonical processed table not found: {path}"
                )

            dataframe = pd.read_csv(
                path
            )

            counts[
                table_name
            ] = load_dataframe(
                connection,
                table_name,
                dataframe,
            )

    connection.execute(
        "ANALYZE"
    )

    return counts