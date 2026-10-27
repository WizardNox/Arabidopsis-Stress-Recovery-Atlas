"""Validate the built analytical pipeline outputs."""

import sqlite3
from pathlib import Path
from typing import Any

import pandas as pd

from arabidopsis_stress_recovery_atlas.database.connection import (
    connect_database,
)
from arabidopsis_stress_recovery_atlas.database.schema import (
    validate_database_integrity,
    validate_schema,
)


class PipelineValidationError(RuntimeError):
    """Raised when a built pipeline artifact fails validation."""


PROCESSED_FILES = {
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


def _load_processed_tables(
    processed_directory: Path,
) -> dict[str, pd.DataFrame]:
    """Load all required canonical processed tables."""

    tables = {}

    for table_name, filename in (
        PROCESSED_FILES.items()
    ):
        path = (
            processed_directory
            / filename
        )

        if not path.is_file():
            raise PipelineValidationError(
                f"Required processed file is missing: {path}"
            )

        tables[
            table_name
        ] = pd.read_csv(
            path
        )

    return tables


def _require_unique_key(
    table: pd.DataFrame,
    columns: list[str],
    table_name: str,
) -> None:
    """Require a unique relational key."""

    if table.duplicated(
        subset=columns
    ).any():
        raise PipelineValidationError(
            f"{table_name} has duplicate key values for: "
            + ", ".join(
                columns
            )
        )


def validate_processed_layer(
    processed_directory: Path,
) -> dict[str, Any]:
    """Validate relational consistency of canonical CSV outputs."""

    processed_directory = Path(
        processed_directory
    )

    tables = (
        _load_processed_tables(
            processed_directory
        )
    )

    sample = tables[
        "sample"
    ]

    gene = tables[
        "gene"
    ]

    response = tables[
        "gene_condition_response"
    ]

    summary = tables[
        "gene_persistence_summary"
    ]

    go = tables[
        "go_term_enrichment"
    ]

    published = tables[
        "published_de_concordance"
    ]

    # -------------------------------------------------------------
    # Relational keys
    # -------------------------------------------------------------

    _require_unique_key(
        sample,
        [
            "sample_id",
        ],
        "sample",
    )

    _require_unique_key(
        gene,
        [
            "gene_id",
        ],
        "gene",
    )

    _require_unique_key(
        response,
        [
            "gene_id",
            "condition",
        ],
        "gene_condition_response",
    )

    _require_unique_key(
        summary,
        [
            "gene_id",
        ],
        "gene_persistence_summary",
    )

    _require_unique_key(
        published,
        [
            "gene_id",
            "condition",
        ],
        "published_de_concordance",
    )

    # GO-Slim legitimately has no GO ID, so normalize
    # missing IDs only for the purpose of checking identity.
    go_identity = (
        go[
            [
                "analysis_type",
                "go_aspect",
                "go_id",
                "term",
            ]
        ]
        .copy()
    )

    go_identity[
        "go_id"
    ] = (
        go_identity[
            "go_id"
        ]
        .fillna("")
    )

    if go_identity.duplicated().any():
        raise PipelineValidationError(
            "go_term_enrichment contains duplicate analytical terms."
        )

    # -------------------------------------------------------------
    # Gene relationships
    # -------------------------------------------------------------

    gene_ids = set(
        gene[
            "gene_id"
        ]
    )

    response_gene_ids = set(
        response[
            "gene_id"
        ]
    )

    summary_gene_ids = set(
        summary[
            "gene_id"
        ]
    )

    published_gene_ids = set(
        published[
            "gene_id"
        ]
    )

    if response_gene_ids != gene_ids:
        raise PipelineValidationError(
            "gene_condition_response does not contain exactly "
            "the canonical gene universe."
        )

    if summary_gene_ids != gene_ids:
        raise PipelineValidationError(
            "gene_persistence_summary does not contain exactly "
            "the canonical gene universe."
        )

    orphan_published = (
        published_gene_ids
        - gene_ids
    )

    if orphan_published:
        raise PipelineValidationError(
            "published_de_concordance contains genes absent "
            "from the gene dimension."
        )

    # -------------------------------------------------------------
    # Condition coverage
    # -------------------------------------------------------------

    stress_conditions = set(
        sample.loc[
            sample[
                "state"
            ]
            != "control",
            "condition",
        ]
    )

    response_conditions = set(
        response[
            "condition"
        ]
    )

    if (
        response_conditions
        != stress_conditions
    ):
        raise PipelineValidationError(
            "Gene-condition table conditions do not match "
            "the experimental design."
        )

    expected_response_rows = (
        len(
            gene_ids
        )
        * len(
            stress_conditions
        )
    )

    if (
        len(
            response
        )
        != expected_response_rows
    ):
        raise PipelineValidationError(
            "Gene-condition response table is incomplete."
        )

    # -------------------------------------------------------------
    # Persistence consistency
    # -------------------------------------------------------------

    if not (
        summary[
            "persistence_breadth"
        ]
        == summary[
            "n_persistent_conditions"
        ]
    ).all():
        raise PipelineValidationError(
            "Persistence breadth disagrees with "
            "n_persistent_conditions."
        )

    if not (
        summary[
            "n_persistent_up"
        ]
        + summary[
            "n_persistent_down"
        ]
        == summary[
            "n_persistent_conditions"
        ]
    ).all():
        raise PipelineValidationError(
            "Persistent up/down counts do not sum to "
            "the total persistent-condition count."
        )

    # -------------------------------------------------------------
    # Return audit summary
    # -------------------------------------------------------------

    return {
        "table_counts": {
            table_name: len(
                table
            )
            for table_name, table in (
                tables.items()
            )
        },
        "gene_count": (
            len(
                gene_ids
            )
        ),
        "stress_condition_count": (
            len(
                stress_conditions
            )
        ),
        "expected_gene_condition_rows": (
            expected_response_rows
        ),
    }


def _get_database_row_counts(
    connection: sqlite3.Connection,
) -> dict[str, int]:
    """Return database row counts for canonical tables."""

    counts = {}

    for table_name in (
        PROCESSED_FILES
    ):
        counts[
            table_name
        ] = (
            connection.execute(
                f"SELECT COUNT(*) FROM {table_name}"
            )
            .fetchone()[0]
        )

    return counts


def validate_built_pipeline(
    processed_directory: Path,
    database_path: Path,
) -> dict[str, Any]:
    """Validate both processed CSVs and the final SQLite database."""

    processed_audit = (
        validate_processed_layer(
            processed_directory
        )
    )

    database_path = Path(
        database_path
    )

    if not database_path.is_file():
        raise PipelineValidationError(
            f"Database does not exist: {database_path}"
        )

    connection = (
        connect_database(
            database_path,
            read_only=True,
        )
    )

    try:
        validate_schema(
            connection
        )

        database_integrity = (
            validate_database_integrity(
                connection
            )
        )

        database_counts = (
            _get_database_row_counts(
                connection
            )
        )

    finally:
        connection.close()

    if (
        database_integrity[
            "integrity_check"
        ]
        != "ok"
    ):
        raise PipelineValidationError(
            "SQLite integrity check failed."
        )

    if (
        database_integrity[
            "foreign_key_violations"
        ]
        != 0
    ):
        raise PipelineValidationError(
            "SQLite foreign-key validation failed."
        )

    processed_counts = (
        processed_audit[
            "table_counts"
        ]
    )

    if (
        database_counts
        != processed_counts
    ):
        raise PipelineValidationError(
            "Database row counts do not match "
            "canonical processed tables."
        )

    return {
        "processed": (
            processed_audit
        ),
        "database_counts": (
            database_counts
        ),
        "database_integrity": (
            database_integrity
        ),
    }