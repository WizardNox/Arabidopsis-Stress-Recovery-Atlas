"""Create and validate the analytical database schema."""

import sqlite3
from pathlib import Path

from arabidopsis_stress_recovery_atlas.config import (
    PROJECT_ROOT,
)


EXPECTED_TABLES = {
    "sample",
    "gene",
    "gene_condition_response",
    "gene_persistence_summary",
    "go_term_enrichment",
    "published_de_concordance",
}


def get_schema_path() -> Path:
    """Return the primary analytical-schema SQL file."""

    return (
        PROJECT_ROOT
        / "sql"
        / "schema"
        / "create_analytical_schema.sql"
    )


def create_schema(
    connection: sqlite3.Connection,
) -> None:
    """Create all analytical tables and indexes."""

    schema_path = (
        get_schema_path()
    )

    sql = schema_path.read_text(
        encoding="utf-8"
    )

    connection.executescript(
        sql
    )


def get_user_tables(
    connection: sqlite3.Connection,
) -> set[str]:
    """Return user-defined SQLite table names."""

    rows = connection.execute(
        """
        SELECT name
        FROM sqlite_master
        WHERE type = 'table'
          AND name NOT LIKE 'sqlite_%'
        """
    ).fetchall()

    return {
        row["name"]
        for row in rows
    }


def validate_schema(
    connection: sqlite3.Connection,
) -> None:
    """Require every expected analytical table to exist."""

    tables = get_user_tables(
        connection
    )

    missing = (
        EXPECTED_TABLES
        - tables
    )

    if missing:
        raise RuntimeError(
            "Database schema is missing tables: "
            + ", ".join(
                sorted(missing)
            )
        )


def validate_database_integrity(
    connection: sqlite3.Connection,
) -> dict[str, object]:
    """Run SQLite integrity and foreign-key checks."""

    integrity_result = (
        connection.execute(
            "PRAGMA integrity_check"
        )
        .fetchone()[0]
    )

    foreign_key_violations = (
        connection.execute(
            "PRAGMA foreign_key_check"
        )
        .fetchall()
    )

    return {
        "integrity_check": (
            integrity_result
        ),
        "foreign_key_violations": (
            len(
                foreign_key_violations
            )
        ),
    }