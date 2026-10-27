"""SQLite database connection utilities."""

import sqlite3
from pathlib import Path

from arabidopsis_stress_recovery_atlas.config import (
    PROJECT_ROOT,
    load_project_config,
)


class DatabaseConnectionError(RuntimeError):
    """Raised when the configured database cannot be opened."""


def get_database_path() -> Path:
    """Return the configured SQLite database path."""

    config = (
        load_project_config()
    )

    database_config = config.get(
        "database",
        {}
    )

    engine = database_config.get(
        "engine"
    )

    if engine != "sqlite":
        raise DatabaseConnectionError(
            "The configured database engine must be 'sqlite'."
        )

    relative_path = database_config.get(
        "path"
    )

    if not relative_path:
        raise DatabaseConnectionError(
            "No database path is configured."
        )

    return (
        PROJECT_ROOT
        / relative_path
    )


def connect_database(
    path: Path | None = None,
    read_only: bool = False,
) -> sqlite3.Connection:
    """Open the analytical SQLite database."""

    if path is None:
        path = (
            get_database_path()
        )

    path = Path(
        path
    )

    if read_only:
        if not path.is_file():
            raise DatabaseConnectionError(
                f"Database does not exist: {path}"
            )

        uri = (
            path.resolve().as_uri()
            + "?mode=ro"
        )

        connection = sqlite3.connect(
            uri,
            uri=True,
        )

    else:
        path.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        connection = sqlite3.connect(
            path
        )

    connection.row_factory = (
        sqlite3.Row
    )

    connection.execute(
        "PRAGMA foreign_keys = ON"
    )

    connection.execute(
        "PRAGMA busy_timeout = 5000"
    )

    return connection