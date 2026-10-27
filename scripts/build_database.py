"""Build the SQLite analytical database from canonical processed tables."""

from arabidopsis_stress_recovery_atlas.database.connection import (
    connect_database,
    get_database_path,
)
from arabidopsis_stress_recovery_atlas.database.load import (
    load_processed_tables,
)
from arabidopsis_stress_recovery_atlas.database.schema import (
    create_schema,
    validate_database_integrity,
    validate_schema,
)


def main() -> None:
    """Create a fresh analytical SQLite database."""

    database_path = (
        get_database_path()
    )

    if database_path.exists():
        database_path.unlink()

    connection = (
        connect_database(
            database_path
        )
    )

    try:
        create_schema(
            connection
        )

        validate_schema(
            connection
        )

        counts = (
            load_processed_tables(
                connection
            )
        )

        audit = (
            validate_database_integrity(
                connection
            )
        )

        if (
            audit[
                "integrity_check"
            ]
            != "ok"
        ):
            raise RuntimeError(
                "SQLite integrity check failed."
            )

        if (
            audit[
                "foreign_key_violations"
            ]
            != 0
        ):
            raise RuntimeError(
                "Database contains foreign-key violations."
            )

        print(
            "Analytical database built successfully."
        )

        print()
        print(
            f"Database: {database_path}"
        )

        print()
        print(
            "Loaded rows:"
        )

        for table_name, count in (
            counts.items()
        ):
            print(
                f"  {table_name}: {count}"
            )

        print()
        print(
            "Integrity check:",
            audit[
                "integrity_check"
            ],
        )

        print(
            "Foreign-key violations:",
            audit[
                "foreign_key_violations"
            ],
        )

    finally:
        connection.close()


if __name__ == "__main__":
    main()