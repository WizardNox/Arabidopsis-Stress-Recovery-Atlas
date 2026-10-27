"""Validate the built analytical pipeline artifacts."""

from arabidopsis_stress_recovery_atlas.config import (
    PROJECT_ROOT,
)
from arabidopsis_stress_recovery_atlas.database.connection import (
    get_database_path,
)
from arabidopsis_stress_recovery_atlas.validation.pipeline import (
    validate_built_pipeline,
)


def main() -> None:
    """Run an independent audit of processed tables and database."""

    processed_directory = (
        PROJECT_ROOT
        / "data"
        / "processed"
    )

    database_path = (
        get_database_path()
    )

    print(
        "Validating analytical pipeline..."
    )

    print()

    audit = (
        validate_built_pipeline(
            processed_directory=(
                processed_directory
            ),
            database_path=(
                database_path
            ),
        )
    )

    print(
        "Processed tables:"
    )

    for table_name, count in (
        audit[
            "processed"
        ][
            "table_counts"
        ].items()
    ):
        print(
            f"  {table_name}: {count}"
        )

    print()

    print(
        "Canonical genes:",
        audit[
            "processed"
        ][
            "gene_count"
        ],
    )

    print(
        "Stress conditions:",
        audit[
            "processed"
        ][
            "stress_condition_count"
        ],
    )

    print(
        "Expected gene-condition rows:",
        audit[
            "processed"
        ][
            "expected_gene_condition_rows"
        ],
    )

    print()

    print(
        "Integrity check:",
        audit[
            "database_integrity"
        ][
            "integrity_check"
        ],
    )

    print(
        "Foreign-key violations:",
        audit[
            "database_integrity"
        ][
            "foreign_key_violations"
        ],
    )

    print()

    print(
        "Pipeline validation passed."
    )


if __name__ == "__main__":
    main()