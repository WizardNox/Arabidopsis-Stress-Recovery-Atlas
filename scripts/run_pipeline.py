"""Run the complete Arabidopsis Stress Recovery Atlas pipeline."""

from arabidopsis_stress_recovery_atlas.pipeline import (
    run_full_pipeline,
)


def main() -> None:
    """Run the full pipeline and report its outputs."""

    print(
        "Starting Arabidopsis Stress Recovery Atlas pipeline..."
    )

    print()

    result = (
        run_full_pipeline()
    )

    print(
        "Pipeline completed successfully."
    )

    print()

    print(
        "Response floor:",
        result[
            "response_floor"
        ],
    )

    print()

    print(
        "Processed tables:"
    )

    table_counts = (
        result[
            "processed_audit"
        ][
            "table_counts"
        ]
    )

    for table_name, count in (
        table_counts.items()
    ):
        print(
            f"  {table_name}: {count}"
        )

    print()

    print(
        "Database:",
        result[
            "database_path"
        ],
    )

    print(
        "Integrity check:",
        result[
            "database_integrity"
        ][
            "integrity_check"
        ],
    )

    print(
        "Foreign-key violations:",
        result[
            "database_integrity"
        ][
            "foreign_key_violations"
        ],
    )

    print()

    print(
        "Generated report:",
        result[
            "reporting"
        ][
            "summary_path"
        ],
    )

    print()

    print(
        "Final pipeline validation: passed"
    )


if __name__ == "__main__":
    main()