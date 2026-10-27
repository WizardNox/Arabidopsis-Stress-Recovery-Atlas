"""Generate reproducible reporting outputs from the analytical database."""

from arabidopsis_stress_recovery_atlas.config import (
    PROJECT_ROOT,
)
from arabidopsis_stress_recovery_atlas.database.connection import (
    get_database_path,
)
from arabidopsis_stress_recovery_atlas.reporting.report import (
    build_reporting_outputs,
)


def main() -> None:
    """Build report tables, figures, and generated scientific summary."""

    database_path = (
        get_database_path()
    )

    output_directory = (
        PROJECT_ROOT
        / "reports"
        / "generated"
    )

    result = (
        build_reporting_outputs(
            database_path=(
                database_path
            ),
            output_directory=(
                output_directory
            ),
        )
    )

    print(
        "Reporting outputs built successfully."
    )

    print()

    print(
        "Summary report:",
        result[
            "summary_path"
        ],
    )

    print()

    print(
        "Reporting tables:"
    )

    for name, path in (
        result[
            "table_paths"
        ].items()
    ):
        print(
            f"  {name}: {path}"
        )

    print()

    print(
        "Figures:"
    )

    for name, path in (
        result[
            "figure_paths"
        ].items()
    ):
        print(
            f"  {name}: {path}"
        )


if __name__ == "__main__":
    main()