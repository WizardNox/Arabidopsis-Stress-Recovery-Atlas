"""Integration tests for reproducible reporting."""

from arabidopsis_stress_recovery_atlas.database.connection import (
    connect_database,
)
from arabidopsis_stress_recovery_atlas.database.load import (
    load_processed_tables,
)
from arabidopsis_stress_recovery_atlas.database.schema import (
    create_schema,
)
from arabidopsis_stress_recovery_atlas.reporting.report import (
    build_reporting_outputs,
)


def test_reporting_outputs_are_generated_from_database(
    tmp_path,
):
    """Reporting should reproduce established database results."""

    database_path = (
        tmp_path
        / "atlas.db"
    )

    report_directory = (
        tmp_path
        / "report"
    )

    connection = (
        connect_database(
            database_path
        )
    )

    try:
        create_schema(
            connection
        )

        load_processed_tables(
            connection
        )

    finally:
        connection.close()

    result = (
        build_reporting_outputs(
            database_path=(
                database_path
            ),
            output_directory=(
                report_directory
            ),
        )
    )

    summary = result[
        "summary"
    ]

    assert summary[
        "sample_count"
    ] == 57

    assert summary[
        "gene_count"
    ] == 26118

    assert summary[
        "stress_condition_count"
    ] == 9

    assert summary[
        "gene_condition_rows"
    ] == 235062

    assert summary[
        "high_breadth_genes"
    ] == 206

    assert summary[
        "significant_go_slim_terms"
    ] == 24

    assert summary[
        "significant_detailed_bp_terms"
    ] == 45

    assert summary[
        "published_shared_records"
    ] == 354

    assert summary[
        "published_direction_agreements"
    ] == 354

    tables = result[
        "tables"
    ]

    assert len(
        tables[
            "condition_recovery_summary"
        ]
    ) == 9

    assert len(
        tables[
            "high_breadth_genes"
        ]
    ) == 206

    assert len(
        tables[
            "significant_go_terms"
        ]
    ) == 69

    assert len(
        tables[
            "published_de_summary"
        ]
    ) == 5

    assert result[
        "summary_path"
    ].is_file()

    for path in (
        result[
            "table_paths"
        ].values()
    ):
        assert path.is_file()

    for path in (
        result[
            "figure_paths"
        ].values()
    ):
        assert path.is_file()

    report_text = (
        result[
            "summary_path"
        ].read_text(
            encoding="utf-8"
        )
    )

    assert (
        "same-study methodological concordance"
        in report_text.lower()
    )

    assert (
        "continued transcription"
        in report_text.lower()
    )