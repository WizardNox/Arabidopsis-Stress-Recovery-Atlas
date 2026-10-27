"""Integration tests for full pipeline orchestration."""

from arabidopsis_stress_recovery_atlas.pipeline import (
    run_full_pipeline,
)
from arabidopsis_stress_recovery_atlas.validation.pipeline import (
    validate_built_pipeline,
)


def test_full_pipeline_builds_complete_analytical_system(
    tmp_path,
):
    """One pipeline run should create and validate all final artifacts."""

    processed_directory = (
        tmp_path
        / "processed"
    )

    database_path = (
        processed_directory
        / "atlas.db"
    )

    report_directory = (
        tmp_path
        / "reports"
    )

    result = (
        run_full_pipeline(
            processed_directory=(
                processed_directory
            ),
            database_path=(
                database_path
            ),
            report_directory=(
                report_directory
            ),
        )
    )

    counts = (
        result[
            "processed_audit"
        ][
            "table_counts"
        ]
    )

    assert counts == {
        "sample": 57,
        "gene": 26118,
        "gene_condition_response": 235062,
        "gene_persistence_summary": 26118,
        "go_term_enrichment": 1036,
        "published_de_concordance": 354,
    }

    assert (
        result[
            "processed_audit"
        ][
            "gene_count"
        ]
        == 26118
    )

    assert (
        result[
            "processed_audit"
        ][
            "stress_condition_count"
        ]
        == 9
    )

    assert (
        result[
            "processed_audit"
        ][
            "expected_gene_condition_rows"
        ]
        == 235062
    )

    assert (
        result[
            "database_integrity"
        ][
            "integrity_check"
        ]
        == "ok"
    )

    assert (
        result[
            "database_integrity"
        ][
            "foreign_key_violations"
        ]
        == 0
    )

    assert database_path.is_file()

    for path in (
        result[
            "processed_paths"
        ].values()
    ):
        assert path.is_file()

    assert (
        result[
            "reporting"
        ][
            "summary_path"
        ].is_file()
    )

    for path in (
        result[
            "reporting"
        ][
            "figure_paths"
        ].values()
    ):
        assert path.is_file()

    independent_audit = (
        validate_built_pipeline(
            processed_directory=(
                processed_directory
            ),
            database_path=(
                database_path
            ),
        )
    )

    assert (
        independent_audit[
            "database_counts"
        ]
        == counts
    )