"""Orchestrate the complete analytical data pipeline."""

from pathlib import Path
from typing import Any

from arabidopsis_stress_recovery_atlas.config import (
    PROJECT_ROOT,
    load_project_config,
)
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
from arabidopsis_stress_recovery_atlas.ingestion.load_expression_data import (
    load_expression_data,
)
from arabidopsis_stress_recovery_atlas.reporting.report import (
    build_reporting_outputs,
)
from arabidopsis_stress_recovery_atlas.transformation.canonical import (
    build_primary_analytical_tables,
    write_primary_analytical_tables,
)
from arabidopsis_stress_recovery_atlas.transformation.canonical_reference import (
    build_reference_analytical_tables,
    write_reference_analytical_tables,
)
from arabidopsis_stress_recovery_atlas.validation.pipeline import (
    validate_built_pipeline,
    validate_processed_layer,
)


def run_full_pipeline(
    processed_directory: Path | None = None,
    database_path: Path | None = None,
    report_directory: Path | None = None,
) -> dict[str, Any]:
    """Run the complete reproducible analytical pipeline."""

    config = (
        load_project_config()
    )

    analysis_config = (
        config[
            "analysis"
        ]
    )

    response_floor_quantile = (
        analysis_config[
            "response_floor"
        ][
            "primary_quantile"
        ]
    )

    high_breadth_minimum = (
        analysis_config[
            "persistence"
        ][
            "high_breadth_min_conditions"
        ]
    )

    enrichment_config = (
        analysis_config[
            "enrichment"
        ]
    )

    enrichment_alpha = (
        enrichment_config[
            "fdr_alpha"
        ]
    )

    minimum_bp_term_size = (
        enrichment_config.get(
            "minimum_bp_term_size",
            10,
        )
    )

    if processed_directory is None:
        processed_directory = (
            PROJECT_ROOT
            / "data"
            / "processed"
        )

    processed_directory = Path(
        processed_directory
    )

    if database_path is None:
        database_path = (
            get_database_path()
        )

    database_path = Path(
        database_path
    )

    if report_directory is None:
        report_directory = (
            PROJECT_ROOT
            / "reports"
            / "generated"
        )

    report_directory = Path(
        report_directory
    )

    expression = (
        load_expression_data()
    )

    primary = (
        build_primary_analytical_tables(
            expression,
            response_floor_quantile=(
                response_floor_quantile
            ),
        )
    )

    reference = (
        build_reference_analytical_tables(
            expression,
            primary[
                "gene_condition_response"
            ],
            primary[
                "gene_persistence_summary"
            ],
            high_breadth_minimum=(
                high_breadth_minimum
            ),
            enrichment_alpha=(
                enrichment_alpha
            ),
            minimum_bp_term_size=(
                minimum_bp_term_size
            ),
        )
    )

    primary_paths = (
        write_primary_analytical_tables(
            primary,
            processed_directory,
        )
    )

    reference_paths = (
        write_reference_analytical_tables(
            reference,
            processed_directory,
        )
    )

    processed_paths = {
        **primary_paths,
        **reference_paths,
    }

    processed_audit = (
        validate_processed_layer(
            processed_directory
        )
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

        database_counts = (
            load_processed_tables(
                connection,
                processed_directory=(
                    processed_directory
                ),
            )
        )

        database_integrity = (
            validate_database_integrity(
                connection
            )
        )

        if (
            database_integrity[
                "integrity_check"
            ]
            != "ok"
        ):
            raise RuntimeError(
                "SQLite integrity check failed."
            )

        if (
            database_integrity[
                "foreign_key_violations"
            ]
            != 0
        ):
            raise RuntimeError(
                "Database contains foreign-key violations."
            )

    finally:
        connection.close()

    final_audit = (
        validate_built_pipeline(
            processed_directory=(
                processed_directory
            ),
            database_path=(
                database_path
            ),
        )
    )

    reporting = (
        build_reporting_outputs(
            database_path=(
                database_path
            ),
            output_directory=(
                report_directory
            ),
        )
    )

    return {
        "response_floor": (
            primary[
                "response_floor"
            ]
        ),
        "processed_paths": (
            processed_paths
        ),
        "processed_audit": (
            processed_audit
        ),
        "database_path": (
            database_path
        ),
        "database_counts": (
            database_counts
        ),
        "database_integrity": (
            database_integrity
        ),
        "final_audit": (
            final_audit
        ),
        "reporting": (
            reporting
        ),
    }