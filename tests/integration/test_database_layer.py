"""Integration tests for the SQLite analytical database."""

import sqlite3

import pytest

from arabidopsis_stress_recovery_atlas.database.connection import (
    connect_database,
)
from arabidopsis_stress_recovery_atlas.database.load import (
    load_processed_tables,
)
from arabidopsis_stress_recovery_atlas.database.schema import (
    create_schema,
    validate_database_integrity,
    validate_schema,
)


@pytest.fixture()
def loaded_database(
    tmp_path,
):
    """Create a temporary database loaded from canonical tables."""

    path = (
        tmp_path
        / "test_atlas.db"
    )

    connection = (
        connect_database(
            path
        )
    )

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

    yield (
        connection,
        counts,
    )

    connection.close()


def test_database_loads_expected_row_counts(
    loaded_database,
):
    """Every canonical processed table should load completely."""

    _, counts = (
        loaded_database
    )

    assert counts == {
        "sample": 57,
        "gene": 26118,
        "gene_condition_response": 235062,
        "gene_persistence_summary": 26118,
        "go_term_enrichment": 1036,
        "published_de_concordance": 354,
    }


def test_database_integrity_and_foreign_keys(
    loaded_database,
):
    """The loaded database should have no relational integrity problems."""

    connection, _ = (
        loaded_database
    )

    audit = (
        validate_database_integrity(
            connection
        )
    )

    assert audit[
        "integrity_check"
    ] == "ok"

    assert audit[
        "foreign_key_violations"
    ] == 0

    with pytest.raises(
        sqlite3.IntegrityError
    ):
        connection.execute(
            """
            INSERT INTO gene_condition_response (
                gene_id,
                condition,
                control_mean,
                endpoint_mean,
                recovery_mean,
                endpoint_shift,
                recovery_shift,
                response_floor,
                endpoint_eligible,
                recovery_above_floor,
                same_direction,
                recovery_outcome,
                persistence_breadth
            )
            VALUES (
                'NOT_A_REAL_GENE',
                'drought',
                0,
                1,
                1,
                1,
                1,
                0.342263,
                1,
                1,
                1,
                'persistent_same_direction',
                1
            )
            """
        )


def test_sql_queries_reproduce_scientific_baseline(
    loaded_database,
):
    """SQL queries should reproduce established project results."""

    connection, _ = (
        loaded_database
    )

    drought_persistent = (
        connection.execute(
            """
            SELECT COUNT(*)
            FROM gene_condition_response
            WHERE condition = 'drought'
              AND recovery_outcome
                  = 'persistent_same_direction'
            """
        )
        .fetchone()[0]
    )

    assert drought_persistent == 1240

    high_breadth = (
        connection.execute(
            """
            SELECT COUNT(*)
            FROM gene_persistence_summary
            WHERE persistence_breadth >= 4
            """
        )
        .fetchone()[0]
    )

    assert high_breadth == 206

    annotated_genes = (
        connection.execute(
            """
            SELECT COUNT(*)
            FROM gene
            WHERE has_araport_annotation = 1
            """
        )
        .fetchone()[0]
    )

    assert annotated_genes == 25971

    significant_go_slim = (
        connection.execute(
            """
            SELECT COUNT(*)
            FROM go_term_enrichment
            WHERE analysis_type = 'go_slim'
              AND is_significant = 1
            """
        )
        .fetchone()[0]
    )

    assert significant_go_slim == 24

    significant_bp = (
        connection.execute(
            """
            SELECT COUNT(*)
            FROM go_term_enrichment
            WHERE analysis_type = 'detailed_bp'
              AND is_significant = 1
            """
        )
        .fetchone()[0]
    )

    assert significant_bp == 45

    published_shared = (
        connection.execute(
            """
            SELECT COUNT(*)
            FROM published_de_concordance
            """
        )
        .fetchone()[0]
    )

    assert published_shared == 354

    direction_agreements = (
        connection.execute(
            """
            SELECT COUNT(*)
            FROM published_de_concordance
            WHERE direction_agrees = 1
            """
        )
        .fetchone()[0]
    )

    assert direction_agreements == 354