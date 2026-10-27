"""Query reproducible scientific results from the analytical database."""

import sqlite3
from typing import Any

import pandas as pd


def _scalar(
    connection: sqlite3.Connection,
    sql: str,
    parameters: tuple = (),
):
    """Return the first value from a one-row SQL query."""

    row = connection.execute(
        sql,
        parameters,
    ).fetchone()

    if row is None:
        raise RuntimeError(
            "Reporting query returned no result."
        )

    return row[0]


def get_project_summary(
    connection: sqlite3.Connection,
) -> dict[str, Any]:
    """Return high-level metrics describing the completed analysis."""

    response_floor_count = _scalar(
        connection,
        """
        SELECT COUNT(DISTINCT response_floor)
        FROM gene_condition_response
        """,
    )

    if response_floor_count != 1:
        raise RuntimeError(
            "Expected exactly one primary response floor "
            "in gene_condition_response."
        )

    response_floor = _scalar(
        connection,
        """
        SELECT MIN(response_floor)
        FROM gene_condition_response
        """,
    )

    published_records = _scalar(
        connection,
        """
        SELECT COUNT(*)
        FROM published_de_concordance
        """,
    )

    direction_agreements = _scalar(
        connection,
        """
        SELECT COUNT(*)
        FROM published_de_concordance
        WHERE direction_agrees = 1
        """,
    )

    return {
        "sample_count": _scalar(
            connection,
            """
            SELECT COUNT(*)
            FROM sample
            """,
        ),
        "gene_count": _scalar(
            connection,
            """
            SELECT COUNT(*)
            FROM gene
            """,
        ),
        "stress_condition_count": _scalar(
            connection,
            """
            SELECT COUNT(DISTINCT condition)
            FROM gene_condition_response
            """,
        ),
        "gene_condition_rows": _scalar(
            connection,
            """
            SELECT COUNT(*)
            FROM gene_condition_response
            """,
        ),
        "response_floor": float(
            response_floor
        ),
        "high_breadth_genes": _scalar(
            connection,
            """
            SELECT COUNT(*)
            FROM gene_persistence_summary
            WHERE persistence_breadth >= 4
            """,
        ),
        "significant_go_slim_terms": _scalar(
            connection,
            """
            SELECT COUNT(*)
            FROM go_term_enrichment
            WHERE analysis_type = 'go_slim'
              AND is_significant = 1
            """,
        ),
        "significant_detailed_bp_terms": _scalar(
            connection,
            """
            SELECT COUNT(*)
            FROM go_term_enrichment
            WHERE analysis_type = 'detailed_bp'
              AND is_significant = 1
            """,
        ),
        "published_shared_records": (
            published_records
        ),
        "published_direction_agreements": (
            direction_agreements
        ),
    }


def get_condition_recovery_summary(
    connection: sqlite3.Connection,
) -> pd.DataFrame:
    """Return recovery-category counts for each stress condition."""

    sql = """
        SELECT
            condition,

            SUM(
                CASE
                    WHEN recovery_outcome
                        = 'persistent_same_direction'
                    THEN 1
                    ELSE 0
                END
            ) AS persistent_genes,

            SUM(
                CASE
                    WHEN recovery_outcome
                        = 'recovered'
                    THEN 1
                    ELSE 0
                END
            ) AS recovered_genes,

            SUM(
                CASE
                    WHEN recovery_outcome
                        = 'substantial_reversal'
                    THEN 1
                    ELSE 0
                END
            ) AS substantial_reversal_genes,

            SUM(
                CASE
                    WHEN recovery_outcome
                        = 'not_endpoint_eligible'
                    THEN 1
                    ELSE 0
                END
            ) AS not_endpoint_eligible_genes

        FROM gene_condition_response

        GROUP BY condition

        ORDER BY
            persistent_genes DESC,
            condition
    """

    return pd.read_sql_query(
        sql,
        connection,
    )


def get_high_breadth_genes(
    connection: sqlite3.Connection,
    minimum_breadth: int = 4,
) -> pd.DataFrame:
    """Return genes persistent across at least a defined number of stresses."""

    sql = """
        SELECT
            s.gene_id,
            s.persistence_breadth,
            s.n_persistent_up,
            s.n_persistent_down,
            g.short_description,
            g.curator_summary,
            g.computational_description

        FROM gene_persistence_summary AS s

        LEFT JOIN gene AS g
            ON s.gene_id = g.gene_id

        WHERE s.persistence_breadth >= ?

        ORDER BY
            s.persistence_breadth DESC,
            s.gene_id
    """

    return pd.read_sql_query(
        sql,
        connection,
        params=(
            minimum_breadth,
        ),
    )


def get_significant_go_terms(
    connection: sqlite3.Connection,
) -> pd.DataFrame:
    """Return all significant GO-Slim and detailed-BP terms."""

    sql = """
        SELECT
            analysis_type,
            go_aspect,
            go_id,
            term,
            foreground_with_term,
            foreground_total,
            background_with_term,
            background_total,
            enrichment_ratio,
            odds_ratio,
            p_value,
            fdr_primary,
            fdr_secondary

        FROM go_term_enrichment

        WHERE is_significant = 1

        ORDER BY
            analysis_type,
            fdr_primary,
            term
    """

    return pd.read_sql_query(
        sql,
        connection,
    )


def get_published_de_summary(
    connection: sqlite3.Connection,
) -> pd.DataFrame:
    """Summarize shared published recovery-DE records by condition."""

    sql = """
        SELECT
            condition,
            COUNT(*) AS shared_gene_records,
            SUM(direction_agrees) AS direction_agreements,
            100.0 * AVG(direction_agrees)
                AS direction_agreement_percent

        FROM published_de_concordance

        GROUP BY condition

        ORDER BY
            shared_gene_records DESC,
            condition
    """

    return pd.read_sql_query(
        sql,
        connection,
    )


def build_reporting_tables(
    connection: sqlite3.Connection,
) -> dict[str, pd.DataFrame]:
    """Build all canonical reporting tables from the database."""

    return {
        "condition_recovery_summary": (
            get_condition_recovery_summary(
                connection
            )
        ),
        "high_breadth_genes": (
            get_high_breadth_genes(
                connection,
                minimum_breadth=4,
            )
        ),
        "significant_go_terms": (
            get_significant_go_terms(
                connection
            )
        ),
        "published_de_summary": (
            get_published_de_summary(
                connection
            )
        ),
    }