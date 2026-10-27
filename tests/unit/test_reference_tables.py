"""Tests for canonical reference and interpretation tables."""

import pandas as pd

from arabidopsis_stress_recovery_atlas.transformation.reference_tables import (
    build_gene_dimension,
    build_go_term_enrichment_table,
    build_published_de_concordance_table,
)


def test_gene_dimension_preserves_unannotated_genes():
    """Missing annotation must not remove informative genes."""

    expression = pd.DataFrame(
        {
            "gene_id": [
                "g1",
                "g2",
                "g3",
            ],
            "sample_1": [
                1.0,
                2.0,
                3.0,
            ],
        }
    )

    annotation = pd.DataFrame(
        {
            "gene_id": [
                "g1",
                "g2",
            ],
            "transcript_count": [
                1,
                2,
            ],
            "short_description": [
                "gene one",
                "gene two",
            ],
        }
    )

    result = build_gene_dimension(
        expression,
        annotation,
    )

    assert len(
        result
    ) == 3

    assert result[
        "gene_id"
    ].nunique() == 3

    g3 = result[
        result[
            "gene_id"
        ]
        == "g3"
    ].iloc[0]

    assert not g3[
        "has_araport_annotation"
    ]

    assert pd.isna(
        g3[
            "short_description"
        ]
    )


def test_go_enrichment_table_preserves_fdr_meaning():
    """GO-Slim and detailed BP should preserve their own FDR columns."""

    go_slim = pd.DataFrame(
        {
            "go_aspect": [
                "P",
            ],
            "go_slim_term": [
                "response to stress",
            ],
            "foreground_with_term": [
                10,
            ],
            "foreground_total": [
                20,
            ],
            "background_with_term": [
                100,
            ],
            "background_total": [
                1000,
            ],
            "foreground_fraction": [
                0.5,
            ],
            "background_fraction": [
                0.1,
            ],
            "enrichment_ratio": [
                5.0,
            ],
            "odds_ratio": [
                9.0,
            ],
            "p_value": [
                0.001,
            ],
            "fdr_global": [
                0.01,
            ],
            "fdr_within_aspect": [
                0.005,
            ],
        }
    )

    bp = pd.DataFrame(
        {
            "go_id": [
                "GO:1234567",
            ],
            "go_term": [
                "example process",
            ],
            "foreground_with_term": [
                8,
            ],
            "foreground_total": [
                20,
            ],
            "background_with_term": [
                80,
            ],
            "background_total": [
                1000,
            ],
            "foreground_fraction": [
                0.4,
            ],
            "background_fraction": [
                0.08,
            ],
            "enrichment_ratio": [
                5.0,
            ],
            "odds_ratio": [
                7.0,
            ],
            "p_value": [
                0.002,
            ],
            "fdr_bh": [
                0.02,
            ],
        }
    )

    result = (
        build_go_term_enrichment_table(
            go_slim,
            bp,
        )
    )

    slim_row = result[
        result[
            "analysis_type"
        ]
        == "go_slim"
    ].iloc[0]

    bp_row = result[
        result[
            "analysis_type"
        ]
        == "detailed_bp"
    ].iloc[0]

    assert slim_row[
        "fdr_primary"
    ] == 0.01

    assert slim_row[
        "fdr_secondary"
    ] == 0.005

    assert bp_row[
        "fdr_primary"
    ] == 0.02

    assert pd.isna(
        bp_row[
            "fdr_secondary"
        ]
    )


def test_published_concordance_preserves_both_contrasts():
    """Original and harmonized published fold changes must both remain."""

    validation = pd.DataFrame(
        {
            "gene_id": [
                "g1",
            ],
            "condition": [
                "drought",
            ],
            "direction": [
                "up",
            ],
            "endpoint_shift": [
                2.0,
            ],
            "recovery_shift": [
                1.0,
            ],
            "persistence_breadth": [
                4,
            ],
            "published_log2fc": [
                -1.5,
            ],
            "published_adj_p": [
                0.01,
            ],
            "published_direction": [
                "down",
            ],
            "published_log2fc_treatment_vs_control": [
                1.5,
            ],
            "published_direction_treatment_vs_control": [
                "up",
            ],
            "direction_agrees": [
                True,
            ],
        }
    )

    result = (
        build_published_de_concordance_table(
            validation
        )
    )

    row = result.iloc[0]

    assert row[
        "published_log2fc_original"
    ] == -1.5

    assert row[
        "published_direction_original"
    ] == "down"

    assert row[
        "published_log2fc_treatment_vs_control"
    ] == 1.5

    assert row[
        "published_direction_treatment_vs_control"
    ] == "up"

    assert row[
        "direction_agrees"
    ]