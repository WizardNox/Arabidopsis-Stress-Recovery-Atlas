"""Regression tests for canonical reference tables."""

from pathlib import Path

import pandas as pd

from arabidopsis_stress_recovery_atlas.ingestion.load_expression_data import (
    load_expression_data,
)
from arabidopsis_stress_recovery_atlas.transformation.canonical import (
    build_primary_analytical_tables,
)
from arabidopsis_stress_recovery_atlas.transformation.canonical_reference import (
    build_reference_analytical_tables,
    write_reference_analytical_tables,
)


def test_reference_analytical_tables_match_stage1(
    tmp_path,
):
    """Canonical reference tables should reproduce Stage 1 results."""

    expression = (
        load_expression_data()
    )

    primary = (
        build_primary_analytical_tables(
            expression,
            response_floor_quantile=0.95,
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
            high_breadth_minimum=4,
            enrichment_alpha=0.05,
            minimum_bp_term_size=10,
        )
    )

    gene = reference[
        "gene"
    ]

    go = reference[
        "go_term_enrichment"
    ]

    concordance = reference[
        "published_de_concordance"
    ]

    # -------------------------------------------------------------
    # Gene dimension
    # -------------------------------------------------------------

    assert len(
        gene
    ) == 26118

    assert gene[
        "gene_id"
    ].nunique() == 26118

    assert not gene[
        "gene_id"
    ].duplicated().any()

    assert int(
        gene[
            "has_araport_annotation"
        ].sum()
    ) == 25971

    assert int(
        (
            ~gene[
                "has_araport_annotation"
            ]
        ).sum()
    ) == 147

    # -------------------------------------------------------------
    # GO enrichment
    # -------------------------------------------------------------

    go_slim = go[
        go[
            "analysis_type"
        ]
        == "go_slim"
    ]

    bp = go[
        go[
            "analysis_type"
        ]
        == "detailed_bp"
    ]

    assert len(
        go_slim
    ) == 98

    assert len(
        bp
    ) == 938

    assert int(
        go_slim[
            "is_significant"
        ].sum()
    ) == 24

    assert int(
        bp[
            "is_significant"
        ].sum()
    ) == 45

    biological_dir = (
        Path(
            "data/interim/biological_interpretation"
        )
    )

    stage1_slim = pd.read_csv(
        biological_dir
        / "high_breadth_go_slim_enrichment_significant.csv"
    )

    stage2_slim_terms = {
        (
            row.go_aspect,
            row.term,
        )
        for row in go_slim[
            go_slim[
                "is_significant"
            ]
        ][
            [
                "go_aspect",
                "term",
            ]
        ].itertuples(
            index=False
        )
    }

    stage1_slim_terms = {
        (
            row.go_aspect,
            row.go_slim_term,
        )
        for row in stage1_slim[
            [
                "go_aspect",
                "go_slim_term",
            ]
        ].itertuples(
            index=False
        )
    }

    assert (
        stage2_slim_terms
        == stage1_slim_terms
    )

    stage1_bp = pd.read_csv(
        biological_dir
        / "high_breadth_bp_enrichment_significant.csv"
    )

    assert set(
        bp.loc[
            bp[
                "is_significant"
            ],
            "go_id",
        ]
    ) == set(
        stage1_bp[
            "go_id"
        ]
    )

    # -------------------------------------------------------------
    # Published-DE concordance
    # -------------------------------------------------------------

    assert len(
        concordance
    ) == 354

    assert concordance[
        "direction_agrees"
    ].all()

    condition_counts = (
        concordance[
            "condition"
        ]
        .value_counts()
        .to_dict()
    )

    assert condition_counts == {
        "drought": 348,
        "cold": 2,
        "heat_drought": 2,
        "high_humidity": 1,
        "high_light": 1,
    }

    stage1_direction = pd.read_csv(
        Path(
            "data/interim/published_de_validation"
        )
        / "published_recovery_direction_validation.csv"
    )

    stage1_pairs = {
        (
            row.gene_id,
            row.condition,
        )
        for row in stage1_direction[
            [
                "gene_id",
                "condition",
            ]
        ].itertuples(
            index=False
        )
    }

    stage2_pairs = {
        (
            row.gene_id,
            row.condition,
        )
        for row in concordance[
            [
                "gene_id",
                "condition",
            ]
        ].itertuples(
            index=False
        )
    }

    assert (
        stage2_pairs
        == stage1_pairs
    )

    # Original published contrast is opposite
    # to this project's treatment-vs-control orientation.
    assert (
        concordance[
            "project_direction"
        ]
        != concordance[
            "published_direction_original"
        ]
    ).all()

    # -------------------------------------------------------------
    # Writing
    # -------------------------------------------------------------

    output_directory = (
        tmp_path
        / "processed"
    )

    paths = (
        write_reference_analytical_tables(
            reference,
            output_directory,
        )
    )

    assert set(
        paths
    ) == {
        "gene",
        "go_term_enrichment",
        "published_de_concordance",
    }

    for path in paths.values():
        assert path.is_file()