"""Regression tests for Stage 1 Gene Ontology enrichment."""

import pandas as pd

from arabidopsis_stress_recovery_atlas.analysis.exploratory import (
    calculate_condition_shifts,
    calculate_control_pairwise_differences,
    calculate_persistence_breadth,
    calculate_response_floor,
    classify_recovery_outcomes,
)
from arabidopsis_stress_recovery_atlas.analysis.go_enrichment import (
    build_significant_bp_membership,
    calculate_bp_enrichment,
    calculate_go_slim_enrichment,
)
from arabidopsis_stress_recovery_atlas.annotation.go import (
    load_go_annotations,
)
from arabidopsis_stress_recovery_atlas.config import (
    PROJECT_ROOT,
    load_project_config,
    load_scientific_baseline,
)
from arabidopsis_stress_recovery_atlas.ingestion.load_expression_data import (
    load_expression_data,
)
from arabidopsis_stress_recovery_atlas.transformation.expression import (
    get_expression_columns,
    get_informative_expression,
    standardize_gene_id_column,
)
from arabidopsis_stress_recovery_atlas.transformation.metadata import (
    build_sample_manifest,
)


def test_go_enrichment_matches_stage1_results():
    """Stage 2 should reproduce Stage 1 GO enrichment."""

    config = load_project_config()
    baseline = load_scientific_baseline()

    expression = get_informative_expression(
        standardize_gene_id_column(
            load_expression_data()
        )
    )

    sample_columns = get_expression_columns(
        expression
    )

    manifest = build_sample_manifest(
        sample_columns
    )

    control_differences = (
        calculate_control_pairwise_differences(
            expression,
            manifest,
        )
    )

    primary_quantile = config[
        "analysis"
    ]["response_floor"]["primary_quantile"]

    response_floor = calculate_response_floor(
        control_differences,
        primary_quantile,
    )

    shifts = calculate_condition_shifts(
        expression,
        manifest,
    )

    classified = classify_recovery_outcomes(
        shifts,
        response_floor,
    )

    breadth = calculate_persistence_breadth(
        classified
    )

    minimum_breadth = config[
        "analysis"
    ]["persistence"][
        "high_breadth_min_conditions"
    ]

    foreground_ids = set(
        breadth.loc[
            breadth["persistence_breadth"]
            >= minimum_breadth,
            "gene_id",
        ]
    )

    background_ids = set(
        expression["gene_id"]
    )

    assert len(foreground_ids) == baseline[
        "functional_analysis"
    ]["high_breadth_foreground_genes"]

    go_annotations = load_go_annotations()

    go_slim_all, go_slim_significant = (
        calculate_go_slim_enrichment(
            go_annotations,
            background_ids,
            foreground_ids,
            alpha=config[
                "analysis"
            ]["enrichment"]["fdr_alpha"],
        )
    )

    (
        bp_all,
        bp_significant,
        bp_background_ids,
        bp_foreground_ids,
    ) = calculate_bp_enrichment(
        go_annotations,
        background_ids,
        foreground_ids,
        minimum_term_size=10,
        root_go_id="GO:0008150",
        alpha=config[
            "analysis"
        ]["enrichment"]["fdr_alpha"],
    )

    membership = (
        build_significant_bp_membership(
            go_annotations,
            bp_significant,
            bp_foreground_ids,
        )
    )

    functional = baseline[
        "functional_analysis"
    ]

    assert len(go_slim_all) == functional[
        "go_slim_terms_tested"
    ]

    assert len(go_slim_significant) == functional[
        "significant_go_slim_terms"
    ]

    assert len(bp_background_ids) == functional[
        "bp_annotated_background_genes"
    ]

    assert len(bp_foreground_ids) == functional[
        "bp_annotated_high_breadth_genes"
    ]

    assert len(bp_all) == functional[
        "detailed_bp_terms_tested"
    ]

    assert len(bp_significant) == functional[
        "significant_detailed_bp_terms"
    ]

    assert len(membership) == functional[
        "bp_term_gene_memberships"
    ]

    assert (
        membership["gene_id"].nunique()
        == functional[
            "bp_unique_foreground_genes"
        ]
    )

    stage1_dir = (
        PROJECT_ROOT
        / "data"
        / "interim"
        / "biological_interpretation"
    )

    stage1_go_slim = pd.read_csv(
        stage1_dir
        / "high_breadth_go_slim_enrichment_significant.csv"
    )

    stage1_bp = pd.read_csv(
        stage1_dir
        / "high_breadth_bp_enrichment_significant.csv"
    )

    assert set(
        go_slim_significant[
            [
                "go_aspect",
                "go_slim_term",
            ]
        ].itertuples(
            index=False,
            name=None,
        )
    ) == set(
        stage1_go_slim[
            [
                "go_aspect",
                "go_slim_term",
            ]
        ].itertuples(
            index=False,
            name=None,
        )
    )

    assert set(
        bp_significant["go_id"]
    ) == set(
        stage1_bp["go_id"]
    )