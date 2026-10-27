"""Regression tests for cross-community BP gene analysis."""

import pandas as pd
from pandas.testing import (
    assert_frame_equal,
)

from arabidopsis_stress_recovery_atlas.analysis.cross_community import (
    build_community_gene_sets,
    build_cross_community_recovery_matrix,
    build_gene_community_membership,
    build_gene_community_profile,
    select_cross_community_genes,
    summarize_cross_community_recovery,
    summarize_gene_community_breadth,
)
from arabidopsis_stress_recovery_atlas.analysis.go_network import (
    build_bp_overlap_graph,
    calculate_bp_term_overlap,
    detect_bp_communities,
    filter_strong_bp_overlap,
    get_largest_component_graph,
)
from arabidopsis_stress_recovery_atlas.config import (
    PROJECT_ROOT,
    load_project_config,
    load_scientific_baseline,
)


def test_cross_community_analysis_matches_stage1():
    """Stage 2 should reproduce the Stage 1 cross-community gene analysis."""

    config = load_project_config()
    baseline = load_scientific_baseline()

    biological_dir = (
        PROJECT_ROOT
        / "data"
        / "interim"
        / "biological_interpretation"
    )

    exploratory_dir = (
        PROJECT_ROOT
        / "data"
        / "interim"
        / "exploratory_analysis"
    )

    significant_enrichment = pd.read_csv(
        biological_dir
        / "high_breadth_bp_enrichment_significant.csv"
    )

    significant_membership = pd.read_csv(
        biological_dir
        / "significant_bp_gene_term_membership.csv"
    )

    overlap = calculate_bp_term_overlap(
        significant_membership
    )

    threshold = config[
        "analysis"
    ]["go_overlap_network"][
        "primary_overlap_threshold"
    ]

    strong = filter_strong_bp_overlap(
        overlap,
        threshold,
    )

    graph = build_bp_overlap_graph(
        significant_enrichment,
        strong,
    )

    largest = get_largest_component_graph(
        graph
    )

    communities, _ = (
        detect_bp_communities(
            largest
        )
    )

    community_gene_sets = (
        build_community_gene_sets(
            significant_membership,
            communities,
        )
    )

    gene_community_membership = (
        build_gene_community_membership(
            community_gene_sets
        )
    )

    breadth = (
        summarize_gene_community_breadth(
            gene_community_membership
        )
    )

    profile = (
        build_gene_community_profile(
            gene_community_membership
        )
    )

    functional = baseline[
        "functional_analysis"
    ]

    assert breadth[
        "gene_id"
    ].nunique() == functional[
        "bp_community_unique_genes"
    ]

    breadth_counts = (
        breadth[
            "n_communities"
        ]
        .value_counts()
        .sort_index()
        .to_dict()
    )

    expected_breadth = {
        int(key): value
        for key, value in functional[
            "bp_gene_community_breadth"
        ].items()
    }

    assert breadth_counts == expected_breadth

    selected = (
        select_cross_community_genes(
            profile,
            minimum_communities=3,
        )
    )

    assert len(selected) == functional[
        "cross_community_genes_min_3"
    ]

    stage1_selected = pd.read_csv(
        biological_dir
        / "high_breadth_cross_community_genes.csv"
    )

    assert set(
        selected["gene_id"]
    ) == set(
        stage1_selected["gene_id"]
    )

    recovery_categories = pd.read_csv(
        exploratory_dir
        / "gene_recovery_categories.csv",
        index_col=0,
    )

    stress_conditions = list(
        config[
            "experimental_design"
        ]["stress_components"].keys()
    )

    recovery_matrix = (
        build_cross_community_recovery_matrix(
            recovery_categories,
            selected["gene_id"],
            stress_conditions,
        )
    )

    stage1_recovery_matrix = (
        pd.read_csv(
            biological_dir
            / "high_breadth_cross_community_recovery_matrix.csv",
            index_col=0,
        )
    )

    recovery_matrix_sorted = (
        recovery_matrix
        .sort_index()
        .sort_index(
            axis=1
        )
    )

    stage1_matrix_sorted = (
        stage1_recovery_matrix
        .sort_index()
        .sort_index(
            axis=1
        )
    )

    assert_frame_equal(
        recovery_matrix_sorted,
        stage1_matrix_sorted,
        check_dtype=False,
    )

    recovery_summary = (
        summarize_cross_community_recovery(
            recovery_matrix
        )
    )

    stage1_recovery_summary = (
        pd.read_csv(
            biological_dir
            / "high_breadth_cross_community_recovery_summary.csv",
            index_col=0,
        )
    )

    expected_columns = [
        "no_qualifying_endpoint_response",
        "endpoint_responders",
        "persistent_same_direction",
        "recovered",
        "substantial_reversal",
        "endpoint_engagement_percent",
        "persistence_among_responders_percent",
        "recovery_among_responders_percent",
        "reversal_among_responders_percent",
    ]

    recovery_summary = (
        recovery_summary[
            expected_columns
        ]
        .sort_index()
    )

    stage1_recovery_summary = (
        stage1_recovery_summary[
            expected_columns
        ]
        .sort_index()
    )

    assert_frame_equal(
        recovery_summary,
        stage1_recovery_summary,
        check_dtype=False,
        check_exact=False,
        rtol=1e-10,
        atol=1e-10,
    )