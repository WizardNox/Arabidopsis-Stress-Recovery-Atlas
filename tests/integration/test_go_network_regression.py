"""Regression tests for the significant-BP overlap network."""

import pandas as pd
import pytest

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


def test_bp_overlap_network_matches_stage1():
    """Stage 2 should reproduce the Stage 1 significant-BP network."""

    baseline = load_scientific_baseline()
    config = load_project_config()

    base = (
        PROJECT_ROOT
        / "data"
        / "interim"
        / "biological_interpretation"
    )

    significant_enrichment = pd.read_csv(
        base
        / "high_breadth_bp_enrichment_significant.csv"
    )

    membership = pd.read_csv(
        base
        / "significant_bp_gene_term_membership.csv"
    )

    overlap = calculate_bp_term_overlap(
        membership
    )

    functional = baseline[
        "functional_analysis"
    ]

    # There are 45 significant BP terms.
    # Every unique pair of terms is compared:
    #
    # 45 × 44 / 2 = 990 pairs.
    assert len(overlap) == functional[
        "significant_bp_term_pairs"
    ]

    # Stage 1 found 562 term pairs with no shared
    # foreground genes.
    assert int(
        (overlap["intersection"] == 0).sum()
    ) == functional[
        "bp_pairs_without_shared_genes"
    ]

    # Therefore 428 pairs shared at least one
    # foreground gene.
    assert int(
        (overlap["intersection"] > 0).sum()
    ) == functional[
        "bp_pairs_with_shared_genes"
    ]

    # Recover the primary overlap-coefficient threshold
    # from project configuration.
    threshold = config[
        "analysis"
    ]["go_overlap_network"][
        "primary_overlap_threshold"
    ]

    strong = filter_strong_bp_overlap(
        overlap,
        threshold,
    )

    # At the primary threshold of 0.50,
    # Stage 1 produced 88 strong term-term links.
    assert len(strong) == functional[
        "primary_bp_strong_edges"
    ]

    # Compare the exact Stage 2 strong relationships
    # with the strong relationships saved in Stage 1.
    stage1_strong = pd.read_csv(
        base
        / "significant_bp_term_overlap_strong.csv"
    )

    stage2_edges = {
        frozenset(
            (
                row.go_id_a,
                row.go_id_b,
            )
        )
        for row in strong.itertuples(
            index=False
        )
    }

    stage1_edges = {
        frozenset(
            (
                row.go_id_a,
                row.go_id_b,
            )
        )
        for row in stage1_strong.itertuples(
            index=False
        )
    }

    assert stage2_edges == stage1_edges

    # Count how many of the 45 significant BP terms
    # participate in at least one strong relationship.
    strong_term_ids = (
        set(strong["go_id_a"])
        | set(strong["go_id_b"])
    )

    assert len(strong_term_ids) == functional[
        "bp_terms_in_strong_edges"
    ]

    # Build the complete graph.
    #
    # All 45 significant BP terms remain as nodes,
    # including terms that have no strong edges.
    graph = build_bp_overlap_graph(
        significant_enrichment,
        strong,
    )

    assert graph.number_of_nodes() == functional[
        "significant_detailed_bp_terms"
    ]

    assert graph.number_of_edges() == functional[
        "primary_bp_strong_edges"
    ]

    # Extract the largest connected component.
    largest = get_largest_component_graph(
        graph
    )

    assert largest.number_of_nodes() == functional[
        "largest_bp_component_terms"
    ]

    assert largest.number_of_edges() == functional[
        "largest_bp_component_edges"
    ]

    # Reproduce the Stage 1 greedy-modularity
    # community analysis on the largest component.
    communities, modularity_score = (
        detect_bp_communities(
            largest
        )
    )

    # Stage 1 divided the 35-term largest component
    # into four descriptive communities.
    assert len(communities) == 4

    # Stage 1 modularity was approximately 0.426.
    assert modularity_score == pytest.approx(
        0.426,
        abs=0.001,
    )