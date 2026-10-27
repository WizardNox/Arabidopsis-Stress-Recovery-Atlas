"""Tests for cross-community BP gene analysis."""

import pandas as pd
import pytest

from arabidopsis_stress_recovery_atlas.analysis.cross_community import (
    build_community_gene_sets,
    build_cross_community_recovery_matrix,
    build_gene_community_membership,
    build_gene_community_profile,
    select_cross_community_genes,
    summarize_cross_community_recovery,
    summarize_gene_community_breadth,
)


def make_membership() -> pd.DataFrame:
    """Create artificial gene-term memberships."""

    return pd.DataFrame(
        {
            "gene_id": [
                "g1",
                "g2",
                "g1",
                "g3",
                "g1",
                "g4",
                "g5",
            ],
            "go_id": [
                "GO:A",
                "GO:A",
                "GO:B",
                "GO:B",
                "GO:C",
                "GO:C",
                "GO:D",
            ],
            "go_term": [
                "term A",
                "term A",
                "term B",
                "term B",
                "term C",
                "term C",
                "term D",
            ],
        }
    )


def test_gene_community_breadth_is_calculated():
    """Genes should be counted once per BP community."""

    communities = [
        {
            "GO:A",
        },
        {
            "GO:B",
        },
        {
            "GO:C",
        },
        {
            "GO:D",
        },
    ]

    gene_sets = build_community_gene_sets(
        make_membership(),
        communities,
    )

    membership = (
        build_gene_community_membership(
            gene_sets
        )
    )

    breadth = (
        summarize_gene_community_breadth(
            membership
        )
    )

    g1 = breadth[
        breadth["gene_id"] == "g1"
    ].iloc[0]

    assert g1["n_communities"] == 3


def test_cross_community_selection_uses_minimum_three():
    """Genes represented in at least three communities should be selected."""

    communities = [
        {
            "GO:A",
        },
        {
            "GO:B",
        },
        {
            "GO:C",
        },
        {
            "GO:D",
        },
    ]

    gene_sets = build_community_gene_sets(
        make_membership(),
        communities,
    )

    membership = (
        build_gene_community_membership(
            gene_sets
        )
    )

    profile = (
        build_gene_community_profile(
            membership
        )
    )

    selected = (
        select_cross_community_genes(
            profile,
            minimum_communities=3,
        )
    )

    assert selected[
        "gene_id"
    ].tolist() == [
        "g1"
    ]

    assert selected.iloc[0][
        "n_communities"
    ] == 3


def test_recovery_summary_uses_correct_denominators():
    """Recovery percentages should use the Stage 1 denominators."""

    recovery = pd.DataFrame(
        {
            "cold": [
                "persistent_same_direction",
                "recovered",
                "no_qualifying_endpoint_response",
            ],
            "drought": [
                "persistent_same_direction",
                "persistent_same_direction",
                "substantial_reversal",
            ],
        },
        index=[
            "g1",
            "g2",
            "g3",
        ],
    )

    matrix = (
        build_cross_community_recovery_matrix(
            recovery,
            selected_gene_ids=[
                "g1",
                "g2",
                "g3",
            ],
            stress_conditions=[
                "cold",
                "drought",
            ],
        )
    )

    summary = (
        summarize_cross_community_recovery(
            matrix
        )
    )

    cold = summary.loc[
        "cold"
    ]

    assert cold[
        "endpoint_responders"
    ] == 2

    assert cold[
        "endpoint_engagement_percent"
    ] == pytest.approx(
        2 / 3 * 100
    )

    assert cold[
        "persistence_among_responders_percent"
    ] == pytest.approx(
        50.0
    )

    assert cold[
        "recovery_among_responders_percent"
    ] == pytest.approx(
        50.0
    )

    drought = summary.loc[
        "drought"
    ]

    assert drought[
        "endpoint_responders"
    ] == 3

    assert drought[
        "persistence_among_responders_percent"
    ] == pytest.approx(
        2 / 3 * 100
    )

    assert drought[
        "reversal_among_responders_percent"
    ] == pytest.approx(
        1 / 3 * 100
    )