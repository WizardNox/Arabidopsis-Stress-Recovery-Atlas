"""Tests for Gene Ontology enrichment calculations."""

import pandas as pd

from arabidopsis_stress_recovery_atlas.analysis.go_enrichment import (
    build_significant_bp_membership,
    calculate_bp_enrichment,
    calculate_go_slim_enrichment,
)


def make_test_go_annotations() -> pd.DataFrame:
    """Create GO annotations for a small artificial enrichment example."""

    records = []

    stress_genes = {
        "gene_a",
        "gene_b",
        "gene_c",
    }

    all_genes = [
        "gene_a",
        "gene_b",
        "gene_c",
        "gene_d",
        "gene_e",
        "gene_f",
        "gene_g",
        "gene_h",
        "gene_i",
        "gene_j",
        "gene_k",
        "gene_l",
    ]

    for gene in all_genes:
        records.append(
            {
                "gene_id": gene,
                "go_id": "GO:0000001",
                "go_term": "common process",
                "aspect": "P",
                "goslim_category": (
                    "other process"
                ),
            }
        )

        if gene in stress_genes:
            records.append(
                {
                    "gene_id": gene,
                    "go_id": "GO:0000002",
                    "go_term": "stress process",
                    "aspect": "P",
                    "goslim_category": (
                        "response to stress"
                    ),
                }
            )

    return pd.DataFrame(records)


def test_go_slim_enrichment_uses_aspect_specific_gene_sets():
    """GO-Slim enrichment should use genes represented in each aspect."""

    annotations = make_test_go_annotations()

    background = {
        f"gene_{letter}"
        for letter in "abcdefghijkl"
    }

    foreground = {
        "gene_a",
        "gene_b",
        "gene_c",
    }

    all_results, _ = (
        calculate_go_slim_enrichment(
            annotations,
            background,
            foreground,
            aspects=("P",),
        )
    )

    process = all_results[
        all_results["go_slim_term"]
        == "response to stress"
    ].iloc[0]

    assert process["foreground_with_term"] == 3
    assert process["foreground_total"] == 3
    assert process["background_with_term"] == 3
    assert process["background_total"] == 12


def test_bp_minimum_term_size_is_applied():
    """Detailed BP terms below the minimum background size should be excluded."""

    annotations = make_test_go_annotations()

    background = {
        f"gene_{letter}"
        for letter in "abcdefghijkl"
    }

    foreground = {
        "gene_a",
        "gene_b",
        "gene_c",
    }

    all_results, _, _, _ = (
        calculate_bp_enrichment(
            annotations,
            background,
            foreground,
            minimum_term_size=4,
        )
    )

    assert set(
        all_results["go_id"]
    ) == {
        "GO:0000001"
    }


def test_bp_root_term_is_excluded():
    """The Biological Process ontology root should not be tested."""

    annotations = make_test_go_annotations()

    root_rows = annotations.iloc[:10].copy()

    root_rows["go_id"] = "GO:0008150"
    root_rows["go_term"] = "biological_process"

    annotations = pd.concat(
        [
            annotations,
            root_rows,
        ],
        ignore_index=True,
    )

    background = {
        f"gene_{letter}"
        for letter in "abcdefghijkl"
    }

    foreground = {
        "gene_a",
        "gene_b",
        "gene_c",
    }

    all_results, _, _, _ = (
        calculate_bp_enrichment(
            annotations,
            background,
            foreground,
            minimum_term_size=4,
        )
    )

    assert "GO:0008150" not in set(
        all_results["go_id"]
    )


def test_significant_membership_keeps_only_foreground_genes():
    """Significant-term membership should be restricted to the foreground."""

    annotations = make_test_go_annotations()

    significant = pd.DataFrame(
        {
            "go_id": [
                "GO:0000002"
            ]
        }
    )

    membership = (
        build_significant_bp_membership(
            annotations,
            significant,
            {
                "gene_a",
                "gene_b",
                "gene_c",
            },
        )
    )

    assert set(
        membership["gene_id"]
    ) == {
        "gene_a",
        "gene_b",
        "gene_c",
    }