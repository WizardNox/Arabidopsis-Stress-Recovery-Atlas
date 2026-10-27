"""Tests for Araport11 gene annotation transformations."""

import pandas as pd

from arabidopsis_stress_recovery_atlas.annotation.genes import (
    get_araport_gene_models,
    summarize_araport_loci,
)


def test_gene_model_suffixes_are_reduced_to_locus_ids():
    """Transcript-model suffixes should map to canonical AGI loci."""

    annotations = pd.DataFrame(
        {
            "name": [
                "AT1G01020.1",
                "AT1G01020.2",
                "NOT_AN_AGI",
            ],
            "gene_model_type": [
                "protein_coding",
                "protein_coding",
                "unknown",
            ],
            "short_description": [
                "protein A",
                "protein A",
                "other",
            ],
            "Curator_summary": [
                None,
                None,
                None,
            ],
            "Computational_description": [
                None,
                None,
                None,
            ],
        }
    )

    models = get_araport_gene_models(
        annotations
    )

    assert models["gene_id"].tolist() == [
        "AT1G01020",
        "AT1G01020",
    ]


def test_locus_summary_does_not_duplicate_gene_ids():
    """Multiple gene models should collapse to one locus row."""

    annotations = pd.DataFrame(
        {
            "name": [
                "AT1G01020.1",
                "AT1G01020.2",
            ],
            "gene_model_type": [
                "protein_coding",
                "protein_coding",
            ],
            "short_description": [
                "Arv1-like protein",
                "Arv1-like protein",
            ],
            "Curator_summary": [
                None,
                None,
            ],
            "Computational_description": [
                None,
                None,
            ],
        }
    )

    models = get_araport_gene_models(
        annotations
    )

    loci = summarize_araport_loci(
        models
    )

    assert len(loci) == 1
    assert loci.iloc[0]["gene_id"] == "AT1G01020"
    assert loci.iloc[0]["gene_model_count"] == 2
    assert (
        loci.iloc[0]["short_description"]
        == "Arv1-like protein"
    )


def test_distinct_model_descriptions_are_preserved():
    """Different non-null model descriptions should not be discarded."""

    annotations = pd.DataFrame(
        {
            "name": [
                "AT1G00010.1",
                "AT1G00010.2",
            ],
            "gene_model_type": [
                "protein_coding",
                "protein_coding",
            ],
            "short_description": [
                "description A",
                "description B",
            ],
            "Curator_summary": [
                None,
                None,
            ],
            "Computational_description": [
                None,
                None,
            ],
        }
    )

    loci = summarize_araport_loci(
        get_araport_gene_models(
            annotations
        )
    )

    description = loci.iloc[0][
        "short_description"
    ]

    assert "description A" in description
    assert "description B" in description