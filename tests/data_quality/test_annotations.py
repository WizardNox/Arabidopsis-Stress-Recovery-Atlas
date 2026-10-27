"""Data-quality tests for Arabidopsis annotation resources."""

from arabidopsis_stress_recovery_atlas.annotation.genes import (
    get_araport_gene_models,
    load_araport_annotations,
    summarize_araport_loci,
)
from arabidopsis_stress_recovery_atlas.annotation.go import (
    get_go_term_memberships,
    get_goslim_memberships,
    load_go_annotations,
)
from arabidopsis_stress_recovery_atlas.ingestion.load_expression_data import (
    load_expression_data,
)
from arabidopsis_stress_recovery_atlas.transformation.expression import (
    get_informative_expression,
    standardize_gene_id_column,
)


def test_araport_annotation_source_matches_expected_structure():
    """Araport11 should retain the audited source dimensions."""

    annotations = load_araport_annotations()

    models = get_araport_gene_models(
        annotations
    )

    loci = summarize_araport_loci(
        models
    )

    assert len(annotations) == 64475
    assert len(models) == 59782
    assert len(loci) == 38208

    assert loci["gene_id"].is_unique


def test_araport_covers_expected_informative_genes():
    """Araport11 should annotate nearly all informative expression loci."""

    expression = get_informative_expression(
        standardize_gene_id_column(
            load_expression_data()
        )
    )

    annotations = summarize_araport_loci(
        get_araport_gene_models(
            load_araport_annotations()
        )
    )

    expression_genes = set(
        expression["gene_id"].astype(str)
    )

    annotation_genes = set(
        annotations["gene_id"].astype(str)
    )

    overlap = (
        expression_genes
        & annotation_genes
    )

    missing = (
        expression_genes
        - annotation_genes
    )

    assert len(expression_genes) == 26118
    assert len(overlap) == 25971
    assert len(missing) == 147


def test_go_annotation_source_matches_expected_structure():
    """The GO source should retain the audited annotation dimensions."""

    annotations = load_go_annotations()

    memberships = get_go_term_memberships(
        annotations
    )

    bp_memberships = get_go_term_memberships(
        annotations,
        aspect="P",
    )

    goslim = get_goslim_memberships(
        annotations
    )

    bp_goslim = get_goslim_memberships(
        annotations,
        aspect="P",
    )

    assert len(annotations) == 337091

    assert (
        annotations["gene_id"].nunique()
        == 31036
    )

    assert (
        annotations["go_id"].nunique()
        == 7265
    )

    assert len(memberships) == 165764
    assert len(bp_memberships) == 55622

    assert len(goslim) == 188875
    assert len(bp_goslim) == 82149

    assert set(
        annotations["aspect"]
        .dropna()
        .unique()
    ) == {
        "C",
        "F",
        "P",
    }


def test_go_coverage_of_informative_expression_genes():
    """GO coverage should match the audited informative-gene universe."""

    expression = get_informative_expression(
        standardize_gene_id_column(
            load_expression_data()
        )
    )

    annotations = load_go_annotations()

    all_go = get_go_term_memberships(
        annotations
    )

    bp_go = get_go_term_memberships(
        annotations,
        aspect="P",
    )

    expression_genes = set(
        expression["gene_id"].astype(str)
    )

    go_genes = set(
        all_go["gene_id"].astype(str)
    )

    bp_genes = set(
        bp_go["gene_id"].astype(str)
    )

    assert (
        len(expression_genes & go_genes)
        == 25918
    )

    assert (
        len(expression_genes & bp_genes)
        == 20524
    )