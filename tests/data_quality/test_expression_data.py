"""Data-quality tests for the GSE226105 expression matrix."""

from arabidopsis_stress_recovery_atlas.config import (
    load_scientific_baseline,
)
from arabidopsis_stress_recovery_atlas.ingestion.load_expression_data import (
    load_expression_data,
)
from arabidopsis_stress_recovery_atlas.transformation.expression import (
    get_all_zero_expression,
    get_expression_columns,
    get_informative_expression,
    standardize_gene_id_column,
)
from arabidopsis_stress_recovery_atlas.validation.expression import (
    validate_expression_matrix,
)

########################################################################################################################

def test_expression_matrix_matches_scientific_baseline():
    """The source expression matrix should match canonical expectations."""

    data = load_expression_data()
    baseline = load_scientific_baseline()

    summary = validate_expression_matrix(data)

    assert (
        summary["gene_rows"]
        == baseline["dataset"]["source_gene_rows"]
    )

    assert (
        summary["expression_samples"]
        == baseline["dataset"]["expression_samples"]
    )

    assert (
        summary["informative_genes"]
        == baseline["dataset"]["informative_genes"]
    )

    assert (
        summary["all_zero_genes"]
        == baseline["dataset"]["all_zero_genes"]
    )

    assert summary["missing_expression_values"] == 0

###########################################################################################################################

def test_canonical_expression_transformation_matches_baseline():
    """Canonical transformation should reproduce Stage 1 gene counts."""

    data = load_expression_data()
    baseline = load_scientific_baseline()

    standardized = standardize_gene_id_column(data)

    informative = get_informative_expression(
        standardized
    )

    all_zero = get_all_zero_expression(
        standardized
    )

    expression_columns = get_expression_columns(
        standardized
    )

    assert "gene_id" in standardized.columns

    assert (
        len(expression_columns)
        == baseline["dataset"]["expression_samples"]
    )

    assert (
        len(informative)
        == baseline["dataset"]["informative_genes"]
    )

    assert (
        len(all_zero)
        == baseline["dataset"]["all_zero_genes"]
    )

    assert (
        len(informative) + len(all_zero)
        == baseline["dataset"]["source_gene_rows"]
    )