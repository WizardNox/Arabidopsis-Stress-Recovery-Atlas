"""Data-quality tests for reconstructed sample metadata."""

from arabidopsis_stress_recovery_atlas.config import (
    load_scientific_baseline,
)
from arabidopsis_stress_recovery_atlas.ingestion.load_expression_data import (
    load_expression_data,
)
from arabidopsis_stress_recovery_atlas.transformation.expression import (
    get_expression_columns,
    standardize_gene_id_column,
)
from arabidopsis_stress_recovery_atlas.transformation.metadata import (
    build_sample_manifest,
)
from arabidopsis_stress_recovery_atlas.validation.samples import (
    validate_sample_manifest,
)


def test_reconstructed_sample_manifest_matches_experimental_design():
    """The expression columns should reproduce the known sample design."""

    data = load_expression_data()
    baseline = load_scientific_baseline()

    standardized = standardize_gene_id_column(data)

    sample_columns = get_expression_columns(
        standardized
    )

    manifest = build_sample_manifest(
        sample_columns
    )

    summary = validate_sample_manifest(
        manifest
    )

    assert (
        manifest["sample_id"].tolist()
        == sample_columns
    )

    assert (
        summary["samples"]
        == baseline["dataset"]["expression_samples"]
    )

    assert (
        summary["control_samples"]
        == baseline["experimental_design"]["control_samples"]
    )

    assert (
        summary["endpoint_samples"]
        == baseline["experimental_design"]["endpoint_samples"]
    )

    assert (
        summary["recovery_samples"]
        == baseline["experimental_design"]["recovery_samples"]
    )

    assert (
        summary["stress_conditions"]
        == baseline["experimental_design"]["stress_conditions"]
    )