"""Build database-ready analytical tables from the primary dataset."""

from pathlib import Path
from typing import Any

import pandas as pd

from arabidopsis_stress_recovery_atlas.analysis.exploratory import (
    calculate_condition_shifts,
    calculate_control_pairwise_differences,
    calculate_persistence_breadth,
    calculate_response_floor,
    classify_recovery_outcomes,
)
from arabidopsis_stress_recovery_atlas.transformation.analytical_tables import (
    build_gene_condition_response,
    build_gene_persistence_summary,
)
from arabidopsis_stress_recovery_atlas.transformation.expression import (
    get_expression_columns,
    get_informative_expression,
)
from arabidopsis_stress_recovery_atlas.transformation.metadata import (
    build_sample_manifest,
)
from arabidopsis_stress_recovery_atlas.validation.samples import (
    validate_sample_manifest,
)


def build_primary_analytical_tables(
    expression: pd.DataFrame,
    response_floor_quantile: float = 0.95,
) -> dict[str, Any]:
    """Build the primary canonical analytical tables.

    The function orchestrates existing tested scientific functions
    and then reshapes their outputs into stable tables.
    """

    informative = (
        get_informative_expression(
            expression
        )
    )

    sample_columns = (
        get_expression_columns(
            informative
        )
    )

    sample_manifest = (
        build_sample_manifest(
            sample_columns
        )
    )

    validate_sample_manifest(
        sample_manifest
    )

    pairwise_differences = (
        calculate_control_pairwise_differences(
            informative,
            sample_manifest,
        )
    )

    response_floor = (
        calculate_response_floor(
            pairwise_differences,
            response_floor_quantile,
        )
    )

    shifts = (
        calculate_condition_shifts(
            informative,
            sample_manifest,
        )
    )

    classified = (
        classify_recovery_outcomes(
            shifts,
            response_floor,
        )
    )

    persistence_breadth = (
        calculate_persistence_breadth(
            classified
        )
    )

    gene_condition_response = (
        build_gene_condition_response(
            classified,
            persistence_breadth,
            response_floor,
        )
    )

    gene_persistence_summary = (
        build_gene_persistence_summary(
            classified,
            persistence_breadth,
        )
    )

    return {
        "sample": (
            sample_manifest.copy()
        ),
        "gene_condition_response": (
            gene_condition_response
        ),
        "gene_persistence_summary": (
            gene_persistence_summary
        ),
        "response_floor": (
            float(
                response_floor
            )
        ),
    }


def write_primary_analytical_tables(
    tables: dict[str, Any],
    output_directory: Path,
) -> dict[str, Path]:
    """Write canonical analytical tables to the processed-data layer."""

    output_directory.mkdir(
        parents=True,
        exist_ok=True,
    )

    required_tables = {
        "sample",
        "gene_condition_response",
        "gene_persistence_summary",
    }

    missing_tables = (
        required_tables
        - set(tables)
    )

    if missing_tables:
        raise ValueError(
            "Canonical table collection is missing: "
            + ", ".join(
                sorted(
                    missing_tables
                )
            )
        )

    output_paths = {
        "sample": (
            output_directory
            / "sample.csv"
        ),
        "gene_condition_response": (
            output_directory
            / "gene_condition_response.csv"
        ),
        "gene_persistence_summary": (
            output_directory
            / "gene_persistence_summary.csv"
        ),
    }

    for table_name, path in (
        output_paths.items()
    ):
        tables[
            table_name
        ].to_csv(
            path,
            index=False,
        )

    return output_paths