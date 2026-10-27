"""Build the remaining canonical reference and interpretation tables."""

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from arabidopsis_stress_recovery_atlas.analysis.cross_study import (
    build_published_de_manifest,
    build_published_recovery_long,
    load_published_de_tables,
    validate_published_directions,
)
from arabidopsis_stress_recovery_atlas.analysis.go_enrichment import (
    calculate_bp_enrichment,
    calculate_go_slim_enrichment,
)
from arabidopsis_stress_recovery_atlas.annotation.genes import (
    get_araport_gene_models,
    load_araport_annotations,
    summarize_araport_loci,
)
from arabidopsis_stress_recovery_atlas.annotation.go import (
    load_go_annotations,
)
from arabidopsis_stress_recovery_atlas.transformation.expression import (
    get_informative_expression,
)
from arabidopsis_stress_recovery_atlas.transformation.reference_tables import (
    build_gene_dimension,
    build_go_term_enrichment_table,
    build_published_de_concordance_table,
)


def build_reference_analytical_tables(
    expression: pd.DataFrame,
    gene_condition_response: pd.DataFrame,
    gene_persistence_summary: pd.DataFrame,
    high_breadth_minimum: int = 4,
    enrichment_alpha: float = 0.05,
    minimum_bp_term_size: int = 10,
) -> dict[str, Any]:
    """Build gene, GO-enrichment, and published-DE canonical tables."""

    informative = (
        get_informative_expression(
            expression
        )
    )

    # -----------------------------------------------------------------
    # Gene dimension
    # -----------------------------------------------------------------

    araport_annotations = (
        load_araport_annotations()
    )

    araport_models = (
        get_araport_gene_models(
            araport_annotations
        )
    )

    araport_loci = (
        summarize_araport_loci(
            araport_models
        )
    )

    gene = build_gene_dimension(
        informative,
        araport_loci,
    )

    # -----------------------------------------------------------------
    # GO enrichment
    # -----------------------------------------------------------------

    background_gene_ids = set(
        gene[
            "gene_id"
        ]
    )

    foreground_gene_ids = set(
        gene_persistence_summary.loc[
            gene_persistence_summary[
                "persistence_breadth"
            ]
            >= high_breadth_minimum,
            "gene_id",
        ]
    )

    go_annotations = (
        load_go_annotations()
    )

    go_slim_all, _ = (
        calculate_go_slim_enrichment(
            go_annotations,
            background_gene_ids,
            foreground_gene_ids,
            alpha=enrichment_alpha,
        )
    )

    bp_all, _, _, _ = (
        calculate_bp_enrichment(
            go_annotations,
            background_gene_ids,
            foreground_gene_ids,
            minimum_term_size=(
                minimum_bp_term_size
            ),
            alpha=enrichment_alpha,
        )
    )

    go_term_enrichment = (
        build_go_term_enrichment_table(
            go_slim_all,
            bp_all,
            alpha=enrichment_alpha,
        )
    )

    # -----------------------------------------------------------------
    # Published recovery-DE concordance
    # -----------------------------------------------------------------

    persistent_gene_table = (
        gene_condition_response[
            gene_condition_response[
                "recovery_outcome"
            ]
            == "persistent_same_direction"
        ][
            [
                "gene_id",
                "condition",
                "endpoint_shift",
                "recovery_shift",
                "persistence_breadth",
            ]
        ]
        .copy()
    )

    persistent_gene_table[
        "direction"
    ] = np.where(
        persistent_gene_table[
            "recovery_shift"
        ]
        > 0,
        "up",
        "down",
    )

    manifest = (
        build_published_de_manifest()
    )

    published_tables = (
        load_published_de_tables(
            manifest=manifest
        )
    )

    published_recovery_long = (
        build_published_recovery_long(
            published_tables,
            manifest,
        )
    )

    direction_validation = (
        validate_published_directions(
            persistent_gene_table,
            published_recovery_long,
        )
    )

    published_de_concordance = (
        build_published_de_concordance_table(
            direction_validation
        )
    )

    return {
        "gene": gene,
        "go_term_enrichment": (
            go_term_enrichment
        ),
        "published_de_concordance": (
            published_de_concordance
        ),
    }


def write_reference_analytical_tables(
    tables: dict[str, Any],
    output_directory: Path,
) -> dict[str, Path]:
    """Write canonical reference tables to the processed-data layer."""

    output_directory.mkdir(
        parents=True,
        exist_ok=True,
    )

    required_tables = {
        "gene",
        "go_term_enrichment",
        "published_de_concordance",
    }

    missing = (
        required_tables
        - set(tables)
    )

    if missing:
        raise ValueError(
            "Reference table collection is missing: "
            + ", ".join(
                sorted(missing)
            )
        )

    paths = {
        "gene": (
            output_directory
            / "gene.csv"
        ),
        "go_term_enrichment": (
            output_directory
            / "go_term_enrichment.csv"
        ),
        "published_de_concordance": (
            output_directory
            / "published_de_concordance.csv"
        ),
    }

    for table_name, path in (
        paths.items()
    ):
        tables[
            table_name
        ].to_csv(
            path,
            index=False,
        )

    return paths