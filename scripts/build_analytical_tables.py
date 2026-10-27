"""Build all canonical analytical tables."""

from arabidopsis_stress_recovery_atlas.config import (
    PROJECT_ROOT,
    load_project_config,
)
from arabidopsis_stress_recovery_atlas.ingestion.load_expression_data import (
    load_expression_data,
)
from arabidopsis_stress_recovery_atlas.transformation.canonical import (
    build_primary_analytical_tables,
    write_primary_analytical_tables,
)
from arabidopsis_stress_recovery_atlas.transformation.canonical_reference import (
    build_reference_analytical_tables,
    write_reference_analytical_tables,
)


def main() -> None:
    """Build and write the complete canonical processed layer."""

    config = (
        load_project_config()
    )

    analysis_config = config[
        "analysis"
    ]

    response_floor_quantile = (
        analysis_config[
            "response_floor"
        ][
            "primary_quantile"
        ]
    )

    high_breadth_minimum = (
        analysis_config[
            "persistence"
        ][
            "high_breadth_min_conditions"
        ]
    )

    enrichment_config = (
        analysis_config[
            "enrichment"
        ]
    )

    enrichment_alpha = (
        enrichment_config[
            "fdr_alpha"
        ]
    )

    minimum_bp_term_size = (
        enrichment_config.get(
            "minimum_bp_term_size",
            10,
        )
    )

    expression = (
        load_expression_data()
    )

    primary = (
        build_primary_analytical_tables(
            expression,
            response_floor_quantile=(
                response_floor_quantile
            ),
        )
    )

    reference = (
        build_reference_analytical_tables(
            expression,
            primary[
                "gene_condition_response"
            ],
            primary[
                "gene_persistence_summary"
            ],
            high_breadth_minimum=(
                high_breadth_minimum
            ),
            enrichment_alpha=(
                enrichment_alpha
            ),
            minimum_bp_term_size=(
                minimum_bp_term_size
            ),
        )
    )

    output_directory = (
        PROJECT_ROOT
        / "data"
        / "processed"
    )

    primary_paths = (
        write_primary_analytical_tables(
            primary,
            output_directory,
        )
    )

    reference_paths = (
        write_reference_analytical_tables(
            reference,
            output_directory,
        )
    )

    all_paths = {
        **primary_paths,
        **reference_paths,
    }

    print(
        "Canonical analytical layer built successfully."
    )

    print()
    print(
        "Response floor:",
        primary[
            "response_floor"
        ],
    )

    print()
    print(
        "Sample rows:",
        len(
            primary[
                "sample"
            ]
        ),
    )

    print(
        "Gene rows:",
        len(
            reference[
                "gene"
            ]
        ),
    )

    print(
        "Gene-condition rows:",
        len(
            primary[
                "gene_condition_response"
            ]
        ),
    )

    print(
        "Gene-summary rows:",
        len(
            primary[
                "gene_persistence_summary"
            ]
        ),
    )

    print(
        "GO-enrichment rows:",
        len(
            reference[
                "go_term_enrichment"
            ]
        ),
    )

    print(
        "Published-DE concordance rows:",
        len(
            reference[
                "published_de_concordance"
            ]
        ),
    )

    print()
    print(
        "Written files:"
    )

    for name, path in (
        all_paths.items()
    ):
        print(
            f"  {name}: {path}"
        )


if __name__ == "__main__":
    main()