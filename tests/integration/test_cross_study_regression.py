"""Regression tests for published recovery-DE concordance."""

import pandas as pd
import pytest
from pandas.testing import (
    assert_frame_equal,
)

from arabidopsis_stress_recovery_atlas.analysis.cross_study import (
    build_persistent_sets_from_recovery_categories,
    build_published_condition_counts,
    build_published_de_manifest,
    build_published_recovery_long,
    build_published_recovery_sets,
    calculate_published_overlap_enrichment,
    compare_persistent_to_published,
    load_published_de_tables,
    summarize_direction_concordance,
    summarize_published_de_workbook,
    validate_published_directions,
)
from arabidopsis_stress_recovery_atlas.config import (
    PROJECT_ROOT,
    load_scientific_baseline,
)


def test_published_de_concordance_matches_stage1():
    """Stage 2 should reproduce the Stage 1 published-DE comparison."""

    baseline = (
        load_scientific_baseline()
    )

    published_baseline = baseline[
        "published_de_concordance"
    ]

    manifest = (
        build_published_de_manifest()
    )

    tables = (
        load_published_de_tables(
            manifest=manifest
        )
    )

    workbook_summary = (
        summarize_published_de_workbook(
            tables,
            manifest,
        )
    )

    assert len(
        workbook_summary
    ) == 14

    assert (
        workbook_summary[
            "n_duplicate_genes"
        ]
        == 0
    ).all()

    assert (
        workbook_summary[
            "n_missing_values"
        ]
        == 0
    ).all()

    published_counts = (
        build_published_condition_counts(
            workbook_summary
        )
    )

    assert int(
        published_counts[
            "recovery_sheet_present"
        ].sum()
    ) == published_baseline[
        "available_recovery_conditions"
    ]

    published_de_dir = (
        PROJECT_ROOT
        / "data"
        / "interim"
        / "published_de_validation"
    )

    stage1_counts = pd.read_csv(
        published_de_dir
        / "published_de_condition_counts.csv"
    )

    assert_frame_equal(
        published_counts[
            stage1_counts.columns
        ]
        .sort_values(
            "condition"
        )
        .reset_index(
            drop=True
        ),
        stage1_counts
        .sort_values(
            "condition"
        )
        .reset_index(
            drop=True
        ),
        check_dtype=False,
    )

    exploratory_dir = (
        PROJECT_ROOT
        / "data"
        / "interim"
        / "exploratory_analysis"
    )

    recovery_categories = (
        pd.read_csv(
            exploratory_dir
            / "gene_recovery_categories.csv"
        )
    )

    persistent_sets = (
        build_persistent_sets_from_recovery_categories(
            recovery_categories
        )
    )

    conditions = [
        column
        for column
        in recovery_categories.columns
        if column
        != "gene_id"
    ]

    published_sets = (
        build_published_recovery_sets(
            tables,
            manifest,
            conditions,
        )
    )

    overlap_summary = (
        compare_persistent_to_published(
            persistent_sets,
            published_sets,
        )
    )

    stage1_overlap_summary = (
        pd.read_csv(
            published_de_dir
            / "persistent_vs_published_recovery_degs.csv"
        )
    )

    assert_frame_equal(
        overlap_summary
        .sort_values(
            "condition"
        )
        .reset_index(
            drop=True
        ),
        stage1_overlap_summary
        .sort_values(
            "condition"
        )
        .reset_index(
            drop=True
        ),
        check_dtype=False,
        check_exact=False,
        rtol=1e-12,
        atol=1e-12,
    )

    informative_gene_ids = set(
        recovery_categories[
            "gene_id"
        ]
        .astype(str)
        .str.strip()
    )

    enrichment = (
        calculate_published_overlap_enrichment(
            persistent_sets,
            published_sets,
            informative_gene_ids,
        )
    )

    stage1_enrichment = (
        pd.read_csv(
            published_de_dir
            / "published_recovery_overlap_enrichment.csv"
        )
    )

    assert set(
        enrichment[
            "condition"
        ]
    ) == set(
        stage1_enrichment[
            "condition"
        ]
    )

    numeric_columns = [
        "our_persistent_genes",
        "published_recovery_degs",
        "shared_genes",
        "expected_overlap_by_chance",
        "overlap_enrichment",
        "odds_ratio",
        "p_value",
        "fdr_bh",
    ]

    enrichment_sorted = (
        enrichment
        .sort_values(
            "condition"
        )
        .reset_index(
            drop=True
        )
    )

    stage1_enrichment_sorted = (
        stage1_enrichment
        .sort_values(
            "condition"
        )
        .reset_index(
            drop=True
        )
    )

    assert_frame_equal(
        enrichment_sorted[
            numeric_columns
        ],
        stage1_enrichment_sorted[
            numeric_columns
        ],
        check_dtype=False,
        check_exact=False,
        rtol=1e-10,
        atol=1e-300,
    )

    drought = (
        enrichment[
            enrichment[
                "condition"
            ]
            == "drought"
        ]
        .iloc[0]
    )

    drought_baseline = (
        published_baseline[
            "drought"
        ]
    )

    assert drought[
        "our_persistent_genes"
    ] == drought_baseline[
        "exploratory_persistent_genes"
    ]

    assert drought[
        "published_recovery_degs"
    ] == drought_baseline[
        "published_recovery_degs"
    ]

    assert drought[
        "shared_genes"
    ] == drought_baseline[
        "shared_genes"
    ]

    assert drought[
        "overlap_enrichment"
    ] == pytest.approx(
        drought_baseline[
            "approximate_overlap_enrichment"
        ],
        rel=0.01,
    )

    assert drought[
        "odds_ratio"
    ] == pytest.approx(
        drought_baseline[
            "approximate_odds_ratio"
        ],
        rel=0.01,
    )

    published_long = (
        build_published_recovery_long(
            tables,
            manifest,
        )
    )

    assert len(
        published_long
    ) == published_baseline[
        "published_recovery_gene_condition_records"
    ]

    persistent_gene_table = (
        pd.read_csv(
            exploratory_dir
            / "persistent_gene_table.csv"
        )
    )

    direction_validation = (
        validate_published_directions(
            persistent_gene_table,
            published_long,
        )
    )

    assert len(
        direction_validation
    ) == published_baseline[
        "shared_gene_condition_records"
    ]

    assert int(
        direction_validation[
            "direction_agrees"
        ].sum()
    ) == published_baseline[
        "direction_agreement_after_harmonization"
    ]

    assert direction_validation[
        "direction_agrees"
    ].all()

    # The source contrast orientation is opposite
    # to this project's treatment-vs-control convention.
    assert (
        direction_validation[
            "direction"
        ]
        != direction_validation[
            "published_direction"
        ]
    ).all()

    direction_summary = (
        summarize_direction_concordance(
            direction_validation
        )
    )

    stage1_direction_summary = (
        pd.read_csv(
            published_de_dir
            / "published_recovery_direction_summary.csv"
        )
    )

    assert_frame_equal(
        direction_summary
        .sort_values(
            "condition"
        )
        .reset_index(
            drop=True
        ),
        stage1_direction_summary
        .sort_values(
            "condition"
        )
        .reset_index(
            drop=True
        ),
        check_dtype=False,
        check_exact=False,
        rtol=1e-12,
        atol=1e-12,
    )

    assert (
        direction_summary[
            "direction_concordance"
        ]
        == 1.0
    ).all()