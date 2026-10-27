"""Tests for published recovery-DE concordance."""

import pandas as pd
import pytest

from arabidopsis_stress_recovery_atlas.analysis.cross_study import (
    build_published_de_manifest,
    build_published_recovery_long,
    calculate_published_overlap_enrichment,
    summarize_direction_concordance,
    validate_published_directions,
)


def test_published_manifest_has_expected_design():
    """The workbook manifest should contain 9 endpoint and 5 recovery sheets."""

    manifest = (
        build_published_de_manifest()
    )

    assert len(
        manifest
    ) == 14

    counts = (
        manifest[
            "state"
        ]
        .value_counts()
        .to_dict()
    )

    assert counts == {
        "endpoint": 9,
        "recovery": 5,
    }

    assert set(
        manifest.loc[
            manifest[
                "state"
            ]
            == "recovery",
            "condition",
        ]
    ) == {
        "cold",
        "drought",
        "heat_drought",
        "high_humidity",
        "high_light",
    }


def test_reoriented_published_direction_preserves_original_value():
    """Published fold change should be preserved while its comparison copy is reversed."""

    manifest = pd.DataFrame(
        [
            (
                "Example rec",
                "recovery",
                "drought",
            )
        ],
        columns=[
            "sheet_name",
            "state",
            "condition",
        ],
    )

    tables = {
        "Example rec": pd.DataFrame(
            {
                "gene_ID": [
                    "g1",
                    "g2",
                ],
                "log2(FC)": [
                    -2.0,
                    1.5,
                ],
                "adj-P": [
                    0.01,
                    0.02,
                ],
            }
        )
    }

    result = (
        build_published_recovery_long(
            tables,
            manifest,
        )
    )

    g1 = result[
        result[
            "gene_id"
        ]
        == "g1"
    ].iloc[0]

    assert g1[
        "published_log2fc"
    ] == -2.0

    assert g1[
        "published_direction"
    ] == "down"

    assert g1[
        "published_log2fc_treatment_vs_control"
    ] == 2.0

    assert g1[
        "published_direction_treatment_vs_control"
    ] == "up"


def test_overlap_enrichment_uses_one_sided_fisher_and_bh():
    """Overlap testing should return Fisher statistics and BH FDR."""

    persistent_sets = {
        "condition_a": {
            "g1",
            "g2",
            "g3",
        },
        "condition_b": {
            "g4",
            "g5",
        },
    }

    published_sets = {
        "condition_a": {
            "g1",
            "g2",
        },
        "condition_b": {
            "g4",
        },
    }

    universe = {
        f"g{i}"
        for i in range(
            1,
            21,
        )
    }

    result = (
        calculate_published_overlap_enrichment(
            persistent_sets,
            published_sets,
            universe,
        )
    )

    assert len(
        result
    ) == 2

    assert "fdr_bh" in (
        result.columns
    )

    assert (
        result[
            "shared_genes"
        ].tolist()
        == [
            2,
            1,
        ]
    )

    assert (
        result[
            "fdr_bh"
        ]
        >= result[
            "p_value"
        ]
    ).all()


def test_direction_concordance_uses_reoriented_sign():
    """Direction validation should compare against the harmonized contrast."""

    persistent = pd.DataFrame(
        {
            "gene_id": [
                "g1",
                "g2",
            ],
            "condition": [
                "drought",
                "drought",
            ],
            "direction": [
                "up",
                "down",
            ],
        }
    )

    published = pd.DataFrame(
        {
            "gene_id": [
                "g1",
                "g2",
            ],
            "condition": [
                "drought",
                "drought",
            ],
            "published_log2fc": [
                -2.0,
                1.5,
            ],
            "published_adj_p": [
                0.01,
                0.02,
            ],
            "published_direction": [
                "down",
                "up",
            ],
            "published_log2fc_treatment_vs_control": [
                2.0,
                -1.5,
            ],
            "published_direction_treatment_vs_control": [
                "up",
                "down",
            ],
        }
    )

    validation = (
        validate_published_directions(
            persistent,
            published,
        )
    )

    assert validation[
        "direction_agrees"
    ].all()

    summary = (
        summarize_direction_concordance(
            validation
        )
    )

    assert summary.iloc[0][
        "shared_genes"
    ] == 2

    assert summary.iloc[0][
        "direction_agreements"
    ] == 2

    assert summary.iloc[0][
        "direction_concordance"
    ] == pytest.approx(
        1.0
    )