"""Tests for sample-metadata reconstruction."""

import pytest

from arabidopsis_stress_recovery_atlas.transformation.metadata import (
    SampleMetadataError,
    build_sample_manifest,
    parse_sample_name,
)


def test_control_sample_is_parsed():
    """Control labels should produce control metadata."""

    metadata = parse_sample_name("CONT 2")

    assert metadata["state"] == "control"
    assert metadata["condition"] == "control"
    assert metadata["replicate"] == 2

    assert metadata["cold"] is False
    assert metadata["heat"] is False
    assert metadata["drought"] is False
    assert metadata["high_light"] is False
    assert metadata["high_humidity"] is False


def test_single_stress_sample_is_parsed():
    """Single-stress labels should map to the canonical condition."""

    metadata = parse_sample_name(
        "endpoint_DROUGHT 3"
    )

    assert metadata["state"] == "endpoint"
    assert metadata["condition"] == "drought"
    assert metadata["replicate"] == 3
    assert metadata["drought"] is True
    assert metadata["heat"] is False


def test_combined_stress_sample_is_parsed():
    """Combined-stress labels should preserve all stress components."""

    metadata = parse_sample_name(
        "recovery_HEAT HL 1"
    )

    assert metadata["state"] == "recovery"
    assert metadata["condition"] == "heat_high_light"
    assert metadata["replicate"] == 1
    assert metadata["heat"] is True
    assert metadata["high_light"] is True
    assert metadata["drought"] is False


def test_triple_stress_t_label_is_parsed():
    """The source T label should map to the triple-stress condition."""

    metadata = parse_sample_name(
        "endpoint_T1"
    )

    assert metadata["state"] == "endpoint"
    assert metadata["replicate"] == 1

    assert (
        metadata["condition"]
        == "heat_drought_high_light"
    )

    assert metadata["heat"] is True
    assert metadata["drought"] is True
    assert metadata["high_light"] is True

    assert (
        metadata["condition"]
        == "heat_drought_high_light"
    )

    assert metadata["heat"] is True
    assert metadata["drought"] is True
    assert metadata["high_light"] is True


def test_unknown_sample_label_raises_error():
    """Unknown source labels should fail rather than be guessed."""

    with pytest.raises(SampleMetadataError):
        parse_sample_name(
            "endpoint_UNKNOWN 1"
        )


def test_manifest_preserves_sample_order():
    """Manifest rows should follow the supplied expression-column order."""

    sample_columns = [
        "CONT 1",
        "endpoint_COLD 1",
        "recovery_HL 2",
    ]

    manifest = build_sample_manifest(
        sample_columns
    )

    assert manifest["sample_id"].tolist() == sample_columns