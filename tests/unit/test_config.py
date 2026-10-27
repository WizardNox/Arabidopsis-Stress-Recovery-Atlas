"""Tests for project configuration loading."""

import pytest

from arabidopsis_stress_recovery_atlas.config import (
    CONFIG_DIR,
    PROJECT_ROOT,
    load_project_config,
    load_scientific_baseline,
    load_sources_config,
)


def test_project_root_contains_config_directory():
    """The package should resolve the repository configuration directory."""

    assert PROJECT_ROOT.is_dir()
    assert CONFIG_DIR.is_dir()


def test_project_config_loads():
    """The main project configuration should load expected settings."""

    config = load_project_config()

    assert config["project"]["primary_dataset"] == "GSE226105"

    assert (
        config["analysis"]["response_floor"]["primary_quantile"]
        == pytest.approx(0.95)
    )

    assert (
        config["analysis"]["persistence"]["high_breadth_min_conditions"]
        == 4
    )

    assert (
        config["experimental_design"]["condition_labels"]["T"]
        == "heat_drought_high_light"
    )

    assert (
        len(config["experimental_design"]["condition_labels"])
        == 9
    )

def test_sources_config_loads():
    """The source registry should contain the required scientific inputs."""

    config = load_sources_config()

    sources = config["sources"]

    assert "expression_matrix" in sources
    assert "published_de_workbook" in sources
    assert "araport11_functional_descriptions" in sources
    assert "ath_go_goslim" in sources

    assert (
        sources["expression_matrix"]["accession"]
        == "GSE226105"
    )


def test_scientific_baseline_loads():
    """The scientific baseline should contain canonical project results."""

    baseline = load_scientific_baseline()

    assert baseline["dataset"]["expression_samples"] == 57
    assert baseline["dataset"]["informative_genes"] == 26118

    assert (
        baseline["primary_response_floor"]["rlog_units"]
        == pytest.approx(0.342263)
    )

    assert (
        baseline["primary_p95_persistence_by_condition"]["drought"]
        == 1240
    )

    assert (
        baseline["persistent_gene_counts"]["at_least_4_stresses"]
        == 206
    )