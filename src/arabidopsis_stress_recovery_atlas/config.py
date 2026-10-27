"""Load and validate project configuration files."""

from pathlib import Path
from typing import Any

import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONFIG_DIR = PROJECT_ROOT / "config"


class ConfigurationError(RuntimeError):
    """Raised when a required project configuration is missing or invalid."""


def _load_yaml(path: Path) -> dict[str, Any]:
    """Load a YAML file and return its contents as a dictionary."""

    if not path.exists():
        raise ConfigurationError(
            f"Configuration file not found: {path}"
        )

    with path.open("r", encoding="utf-8") as file:
        data = yaml.safe_load(file)

    if data is None:
        raise ConfigurationError(
            f"Configuration file is empty: {path}"
        )

    if not isinstance(data, dict):
        raise ConfigurationError(
            f"Configuration file must contain a mapping: {path}"
        )

    return data


def _require_sections(
    config: dict[str, Any],
    required_sections: set[str],
    config_name: str,
) -> None:
    """Check that a configuration contains all required top-level sections."""

    missing_sections = required_sections - config.keys()

    if missing_sections:
        missing_text = ", ".join(sorted(missing_sections))

        raise ConfigurationError(
            f"{config_name} is missing required sections: {missing_text}"
        )


def load_project_config() -> dict[str, Any]:
    """Load the main project configuration."""

    config = _load_yaml(CONFIG_DIR / "project.yaml")

    _require_sections(
        config,
        {
            "project",
            "experimental_design",
            "paths",
            "analysis",
            "robustness",
        },
        "project.yaml",
)

    return config


def load_sources_config() -> dict[str, Any]:
    """Load the external-source registry."""

    config = _load_yaml(CONFIG_DIR / "sources.yaml")

    _require_sections(
        config,
        {"sources"},
        "sources.yaml",
    )

    return config


def load_scientific_baseline() -> dict[str, Any]:
    """Load the expected scientific results used for regression testing."""

    config = _load_yaml(CONFIG_DIR / "scientific_baseline.yaml")

    _require_sections(
        config,
        {
            "project",
            "dataset",
            "experimental_design",
            "primary_response_floor",
            "persistent_gene_counts",
            "primary_p95_persistence_by_condition",
            "directional_consistency",
            "functional_analysis",
            "published_de_concordance",
            "robustness",
        },
        "scientific_baseline.yaml",
    )

    return config