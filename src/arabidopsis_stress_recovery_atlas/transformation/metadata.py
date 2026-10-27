"""Reconstruct canonical sample metadata from expression-column names."""

import re
from typing import Any

import pandas as pd

from arabidopsis_stress_recovery_atlas.config import (
    load_project_config,
)


class SampleMetadataError(ValueError):
    """Raised when a sample label cannot be interpreted."""


COMPONENT_COLUMNS = (
    "cold",
    "heat",
    "drought",
    "high_light",
    "high_humidity",
)


CONTROL_PATTERN = re.compile(
    r"^CONT\s+(?P<replicate>\d+)$"
)


TREATMENT_PATTERN = re.compile(
    r"^(?P<state>endpoint|recovery)_"
    r"(?P<condition>.+?)\s*"
    r"(?P<replicate>\d+)$"
)


def parse_sample_name(
    sample_name: str,
) -> dict[str, Any]:
    """Convert one source sample name into canonical biological metadata."""

    design = load_project_config()["experimental_design"]

    control_match = CONTROL_PATTERN.fullmatch(sample_name)

    if control_match:
        metadata = {
            "sample_id": sample_name,
            "state": "control",
            "condition": "control",
            "source_condition_label": design["control_label"],
            "replicate": int(
                control_match.group("replicate")
            ),
        }

        for component in COMPONENT_COLUMNS:
            metadata[component] = False

        return metadata

    treatment_match = TREATMENT_PATTERN.fullmatch(
        sample_name
    )

    if treatment_match is None:
        raise SampleMetadataError(
            f"Unrecognized sample name: {sample_name}"
        )

    state = treatment_match.group("state")
    source_condition = treatment_match.group(
        "condition"
    ).strip()
    
    replicate = int(
        treatment_match.group("replicate")
    )

    condition_labels = design["condition_labels"]

    if source_condition not in condition_labels:
        raise SampleMetadataError(
            "Unrecognized stress-condition label "
            f"'{source_condition}' in sample '{sample_name}'."
        )

    condition = condition_labels[source_condition]

    stress_components = set(
        design["stress_components"][condition]
    )

    metadata = {
        "sample_id": sample_name,
        "state": state,
        "condition": condition,
        "source_condition_label": source_condition,
        "replicate": replicate,
    }

    for component in COMPONENT_COLUMNS:
        metadata[component] = (
            component in stress_components
        )

    return metadata


def build_sample_manifest(
    sample_columns: list[str],
) -> pd.DataFrame:
    """Build the canonical sample manifest from expression-column names."""

    records = [
        parse_sample_name(sample_name)
        for sample_name in sample_columns
    ]

    manifest = pd.DataFrame(records)

    return manifest