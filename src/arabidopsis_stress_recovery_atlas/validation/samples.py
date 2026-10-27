"""Validate the reconstructed experimental sample manifest."""

from typing import Any

import pandas as pd

from arabidopsis_stress_recovery_atlas.config import (
    load_project_config,
    load_scientific_baseline,
)
from arabidopsis_stress_recovery_atlas.transformation.metadata import (
    COMPONENT_COLUMNS,
)


class SampleValidationError(ValueError):
    """Raised when the sample manifest violates experimental-design rules."""


REQUIRED_MANIFEST_COLUMNS = {
    "sample_id",
    "state",
    "condition",
    "source_condition_label",
    "replicate",
    *COMPONENT_COLUMNS,
}


def validate_sample_manifest(
    manifest: pd.DataFrame,
) -> dict[str, Any]:
    """Validate the sample manifest against the known experimental design."""

    baseline = load_scientific_baseline()
    design = load_project_config()["experimental_design"]

    missing_columns = (
        REQUIRED_MANIFEST_COLUMNS - set(manifest.columns)
    )

    if missing_columns:
        raise SampleValidationError(
            "Sample manifest is missing required columns: "
            + ", ".join(sorted(missing_columns))
        )

    if manifest["sample_id"].duplicated().any():
        duplicate_count = int(
            manifest["sample_id"].duplicated().sum()
        )

        raise SampleValidationError(
            "Sample manifest contains "
            f"{duplicate_count} duplicate sample identifiers."
        )

    allowed_states = set(design["states"])

    observed_states = set(
        manifest["state"].unique()
    )

    unexpected_states = observed_states - allowed_states

    if unexpected_states:
        raise SampleValidationError(
            "Unexpected experimental states: "
            + ", ".join(sorted(unexpected_states))
        )

    expected_state_counts = {
        "control": baseline["experimental_design"][
            "control_samples"
        ],
        "endpoint": baseline["experimental_design"][
            "endpoint_samples"
        ],
        "recovery": baseline["experimental_design"][
            "recovery_samples"
        ],
    }

    observed_state_counts = (
        manifest["state"]
        .value_counts()
        .to_dict()
    )

    if observed_state_counts != expected_state_counts:
        raise SampleValidationError(
            "Observed sample-state counts do not match "
            "the scientific baseline. "
            f"Expected {expected_state_counts}, "
            f"found {observed_state_counts}."
        )

    stress_rows = manifest[
        manifest["state"] != "control"
    ]

    observed_conditions = set(
        stress_rows["condition"].unique()
    )

    expected_conditions = set(
        design["condition_labels"].values()
    )

    if observed_conditions != expected_conditions:
        raise SampleValidationError(
            "Observed stress conditions do not match "
            "the configured experimental design."
        )

    expected_condition_count = baseline[
        "experimental_design"
    ]["stress_conditions"]

    if len(observed_conditions) != expected_condition_count:
        raise SampleValidationError(
            "Unexpected number of stress conditions. "
            f"Expected {expected_condition_count}, "
            f"found {len(observed_conditions)}."
        )

    expected_replicates = set(
        range(
            1,
            baseline["experimental_design"][
                "biological_replicates_per_stress_state"
            ]
            + 1,
        )
    )

    for state in ("endpoint", "recovery"):
        for condition in sorted(expected_conditions):
            group = manifest[
                (manifest["state"] == state)
                & (manifest["condition"] == condition)
            ]

            observed_replicates = set(
                group["replicate"].tolist()
            )

            if observed_replicates != expected_replicates:
                raise SampleValidationError(
                    "Unexpected replicate structure for "
                    f"{state} / {condition}. "
                    f"Expected {sorted(expected_replicates)}, "
                    f"found {sorted(observed_replicates)}."
                )

    control_rows = manifest[
        manifest["state"] == "control"
    ]

    control_replicates = set(
        control_rows["replicate"].tolist()
    )

    if control_replicates != expected_replicates:
        raise SampleValidationError(
            "Unexpected control replicate structure."
        )

    if not (
        control_rows[list(COMPONENT_COLUMNS)]
        .eq(False)
        .all()
        .all()
    ):
        raise SampleValidationError(
            "Control samples must not contain stress-component flags."
        )

    for condition, components in design[
        "stress_components"
    ].items():
        expected_components = set(components)

        condition_rows = stress_rows[
            stress_rows["condition"] == condition
        ]

        for component in COMPONENT_COLUMNS:
            expected_value = (
                component in expected_components
            )

            if not (
                condition_rows[component]
                == expected_value
            ).all():
                raise SampleValidationError(
                    "Stress-component flags are inconsistent "
                    f"for condition '{condition}'."
                )

    return {
        "samples": len(manifest),
        "control_samples": observed_state_counts[
            "control"
        ],
        "endpoint_samples": observed_state_counts[
            "endpoint"
        ],
        "recovery_samples": observed_state_counts[
            "recovery"
        ],
        "stress_conditions": len(
            observed_conditions
        ),
    }