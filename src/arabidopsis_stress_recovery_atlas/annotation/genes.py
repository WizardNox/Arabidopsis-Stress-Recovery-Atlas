"""Load and standardize Araport11 functional gene annotations."""

from pathlib import Path

import pandas as pd

from arabidopsis_stress_recovery_atlas.config import (
    PROJECT_ROOT,
    load_sources_config,
)


class GeneAnnotationError(ValueError):
    """Raised when gene annotations cannot be loaded or standardized."""


ARAPORT_DESCRIPTION_COLUMNS = (
    "short_description",
    "Curator_summary",
    "Computational_description",
)


def get_araport_path() -> Path:
    """Return the configured path to the Araport11 annotation file."""

    sources = load_sources_config()

    try:
        relative_path = sources["sources"][
            "araport11_functional_descriptions"
        ]["local_path"]
    except KeyError as exc:
        raise GeneAnnotationError(
            "Araport11 source is not configured correctly."
        ) from exc

    return PROJECT_ROOT / relative_path


def load_araport_annotations() -> pd.DataFrame:
    """Load the Araport11 functional-description source table."""

    path = get_araport_path()

    if not path.is_file():
        raise GeneAnnotationError(
            f"Araport11 annotation file not found: {path}"
        )

    try:
        data = pd.read_csv(
            path,
            sep="\t",
            encoding="cp1252",
            dtype=str,
            na_values=["NULL"],
        )
    except (
        pd.errors.EmptyDataError,
        pd.errors.ParserError,
        UnicodeDecodeError,
        OSError,
    ) as exc:
        raise GeneAnnotationError(
            f"Could not load Araport11 annotations: {path}"
        ) from exc

    required_columns = {
        "name",
        "gene_model_type",
        "short_description",
        "Curator_summary",
        "Computational_description",
    }

    missing_columns = (
        required_columns - set(data.columns)
    )

    if missing_columns:
        raise GeneAnnotationError(
            "Araport11 table is missing required columns: "
            + ", ".join(sorted(missing_columns))
        )

    return data


def get_araport_gene_models(
    annotations: pd.DataFrame,
) -> pd.DataFrame:
    """Return AGI gene-model rows with a canonical locus-level gene_id."""

    models = annotations.copy()

    models["gene_id"] = models["name"].str.extract(
        r"^(AT[1-5MC]G\d+)(?:\.\d+)?$",
        expand=False,
    )

    models = models[
        models["gene_id"].notna()
    ].copy()

    return models


def _combine_unique_text(
    values: pd.Series,
) -> str | None:
    """Combine distinct non-empty annotation texts without losing variants."""

    cleaned = (
        values
        .dropna()
        .astype(str)
        .str.strip()
    )

    cleaned = cleaned[
        cleaned != ""
    ]

    unique_values = list(
        dict.fromkeys(cleaned.tolist())
    )

    if not unique_values:
        return None

    return " || ".join(unique_values)


def summarize_araport_loci(
    gene_models: pd.DataFrame,
) -> pd.DataFrame:
    """Collapse gene-model annotations to one row per Arabidopsis locus."""

    if "gene_id" not in gene_models.columns:
        raise GeneAnnotationError(
            "Gene-model table must contain canonical column 'gene_id'."
        )

    grouped = (
        gene_models
        .groupby("gene_id", sort=True)
        .agg(
            gene_model_count=(
                "name",
                "nunique",
            ),
            gene_models=(
                "name",
                lambda values: " | ".join(
                    sorted(
                        values.dropna().unique()
                    )
                ),
            ),
            gene_model_types=(
                "gene_model_type",
                _combine_unique_text,
            ),
            short_description=(
                "short_description",
                _combine_unique_text,
            ),
            curator_summary=(
                "Curator_summary",
                _combine_unique_text,
            ),
            computational_description=(
                "Computational_description",
                _combine_unique_text,
            ),
        )
        .reset_index()
    )

    return grouped