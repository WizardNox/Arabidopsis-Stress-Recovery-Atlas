"""Load the processed GSE226105 expression matrix."""

from pathlib import Path

import pandas as pd

from arabidopsis_stress_recovery_atlas.config import (
    PROJECT_ROOT,
    load_sources_config,
)


class ExpressionDataLoadError(RuntimeError):
    """Raised when the expression source cannot be loaded."""


def get_expression_matrix_path() -> Path:
    """Return the configured path to the expression matrix."""

    sources = load_sources_config()

    try:
        relative_path = sources["sources"]["expression_matrix"]["local_path"]
    except KeyError as exc:
        raise ExpressionDataLoadError(
            "Expression-matrix source is not configured correctly."
        ) from exc

    return PROJECT_ROOT / relative_path


def load_expression_data() -> pd.DataFrame:
    """Load the processed expression matrix without modifying it."""

    path = get_expression_matrix_path()

    if not path.is_file():
        raise ExpressionDataLoadError(
            f"Expression matrix not found: {path}"
        )

    try:
        data = pd.read_csv(path)
    except (
        pd.errors.EmptyDataError,
        pd.errors.ParserError,
        UnicodeDecodeError,
        OSError,
    ) as exc:
        raise ExpressionDataLoadError(
            f"Could not load expression matrix: {path}"
        ) from exc

    return data