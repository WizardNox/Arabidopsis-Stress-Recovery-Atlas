"""Load and standardize Arabidopsis Gene Ontology annotations."""

import io
import zipfile
from pathlib import Path

import pandas as pd

from arabidopsis_stress_recovery_atlas.config import (
    PROJECT_ROOT,
    load_sources_config,
)


class GOAnnotationError(ValueError):
    """Raised when Gene Ontology annotations cannot be loaded or standardized."""


GO_COLUMNS = [
    "gene_id",
    "tair_object",
    "gene_symbol",
    "relationship",
    "go_term",
    "go_id",
    "tair_keyword_id",
    "aspect",
    "goslim_category",
    "evidence_code",
    "evidence_description",
    "evidence_with",
    "reference",
    "annotator",
    "date_modified",
]


def get_go_annotation_path() -> Path:
    """Return the configured path to the Arabidopsis GO annotation source."""

    sources = load_sources_config()

    try:
        relative_path = sources[
            "sources"
        ][
            "ath_go_goslim"
        ][
            "local_path"
        ]

    except KeyError as exc:
        raise GOAnnotationError(
            "GO annotation source is not configured correctly."
        ) from exc

    return (
        PROJECT_ROOT
        / relative_path
    )


def _count_leading_metadata_lines(
    path: Path,
) -> int:
    """Count leading metadata lines in an uncompressed text source."""

    count = 0

    with path.open(
        "r",
        encoding="cp1252",
    ) as file:

        for line in file:

            if line.startswith(
                "!"
            ):
                count += 1

            else:
                break

    return count


def _load_go_text_file(
    path: Path,
) -> pd.DataFrame:
    """Load GO annotations from the uncompressed text representation."""

    metadata_lines = (
        _count_leading_metadata_lines(
            path
        )
    )

    return pd.read_csv(
        path,
        sep="\t",
        encoding="cp1252",
        header=None,
        names=GO_COLUMNS,
        skiprows=metadata_lines,
        dtype=str,
        keep_default_na=True,
    )


def _load_go_archive(
    path: Path,
) -> pd.DataFrame:
    """Load GO annotations directly from the TAIR archive."""

    try:
        with zipfile.ZipFile(
            path,
            "r",
        ) as archive:

            members = [
                name
                for name in archive.namelist()
                if (
                    not name.endswith("/")
                    and Path(
                        name
                    ).name
                    == "ATH_GO_GOSLIM.txt"
                )
            ]

            if len(
                members
            ) != 1:
                raise GOAnnotationError(
                    "GO annotation archive must contain exactly one "
                    "ATH_GO_GOSLIM.txt file."
                )

            with archive.open(
                members[0],
                "r",
            ) as binary_file:

                with io.TextIOWrapper(
                    binary_file,
                    encoding="cp1252",
                ) as text_file:

                    return pd.read_csv(
                        text_file,
                        sep="\t",
                        comment="!",
                        header=None,
                        names=GO_COLUMNS,
                        dtype=str,
                        keep_default_na=True,
                    )

    except (
        zipfile.BadZipFile,
        KeyError,
        OSError,
    ) as exc:
        raise GOAnnotationError(
            f"Could not load GO annotation archive: {path}"
        ) from exc


def load_go_annotations() -> pd.DataFrame:
    """Load the 15-column Arabidopsis GO/GO-Slim annotation table.

    Both the uncompressed text representation and the original
    TAIR archive are supported.
    """

    path = (
        get_go_annotation_path()
    )

    if not path.is_file():
        raise GOAnnotationError(
            f"GO annotation source not found: {path}"
        )

    try:

        if zipfile.is_zipfile(
            path
        ):
            data = (
                _load_go_archive(
                    path
                )
            )

        else:
            data = (
                _load_go_text_file(
                    path
                )
            )

    except (
        pd.errors.EmptyDataError,
        pd.errors.ParserError,
        UnicodeDecodeError,
        OSError,
    ) as exc:
        raise GOAnnotationError(
            f"Could not load GO annotations: {path}"
        ) from exc

    if len(
        data.columns
    ) != len(
        GO_COLUMNS
    ):
        raise GOAnnotationError(
            "GO annotation table does not contain "
            "the expected 15 fields."
        )

    return data


def get_go_term_memberships(
    annotations: pd.DataFrame,
    aspect: str | None = None,
) -> pd.DataFrame:
    """Return unique gene-to-GO-term memberships."""

    required_columns = {
        "gene_id",
        "go_id",
        "go_term",
        "aspect",
    }

    missing_columns = (
        required_columns
        - set(
            annotations.columns
        )
    )

    if missing_columns:
        raise GOAnnotationError(
            "GO annotation table is missing required columns: "
            + ", ".join(
                sorted(
                    missing_columns
                )
            )
        )

    memberships = annotations[
        [
            "gene_id",
            "go_id",
            "go_term",
            "aspect",
        ]
    ].copy()

    memberships = (
        memberships.dropna(
            subset=[
                "gene_id",
                "go_id",
                "go_term",
                "aspect",
            ]
        )
    )

    if aspect is not None:
        memberships = memberships[
            memberships[
                "aspect"
            ]
            == aspect
        ]

    return (
        memberships
        .drop_duplicates()
        .reset_index(
            drop=True
        )
    )


def get_goslim_memberships(
    annotations: pd.DataFrame,
    aspect: str | None = None,
) -> pd.DataFrame:
    """Return unique gene-to-GO-Slim-category memberships."""

    required_columns = {
        "gene_id",
        "aspect",
        "goslim_category",
    }

    missing_columns = (
        required_columns
        - set(
            annotations.columns
        )
    )

    if missing_columns:
        raise GOAnnotationError(
            "GO annotation table is missing required columns: "
            + ", ".join(
                sorted(
                    missing_columns
                )
            )
        )

    memberships = annotations[
        [
            "gene_id",
            "aspect",
            "goslim_category",
        ]
    ].copy()

    memberships = (
        memberships.dropna(
            subset=[
                "gene_id",
                "aspect",
                "goslim_category",
            ]
        )
    )

    memberships[
        "goslim_category"
    ] = (
        memberships[
            "goslim_category"
        ]
        .str.strip()
    )

    memberships = memberships[
        memberships[
            "goslim_category"
        ]
        != ""
    ]

    if aspect is not None:
        memberships = memberships[
            memberships[
                "aspect"
            ]
            == aspect
        ]

    return (
        memberships
        .drop_duplicates()
        .reset_index(
            drop=True
        )
    )