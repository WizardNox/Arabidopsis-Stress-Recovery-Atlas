"""Build overlap networks among significant Biological Process terms."""

from itertools import combinations

import networkx as nx
import pandas as pd
from networkx.algorithms.community import (
    greedy_modularity_communities,
    modularity,
)


class GONetworkError(ValueError):
    """Raised when a GO-term overlap network cannot be constructed."""


def calculate_bp_term_overlap(
    significant_membership: pd.DataFrame,
) -> pd.DataFrame:
    """Calculate pairwise gene overlap among significant BP terms."""

    required_columns = {
        "gene_id",
        "go_id",
        "go_term",
    }

    missing_columns = (
        required_columns
        - set(significant_membership.columns)
    )

    if missing_columns:
        raise GONetworkError(
            "Significant BP membership is missing required columns: "
            + ", ".join(sorted(missing_columns))
        )

    term_gene_sets = (
        significant_membership
        .groupby(
            [
                "go_id",
                "go_term",
            ]
        )["gene_id"]
        .apply(set)
    )

    terms = list(
        term_gene_sets.index
    )

    records = []

    for (
        go_id_a,
        go_term_a,
    ), (
        go_id_b,
        go_term_b,
    ) in combinations(
        terms,
        2,
    ):
        genes_a = term_gene_sets.loc[
            (go_id_a, go_term_a)
        ]

        genes_b = term_gene_sets.loc[
            (go_id_b, go_term_b)
        ]

        intersection = len(
            genes_a & genes_b
        )

        union = len(
            genes_a | genes_b
        )

        jaccard = (
            intersection / union
            if union > 0
            else 0.0
        )

        smaller_set = min(
            len(genes_a),
            len(genes_b),
        )

        overlap_coefficient = (
            intersection / smaller_set
            if smaller_set > 0
            else 0.0
        )

        records.append(
            {
                "go_id_a": go_id_a,
                "go_term_a": go_term_a,
                "go_id_b": go_id_b,
                "go_term_b": go_term_b,
                "genes_a": len(genes_a),
                "genes_b": len(genes_b),
                "intersection": intersection,
                "union": union,
                "jaccard": jaccard,
                "overlap_coefficient": (
                    overlap_coefficient
                ),
            }
        )

    return pd.DataFrame(records)


def filter_strong_bp_overlap(
    overlap: pd.DataFrame,
    threshold: float,
) -> pd.DataFrame:
    """Keep term pairs whose overlap coefficient reaches the threshold."""

    if not 0 <= threshold <= 1:
        raise GONetworkError(
            "Overlap threshold must be between 0 and 1."
        )

    return (
        overlap[
            overlap["overlap_coefficient"]
            >= threshold
        ]
        .copy()
        .reset_index(drop=True)
    )


def build_bp_overlap_graph(
    significant_enrichment: pd.DataFrame,
    strong_overlap: pd.DataFrame,
) -> nx.Graph:
    """Build the significant-BP overlap graph."""

    graph = nx.Graph()

    for row in significant_enrichment.itertuples(
        index=False
    ):
        graph.add_node(
            row.go_id,
            go_term=row.go_term,
        )

    for row in strong_overlap.itertuples(
        index=False
    ):
        graph.add_edge(
            row.go_id_a,
            row.go_id_b,
            overlap_coefficient=(
                row.overlap_coefficient
            ),
            jaccard=row.jaccard,
            intersection=row.intersection,
        )

    return graph


def summarize_components(
    graph: nx.Graph,
    significant_enrichment: pd.DataFrame,
) -> pd.DataFrame:
    """Summarize connected components of the BP overlap graph."""

    components = sorted(
        nx.connected_components(graph),
        key=len,
        reverse=True,
    )

    records = []

    for component_id, component in enumerate(
        components,
        start=1,
    ):
        for go_id in component:
            term_row = (
                significant_enrichment[
                    significant_enrichment[
                        "go_id"
                    ]
                    == go_id
                ]
                .iloc[0]
            )

            records.append(
                {
                    "component_id": (
                        component_id
                    ),
                    "component_size": (
                        len(component)
                    ),
                    "go_id": go_id,
                    "go_term": (
                        term_row["go_term"]
                    ),
                    "n_strong_links": (
                        graph.degree(go_id)
                    ),
                    "foreground_with_term": (
                        term_row[
                            "foreground_with_term"
                        ]
                    ),
                    "enrichment_ratio": (
                        term_row[
                            "enrichment_ratio"
                        ]
                    ),
                    "fdr_bh": (
                        term_row["fdr_bh"]
                    ),
                }
            )

    return (
        pd.DataFrame(records)
        .sort_values(
            [
                "component_id",
                "n_strong_links",
                "fdr_bh",
            ],
            ascending=[
                True,
                False,
                True,
            ],
        )
        .reset_index(drop=True)
    )


def get_largest_component_graph(
    graph: nx.Graph,
) -> nx.Graph:
    """Return a copy of the graph's largest connected component."""

    components = sorted(
        nx.connected_components(graph),
        key=len,
        reverse=True,
    )

    if not components:
        raise GONetworkError(
            "GO overlap graph contains no components."
        )

    return graph.subgraph(
        components[0]
    ).copy()


def detect_bp_communities(
    graph: nx.Graph,
) -> tuple[list[set[str]], float]:
    """Detect communities in a connected BP graph."""

    if graph.number_of_nodes() == 0:
        raise GONetworkError(
            "Cannot detect communities in an empty graph."
        )

    communities = list(
        greedy_modularity_communities(
            graph,
            weight="overlap_coefficient",
        )
    )

    community_sets = [
        set(community)
        for community in communities
    ]

    score = modularity(
        graph,
        communities,
        weight="overlap_coefficient",
    )

    return community_sets, float(score)