"""Robustness analysis for GO-term overlap-network thresholds."""

from itertools import combinations
from collections.abc import Iterable

import networkx as nx
import pandas as pd

from arabidopsis_stress_recovery_atlas.analysis.go_network import (
    build_bp_overlap_graph,
    detect_bp_communities,
    filter_strong_bp_overlap,
    get_largest_component_graph,
)


def summarize_network_thresholds(
    significant_enrichment: pd.DataFrame,
    overlap: pd.DataFrame,
    thresholds: Iterable[float],
) -> tuple[
    pd.DataFrame,
    dict[float, set[str]],
]:
    """Summarize BP overlap networks across multiple thresholds."""

    summary_records = {}
    largest_components = {}

    for threshold in thresholds:
        strong = filter_strong_bp_overlap(
            overlap,
            threshold,
        )

        graph = build_bp_overlap_graph(
            significant_enrichment,
            strong,
        )

        components = sorted(
            nx.connected_components(graph),
            key=len,
            reverse=True,
        )

        largest = get_largest_component_graph(
            graph
        )

        communities, modularity_score = (
            detect_bp_communities(
                largest
            )
        )

        largest_components[
            float(threshold)
        ] = set(
            largest.nodes()
        )

        summary_records[
            float(threshold)
        ] = {
            "overlap_threshold": (
                float(threshold)
            ),
            "n_nodes": (
                graph.number_of_nodes()
            ),
            "n_edges": (
                graph.number_of_edges()
            ),
            "n_components": (
                len(components)
            ),
            "component_sizes": sorted(
                [
                    len(component)
                    for component in components
                ],
                reverse=True,
            ),
            "largest_component_size": (
                largest.number_of_nodes()
            ),
            "largest_component_edges": (
                largest.number_of_edges()
            ),
            "n_isolated_terms": (
                len(
                    list(
                        nx.isolates(graph)
                    )
                )
            ),
            "n_communities": (
                len(communities)
            ),
            "community_sizes": sorted(
                [
                    len(community)
                    for community
                    in communities
                ],
                reverse=True,
            ),
            "modularity": (
                float(modularity_score)
            ),
        }

    summary = pd.DataFrame(
        summary_records.values()
    ).sort_values(
        "overlap_threshold"
    )

    return (
        summary.reset_index(drop=True),
        largest_components,
    )


def compare_largest_components(
    largest_components: dict[
        float,
        set[str],
    ],
) -> pd.DataFrame:
    """Compare largest-component membership across thresholds."""

    records = []

    thresholds = sorted(
        largest_components
    )

    for threshold_a, threshold_b in combinations(
        thresholds,
        2,
    ):
        terms_a = largest_components[
            threshold_a
        ]

        terms_b = largest_components[
            threshold_b
        ]

        shared_terms = len(
            terms_a & terms_b
        )

        union_terms = len(
            terms_a | terms_b
        )

        jaccard = (
            shared_terms / union_terms
            if union_terms > 0
            else 0.0
        )

        records.append(
            {
                "threshold_a": threshold_a,
                "threshold_b": threshold_b,
                "size_a": len(terms_a),
                "size_b": len(terms_b),
                "shared_terms": shared_terms,
                "union_terms": union_terms,
                "jaccard": jaccard,
            }
        )

    return pd.DataFrame(records)