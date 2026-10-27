"""Regression tests for GO-network threshold robustness."""

import pandas as pd
import pytest

from arabidopsis_stress_recovery_atlas.analysis.go_network import (
    calculate_bp_term_overlap,
)
from arabidopsis_stress_recovery_atlas.analysis.go_network_robustness import (
    compare_largest_components,
    summarize_network_thresholds,
)
from arabidopsis_stress_recovery_atlas.config import (
    PROJECT_ROOT,
    load_project_config,
    load_scientific_baseline,
)


def test_go_network_threshold_robustness_matches_stage1():
    """Stage 2 should reproduce Stage 1 network-threshold robustness."""

    baseline = load_scientific_baseline()
    config = load_project_config()

    base = (
        PROJECT_ROOT
        / "data"
        / "interim"
        / "biological_interpretation"
    )

    significant_enrichment = pd.read_csv(
        base
        / "high_breadth_bp_enrichment_significant.csv"
    )

    membership = pd.read_csv(
        base
        / "significant_bp_gene_term_membership.csv"
    )

    overlap = calculate_bp_term_overlap(
        membership
    )

    thresholds = config[
        "robustness"
    ]["go_overlap_thresholds"]

    summary, largest_components = (
        summarize_network_thresholds(
            significant_enrichment,
            overlap,
            thresholds,
        )
    )

    comparison = compare_largest_components(
        largest_components
    )

    expected = baseline[
        "robustness"
    ]["go_overlap_network"]

    key_by_threshold = {
        0.40: "t040",
        0.50: "t050",
        0.60: "t060",
    }

    for threshold in thresholds:
        threshold = float(
            threshold
        )

        key = key_by_threshold[
            threshold
        ]

        expected_row = expected[
            key
        ]

        actual = summary[
            summary[
                "overlap_threshold"
            ]
            == threshold
        ].iloc[0]

        assert actual[
            "n_nodes"
        ] == expected_row[
            "nodes"
        ]

        assert actual[
            "n_edges"
        ] == expected_row[
            "edges"
        ]

        assert actual[
            "n_components"
        ] == expected_row[
            "components"
        ]

        assert actual[
            "component_sizes"
        ] == expected_row[
            "component_sizes"
        ]

        assert actual[
            "largest_component_size"
        ] == expected_row[
            "largest_component_size"
        ]

        assert actual[
            "largest_component_edges"
        ] == expected_row[
            "largest_component_edges"
        ]

        assert actual[
            "n_isolated_terms"
        ] == expected_row[
            "isolated_terms"
        ]

        assert actual[
            "n_communities"
        ] == expected_row[
            "communities"
        ]

        assert actual[
            "community_sizes"
        ] == expected_row[
            "community_sizes"
        ]

        assert actual[
            "modularity"
        ] == pytest.approx(
            expected_row[
                "modularity"
            ],
            abs=1e-6,
        )

    expected_similarity = expected[
        "largest_component_similarity"
    ]

    pair_keys = {
        (
            0.40,
            0.50,
        ): "t040_vs_t050",
        (
            0.50,
            0.60,
        ): "t050_vs_t060",
        (
            0.40,
            0.60,
        ): "t040_vs_t060",
    }

    for thresholds_pair, key in (
        pair_keys.items()
    ):
        threshold_a, threshold_b = (
            thresholds_pair
        )

        actual = comparison[
            (
                comparison["threshold_a"]
                == threshold_a
            )
            & (
                comparison["threshold_b"]
                == threshold_b
            )
        ].iloc[0]

        expected_pair = (
            expected_similarity[key]
        )

        assert actual[
            "shared_terms"
        ] == expected_pair[
            "shared_terms"
        ]

        assert actual[
            "jaccard"
        ] == pytest.approx(
            expected_pair[
                "jaccard"
            ],
            abs=1e-6,
        )

    # The stricter largest components should
    # remain nested within the looser ones.
    assert largest_components[
        0.60
    ].issubset(
        largest_components[
            0.50
        ]
    )

    assert largest_components[
        0.50
    ].issubset(
        largest_components[
            0.40
        ]
    )