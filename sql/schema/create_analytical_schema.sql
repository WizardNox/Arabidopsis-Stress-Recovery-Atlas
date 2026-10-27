PRAGMA foreign_keys = ON;


CREATE TABLE sample (
    sample_id TEXT PRIMARY KEY,

    state TEXT NOT NULL
        CHECK (
            state IN (
                'control',
                'endpoint',
                'recovery'
            )
        ),

    condition TEXT NOT NULL
        CHECK (
            condition IN (
                'control',
                'cold',
                'cold_high_light',
                'drought',
                'heat',
                'heat_drought',
                'heat_drought_high_light',
                'heat_high_light',
                'high_humidity',
                'high_light'
            )
        ),

    source_condition_label TEXT NOT NULL,

    replicate INTEGER NOT NULL
        CHECK (replicate >= 1),

    cold INTEGER NOT NULL
        CHECK (cold IN (0, 1)),

    heat INTEGER NOT NULL
        CHECK (heat IN (0, 1)),

    drought INTEGER NOT NULL
        CHECK (drought IN (0, 1)),

    high_light INTEGER NOT NULL
        CHECK (high_light IN (0, 1)),

    high_humidity INTEGER NOT NULL
        CHECK (high_humidity IN (0, 1))
);


CREATE TABLE gene (
    gene_id TEXT PRIMARY KEY,

    has_araport_annotation INTEGER NOT NULL
        CHECK (
            has_araport_annotation
            IN (0, 1)
        ),

    gene_model_count INTEGER
        CHECK (
            gene_model_count IS NULL
            OR gene_model_count >= 1
        ),

    gene_models TEXT,

    gene_model_types TEXT,

    short_description TEXT,

    curator_summary TEXT,

    computational_description TEXT
);


CREATE TABLE gene_condition_response (
    gene_id TEXT NOT NULL,

    condition TEXT NOT NULL
        CHECK (
            condition IN (
                'cold',
                'cold_high_light',
                'drought',
                'heat',
                'heat_drought',
                'heat_drought_high_light',
                'heat_high_light',
                'high_humidity',
                'high_light'
            )
        ),

    control_mean REAL NOT NULL,

    endpoint_mean REAL NOT NULL,

    recovery_mean REAL NOT NULL,

    endpoint_shift REAL NOT NULL,

    recovery_shift REAL NOT NULL,

    response_floor REAL NOT NULL
        CHECK (response_floor >= 0),

    endpoint_eligible INTEGER NOT NULL
        CHECK (
            endpoint_eligible
            IN (0, 1)
        ),

    recovery_above_floor INTEGER NOT NULL
        CHECK (
            recovery_above_floor
            IN (0, 1)
        ),

    same_direction INTEGER NOT NULL
        CHECK (
            same_direction
            IN (0, 1)
        ),

    recovery_outcome TEXT NOT NULL
        CHECK (
            recovery_outcome IN (
                'not_endpoint_eligible',
                'persistent_same_direction',
                'recovered',
                'substantial_reversal'
            )
        ),

    persistence_breadth INTEGER NOT NULL
        CHECK (
            persistence_breadth
            BETWEEN 0 AND 9
        ),

    PRIMARY KEY (
        gene_id,
        condition
    ),

    FOREIGN KEY (
        gene_id
    )
    REFERENCES gene (
        gene_id
    )
    ON DELETE CASCADE
);


CREATE TABLE gene_persistence_summary (
    gene_id TEXT PRIMARY KEY,

    persistence_breadth INTEGER NOT NULL
        CHECK (
            persistence_breadth
            BETWEEN 0 AND 9
        ),

    n_persistent_conditions INTEGER NOT NULL
        CHECK (
            n_persistent_conditions
            BETWEEN 0 AND 9
        ),

    n_persistent_up INTEGER NOT NULL
        CHECK (
            n_persistent_up
            BETWEEN 0 AND 9
        ),

    n_persistent_down INTEGER NOT NULL
        CHECK (
            n_persistent_down
            BETWEEN 0 AND 9
        ),

    FOREIGN KEY (
        gene_id
    )
    REFERENCES gene (
        gene_id
    )
    ON DELETE CASCADE,

    CHECK (
        persistence_breadth
        = n_persistent_conditions
    ),

    CHECK (
        n_persistent_up
        + n_persistent_down
        = n_persistent_conditions
    )
);


CREATE TABLE go_term_enrichment (
    enrichment_id INTEGER PRIMARY KEY AUTOINCREMENT,

    analysis_type TEXT NOT NULL
        CHECK (
            analysis_type IN (
                'go_slim',
                'detailed_bp'
            )
        ),

    go_aspect TEXT NOT NULL
        CHECK (
            go_aspect IN (
                'P',
                'F',
                'C'
            )
        ),

    go_id TEXT,

    term TEXT NOT NULL,

    foreground_with_term INTEGER NOT NULL
        CHECK (
            foreground_with_term >= 0
        ),

    foreground_total INTEGER NOT NULL
        CHECK (
            foreground_total > 0
        ),

    background_with_term INTEGER NOT NULL
        CHECK (
            background_with_term >= 0
        ),

    background_total INTEGER NOT NULL
        CHECK (
            background_total > 0
        ),

    foreground_fraction REAL NOT NULL
        CHECK (
            foreground_fraction
            BETWEEN 0 AND 1
        ),

    background_fraction REAL NOT NULL
        CHECK (
            background_fraction
            BETWEEN 0 AND 1
        ),

    enrichment_ratio REAL NOT NULL,

    odds_ratio REAL NOT NULL,

    p_value REAL NOT NULL
        CHECK (
            p_value
            BETWEEN 0 AND 1
        ),

    fdr_primary REAL NOT NULL
        CHECK (
            fdr_primary
            BETWEEN 0 AND 1
        ),

    fdr_secondary REAL
        CHECK (
            fdr_secondary IS NULL
            OR fdr_secondary
            BETWEEN 0 AND 1
        ),

    is_significant INTEGER NOT NULL
        CHECK (
            is_significant
            IN (0, 1)
        ),

    CHECK (
        (
            analysis_type = 'go_slim'
            AND go_id IS NULL
            AND fdr_secondary IS NOT NULL
        )
        OR
        (
            analysis_type = 'detailed_bp'
            AND go_aspect = 'P'
            AND go_id IS NOT NULL
            AND fdr_secondary IS NULL
        )
    )
);


CREATE TABLE published_de_concordance (
    gene_id TEXT NOT NULL,

    condition TEXT NOT NULL,

    project_direction TEXT NOT NULL
        CHECK (
            project_direction IN (
                'up',
                'down'
            )
        ),

    endpoint_shift REAL NOT NULL,

    recovery_shift REAL NOT NULL,

    persistence_breadth INTEGER NOT NULL
        CHECK (
            persistence_breadth
            BETWEEN 0 AND 9
        ),

    published_log2fc_original REAL NOT NULL,

    published_adj_p REAL NOT NULL
        CHECK (
            published_adj_p
            BETWEEN 0 AND 1
        ),

    published_direction_original TEXT NOT NULL
        CHECK (
            published_direction_original IN (
                'up',
                'down'
            )
        ),

    published_log2fc_treatment_vs_control REAL NOT NULL,

    published_direction_treatment_vs_control TEXT NOT NULL
        CHECK (
            published_direction_treatment_vs_control IN (
                'up',
                'down'
            )
        ),

    direction_agrees INTEGER NOT NULL
        CHECK (
            direction_agrees
            IN (0, 1)
        ),

    PRIMARY KEY (
        gene_id,
        condition
    ),

    FOREIGN KEY (
        gene_id
    )
    REFERENCES gene (
        gene_id
    )
    ON DELETE CASCADE,

    FOREIGN KEY (
        gene_id,
        condition
    )
    REFERENCES gene_condition_response (
        gene_id,
        condition
    )
    ON DELETE CASCADE
);


CREATE INDEX idx_sample_state_condition
ON sample (
    state,
    condition
);


CREATE INDEX idx_gene_araport
ON gene (
    has_araport_annotation
);


CREATE INDEX idx_response_condition_outcome
ON gene_condition_response (
    condition,
    recovery_outcome
);


CREATE INDEX idx_response_breadth
ON gene_condition_response (
    persistence_breadth
);


CREATE INDEX idx_persistence_breadth
ON gene_persistence_summary (
    persistence_breadth
);


CREATE INDEX idx_go_significant
ON go_term_enrichment (
    analysis_type,
    is_significant
);


CREATE INDEX idx_go_id
ON go_term_enrichment (
    go_id
);


CREATE UNIQUE INDEX idx_go_identity
ON go_term_enrichment (
    analysis_type,
    go_aspect,
    IFNULL(
        go_id,
        ''
    ),
    term
);


CREATE INDEX idx_published_condition
ON published_de_concordance (
    condition
);