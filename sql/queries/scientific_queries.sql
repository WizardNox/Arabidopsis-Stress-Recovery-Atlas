-- ============================================================
-- 1. Genes persistent across the largest number of stresses
-- ============================================================

SELECT
    s.gene_id,
    s.persistence_breadth,
    s.n_persistent_up,
    s.n_persistent_down,
    g.short_description,
    g.curator_summary
FROM gene_persistence_summary AS s
LEFT JOIN gene AS g
    ON s.gene_id = g.gene_id
WHERE s.persistence_breadth >= 4
ORDER BY
    s.persistence_breadth DESC,
    s.gene_id;


-- ============================================================
-- 2. Recovery-outcome counts by stress condition
-- ============================================================

SELECT
    condition,
    recovery_outcome,
    COUNT(*) AS gene_count
FROM gene_condition_response
GROUP BY
    condition,
    recovery_outcome
ORDER BY
    condition,
    recovery_outcome;


-- ============================================================
-- 3. Drought-persistent genes
-- ============================================================

SELECT
    r.gene_id,
    r.endpoint_shift,
    r.recovery_shift,
    r.persistence_breadth,
    g.short_description
FROM gene_condition_response AS r
LEFT JOIN gene AS g
    ON r.gene_id = g.gene_id
WHERE
    r.condition = 'drought'
    AND r.recovery_outcome
        = 'persistent_same_direction'
ORDER BY
    r.persistence_breadth DESC,
    ABS(r.recovery_shift) DESC;


-- ============================================================
-- 4. Significant detailed BP enrichment terms
-- ============================================================

SELECT
    go_id,
    term,
    foreground_with_term,
    enrichment_ratio,
    odds_ratio,
    fdr_primary
FROM go_term_enrichment
WHERE
    analysis_type = 'detailed_bp'
    AND is_significant = 1
ORDER BY
    fdr_primary,
    enrichment_ratio DESC;


-- ============================================================
-- 5. Published recovery-DE concordance
-- ============================================================

SELECT
    p.gene_id,
    p.condition,
    p.project_direction,
    p.published_direction_original,
    p.published_direction_treatment_vs_control,
    p.direction_agrees,
    p.recovery_shift,
    p.published_log2fc_treatment_vs_control
FROM published_de_concordance AS p
ORDER BY
    p.condition,
    p.gene_id;