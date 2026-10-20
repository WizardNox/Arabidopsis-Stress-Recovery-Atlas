# Arabidopsis Stress Recovery Atlas

**A computational analysis of context-dependent transcriptomic recovery after abiotic stress in _Arabidopsis thaliana_.**

**Author:** Swapnanil Sinha  
**Primary dataset:** NCBI GEO `GSE226105` · BioProject `PRJNA938512`  
**Scientific phase:** Complete  
**Next stage:** Reproducible analytical data platform

> **What remains after the stress is removed?**

Most stress-transcriptomics analyses focus on the acute response: which genes change under stress, how strongly they change, and which biological functions are enriched. This project asks the complementary question of what remains after plants return to standard conditions.

Using the public _Arabidopsis thaliana_ RNA-seq experiment `GSE226105`, the analysis compares stress endpoint and recovery states across nine single and combined abiotic-stress conditions. It examines transcriptome-level displacement from control, classifies gene-level recovery behavior, identifies persistent responses that recur across stresses, interprets those genes through functional annotation and Gene Ontology analysis, compares the exploratory results with the source study's published DESeq2 recovery results, and tests whether the main conclusions survive substantial changes in the empirical response threshold.

The main result is not that the transcriptome fails to recover. The data are broadly consistent with recovery after stress. Instead, the analysis shows that broad recovery can coexist with smaller residual expression patterns whose size, direction, recurrence, and functional composition differ strongly among stress conditions.

> **Scope note:** this repository analyzes the processed rlog expression matrix and published supplementary differential-expression results from the source study. It does not process FASTQ files or replace formal count-based differential-expression analysis.

---

## Contents

- [Research question](#research-question)
- [Dataset and experimental design](#dataset-and-experimental-design)
- [Analysis strategy](#analysis-strategy)
- [Key analytical definitions](#key-analytical-definitions)
- [Main findings](#main-findings)
- [Functional interpretation](#functional-interpretation)
- [GO redundancy and overlap network](#go-redundancy-and-overlap-network)
- [Published DESeq2 concordance](#published-deseq2-concordance)
- [Robustness analysis](#robustness-analysis)
- [Interpretation guardrails](#interpretation-guardrails)
- [Limitations](#limitations)
- [Scientific notebooks](#scientific-notebooks)
- [Repository structure](#repository-structure)
- [Reproducing the scientific analysis](#reproducing-the-scientific-analysis)
- [Generated outputs](#generated-outputs)
- [Next stage: analytical data platform](#next-stage-analytical-data-platform)
- [Future scientific work](#future-scientific-work)
- [Data provenance and attribution](#data-provenance-and-attribution)
- [Citation and license](#citation-and-license)
- [Final perspective](#final-perspective)

---

## Research question

The project is built around one main question:

> **How does the _Arabidopsis thaliana_ transcriptome differ between stress endpoint and recovery states, and which residual expression patterns are shared versus context-dependent across stress conditions?**

This is broken into six connected questions:

1. Is the processed expression matrix internally consistent enough for transcriptome-level exploratory analysis?
2. How far is each stress condition displaced from control at endpoint, and how much displacement remains in the recovery cohort?
3. Which qualifying gene-level endpoint responses recover, remain persistently displaced in the same direction, or substantially reverse?
4. Which persistent genes are stress-specific, and which recur across multiple stress conditions?
5. Are broadly persistent genes associated with coherent biological functions?
6. Are the main patterns consistent with the available published DESeq2 recovery results, and do they remain recognizable when the empirical response threshold is changed?

The working expectation was that different stresses would show different degrees of recovery rather than a simple relationship in which a larger endpoint response automatically produces a larger residual state.

---

## Dataset and experimental design

### Primary dataset

| Item | Value |
|---|---|
| Organism | _Arabidopsis thaliana_ |
| Genotype | wild-type Col-0 |
| Tissue | rosette leaves |
| GEO accession | `GSE226105` |
| BioProject | `PRJNA938512` |
| Processed matrix | `GSE226105_rLog-Normalized_counts.csv` |
| Gene rows | 28,775 |
| Expression samples | 57 |
| Informative genes used for variation-based analyses | 26,118 |
| All-zero genes excluded from variation-based analyses | 2,657 |

The processed matrix contains one gene identifier column and 57 rlog-transformed expression columns.

### Experimental states

The experiment contains three sample states:

- **Control:** 3 samples
- **Endpoint:** 27 samples
- **Recovery:** 27 samples

There are three biological replicates for each stress-state combination.

The nine standardized stress conditions are:

```text
cold
cold_high_light
drought
heat
heat_drought
heat_drought_high_light
heat_high_light
high_humidity
high_light
```

The abbreviated raw condition label `T` corresponds to the triple treatment:

```text
heat + drought + high light
```

### Endpoint and recovery are separate cohorts

The experiment does **not** repeatedly measure the same plants.

Endpoint and recovery samples are parallel age-matched day-18 cohorts:

- endpoint plants experienced stress during the final three days before harvest;
- recovery plants experienced the corresponding stress earlier, returned to standard conditions for three days, and were then harvested at the same age.

The analysis therefore compares **group-level experimental states**, not within-plant longitudinal trajectories.

---

## Analysis strategy

The scientific phase follows a deliberate sequence:

```text
processed expression matrix
        |
        v
experimental-design reconstruction
        |
        v
expression-matrix audit
        |
        v
sample-level quality control
        |
        v
endpoint and recovery displacement from control
        |
        v
gene-level recovery classification
        |
        v
cross-stress persistence and persistence breadth
        |
        v
functional annotation and GO enrichment
        |
        v
GO-term overlap and descriptive network analysis
        |
        v
comparison with published DESeq2 recovery results
        |
        v
response-floor and functional robustness checks
```

Each later analysis depends on assumptions checked earlier. The processed matrix is not treated as automatically correct simply because it is publicly available.

### Expression-matrix audit

The initial notebook verifies:

- matrix dimensions;
- gene identifiers;
- numeric expression columns;
- missing and infinite values;
- all-zero rows;
- repeated numerical profiles;
- sample labels;
- experimental state;
- stress condition;
- replicate number;
- component stressors.

The resulting design contains exactly 3 control, 27 endpoint, and 27 recovery samples.

### Sample-level quality control

Quality control is performed on the **processed rlog matrix**, not on FASTQ reads.

The QC analysis includes:

- sample-wise expression distributions;
- transcriptome-wide Pearson correlations;
- within-group and between-group similarity;
- nearest-neighbor relationships;
- PCA using the 500 most variable genes;
- hierarchical clustering;
- a second PCA using all 26,118 informative genes.

Two endpoint samples received particular attention:

- `endpoint_T 2`, an atypical triple-stress endpoint replicate;
- `endpoint_HEAT DROUGHT 3`, which showed greater within-treatment variation among highly variable genes while retaining the expected treatment identity.

Neither sample was removed because the available evidence did not establish technical failure or clear misassignment.

> **Sample-to-sample expression correlations**  
>![Sample-to-sample expression correlation](visualizations/sample_to_sample_expression_correlation_heatmap.png)  
> Shows transcriptome-wide Pearson correlation structure across all 57 samples and helps identify samples with unusual similarity patterns.

> **PCA of the 500 most variable genes by stress condition**  
>![PCA of the 500 most variable genes by stress condition](visualizations/pca_of_the_500_most_variable_genes_by_stress_conditions.png)
> Shows the dominant sample structure across stress conditions and highlights the atypical position of `endpoint_T 2`.

---

## Key analytical definitions

### Informative gene

A gene with at least one non-zero rlog value across the 57 samples.

Genes equal to zero in every sample are excluded from variation-based analyses.

### Endpoint shift

For gene `g` under stress condition `s`:

```text
endpoint_shift(g, s) = mean(endpoint_s) - mean(control)
```

This measures how far the endpoint group lies from the untreated control baseline for that gene.

### Recovery shift

```text
recovery_shift(g, s) = mean(recovery_s) - mean(control)
```

This measures how far the recovery group remains from the same control baseline.

Signed shifts preserve direction relative to control. Absolute shifts are used when the question concerns magnitude rather than direction.

### Empirical response floor

Small differences in transformed expression should not automatically be treated as meaningful biological responses.

The primary exploratory response floor is derived from absolute pairwise differences among the three control replicates across all 26,118 informative genes:

```text
26,118 genes x 3 control-control pairs = 78,354 control-derived differences
```

The primary threshold is the 95th percentile of that distribution:

```text
p95 response floor = 0.342263 rlog units
```

A gene-condition endpoint response must satisfy:

```text
|endpoint_shift| >= response_floor
```

before it enters the recovery classification.

This response floor is an **exploratory magnitude threshold**. It is not a p-value, not a formal significance cutoff, not a gene-specific variance model, and not a replacement for DESeq2 or another count-based differential-expression method.

### Comparative Reversibility Index

For qualifying endpoint responses:

```text
RI = 1 - |Recovery - Control| / |Endpoint - Control|
```

Interpretation:

- `RI ~ 1`: recovery is close to control relative to the endpoint state;
- `RI ~ 0`: recovery displacement is similar in magnitude to endpoint displacement;
- `RI < 0`: recovery is farther from control than the endpoint state.

Because endpoint and recovery samples are separate cohorts, this is a **comparative state index**, not a within-plant trajectory.

### Recovery categories

For genes with a qualifying endpoint response:

**Recovered**

```text
|recovery_shift| < response_floor
```

**Persistent, same direction**

```text
|recovery_shift| >= response_floor
and
sign(recovery_shift) == sign(endpoint_shift)
```

**Substantial reversal**

```text
|recovery_shift| >= response_floor
and
sign(recovery_shift) != sign(endpoint_shift)
```

Genes without a qualifying endpoint response are kept separate rather than being mislabeled as recovered.

### Persistence breadth

Persistence breadth is the number of stress conditions in which a gene is classified as persistent in the same direction.

The primary high-breadth definition is:

```text
persistent in at least 4 of 9 stress conditions
```

This is an operational prioritization rule, not a natural biological boundary. The robustness analysis also examines nearby breadth definitions and alternative response floors.

---

## Main findings

### At a glance

| Result | Primary analysis |
|---|---:|
| Informative genes | 26,118 |
| Primary response floor | 0.342263 rlog |
| Genes persistent in at least one stress | 1,609 |
| Genes persistent in at least two stresses | 726 |
| High-breadth genes, persistence in at least four stresses | 206 |
| Directionally consistent among multi-stress persistent genes | 720 / 726 |
| Consistently above control among directionally consistent genes | 711 / 720 |
| Significant GO-Slim terms for the 206-gene set | 24 |
| Significant detailed Biological Process terms | 45 |
| Shared drought genes with published recovery DEGs | 348 |
| Shared published/exploratory gene-condition records agreeing in direction after contrast harmonization | 354 / 354 |

### 1. The transcriptome is broadly recovery-oriented, but the degree of residual expression differs strongly by stress

Condition-level mean absolute displacement from control shows that endpoint response and recovery-state residual are not the same property.

| Stress | Endpoint mean absolute shift | Recovery mean absolute shift | Residual fraction |
|---|---:|---:|---:|
| heat+drought | 0.654 | 0.063 | 0.096 |
| drought | 0.458 | 0.120 | 0.262 |
| cold+high light | 0.391 | 0.097 | 0.249 |
| heat+drought+high light | 0.290 | 0.062 | 0.214 |
| cold | 0.275 | 0.050 | 0.184 |
| high light | 0.130 | 0.100 | 0.770 |
| high humidity | 0.116 | 0.051 | 0.443 |
| heat | 0.107 | 0.049 | 0.460 |
| heat+high light | 0.095 | 0.060 | 0.627 |

The residual fraction is a descriptive ratio and is interpreted alongside the absolute endpoint and recovery shifts.

Two conditions illustrate the contrast particularly clearly:

- **heat+drought** has the largest endpoint displacement but a small residual fraction after recovery;
- **high light** has a much smaller endpoint displacement but a substantially larger residual fraction.

The appropriate conclusion is therefore not that endpoint magnitude formally predicts persistence, but that **larger endpoint displacement does not consistently correspond to greater persistence after recovery**.

> **Transcriptomic displacement from control at endpoint and recovery**  
>![Transcriptomic displacement from control at endpoint and recovery](visualizations/transcriptomic_displacement_from_control_at_endpoint_and_recovery.png)  
> Shows condition-level mean absolute rlog displacement from control at endpoint and recovery.

> **Stress response magnitude versus residual transcriptomic fraction**  
>![Stress response magnitude versus residual transcriptomic fraction](visualizations/stress_response_magnitude_vs_residual_transcriptomic_fraction.png)  
> Shows why a larger endpoint displacement does not necessarily correspond to a larger proportional residual after recovery.

### 2. Gene-level recovery is strongly stress-dependent

At the primary p95 response floor:

| Condition | Endpoint eligible | Recovered | Persistent same direction | Substantial reversal |
|---|---:|---:|---:|---:|
| cold | 8,915 | 8,838 | 71 | 6 |
| cold+high light | 12,210 | 11,677 | 187 | 346 |
| drought | 11,942 | 10,613 | 1,240 | 89 |
| heat | 1,310 | 1,226 | 70 | 14 |
| heat+drought | 14,088 | 13,690 | 344 | 54 |
| heat+drought+high light | 8,609 | 8,323 | 255 | 31 |
| heat+high light | 814 | 729 | 72 | 13 |
| high humidity | 1,712 | 1,661 | 28 | 23 |
| high light | 1,950 | 1,154 | 779 | 17 |

For every condition, the three recovery categories sum exactly to the number of qualifying endpoint responses.

Persistent fractions among qualifying endpoint responders vary sharply. Selected examples:

```text
high light           39.95%
drought              10.38%
heat + high light     8.85%
heat                   5.34%
heat + drought         2.44%
cold                   0.80%
```

A condition can therefore have many endpoint-responsive genes yet retain only a small persistent fraction.

> **Recovery outcomes of endpoint-responsive genes**  
>![Recovery outcomes of endpoint-responsive genes](visualizations/recovery_trajectories_of_endpoint_responsive_genes.png)
> Shows the recovered, persistent-same-direction, and substantial-reversal fractions among qualifying endpoint responders for each stress condition.

### 3. Most persistent genes are stress-specific, but a recurrent multi-stress component exists

Across the 26,118 informative genes:

```text
24,509 genes were never persistent
1,609 genes were persistent in at least one stress
```

#### Persistence breadth at the primary p95 floor

| Persistence breadth | Number of genes |
|---|---:|
| 1 stress | 883 |
| 2 stresses | 350 |
| 3 stresses | 170 |
| 4 stresses | 126 |
| 5 stresses | 48 |
| 6 stresses | 17 |
| 7 stresses | 13 |
| 8 stresses | 2 |
| 9 stresses | 0 |

Thus, 883 of 1,609 persistent genes, or 54.9%, are stress-specific under the primary definition.

At the same time, 726 genes are persistent in at least two stress conditions, and 206 satisfy the primary high-breadth definition of persistence in at least four.

### 4. Multi-stress persistence is exceptionally directionally coherent

Among the 726 genes persistent in at least two stresses:

```text
720 / 726 = 99.17%
```

show the same direction relative to control across every condition in which they are persistent.

Among those 720 directionally consistent genes:

```text
711 / 720 = 98.75%
```

remain **above control**.

This is one of the most stable findings in the project. Recurrent multi-stress persistence is overwhelmingly associated with residual elevation of RNA abundance rather than a mixture of positive and negative residual states.

### 5. Persistent-gene sharing is context-dependent

Stress-specific persistent sets vary considerably in how much of their membership is shared with other conditions.

Examples of persistent genes shared with at least one other stress:

```text
heat+drought+high light    97.25%
heat+high light            97.22%
heat+drought               97.09%
drought                    53.63%
cold+high light            39.04%
```

Notable pairwise Jaccard similarities include:

```text
drought <-> high light                     ~0.410
heat+drought <-> heat+drought+high light   ~0.400
heat+drought <-> high light                ~0.351
```

The data therefore do not support either one universal residual stress response or nine completely separate responses. Instead, there is a mixture of stress-specific persistence and partially shared cross-stress recurrence.

> **Similarity of persistent gene sets across stress conditions**  
>![Similarity of persistent gene sets across stress conditions](visualizations/similarity_of_persistent_gene_programs_across_stress_conditions.png)  
> Shows pairwise Jaccard similarity among the condition-specific persistent-gene sets.

---

## Functional interpretation

### Araport11 annotation

Persistent genes were mapped to Araport11 functional descriptions at the locus level.

All 1,609 persistent genes received matched locus-level annotations. The primary high-breadth set contains 206 genes.

### GO-Slim enrichment

GO-Slim enrichment uses the 26,118 informative genes as the expression-analysis universe, with aspect-specific restriction to genes that are actually annotatable in each GO branch.

The 206-gene high-breadth foreground produces:

```text
24 significant GO-Slim terms at BH FDR < 0.05
```

Selected enriched themes include:

| GO-Slim theme | Approx. enrichment |
|---|---:|
| response to biotic stimulus | 6.87x |
| response to external stimulus | 6.41x |
| response to abiotic stimulus | 3.93x |
| response to chemical | 3.86x |
| response to stress | 3.59x |
| response to endogenous stimulus | 3.27x |
| cell wall | 3.26x |
| signal transduction | 2.48x |
| response to light stimulus | 2.36x |
| kinase activity | 2.26x |
| plasma membrane | 2.14x |
| extracellular region | 1.91x |

The enrichment of biotic- and defense-related terms does **not** imply that the experiment applied biotic stress. Plant response networks overlap extensively, and genes annotated through defense, hormone, redox, and signaling biology can also participate in abiotic-stress responses.

### Detailed Biological Process enrichment

Detailed GO Biological Process testing used 938 eligible terms after removing very small background categories and the generic Biological Process root.

The analysis identified:

```text
45 significant Biological Process terms at BH FDR < 0.05
```

Prominent themes include responses to salicylic acid, jasmonic acid, abscisic acid, water deprivation, osmotic stress, salt stress, cold, heat acclimation, oxidative stress, reactive oxygen species, hypoxia, bacterial and fungal stimuli, systemic acquired resistance, leaf senescence, and camalexin biosynthesis.

These terms are highly overlapping. They should not be read as 45 independent pathways.

---

## GO redundancy and overlap network

To avoid treating overlapping GO terms as independent biological discoveries, the project reconstructs the high-breadth genes supporting each significant Biological Process term and directly measures term overlap.

The 45 significant terms contain:

```text
388 term-gene memberships
120 unique high-breadth genes represented across the significant terms
```

All 990 pairwise term comparisons were evaluated:

```text
562 / 990 shared no foreground genes                 56.8%
428 / 990 shared at least one foreground gene       43.2%
88 / 990 had overlap coefficient >= 0.50             8.9%
17 / 990 had overlap coefficient >= 0.75             1.7%
9 / 990 had overlap coefficient = 1.00               0.9%
```

This shows that GO redundancy is strong in some local groups but is not global across the entire significant-term set.

### Strong-overlap network

A descriptive network connects two GO terms when:

```text
overlap coefficient >= 0.50
```

At the primary threshold:

```text
45 total significant terms
88 strong-overlap edges
39 terms participate in at least one strong edge
largest connected component = 35 terms
86 / 88 edges lie within that largest component
```

Weighted greedy modularity separates the connected structure into four descriptive functional neighborhoods:

1. pathogen / defense-response;
2. defense-signaling / hormone / redox / damage;
3. abiotic stress / acclimation;
4. specialized defense.

Primary modularity:

```text
Q = 0.426
```

This is a network structure score, not a p-value and not evidence that these neighborhoods are statistically proven biological modules.

The network was also examined at overlap thresholds of 0.40, 0.50, and 0.60. The central connected structure remains recognizable while peripheral terms change as the criterion becomes stricter.

> **Strong-overlap network of the largest significant Biological Process component**  
>![Strong-overlap network of the largest significant Biological Process component](visualizations/strong_overlap_network_of_the_largest_significant_bp_component.png)  
> Shows the locally connected GO terms produced by shared high-breadth foreground genes. The network summarizes annotation overlap and should not be interpreted as a causal pathway map.

### Genes spanning multiple functional neighborhoods

The four main GO-overlap neighborhoods contain 102 unique high-breadth foreground genes.

Community breadth:

```text
55 genes occur in 1 community
34 genes occur in 2 communities
11 genes occur in 3 communities
2 genes occur in all 4 communities
```

Thirteen genes occur in at least three communities while also showing broad experimental persistence:

```text
AT3G56400
AT4G02380
AT1G21250
AT3G26830
AT3G49120
AT2G40750
AT3G01420
AT1G73260
AT2G38470
AT1G15520
AT1G19020
AT1G80840
AT5G45630
```

This set includes genes annotated as WRKY transcription factors, a senescence-associated gene, a cell-wall-associated kinase, cytochrome P450 proteins, peroxidases, and other stress-associated proteins.

These genes are useful because they combine two different forms of breadth:

- persistence across several experimental stress conditions;
- representation across several GO-overlap neighborhoods.

They are **not** automatically master regulators or causal controllers. The 13-gene set is deliberately selected for these properties and must not be treated as representative of the transcriptome.

> **Recovery-state profiles of cross-community persistent genes**  
>![Recovery-state profiles of cross-community persistent genes](visualizations/recovery_state_profiles_of_cross_community_persistent_genes.png)  
> Shows the gene-by-stress recovery classifications of the selected 13 genes and makes their strong context dependence visible.

---

## Published DESeq2 concordance

The rlog-based recovery classification is descriptive. To ask whether it captures a residual expression signal also visible under formal count-based analysis, the project compares it with the source study's published recovery DEG tables.

### Available recovery DEG lists

Exact gene-level recovery DEG lists are available in the supplementary workbook for five conditions:

| Condition | Published recovery DEG records |
|---|---:|
| cold | 5 |
| drought | 567 |
| heat+drought | 3 |
| high humidity | 2 |
| high light | 2 |

Total:

```text
579 gene-condition records
```

All 579 occur within the project's 26,118-gene informative universe.

Exact recovery gene lists are not available in the workbook for:

```text
cold+high light
heat
heat+drought+high light
heat+high light
```

These conditions are treated as **unavailable**, not as zero-DEG conditions.

### Drought provides the strongest overlap test

For drought:

```text
Exploratory persistent genes:                1,240
Published recovery DEGs:                       567
Shared genes:                                  348
Expected overlap by chance:                  ~26.9
Overlap enrichment:                          ~12.9x
Odds ratio:                                   ~43.9
BH-adjusted p-value:                       ~2.47e-323
```

Approximately 61.4% of the published drought recovery DEG set overlaps the exploratory drought persistent set, while 28.1% of the exploratory persistent set also meets the published DEG criteria.

The two methods are not expected to produce identical gene sets because they use different definitions: one is based on descriptive rlog displacement relative to an empirical response floor, while the other comes from formal count-based differential-expression analysis.

### Contrast orientation and direction agreement

The supplementary workbook uses the opposite fold-change sign convention from the treatment-minus-control orientation used in this project.

The original published values are preserved unchanged. A separate sign-reoriented column is created only for direction comparisons.

After harmonizing contrast orientation:

```text
354 / 354 shared gene-condition records agree in direction
```

The reoriented published fold changes show a moderate positive correlation with the exploratory recovery shifts:

```text
r ~ 0.52
```

This is a **methodological concordance check within the same experiment**, not independent biological validation.

The published workbook contains filtered significant DEG lists rather than complete model output for every tested gene, so the comparison is necessarily limited to the reported records.

---

## Robustness analysis

The primary p95 response floor is operational rather than biologically unique. The robustness notebook therefore repeats the classification with more permissive and more stringent control-derived thresholds:

```text
p90 = 0.272532 rlog units
p95 = 0.342263 rlog units   <- primary analysis
p99 = 0.554173 rlog units
```

Before any sensitivity comparison, the p95 implementation reproduced the original recovery classification exactly.

### Condition-level persistence remains recognizable

Persistent fractions among qualifying endpoint responders:

| Condition | p90 | p95 | p99 |
|---|---:|---:|---:|
| cold | 1.37% | 0.80% | 0.60% |
| cold+high light | 2.98% | 1.53% | 0.56% |
| drought | 14.49% | 10.38% | 5.96% |
| heat | 5.20% | 5.34% | 5.91% |
| heat+drought | 3.18% | 2.44% | 1.14% |
| heat+drought+high light | 3.92% | 2.96% | 1.93% |
| heat+high light | 7.82% | 8.85% | 11.48% |
| high humidity | 2.24% | 1.64% | 1.49% |
| high light | 33.30% | 39.95% | 62.24% |

High light has the highest persistent fraction at all three thresholds.

The number of persistent genes changes substantially with the threshold, as expected. The relative differences among stress conditions are much more stable.

### Stress-condition ranking is highly stable

Spearman rank correlations of persistent fractions are:

```text
p90 vs p95: rho = 0.983
p95 vs p99: rho = 0.950
p90 vs p99: rho = 0.900
```

Thus, even the p90-to-p99 comparison preserves a strong relative ordering of conditions.

### Stricter thresholds narrow rather than replace persistent sets

For every stress condition:

```text
p99 persistent genes subset of p95 persistent genes subset of p90 persistent genes
```

This nesting follows from the classification rule itself and is therefore mainly a consistency check.

Retention quantifies how strongly the sets shrink. For example:

```text
p90 -> p95 retention
drought          64.1%
high light       68.7%
heat+drought     71.7%

p95 -> p99 retention
drought          40.1%
high light       46.3%
```

The stricter thresholds progressively retain genes that continue to satisfy the stronger endpoint and recovery displacement criteria rather than replacing the persistent set with an unrelated collection.

### A recurrent multi-stress component remains under alternative floors

Number of genes persistent in at least four stresses:

```text
p90: 344
p95: 206
p99:  62
```

At p95:

```text
persistent in at least 3 stresses: 376
persistent in at least 4 stresses: 206
persistent in at least 5 stresses:  80
```

This supports using `>=4` as an intermediate operational prioritization rule rather than a biological discontinuity.

The maximum breadth is threshold-dependent:

```text
p90 maximum breadth: 9 stresses
p95 maximum breadth: 8 stresses
p99 maximum breadth: 7 stresses
```

Therefore, exact statements about the largest observed breadth depend on the response floor, while the existence of recurrent multi-stress persistence is much more stable.

### Directional consistency is highly robust

Among genes persistent in at least two stresses:

| Floor | Multi-stress genes | Directionally consistent | Of consistent genes, above control |
|---|---:|---:|---:|
| p90 | 1,124 | 98.04% | 98.19% |
| p95 | 726 | 99.17% | 98.75% |
| p99 | 305 | 99.67% | 99.67% |

The dominant above-control direction is therefore not specific to the p95 threshold.

### Broad functional interpretation also survives threshold changes

For genes persistent in at least four stresses:

```text
p90 foreground: 344 genes -> 27 significant GO-Slim terms
p95 foreground: 206 genes -> 24 significant GO-Slim terms
p99 foreground:  62 genes -> 14 significant GO-Slim terms
```

Of the 24 GO-Slim terms significant in the primary p95 analysis:

```text
23 remain enriched above background at all three floors
12 remain at least 2-fold enriched at all three floors
12 remain FDR-significant at all three floors
```

The most stable themes include stress and external-stimulus responses, abiotic and chemical response, endogenous-stimulus response, signal transduction, plasma membrane, extracellular region, cell wall, and secondary metabolism.

Not every individual GO-Slim term is equally robust. Some terms lose formal FDR significance as the p99 foreground shrinks to 62 genes. This is why functional robustness is interpreted using both effect size and adjusted significance rather than significance status alone.

### Robustness conclusion

The sensitivity analysis does **not** show that the exact number or identity of persistent genes is universal.

It shows that:

> **The exact membership and size of persistent sets depend on the operational response floor, while the main stress-specific ordering, recurrent multi-stress persistence, directional pattern, and broad functional interpretation are considerably more stable.**

---

## Interpretation guardrails

The project deliberately uses conservative biological language.

### Persistent RNA abundance is not stress memory

A residual expression difference after recovery does not establish:

- physiological memory;
- priming;
- chromatin memory;
- altered future responsiveness;
- causal adaptation.

Those claims require additional experiments.

### RNA abundance is not continued transcription

Steady-state RNA abundance reflects both RNA production and RNA turnover. A persistent rlog difference cannot show that transcription itself remains continuously active.

### RNA abundance is not protein activity

Transcript abundance does not guarantee corresponding changes in protein abundance, post-translational state, enzyme activity, metabolite concentration, or phenotype.

### GO enrichment does not prove a pathway mechanism

GO terms are hierarchical and overlapping. Enrichment identifies overrepresented annotations, not experimentally demonstrated causal pathways.

### GO-overlap communities are descriptive

The network communities summarize patterns of shared foreground genes among enriched GO terms. No randomized-network null model was used, so they are not treated as statistically proven biological modules.

### Cross-community genes are not automatically master regulators

A gene appearing across several annotation neighborhoods is broadly represented in the enrichment structure. That does not prove that it coordinates or causally controls those functions.

### Published DESeq2 agreement is not external validation

The exploratory analysis and the published DESeq2 results come from the same underlying experiment. Their agreement supports methodological consistency, but independent datasets are still required to assess generalizability.

---

## Limitations

### Processed rlog matrix rather than raw counts

The primary analysis uses the source study's processed rlog matrix. This is appropriate for PCA, clustering, correlations, transformed-expression comparisons, and descriptive recovery analysis, but it is not a replacement for formal count-based differential-expression modeling.

### No raw-read reprocessing

FASTQ-level QC was not recomputed. This project therefore does not independently reassess per-base read quality, adapter contamination, mapping rate, library complexity, or gene-assignment quality.

### Endpoint and recovery are cross-sectional groups

The analysis compares separate day-18 cohorts. Recovery categories and the reversibility index do not track the same individual plants through time.

### Atypical triple-stress endpoint replicate retained

`endpoint_T 2` is unusual in correlation and PCA analyses but was retained because the available evidence did not justify deletion. The complete downstream analysis was not rerun after removing this sample, so conclusions involving the triple-stress endpoint group retain that sensitivity caveat.

### The empirical response floor is not a formal uncertainty model

The pooled control-derived threshold is intentionally simple and transparent. It does not model gene-specific dispersion or condition-specific replicate uncertainty.

### GO results depend on annotation

Functional enrichment depends on annotation completeness, specificity, ontology structure, and existing biological knowledge.

### The selected 13-gene set is biased by design

Those genes were chosen because they already combine broad experimental persistence with broad functional-community representation. Their behavior cannot be generalized to the whole transcriptome.

### Published recovery DEG information is incomplete

The supplementary workbook contains filtered significant DEG lists rather than complete DESeq2 results for all genes. Exact gene-level recovery lists are available for only five of the nine stress conditions.

### Independent biological replication is still needed

The strongest next scientific validation would test the findings in independent, suitably comparable Arabidopsis recovery experiments while accounting for tissue, developmental stage, stress intensity, exposure duration, recovery duration, sequencing protocol, and expression-quantification method.

---

## Scientific notebooks

The scientific phase is organized into six notebooks. The filenames intentionally do not use numeric prefixes.

| Notebook | Role | Main question |
|---|---|---|
| [`data_exploration.ipynb`](notebooks/data_exploration.ipynb) | Matrix audit and experimental-design reconstruction | What data are present, and do they match the experiment? |
| [`quality_control.ipynb`](notebooks/quality_control.ipynb) | Sample-level QC, correlations, PCA, clustering, replicate review | Can the processed expression structure support downstream comparison? |
| [`exploratory_analysis.ipynb`](notebooks/exploratory_analysis.ipynb) | Endpoint/recovery shifts, response floor, recovery categories, persistence breadth | How do recovery and persistence differ among stresses? |
| [`biological_interpretation.ipynb`](notebooks/biological_interpretation.ipynb) | Araport11 annotation, GO enrichment, GO overlap, descriptive communities, candidate genes | What functions characterize recurrent persistence? |
| [`published_de_validation.ipynb`](notebooks/published_de_validation.ipynb) | Audit and comparison with published DESeq2 recovery results | Is the exploratory residual signal consistent with the source study's formal analysis? |
| [`robustness_check.ipynb`](notebooks/robustness_check.ipynb) | p90/p95/p99 sensitivity, ranking, breadth, direction, GO-Slim robustness | Do the main conclusions depend strongly on the chosen response floor? |

Later notebooks consume explicit files generated by earlier analyses rather than depending on hidden interactive kernel state.

---

## Repository structure

The tree below is intended to let a reader see where the scientific inputs, notebooks, derived results, figures, and future engineering components belong. It separates the **current scientific-phase structure** from areas that are explicitly marked as **planned for Stage 2**.

<details open>
<summary><strong>Full repository tree</strong></summary>

```text
Arabidopsis-Stress-Recovery-Atlas/
├── README.md
├── LICENSE
├── .gitignore
├── pyproject.toml
├── requirements.txt
│
├── data/
│   ├── raw/
│   │   └── GSE226105_rLog-Normalized_counts.csv
│   │
│   ├── external/
│   │   ├── GSE226105_supplementary/
│   │   │   └── mmc16.xlsx
│   │   └── arabidopsis_annotation/
│   │       ├── Araport11_functional_descriptions_20240931.txt
│   │       └── ATH_GO_GOSLIM.txt
│   │
│   ├── interim/
│   │   ├── sample_qc_flags.csv
│   │   │
│   │   ├── exploratory_analysis/
│   │   │   ├── condition_response_magnitude.csv
│   │   │   ├── recovery_category_counts.csv
│   │   │   ├── recovery_category_percent.csv
│   │   │   ├── gene_recovery_categories.csv
│   │   │   ├── persistent_gene_table.csv
│   │   │   ├── high_breadth_persistent_genes.csv
│   │   │   ├── condition_persistence_summary.csv
│   │   │   ├── persistent_jaccard_similarity.csv
│   │   │   └── persistent_gene_intersections.csv
│   │   │
│   │   ├── biological_interpretation/
│   │   │   ├── persistent_gene_annotation.csv
│   │   │   ├── high_breadth_go_slim_enrichment_all.csv
│   │   │   ├── high_breadth_go_slim_enrichment_significant.csv
│   │   │   ├── high_breadth_bp_enrichment_all.csv
│   │   │   ├── high_breadth_bp_enrichment_significant.csv
│   │   │   ├── significant_bp_gene_term_membership.csv
│   │   │   ├── significant_bp_term_overlap_all.csv
│   │   │   ├── significant_bp_term_overlap_strong.csv
│   │   │   ├── significant_bp_network_components.csv
│   │   │   ├── significant_bp_network_communities.csv
│   │   │   ├── significant_bp_community_density.csv
│   │   │   ├── significant_bp_between_community_density.csv
│   │   │   ├── significant_bp_community_gene_overlap.csv
│   │   │   ├── significant_bp_gene_community_breadth.csv
│   │   │   ├── high_breadth_cross_community_genes.csv
│   │   │   ├── high_breadth_cross_community_recovery_matrix.csv
│   │   │   └── high_breadth_cross_community_recovery_summary.csv
│   │   │
│   │   ├── published_de_validation/
│   │   │   ├── published_de_workbook_summary.csv
│   │   │   ├── published_de_condition_counts.csv
│   │   │   ├── persistent_vs_published_recovery_degs.csv
│   │   │   ├── published_recovery_overlap_enrichment.csv
│   │   │   ├── published_recovery_direction_validation.csv
│   │   │   └── published_recovery_direction_summary.csv
│   │   │
│   │   └── robustness_checks/
│   │       ├── response_floor_sensitivity.csv
│   │       ├── persistent_count_sensitivity.csv
│   │       ├── persistent_fraction_sensitivity.csv
│   │       ├── recovered_fraction_sensitivity.csv
│   │       ├── persistent_fraction_rank_stability.csv
│   │       ├── persistent_fraction_ranks.csv
│   │       ├── persistent_set_nesting.csv
│   │       ├── persistent_gene_retention.csv
│   │       ├── persistence_breadth_sensitivity.csv
│   │       ├── directional_consistency_sensitivity.csv
│   │       ├── go_slim_significant_counts_by_floor.csv
│   │       ├── go_slim_significant_term_overlap.csv
│   │       ├── primary_go_slim_term_sensitivity.csv
│   │       ├── go_slim_enrichment_p90_all.csv
│   │       ├── go_slim_enrichment_p95_all.csv
│   │       ├── go_slim_enrichment_p99_all.csv
│   │       ├── go_slim_enrichment_p90_significant.csv
│   │       ├── go_slim_enrichment_p95_significant.csv
│   │       └── go_slim_enrichment_p99_significant.csv
│   │
│   └── processed/
│
├── notebooks/
│   ├── data_exploration.ipynb
│   ├── quality_control.ipynb
│   ├── exploratory_analysis.ipynb
│   ├── biological_interpretation.ipynb
│   ├── published_de_validation.ipynb
│   └── robustness_check.ipynb
│
├── visualizations/
├── distribution_of_expression_values_across_all_genes_and_samples.png
├── endpoint_samples_in_pca_space_by_stress_condition.png
├── expression_distributions_across_individual_samples.png
├── hierarchial_clustering_of_all_samples.png
├── hierarchial_clustering_of_endpoint_samples.png
├── hierarchial_clustering_of_persistent_transcriptional_programs.png
├── hierarchial_clustering_of_significant_bp_terms_based_on_foreground_gene_overlap.png
├── pca_of_the_500_most_variable_genes_by_stress_conditions.png
├── pca_of_the_500_most_variable_genes_with_endpoint_t2.png
├── pca_using_all_26118_informative_genes.png
├── recovery_outcomes_of_functionally_broad_persistent_genes_across_stress_conditions.png
├── recovery_state_profiles_of_cross_community_persistent_genes.png
├── recovery_trajectories_of_endpoint_responsive_genes.png
├── sample_to_sample_expression_correlation_heatmap.png
├── similarity_of_persistent_gene_programs_across_stress_conditions.png
├── stress_response_magnitude_vs_residual_transcriptomic_fraction.png
├── strong_overlap_network_of_the_largest_significant_bp_component.png
└── transcriptomic_displacement_from_control_at_endpoint_and_recovery.png
│
├── src/
│   └── arabidopsis_stress_recovery_atlas/          # Stage 2 planned
│       ├── ingestion/
│       ├── validation/
│       ├── transformation/
│       ├── annotation/
│       ├── analysis/
│       ├── database/
│       └── reporting/
│
├── tests/                                           # Stage 2 planned
│   ├── unit/
│   ├── integration/
│   └── data_quality/
│
├── sql/                                             # Stage 2 planned
│   ├── schema/
│   ├── transformations/
│   └── queries/
│
├── docs/                                            # Stage 2 / extended documentation
│   ├── biology/
│   ├── data_dictionary/
│   ├── architecture/
│   ├── database/
│   └── decisions/
│
├── reports/
│   ├── scientific/
│   ├── data_audit/
│   └── technical/
│
└── scripts/                                         # Stage 2 planned
    ├── download_data.py
    ├── run_validation.py
    └── run_pipeline.py
```

## Reproducing the scientific analysis

### Environment

The completed scientific phase uses Python 3.10.6.

Clone the repository:

```bash
git clone https://github.com/WizardNox/Arabidopsis-Stress-Recovery-Atlas.git
cd Arabidopsis-Stress-Recovery-Atlas
```

Create and activate a virtual environment.

Windows:

```powershell
python -m venv .venv
.venv\Scripts\activate
```

macOS/Linux:

```bash
python -m venv .venv
source .venv/bin/activate
```

Install the repository dependencies:

```bash
python -m pip install --upgrade pip
pip install -r requirements.txt
```

Launch JupyterLab:

```bash
jupyter lab
```

The environment includes the scientific Python stack used throughout the project, including pandas, NumPy, SciPy, statsmodels, scikit-learn, matplotlib, NetworkX, openpyxl, JupyterLab, and related utilities. `requirements.txt` should be treated as the package source of truth for the current repository.

### Required source inputs

The scientific notebooks expect the processed expression matrix at:

```text
data/raw/GSE226105_rLog-Normalized_counts.csv
```

Functional interpretation uses Arabidopsis annotation resources under:

```text
data/external/arabidopsis_annotation/
```

including the ATH GO/GO-Slim annotation file:

```text
data/external/arabidopsis_annotation/ATH_GO_GOSLIM.txt
```

and the Araport11 functional-description resource used by the biological-interpretation notebook.

Published differential-expression concordance uses the source-study supplementary workbook at:

```text
data/external/GSE226105_supplementary/mmc16.xlsx
```

These public/source-study files retain their original provenance and usage conditions. Automated source acquisition and provenance capture are part of the planned engineering stage; the current scientific workflow expects the required source files to be placed locally.

### Execution order

Run the notebooks in this order:

```text
data_exploration.ipynb
    ->
quality_control.ipynb
    ->
exploratory_analysis.ipynb
    ->
biological_interpretation.ipynb
    ->
published_de_validation.ipynb
    ->
robustness_check.ipynb
```

### Clean-run standard

Before a scientific notebook is considered complete:

1. restart the kernel;
2. run every cell from top to bottom;
3. verify that the expected outputs are regenerated;
4. confirm that no hidden notebook state is required;
5. commit the notebook together with the intended derived outputs.

This standard was applied to the completed scientific notebooks.

---

## Generated outputs

Derived results are saved under `data/interim/` so that later analyses can use explicit artifacts rather than values left in notebook memory.

The main output areas are:

```text
data/interim/
├── exploratory_analysis/
├── biological_interpretation/
├── published_de_validation/
├── robustness_checks/
├── sample_manifest.csv
└── sample_qc_flags.csv
```

These directories contain the gene-level recovery classifications, persistence summaries, overlap tables, enrichment results, GO term-gene membership, network/community tables, published-DE concordance results, and threshold-sensitivity outputs used in the scientific interpretation.

The full file-level layout is shown in [Repository structure](#repository-structure). The output files are deliberately more detailed than the README: the README summarizes the scientific story, while the notebooks and exported tables preserve the complete computational record.

---

## Next stage: analytical data platform

The scientific phase establishes what the project needs to compute and which results are biologically meaningful.

The next stage will move the stable parts of that workflow out of notebooks and into a reusable analytical data system. The goal is not to add infrastructure for its own sake. It is to make the current analysis easier to rerun, validate, query, and maintain.

Planned priorities are:

- controlled source ingestion and provenance;
- immutable raw inputs;
- schema and biological-design validation;
- deterministic sample and identifier normalization;
- reusable recovery-classification and enrichment functions;
- stable analytical tables;
- automated data-quality tests;
- queryable storage for samples, genes, recovery states, enrichment results, and published-DE comparisons;
- reproducible reporting from the stored analytical results.

A future API or interactive interface will be added only if it solves a real access need.

The key engineering requirement is simple: the new pipeline should reproduce the already established scientific results from a clean environment rather than silently changing the scientific analysis.

---

## Future scientific work

Potential follow-up work includes:

- testing the persistent-expression patterns in independent Arabidopsis stress-recovery datasets;
- comparing with additional transcriptomic studies when tissue, timing, and stress design are sufficiently compatible;
- performing formal analysis from raw counts where suitable data are available;
- investigating high-breadth persistent candidates experimentally or with independent molecular evidence;
- integrating protein, metabolite, or phenotype data when biologically justified;
- testing stress memory or priming only in experiments explicitly designed to measure altered future responsiveness.

---

## Data provenance and attribution

This project builds on public biological data and annotation resources that remain attributable to their original authors and repositories.

### Primary study

Garcia-Molina A, Pastor V. **Systemic analysis of metabolome reconfiguration in Arabidopsis after abiotic stressors uncovers metabolites that modulate defense against pathogens.** _Plant Communications_. 2024;5:100645. DOI: `10.1016/j.xplc.2023.100645`.

### Public data

- NCBI GEO: `GSE226105`
- BioProject: `PRJNA938512`

### RNA-seq processing reported by the source study

The original study performed the upstream RNA-seq processing and formal differential-expression analysis using tools including:

- Trimmomatic;
- STAR;
- featureCounts;
- DESeq2.

This repository begins from the public processed rlog matrix and published supplementary results rather than claiming those upstream steps as part of the present reanalysis.

### Functional annotation

The project uses:

- Araport11 functional descriptions;
- Gene Ontology / `ATH_GO_GOSLIM`.

External datasets, supplementary tables, and annotation resources retain their original provenance and licensing or usage conditions.

---

## Citation and license

If you use this repository, please cite the original source study and this project where relevant.

Suggested project citation until a fixed archival release is available:

```text
Sinha, S. Arabidopsis Stress Recovery Atlas:
computational analysis of context-dependent transcriptomic recovery
and persistent expression in Arabidopsis thaliana.
GitHub repository.
https://github.com/WizardNox/Arabidopsis-Stress-Recovery-Atlas
```

See the repository's `LICENSE` file for the license governing original project code and documentation.

Public datasets, source-study supplementary files, and external annotation resources remain subject to their original terms and do not automatically inherit the repository's code license.

---

## Final perspective

The **Arabidopsis Stress Recovery Atlas** follows one continuous scientific question from data audit to biological interpretation:

```text
public expression data
    ->
experimental design
    ->
quality control
    ->
endpoint and recovery comparison
    ->
gene-level recovery states
    ->
cross-stress persistence
    ->
functional interpretation
    ->
GO redundancy analysis
    ->
published-method concordance
    ->
robustness testing
```

The central conclusion is deliberately conservative:

> **The GSE226105 transcriptome is broadly consistent with recovery after abiotic stress, but that broad recovery contains smaller, structured, predominantly above-control residual expression patterns whose recurrence depends strongly on stress context.**

These patterns are recurrent, functionally structured, consistent with the available published recovery-DEG results, and robust to substantial changes in the empirical response threshold. They are interpreted as **persistent RNA-abundance states**, not as proof of physiological stress memory, continued transcription, or causal regulatory mechanisms.
