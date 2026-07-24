---
title: "Dengue Lipidomics — Final Summary Report"
subtitle: "Plasma, PBMC, and Oxylipin Results Across the Full Analysis"
date: "`r Sys.Date()`"
output:
  html_document:
    toc: true
    toc_depth: 3
    toc_float: true
    number_sections: true
    theme: flatly
    highlight: tango
---

**Scope:** synthesizes results across all four analysis arms of this project — plasma cross-sectional, plasma longitudinal, PBMC, and oxylipins. Detail, methods, and full result tables live in each arm's own report (linked below); this document is the cross-arm summary and discussion.

---

# Objective

Identify lipid biomarkers of dengue severity (healthy / mild / severe) using plasma and PBMC lipidomics, plus a targeted plasma oxylipin (eicosanoid) panel, cross-sectionally and — where repeated timepoints exist — longitudinally (D0→D3→D10→D60).

# At a glance

| Arm | Compartment | Design | Headline result |
|---|---|---|---|
| Plasma (cross-sectional) | Circulating plasma | Healthy/Mild/Severe, D0 | **Healthy cleanly separates from dengue** in PCA/t-SNE/UMAP; UMAP alone also separates mild from severe |
| Plasma (longitudinal) | Circulating plasma | Same 43 mild/severe patients, D0→D60, LMM | **97/119 lipids** change over time; **23/119** show a severity-dependent trajectory, diverging mainly by **D3** |
| PBMC | Immune cell membranes | Healthy/Mild/Severe/Diseased(pooled), D0 | **No individual lipid reaches significance** in any comparison; total (non-normalized) abundance does differ (mild vs healthy, p=0.021) |
| Oxylipins | Circulating plasma (targeted eicosanoid panel) | Mild/Severe only, D0/D3/D10/D60, cross-sectional + LMM | **No lipid or class reaches significance**; one recurring nominal lead (LA-derived DiHOME/TriHOME) |

# 1. Plasma lipidome — the study's strongest signal

## Cross-sectional (D0, healthy/mild/severe)

Cross-tabulating each dimensionality-reduction method's unsupervised k-means clusters against the true clinical labels (not just eyeballing the plot) confirms a real, structured separation:

- **PCA and t-SNE**: healthy controls split into two subgroups, but every mild and every severe patient falls into one shared cluster — a clean **healthy-vs-dengue** split. Mild and severe are *not* distinguished from each other by either method.
- **UMAP**: achieves a **perfect three-way match** — cluster 1 = 100% mild, cluster 2 = 100% severe, cluster 3 = 100% healthy, with zero cross-contamination. UMAP is the only method of the three that separates severity level, not just disease status.

## Longitudinal (D0→D3→D10→D60, 43 matched mild/severe patients)

A per-lipid linear mixed model (`value ~ timepoint * severity + day_of_fever + (1|patient_id)`) avoids the pseudoreplication of treating repeated draws as independent samples:

- **97 of 119 tested lipids** show a significant time effect (FDR<0.05) — mostly a sharp drop from D0 to D10 with partial rebound by D60.
- **23 of 119** show a significant timepoint×severity interaction — a trajectory *shape* that differs by severity, not just level. These are dominated by **DG, TG, PI, and LPC** species (corroborated independently at the class level).
- The divergence is concentrated at **D3**: mild patients typically show a sharp early shift that severe patients don't make on the same schedule (or reverse later) — the strongest lead for a clinically useful early-triage signal in this dataset, since it's present well before D10/D60.
- This is discovery-stage, not a validated biomarker: no independent cohort, no discriminative-performance testing (ROC/AUC) has been done yet.

# 2. PBMC lipidome — consistently null at the individual-lipid level

Every PBMC comparison run this session came back with **zero individually significant lipids**:

| Comparison | n | Nominally-significant lipids (raw p) | FDR-significant |
|---|---|---|---|
| Healthy vs Mild | 26 | 18 | 0 |
| Healthy vs Severe | ~20 | 1 (`FA 18:1`) | 0 |
| Mild vs Severe | ~19 | 11 | 0 |
| Healthy vs Diseased (mild+severe pooled) | 34 (15H/19D) | 8 | 0 |
| Three-group ANOVA | 34 | — | 0/121 |

Two exceptions worth flagging:
- **Total (non-normalized) lipid abundance differs by group** (Kruskal-Wallis p=0.021 across 3 groups; mild vs healthy pairwise p.adj=0.021) — a signal on the *overall lipid burden* axis that Percentage normalization, used for every individual-lipid test above, specifically removes. This is the one PBMC result that clears significance anywhere in this arm.
- Unsupervised clustering (PCA/t-SNE/UMAP) shows **no correspondence to severity at all** in PBMC, unlike plasma's clean structure — cluster composition is a mix of all three groups regardless of method.

# 3. Oxylipins — null, with one recurring nominal lead

Rebuilt from the platform's raw Results sheet (47 compounds; 20 reliably detected at a 70% presence threshold — resolvins, protectins, maresin, and most minor prostaglandins/leukotrienes are below detection in this cohort and untestable regardless of sample size).

- **Cross-sectional** (D0/D3/D10/D60, mild vs severe — no healthy arm exists for oxylipins): 0/19-21 significant at every timepoint.
- **Longitudinal LMM** (same design as the plasma longitudinal analysis, which did find signal): 0/20 species-level, 0/4 family-level, 0/14 class-level significant time or interaction effects.
- **One consistent lead**: `12,13-DiHOME`/`9,10-DiHOME` and `9,12,13-TriHOME` (all LA-derived) recur as the closest-to-significant hits across three independent cross-sectional timepoints *and* the longitudinal severity test — never FDR-significant, but a repeated pattern is more credible than an isolated nominal p-value.

# 4. Why the three arms disagree

**Plasma vs. PBMC** is most plausibly a **compartment** difference, not (only) a power problem: plasma directly reflects systemic circulating lipid metabolism — exactly what dengue's vascular leakage and hepatic dysfunction perturb — while PBMC lipidome is intracellular membrane composition, a pool cells actively buffer against external change, in a mixed cell population that can dilute a real single-subtype signal. This is compounded by real, measurable power/quality differences: PBMC's matched sample count is smaller than the nominal cohort (34 of 40 — 5 healthy and 1 severe patient had insufficient material), one sample (`BS-007-RAN-02-M3`) is missing ~24% of its panel, and PBMC is explicitly the newer, less mature dataset in this project.

**Plasma vs. oxylipins** needs a different explanation, since oxylipins are plasma-derived too. The most likely drivers: (1) oxylipins have no healthy arm, so they can only ever test the *harder* mild-vs-severe contrast — even the main plasma panel's severity signal only emerged via the longitudinal design, and running that identical design on oxylipins still found nothing, arguing against "wrong test" and toward a genuinely weaker effect; (2) an order-of-magnitude smaller reliably-measured panel (20 vs. 119+); (3) eicosanoids are acutely-produced, rapidly-cleared signaling molecules — a single blood draw is more likely to catch a noisy transient moment than the stable pool differences the main lipid panel measures; (4) the molecules most mechanistically relevant to an "inflammation resolution" hypothesis (resolvins, protectins, maresin) were never even measurable on this assay.

# 5. Recommended next steps

- **Plasma**: prioritize the D3 DG/TG/PI/LPC divergence for validation — build a D0-only supervised classifier to test standalone discriminative value (ROC/AUC), and seek an independent cohort.
- **PBMC**: if this arm remains a priority, a larger cohort and/or cell-subtype-resolved (e.g. sorted monocyte) lipidomics would directly test whether bulk-PBMC dilution is masking a real subtype-specific signal; the total-abundance finding (not the individual-lipid tests) is the one positive lead worth following up here.
- **Oxylipins**: the DiHOME/TriHOME lead is the best candidate for a targeted, single-hypothesis follow-up test (no multiple-testing burden from the rest of the panel). Resolvin/protectin/maresin biology would need a dedicated SPM-targeted assay, not reanalysis of this panel.
- **Across arms**: none of the current findings constitute a validated biomarker — everything here is a within-cohort association. The gap to a clinically usable triage biomarker is discriminative-performance testing and independent validation, neither done yet for any arm.

## Related reports
- [Lipidomics_Analysis_Report.Rmd](Lipidomics_Analysis_Report.Rmd) — plasma cross-sectional, full detail
- [Longitudinal/Longitudinal_Analysis_Report.md](Longitudinal/Longitudinal_Analysis_Report.md) — plasma longitudinal, full detail
- [PBMC_Lipidomics_Analysis_Report.Rmd](PBMC_Lipidomics_Analysis_Report.Rmd) — PBMC, full detail
- [Oxylipids/Oxylipid_Analysis_Report.md](Oxylipids/Oxylipid_Analysis_Report.md) — oxylipins, full detail
