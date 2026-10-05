# Theory verification scripts

Monte Carlo and exact checks for the identities in the theory supplement
(`papers/theory/theory-supplement.md` in the project files). Python 3 with numpy.

| Script | Supplement section | Checks |
|---|---|---|
| `checks/check_master_identity.py` | S1.1–S1.4 | E[z_ka z_lb] = c_ab r_kl + t_ab (R Σ̃ R)_kl from individual-level data with LD, sample overlap and non-normal effects |
| `checks/check_factor2.py` | S1.5 | 2022 cross-block toy simulation: slope = ρ with pairs counted once, ρ/2 with the chapter-6 form |
| `checks/check_lrcq.py` | S1.6, S3 | fast OLS/WLS = brute force; bias D⁻¹diag(RCR) from local ρ; s²-weighted estimand; intercept non-identifiability |
| `checks/check_lrcp.py` | S4 | cross-block regression and moment estimators of ρ√(w_i w_j), sandwich SEs, within-block Σ̃ estimator |
| `checks/check_crosstrait.py` | S3.10 | exact variance of LRCQ/LRCP estimators with correlated, overlapping traits |
| `checks/check_distal_gls.py` | S4.6 | collapsed GLS for distal SNP sets: identity, projection bias, prescreen calibration |
| `checks/check_tight_ld.py` | S3.11, S4.7 | ridge vs tag-set LRCQ under tight LD; tag-projected LRCP correlations |
| `checks/check_category.py` | S7 | category enrichment, annotation regression, phenotype-category contrasts |
| `checks/check_category_tight_ld.py` | S7.2a | categories that split tight LD clumps: clump-level vs annotation-regression vs ridge estimands |
| `checks/check_contrast_plugin.py` | S7.4a | plug-in variance of phenotype-category contrasts under sparse w: naive, debiased, bootstrap (`--boot`) |
| `checks/check_gene_level.py` | S8 | gene-level burden covariance; conditional Gaussian test and max-T FWER under trait dependence |
| `checks/check_moments.py` | S2.4, S3.4, S5 | Isserlis fourth moments, point-normal excess kurtosis, one-mediator second/fourth moments |

| `power/power_tables.py` | S9 | analytic power tables (CSV): per-clump and category LRCQ, pair LRCP, q_eff under overlap |
| `power/spotcheck_power.py` | S9 | simulation check of (9.1) with distinct phenotypic and genetic q_eff |

Run: `cd theory/checks && python3 <script>.py`; `python3 theory/power/power_tables.py <outdir>`.
