# lrcpq

R package for the unified matrix-normal framework of thesis chapter 6:

- **LRCQ**: per-SNP phenome-wide genetic-effect enrichment `w` from squared Z-scores of many GWAS.
- **LRCP**: genetic-effect correlation `Rβ` between SNPs (or genes) from Z-score cross-products.

Model: `Y = XB + E`, `B ~ MN(0, U, V)` with `U = W^{1/2} Rβ W^{1/2}` and `V = H_g / M`; GWAS Z-scores follow `Z = R B diag(√n) + E`, `E ~ MN(0, R, C)`, where `C` holds the inter-GWAS intercepts (`r_ab · o_ab`).

## Install

```r
# from a clone of unfated/effect-size-correlation
install.packages("lrcpq", repos = NULL, type = "source")
```

## Four-stage estimator

| Stage | Parameter | Function |
|---|---|---|
| 1 | `h²`, `h_ab` | `ldsc_h2()`, `ldsc_gcov()`, `ldsc_matrix()` |
| 2 | intercepts `C` | `ldsc_matrix()` (cross-trait intercepts), `make_intercept()` for known overlap |
| 3 | `w` (LRCQ) | `lrcq()` / `lrcq_window()`; `method = "ols" / "wls" / "irls" / "equal"`; `rectify_w()` methods A, B, C |
| 4 | `Rβ` (LRCP) | `lrcp()` (in progress) |

Stage 3 uses the fast OLS/WLS solution: stacking q traits gives the design `s ⊗ D` (`D = R∘R`), so the regression collapses to one m-by-m solve per LD window, `(D Ω D) w = D u`.

## Simulation in five lines

```r
library(lrcpq)
blocks <- make_ld_ar1(rep(50, 10), rho = 0.5)          # toy LD, or real LD blocks
ld     <- ld_from_blocks(blocks)
sim    <- simulate_lrcpq(blocks, n = rep(2e4, 30), prop_nonzero = 0.1)
fit    <- lrcq(sim$Z, ld, sim$n, sim$gcov, M = ld$m, intercept = sim$intercept)
cor(fit$w_raw, sim$w)
```

## Unsettled maths, kept switchable

- `pair_sum = "once"` (default, confirmed by the theory supplement) counts each SNP pair once in the `Rβ` cross terms; `"chapter"` reproduces the thesis equations, which count each pair twice.
- `lrcq(..., cross_terms = FALSE)` is the chi-square-only LRCQ of the abstract; `TRUE` adds the cross terms for a given `Rβ`.

Intercepts are taken from stage 2 and held fixed: with a free per-trait intercept inside a window, `w` is not identifiable.
