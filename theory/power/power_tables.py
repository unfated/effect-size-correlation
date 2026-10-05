"""Analytic power and sample-size tables for LRCQ and LRCP (supplement §S9).

Model per unit (tag/clump, LD resolved, r^2<0.5 tags so (D^-1)_ii ~ kappa_LD):
  s = n h2 / m per trait, ||s||^2 = q_eff * s^2 (q_eff = effective number of traits, §S3.10).
LRCQ per clump, H0: w=0 vs H1: w:
  Var(w_hat | w) = kappa_LD * 2 * (1 + s w)^2 / (q_eff s^2)
LRCQ category (|c|_eff effectively independent units, mean enrichment E vs baseline 1):
  Var(E_hat | E) = 2 * (1 + s E)^2 / (q_eff s^2 |c|_eff)
LRCP pair (tags, distal), H0: rho=0 vs H1: rho, enrichments w_i, w_j:
  Var(C_hat) = [(1 + s w_i)(1 + s w_j) + s^2 C^2] / (q_eff s^2),  C = rho sqrt(w_i w_j)
  test on C_hat; power = P(|C_hat| > z_{a/2} SE0) under H1.
Writes CSVs to the output directory given as argv[1].
"""
import csv, itertools, sys
from math import sqrt
from statistics import NormalDist

N = NormalDist()
out = sys.argv[1] if len(sys.argv) > 1 else "."
m = 1.2e6           # HapMap3-scale number of SNPs defining per-SNP heritability scale
kappa_LD = 1.4      # median (D^-1)_ii after r^2<0.5 pruning (LRCQ thread, chr22 UKB window)


def pow_one_sided(mu, se0, se1, alpha):
    return 1 - N.cdf((N.inv_cdf(1 - alpha) * se0 - mu) / se1)


def pow_two_sided(mu, se0, se1, alpha):
    z = N.inv_cdf(1 - alpha / 2)
    return (1 - N.cdf((z * se0 - mu) / se1)) + N.cdf((-z * se0 - mu) / se1)


rows = []
for q_eff, n, h2, w in itertools.product([50, 100, 300, 1000], [1e5, 4e5, 1e6], [0.05, 0.1, 0.2], [2, 5, 10, 20, 50]):
    s = n * h2 / m
    se0 = sqrt(kappa_LD * 2 / (q_eff * s * s)); se1 = sqrt(kappa_LD * 2 * (1 + s * w) ** 2 / (q_eff * s * s))
    rows.append(dict(q_eff=q_eff, n=int(n), h2=h2, s=round(s, 4), w=w, se_null=round(se0, 3),
                     power_gw=round(pow_one_sided(w, se0, se1, 5e-7), 3),
                     power_nominal=round(pow_one_sided(w, se0, se1, 0.05), 3)))
with open(f"{out}/power_lrcq_clump.csv", "w", newline="") as f:
    wr = csv.DictWriter(f, fieldnames=rows[0].keys()); wr.writeheader(); wr.writerows(rows)

rows = []
for q_eff, n, h2, c_eff, E in itertools.product([50, 100, 300], [1e5, 4e5], [0.05, 0.1, 0.2], [100, 1000, 10000], [1.25, 1.5, 2, 3]):
    s = n * h2 / m
    se0 = sqrt(2 * (1 + s) ** 2 / (q_eff * s * s * c_eff)); se1 = sqrt(2 * (1 + s * E) ** 2 / (q_eff * s * s * c_eff))
    rows.append(dict(q_eff=q_eff, n=int(n), h2=h2, c_eff=c_eff, enrichment=E, se_null=round(se0, 4),
                     power_nominal=round(pow_one_sided(E - 1, se0, se1, 0.05), 3),
                     power_bonf100=round(pow_one_sided(E - 1, se0, se1, 0.0005), 3)))
with open(f"{out}/power_lrcq_category.csv", "w", newline="") as f:
    wr = csv.DictWriter(f, fieldnames=rows[0].keys()); wr.writeheader(); wr.writerows(rows)

rows = []
for q_eff, n, h2, w, rho, K in itertools.product([100, 300, 1000], [1e5, 4e5, 1e6], [0.1, 0.2], [5, 10, 20, 50], [0.2, 0.5], [100, 1000]):
    s = n * h2 / m
    C = rho * w
    se0 = sqrt((1 + s * w) ** 2 / (q_eff * s * s)); se1 = sqrt(((1 + s * w) ** 2 + s * s * C * C) / (q_eff * s * s))
    alpha = 0.05 / (K * (K - 1) / 2)
    rows.append(dict(q_eff=q_eff, n=int(n), h2=h2, w_i=w, w_j=w, rho=rho, K_screened=K, alpha=f"{alpha:.1e}",
                     se_rho_null=round(se0 / w, 3), power=round(pow_two_sided(C, se0, se1, alpha), 3)))
with open(f"{out}/power_lrcp_pair.csv", "w", newline="") as f:
    wr = csv.DictWriter(f, fieldnames=rows[0].keys()); wr.writeheader(); wr.writerows(rows)
print("written to", out)
