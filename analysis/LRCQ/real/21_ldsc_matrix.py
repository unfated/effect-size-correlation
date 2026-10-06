"""Stages 1-2 for all traits and trait pairs at once ("matrix LDSC").

For common regression weights v_k (LDSC-style 1/l_k, l floored at 1) the
cross-trait LDSC regression of z_ka z_kb on [1, l_k] for every pair (a, b)
needs only two q x q matrices, Z' V Z and Z' V L Z, so all q(q+1)/2
regressions are solved together. SNP-block jackknife (200 blocks) gives SEs.
Missing Z are set to 0 and the regression uses SNPs observed in >= 99% of
traits (mean imputation of a few missing values attenuates the
corresponding products negligibly).

Outputs (out_prefix): .intercept.tsv (q x q, c_ab), .gcov.tsv (q x q, h_ab on the
observed scale with M = M_5_50), .intercept_se.tsv, .gcov_se.tsv,
.h2.tsv (per-trait h2, intercept, SEs, plus Pan-UKB manifest values).
Usage: 21_ldsc_matrix.py <zprefix> <ldscore.gz> <M> <out_prefix> [--chisq-max 80]
"""
import sys
import numpy as np
import pandas as pd

zp, ldf, M, out = sys.argv[1], sys.argv[2], float(sys.argv[3]), sys.argv[4]
chisq_max = 80.0
Z = np.load(zp + ".Z.npy", mmap_mode="r")
traits = pd.read_csv(zp + ".traits.tsv", sep="\t")
q = Z.shape[1]
ld = pd.read_csv(ldf, sep="\t")
uni = pd.read_csv("/home/user/data/ref/snp_universe.tsv", sep="\t")
l = ld.set_index("SNP").reindex(uni.ID)["L2"].values
obs = np.mean(~np.isnan(Z), axis=1)
ok = (obs >= 0.99) & np.isfinite(l)
# exclude MHC and SNPs with any extreme chi-square (LDSC convention, chi2 > max(80, 0.001 N))
mhc = (uni.CHR.values == 6) & (uni.BP.values > 25e6) & (uni.BP.values < 34e6)
ok &= ~mhc
idx = np.where(ok)[0]
N = traits.N.values.astype(float)
n_blocks = 200
blocks = np.array_split(idx, n_blocks)

def stats(rows):
    z = np.nan_to_num(np.asarray(Z[rows], dtype=np.float64))
    big = (z ** 2 > chisq_max).any(1)
    z = z[~big]; lk = np.maximum(l[rows][~big], 1.0)
    v = 1.0 / lk
    S0 = v.sum(); S1 = (v * lk).sum(); S2 = (v * lk * lk).sum()
    A = (z * v[:, None]).T @ z
    B = (z * (v * lk)[:, None]).T @ z
    return np.array([S0, S1, S2]), A, B

parts = [stats(b) for b in blocks]
S = sum(p[0] for p in parts); A = sum(p[1] for p in parts); B = sum(p[2] for p in parts)

def solve(S, A, B):
    det = S[0] * S[2] - S[1] ** 2
    icpt = (S[2] * A - S[1] * B) / det
    slope = (S[0] * B - S[1] * A) / det
    return icpt, slope

icpt, slope = solve(S, A, B)
# slope = sqrt(n_a n_b) h_ab / M  ->  h_ab
sn = np.sqrt(np.outer(N, N))
gcov = slope * M / sn
# delete-one-block jackknife
ji, jg = [], []
for p in parts:
    i_, s_ = solve(S - p[0], A - p[1], B - p[2])
    ji.append(i_); jg.append(s_ * M / sn)
ji = np.array(ji); jg = np.array(jg)
nb = len(parts)
se = lambda x: np.sqrt((nb - 1) / nb * ((x - x.mean(0)) ** 2).sum(0))
icpt_se, gcov_se = se(ji), se(jg)
ids = traits.trait_id.values
for name, mat in [("intercept", icpt), ("gcov", gcov), ("intercept_se", icpt_se), ("gcov_se", gcov_se)]:
    pd.DataFrame(mat, index=ids, columns=ids).to_csv(f"{out}.{name}.tsv", sep="\t")
# Univariate LDSC per trait with LDSC's heteroscedasticity weights
# 1 / (l * 2 (c + N h2 l / M)^2), two reweighting steps from the matrix fit.
zz = np.nan_to_num(np.asarray(Z[idx], dtype=np.float64)) ** 2
lk = np.maximum(l[idx], 1.0)
uh2, uint = np.diag(gcov).copy(), np.diag(icpt).copy()
for a in range(q):
    y = zz[:, a]; keep = y < max(chisq_max, 0.001 * N[a])
    y = y[keep]; x = lk[keep]
    c, h2 = uint[a], max(uh2[a], 1e-4)
    for _ in range(2):
        v = 1.0 / (x * 2 * (c + N[a] * h2 * x / M) ** 2)
        X = np.column_stack([np.ones_like(x), x])
        c, b = np.linalg.solve(X.T @ (X * v[:, None]), X.T @ (y * v))
        h2 = max(b * M / N[a], 1e-4)
    uh2[a], uint[a] = b * M / N[a], c
h = pd.DataFrame({"trait_id": ids, "N": N, "h2": uh2, "h2_matrix": np.diag(gcov), "h2_se": np.diag(gcov_se),
                  "intercept": uint, "intercept_matrix": np.diag(icpt), "intercept_se": np.diag(icpt_se),
                  "h2_panukb_ldsc": traits.h2_ldsc_obs.values,
                  "intercept_panukb": traits.ldsc_intercept.values})
h.to_csv(f"{out}.h2.tsv", sep="\t", index=False)
d = np.sqrt(np.clip(np.diag(gcov), 1e-6, None))
rg = gcov / np.outer(d, d)
pd.DataFrame(rg, index=ids, columns=ids).to_csv(f"{out}.rg.tsv", sep="\t")
print(f"{len(idx)} SNPs used; h2 cor with Pan-UKB LDSC: {np.corrcoef(h.h2, h.h2_panukb_ldsc)[0,1]:.3f}; "
      f"intercept cor: {np.corrcoef(h.intercept, h.intercept_panukb)[0,1]:.3f}")
print(h[["trait_id", "h2", "h2_panukb_ldsc", "intercept", "intercept_panukb"]].head(10).to_string())
