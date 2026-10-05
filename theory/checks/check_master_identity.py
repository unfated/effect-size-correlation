"""Master moment identity from individual-level data (supplement §S1.1-S1.4).

Simulates genotypes with LD, two traits with sample overlap, phenotypic
(residual) correlation, genetic covariance and *non-normal* (point-normal)
effects with a non-trivial SNP-wise covariance Sigma~. GWAS Z-scores are
computed per SNP by marginal regression with estimated standard errors.

Checks   E[z_ka z_lb] = c_ab r_kl + t_ab (R Sigma~ R)_kl
with     t_ab = sqrt(n_a n_b) h_ab / m,   c_ab = n_ab rho_P / sqrt(n_a n_b), c_aa = 1.
"""
import numpy as np
from common import ar1_ld

rng = np.random.default_rng(11)
m = 6
n_a, n_b, n_ab = 3000, 2000, 1200
R = ar1_ld(m, 0.6)
Lr = np.linalg.cholesky(R)
# SNP-wise covariance Sigma~ (trace m): unequal w and one correlated pair (0,4)
w = np.array([2.0, 0.5, 0.0, 1.5, 1.5, 0.5])
Rb = np.eye(m); Rb[0, 4] = Rb[4, 0] = 0.7
Sig = np.sqrt(w)[:, None] * Rb * np.sqrt(w)[None, :]
# trait genetic covariance: per-SNP-scale h2 and h_ab (sum_k E beta_k^2 = h2)
h2 = np.array([0.04, 0.03]); rg = 0.5
hab = rg * np.sqrt(h2[0] * h2[1])
V = np.array([[h2[0], hab], [hab, h2[1]]]) / m
# non-normal effects with Cov(vec B) = V ⊗ Sig exactly: B = Sig^{1/2} G V^{1/2}', G iid point-normal
pi = 0.3
Sh = np.linalg.cholesky(Sig + 1e-12 * np.eye(m))
Vh = np.linalg.cholesky(V)
rho_P_resid = 0.4

n_tot = n_a + n_b - n_ab
reps = 6000
acc = np.zeros((2 * m, 2 * m)); acc2 = np.zeros((2 * m, 2 * m))
for _ in range(reps):
    G = rng.standard_normal((m, 2)) * (rng.random((m, 2)) < pi) / np.sqrt(pi)   # point-normal, unit variance
    B = Sh @ G @ Vh.T
    X = rng.standard_normal((n_tot, m)) @ Lr.T
    X = (X - X.mean(0)) / X.std(0)
    g = X @ B
    e = rng.standard_normal((n_tot, 2))
    e[:, 1] = rho_P_resid * e[:, 0] + np.sqrt(1 - rho_P_resid ** 2) * e[:, 1]
    sd_e = np.sqrt(1 - h2)                           # residual scale so Var(y) ~ 1
    Y = g + e * sd_e
    ia = np.arange(0, n_a)                            # trait a: first n_a people
    ib = np.arange(n_a - n_ab, n_tot)                 # trait b: last n_b people; overlap n_ab
    zs = []
    for idx, col in ((ia, 0), (ib, 1)):
        Xs = X[idx]; Xs = (Xs - Xs.mean(0)) / Xs.std(0)
        y = Y[idx, col]; y = (y - y.mean()) / y.std()
        nn = len(idx)
        bhat = Xs.T @ y / nn
        se = np.sqrt((1 - bhat ** 2) / (nn - 2))
        zs.append(bhat / se)
    z = np.concatenate(zs)
    acc += np.outer(z, z); acc2 += np.outer(z, z) ** 2
emp = acc / reps

RSR = R @ Sig @ R
rhoP = rho_P_resid * np.sqrt((1 - h2[0]) * (1 - h2[1])) + hab   # total phenotypic covariance (Var(y)≈1)
t = {(0, 0): n_a * h2[0] / m, (1, 1): n_b * h2[1] / m, (0, 1): np.sqrt(n_a * n_b) * hab / m}
c = {(0, 0): 1.0, (1, 1): 1.0, (0, 1): n_ab * rhoP / np.sqrt(n_a * n_b)}
theo = np.zeros_like(emp)
for (a, b), tv in t.items():
    blk = c[(a, b)] * R + tv * RSR
    theo[a * m:(a + 1) * m, b * m:(b + 1) * m] = blk
    theo[b * m:(b + 1) * m, a * m:(a + 1) * m] = blk.T
mc_se = np.sqrt((acc2 / reps - emp ** 2) / reps)          # empirical MC se (valid without normality)
zdev = (emp - theo) / mc_se
print("max |emp - theory| / MC se over all 78 distinct moments:", np.abs(zdev).max().round(2))
print("mean (emp - theory)/MC se:", zdev.mean().round(3))
signal = theo - np.kron(np.array([[1, c[(0, 1)]], [c[(0, 1)], 1]]), R)
print("relative error of the genetic part (median over moments with signal>1): %.3f" %
      np.median(np.abs((emp - theo)[signal > 1] / signal[signal > 1])))
print("example: E z_a0 z_b4  emp=%.3f theory=%.3f" % (emp[0, m + 4], theo[0, m + 4]))
print("example: E z_a0^2     emp=%.3f theory=%.3f" % (emp[0, 0], theo[0, 0]))
# the chapter's double-counted form for comparison
ch = R @ np.diag(w) @ R + 2 * (RSR - R @ np.diag(w) @ R)
print("chapter form E z_a0^2 would be %.3f" % (1 + t[(0, 0)] * ch[0, 0]))
