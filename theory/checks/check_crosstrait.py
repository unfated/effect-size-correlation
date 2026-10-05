"""Cross-trait dependence: exact variance of LRCQ / LRCP linear estimators (supplement §S3.10).

Z has Cov(vec Z) = Gamma ⊗ R + T ⊗ G with Gamma = intercept matrix (c_ab),
T = (t_ab), G = R Sigma~ R. For an estimator theta = sum_a s_a z_a' A z_a / K
(A symmetric), Gaussian Z gives

  Var(theta) = (2/K^2) [ tau_RR s'(Gamma∘Gamma)s + 2 tau_RG s'(Gamma∘T)s + tau_GG s'(T∘T)s ],
  tau_XY = tr(A X A Y).

Compared with the independence formula (Gamma, T diagonal).
"""
import numpy as np
from common import ar1_ld, block_diag

rng = np.random.default_rng(17)
b, q = 25, 120
R1, R2 = ar1_ld(b, 0.75), ar1_ld(b, 0.6)
R = block_diag(R1, R2); m = 2 * b
w = np.ones(m); i, j = 8, b + 6; w[i] = w[j] = 20.0
Sig = np.diag(w); Sig[i, j] = Sig[j, i] = 0.4 * 20
G = R @ Sig @ R

# traits: 4 clusters, r_g 0.6 within, 0.1 between; full overlap (UKB-like) with phenotypic corr 0.5*r_g + 0.1
cl = np.repeat(np.arange(4), q // 4)
rg = np.where(cl[:, None] == cl[None, :], 0.6, 0.1); np.fill_diagonal(rg, 1)
n = np.full(q, 3.5e5); h2 = rng.uniform(0.05, 0.3, q); M = 1.2e6
s = n * h2 / M
T = np.sqrt(np.outer(n, n)) * rg * np.sqrt(np.outer(h2, h2)) / M
Gam = 0.5 * rg + 0.1 * (1 - np.eye(q)); np.fill_diagonal(Gam, 1.0)

LR, LG = np.linalg.cholesky(R), np.linalg.cholesky(G + 1e-9 * np.eye(m))
LGa, LT = np.linalg.cholesky(Gam), np.linalg.cholesky(T)

D = R ** 2
d = np.linalg.solve(D, np.eye(m)[i])
A_q = np.diag(d); K_q = s @ s                                   # LRCQ OLS w_i
a = np.zeros(m); a[:b] = R1[:, i]; bb = np.zeros(m); bb[b:] = R2[:, j - b]
A_p = 0.5 * (np.outer(a, bb) + np.outer(bb, a))                 # LRCP REG c_ij
K_p = (s @ s) * (R1 ** 2).sum(0)[i] * (R2 ** 2).sum(0)[j - b]


def var_exact(A, K, Gm, Tm):
    tau = lambda X, Y: np.trace(A @ X @ A @ Y)
    return 2 / K ** 2 * (tau(R, R) * s @ (Gm * Gm) @ s + 2 * tau(R, G) * s @ (Gm * Tm) @ s
                         + tau(G, G) * s @ (Tm * Tm) @ s)


reps = 6000
est_q, est_p = [], []
for _ in range(reps):
    Z = LR @ rng.standard_normal((m, q)) @ LGa.T + LG @ rng.standard_normal((m, q)) @ LT.T
    est_q.append(((Z ** 2 - 1) @ s @ d) / K_q)
    est_p.append((a @ Z) * (bb @ Z) @ s / K_p)
for name, est, A, K, truth in (("LRCQ w_i", est_q, A_q, K_q, w[i]), ("LRCP c_ij", est_p, A_p, K_p, Sig[i, j])):
    ex = np.sqrt(var_exact(A, K, Gam, T))
    ind = np.sqrt(var_exact(A, K, np.diag(np.diag(Gam)), np.diag(np.diag(T))))
    print(f"{name:10s} mean {np.mean(est):7.2f} (true {truth:.1f})  empirical SD {np.std(est):.3f}  "
          f"exact cross-trait SD {ex:.3f}  independence SD {ind:.3f}  (q_eff ≈ {q * (ind / ex) ** 2:.0f})")
