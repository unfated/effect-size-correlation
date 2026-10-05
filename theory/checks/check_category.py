"""Category-level LRCQ (supplement §S7).

1. Set-average enrichment E_c = 1_c' w / |c| estimated by E_hat = 1_c' w_hat_OLS / |c|:
   unbiased; exact SD from the weighted Theorem 3.10 (trait weights u, matrix A = Diag(D^-1 1_c)/|c|).
2. Annotation regression w = A tau: tau_hat = (A'D^2A)^{-1} A'D^2 w_hat_OLS (= stacked OLS) is unbiased for
   the D^2-weighted projection (A'D^2A)^{-1}A'D^2 w (= tau when w lies in the annotation span).
3. Phenotype-category contrast Delta_c = E_c^(g1) - E_c^(g2): null calibration of the Wald test using
   the weighted variance formula, with correlated, fully overlapping traits.
"""
import numpy as np
from common import ar1_ld, block_diag

rng = np.random.default_rng(31)
nb, bs = 8, 30
R = block_diag(*[ar1_ld(bs, rng.uniform(0.5, 0.85)) for _ in range(nb)]); m = nb * bs
D = R ** 2
annot = rng.random(m) < 0.2
w = np.where(annot, 3.0, 1.0) * np.where(rng.random(m) < 0.3, 1.0, 0.0)   # sparse, enriched in annotation
w = w / w.mean()
G = R @ np.diag(w) @ R

q = 120
cl = np.repeat(np.arange(4), q // 4)
rg = np.where(cl[:, None] == cl[None, :], 0.6, 0.1); np.fill_diagonal(rg, 1)
n = np.full(q, 3.5e5); h2 = rng.uniform(0.05, 0.3, q); M = 2.0e5
s = n * h2 / M
T = np.sqrt(np.outer(n, n)) * rg * np.sqrt(np.outer(h2, h2)) / M
Gam = 0.5 * rg + 0.1 * (1 - np.eye(q)); np.fill_diagonal(Gam, 1.0)
LR, LG = np.linalg.cholesky(R), np.linalg.cholesky(G + 1e-9 * np.eye(m))
LGa, LT = np.linalg.cholesky(Gam), np.linalg.cholesky(T)
Dinv = np.linalg.inv(D)


def var_w(Amat, u):
    """Var of sum_a u_a z_a' A z_a for Cov(vec Z) = Gam ⊗ R + T ⊗ G (Gaussian)."""
    tau = lambda X, Y: np.trace(Amat @ X @ Amat @ Y)
    return 2 * (tau(R, R) * u @ (Gam * Gam) @ u + 2 * tau(R, G) * u @ (Gam * T) @ u + tau(G, G) * u @ (T * T) @ u)


c1 = annot.astype(float)
A_set = np.diag(Dinv @ c1) / c1.sum()
g1, g2 = cl < 2, cl >= 2
u_all = s / (s @ s)
u_con = np.where(g1, s / (s[g1] @ s[g1]), -s / (s[g2] @ s[g2]))
Ann = np.column_stack([np.ones(m), c1])
P = np.linalg.solve(Ann.T @ D @ D @ Ann, Ann.T @ D @ D)
tau_true = np.linalg.lstsq(Ann, w, rcond=None)[0]

reps = 4000
E, Dl, TA = [], [], []
for _ in range(reps):
    Z = LR @ rng.standard_normal((m, q)) @ LGa.T + LG @ rng.standard_normal((m, q)) @ LT.T
    Y2 = Z ** 2 - 1
    E.append(c1 @ Dinv @ (Y2 @ u_all) / c1.sum())
    Dl.append(c1 @ Dinv @ (Y2 @ u_con) / c1.sum())
    TA.append(P @ Dinv @ (Y2 @ u_all))
E, Dl, TA = np.array(E), np.array(Dl), np.array(TA)
sd_E = np.sqrt(var_w(A_set, u_all)); sd_D = np.sqrt(var_w(A_set, u_con))
indep = lambda Amat, u: np.sqrt(2 * (np.trace(Amat @ R @ Amat @ R) * u @ np.diag(np.diag(Gam ** 2)) @ u
                                   + 2 * np.trace(Amat @ R @ Amat @ G) * u @ np.diag(np.diag(Gam * T)) @ u
                                   + np.trace(Amat @ G @ Amat @ G) * u @ np.diag(np.diag(T * T)) @ u))
print(f"1. E_c true {c1 @ w / c1.sum():.3f}; mean est {E.mean():.3f}; SD emp {E.std():.3f} vs formula {sd_E:.3f} "
      f"(independence {indep(A_set, u_all):.3f})")
print(f"   per-SNP SD (median over annotated SNPs) for comparison: "
      f"{np.median([np.sqrt(var_w(np.diag(Dinv[:, i]), u_all)) for i in np.flatnonzero(annot)]):.3f}")
print(f"2. tau target (A'D^2A)^-1 A'D^2 w = {(P @ w).round(3)}; mean tau_hat {TA.mean(0).round(3)}; "
      f"ordinary LS fit of w would be {tau_true.round(3)} (differs: w is not in the annotation span)")
print(f"3. contrast null: mean {Dl.mean():.3f}; SD emp {Dl.std():.3f} vs formula {sd_D:.3f}; "
      f"type-I (|Z|>1.96) {np.mean(np.abs(Dl / sd_D) > 1.96):.3f}; with independence SD "
      f"{np.mean(np.abs(Dl / indep(A_set, u_con)) > 1.96):.3f}")
