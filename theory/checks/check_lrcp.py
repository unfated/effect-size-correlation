"""LRCP checks (supplement §S4).

Cross-block pair of LD blocks (R1, R2), q traits, small per-trait signal s_a.
Estimators of c_ij = rho_ij sqrt(w_i w_j):
  REG  (prescreened single-pair OLS, the 2022 regression)
       c_hat = [R1 Z1 S Z2' R2]_ij / (||s||^2 l_i l_j),  l = LD score
  MOM  (full cross-block moment inversion)
       C_hat = R1^{-1} Z1 S Z2' R2^{-1} / ||s||^2
Predicted null SEs (c = 1 intercepts, s small):
  SE_REG = sqrt((R1^3)_ii (R2^3)_jj) / (||s|| l_i l_j)
  SE_MOM = sqrt((R1^-1)_ii (R2^-1)_jj) / ||s||
Also checks the within-block full Sigma~ estimator
  Sigma_hat = R^{-1} (Z S Z' - (sum_a c_a s_a) R) R^{-1} / ||s||^2  (unbiased).
"""
import numpy as np
from common import ar1_ld, mvn, block_diag

rng = np.random.default_rng(3)
b, q = 30, 300
R1, R2 = ar1_ld(b, 0.8), ar1_ld(b, 0.6)
R = block_diag(R1, R2)
s = rng.uniform(0.02, 0.15, q)            # n_a h2_a / m for realistic GWAS (e.g. 3e5*0.2/1e6=0.06)
w = np.ones(2 * b); i, j = 10, b + 12
w[i] = w[j] = 50.0                        # strongly enriched pleiotropic SNPs
rho = 0.5
Rb = np.eye(2 * b); Rb[i, j] = Rb[j, i] = rho
Sig = np.sqrt(w)[:, None] * Rb * np.sqrt(w)[None, :]
c_true = rho * np.sqrt(w[i] * w[j])

reps = 2000
reg, mom, diag_w = [], [], []
LR = np.linalg.cholesky(R)
LS = np.linalg.cholesky(Sig)
for _ in range(reps):
    beta = (LS @ rng.standard_normal((2 * b, q))) * np.sqrt(s)       # per-trait scale s_a (n_a absorbed)
    Z = R @ beta + LR @ rng.standard_normal((2 * b, q))             # z_a = sqrt(n_a) R beta_a + u_a
    Z1, Z2 = Z[:b], Z[b:]
    M = (Z1 * s) @ Z2.T
    l1, l2 = (R1 ** 2).sum(0), (R2 ** 2).sum(0)
    reg.append((R1 @ M @ R2)[i, j - b] / ((s @ s) * l1[i] * l2[j - b]))
    mom.append(np.linalg.solve(R1, np.linalg.solve(R2, M.T).T)[i, j - b] / (s @ s))
    Sh = np.linalg.solve(R, np.linalg.solve(R, (Z * s) @ Z.T - s.sum() * R).T) / (s @ s)
    diag_w.append(np.diag(Sh))
reg, mom = np.array(reg), np.array(mom)
l1, l2 = (R1 ** 2).sum(0), (R2 ** 2).sum(0)
R1c, R2c = R1 @ R1 @ R1, R2 @ R2 @ R2
se_reg = np.sqrt(R1c[i, i] * R2c[j - b, j - b]) / (np.sqrt(s @ s) * l1[i] * l2[j - b])
se_mom = np.sqrt(np.linalg.inv(R1)[i, i] * np.linalg.inv(R2)[j - b, j - b]) / np.sqrt(s @ s)
print(f"true c_ij = {c_true:.2f}")
print(f"REG: mean {reg.mean():.2f}  sd {reg.std():.2f}  predicted null SE {se_reg:.2f}")
print(f"MOM: mean {mom.mean():.2f}  sd {mom.std():.2f}  predicted null SE {se_mom:.2f}")
dw = np.mean(diag_w, 0)
print("within-block Sigma~ estimator: mean w_hat at i, j =", dw[i].round(1), dw[j].round(1),
      " mean of the rest =", np.delete(dw, [i, j]).mean().round(2))

# general (non-null) sandwich SE for REG and MOM: Var(sum_a s_a x_a y_a), x=(A1 z1a)_i, y=(A2 z2a)_j
def sandwich(A1, A2):
    S11, S22, S12 = Sig[:b, :b], Sig[b:, b:], Sig[:b, b:]
    tot = 0.0
    for sa in s:
        V1 = R1 + sa * R1 @ S11 @ R1
        V2 = R2 + sa * R2 @ S22 @ R2
        V12 = sa * R1 @ S12 @ R2
        vx = (A1 @ V1 @ A1.T)[i, i]; vy = (A2 @ V2 @ A2.T)[j - b, j - b]; cxy = (A1 @ V12 @ A2.T)[i, j - b]
        tot += sa ** 2 * (vx * vy + cxy ** 2)
    return np.sqrt(tot)
print("non-null sandwich SE: REG %.2f   MOM %.2f" % (
    sandwich(R1, R2) / ((s @ s) * l1[i] * l2[j - b]),
    sandwich(np.linalg.inv(R1), np.linalg.inv(R2)) / (s @ s)))

# null check: rho = 0, w = 1
reg0, mom0 = [], []
for _ in range(reps):
    Z = R @ (rng.standard_normal((2 * b, q)) * np.sqrt(s)) + LR @ rng.standard_normal((2 * b, q))
    M = (Z[:b] * s) @ Z[b:].T
    reg0.append((R1 @ M @ R2)[i, j - b] / ((s @ s) * l1[i] * l2[j - b]))
    mom0.append(np.linalg.solve(R1, np.linalg.solve(R2, M.T).T)[i, j - b] / (s @ s))
print("null (w=1, rho=0): REG mean %.3f sd %.3f | MOM mean %.3f sd %.3f" % (
    np.mean(reg0), np.std(reg0), np.mean(mom0), np.std(mom0)))
Sig = np.eye(2 * b)
print("null sandwich SE: REG %.3f   MOM %.3f  (the small-s formulas above ignore polygenic background)" % (
    sandwich(R1, R2) / ((s @ s) * l1[i] * l2[j - b]),
    sandwich(np.linalg.inv(R1), np.linalg.inv(R2)) / (s @ s)))
