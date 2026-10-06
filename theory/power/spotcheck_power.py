"""Simulation spot-check of the §S9 analytic SEs (independent SNPs, R = I, D = I so kappa_LD = 1).
Traits: q in K clusters with within-cluster phenotypic (overlap) correlation r_w; equal s.
Checks SD(w_hat_i) and SD(C_hat_ij) against the formulas with q_eff = q^2 / sum Gamma^2."""
import numpy as np
rng = np.random.default_rng(3)
q, K, rw, s, m, reps = 200, 10, 0.4, 0.05, 40, 4000
G = np.kron(np.eye(K), np.full((q // K, q // K), rw)); np.fill_diagonal(G, 1)
qeff = q * q / (G ** 2).sum()
Lg = np.linalg.cholesky(G)
w = np.ones(m); w[:10] = 20.0; rho = 0.5                       # SNPs 0,1 enriched and correlated
Sig = np.diag(w); Sig[0, 1] = Sig[1, 0] = rho * 20.0
Ls = np.linalg.cholesky(Sig)
wh, ch = [], []
for _ in range(reps):
    Z = rng.standard_normal((m, q)) @ Lg.T + np.sqrt(s) * Ls @ rng.standard_normal((m, q))   # genetic part indep across traits
    y = (Z ** 2 - 1).sum(1) / (q * s)
    wh.append(y); ch.append((Z[0] * Z[1]).sum() / (q * s))
wh, ch = np.array(wh), np.array(ch)
qP = q * q / (G ** 2).sum(); qPG = q * q / (G * np.eye(q)).sum(); qG = q                 # r_g = I here
pred_w = lambda wi: np.sqrt(2 * (1 / qP + 2 * s * wi / qPG + (s * wi) ** 2 / qG) / s ** 2)
naive = lambda wi: np.sqrt(2 * (1 + s * wi) ** 2 / (qP * s * s))
print(f"q_P = {qP:.1f}, q_G = {qG}")
print(f"SD w_hat (w=1):  sim {wh[:, 20:].std(0).mean():.2f}  formula {pred_w(1):.2f}  (single q_eff {naive(1):.2f})")
print(f"SD w_hat (w=20): sim {wh[:, 2:10].std(0).mean():.2f}  formula {pred_w(20):.2f}  (single q_eff {naive(20):.2f})")
C = rho * 20
pred_c = np.sqrt((1 / qP + 2 * s * 20 / qPG + (s * 20) ** 2 / qG + (s * C) ** 2 / qG) / s ** 2)
print(f"C_hat mean {ch.mean():.2f} (true {C:.1f}); SD sim {ch.std():.2f}  formula {pred_c:.2f}")
