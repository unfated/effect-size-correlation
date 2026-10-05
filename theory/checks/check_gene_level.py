"""Gene-level LRCP and multiple testing (supplement §S8).

Two distal windows A, B; tags grouped into 3 genes per window. Gene-level burden covariance
  C_GH = v_G' C v_H   (v = indicator of the gene's tags, alleles as coded)
estimated by  C_hat_GH = sum_a s_a (g_G' z_Aa)(h_H' z_Ba) / ||s||^2,  g_G = R_A^-1 v_G, h_H = R_B^-1 v_H.
Conditional test: given Z_B, C_hat_GH is Gaussian (Gaussian Z_A) with
  Cov(C_hat_GH, C_hat_G'H' | Z_B) = sum_ab s_a s_b y_aH y_bH' [Gam_ab g_G'R_A g_G' + T_ab g_G'G_A g_G'] / ||s||^4,
valid for ANY cross-trait dependence. Checks per-pair type-I error and the FWER of a max-T test whose
critical value is simulated from that conditional Gaussian, with correlated, overlapping traits.
Also checks unbiasedness under the alternative.
"""
import numpy as np
from common import ar1_ld

rng = np.random.default_rng(53)
b = 9
RA, RB = ar1_ld(b, 0.5), ar1_ld(b, 0.5)
genes = [np.arange(0, 3), np.arange(3, 6), np.arange(6, 9)]
V = np.zeros((b, 3))
for k, g in enumerate(genes):
    V[g, k] = 1.0
GA_ = np.linalg.solve(RA, V); HB_ = np.linalg.solve(RB, V)
q = 120
cl = np.repeat(np.arange(4), q // 4)
rg = np.where(cl[:, None] == cl[None, :], 0.6, 0.1); np.fill_diagonal(rg, 1)
s = rng.uniform(0.05, 0.2, q)
T = np.sqrt(np.outer(s, s)) * rg
Gam = 0.5 * rg + 0.1 * (1 - np.eye(q)); np.fill_diagonal(Gam, 1.0)
LGa, LT = np.linalg.cholesky(Gam), np.linalg.cholesky(T)
wA = np.ones(b); wA[1] = 8; wA[4] = 8; wB = np.ones(b); wB[7] = 8; wB[2] = 8


def draw(CAB):
    Sig = np.block([[np.diag(wA), CAB], [CAB.T, np.diag(wB)]])
    R = np.block([[RA, np.zeros((b, b))], [np.zeros((b, b)), RB]])
    G = R @ Sig @ R
    LR, LG = np.linalg.cholesky(R), np.linalg.cholesky(G + 1e-9 * np.eye(2 * b))
    Z = LR @ rng.standard_normal((2 * b, q)) @ LGa.T + LG @ rng.standard_normal((2 * b, q)) @ LT.T
    return Z[:b], Z[b:]


GA = RA @ np.diag(wA) @ RA
gRg = GA_.T @ RA @ GA_; gGg = GA_.T @ GA @ GA_          # 3x3 (gene x gene) in window A
ss = s @ s


def analyse(ZA, ZB, nsim=4000):
    X = GA_.T @ ZA                                       # 3 x q
    Y = HB_.T @ ZB                                       # 3 x q
    C = (X * s) @ Y.T / ss                               # 3 x 3 gene pairs
    # conditional covariance over the 9 gene pairs (G, H) given Y
    Ys = Y * s
    K_R = Ys @ Gam @ Ys.T; K_T = Ys @ T @ Ys.T          # (H,H') q-forms
    Cov = (np.kron(gRg, K_R) + np.kron(gGg, K_T)) / ss ** 2   # index (G,H) row-major
    sd = np.sqrt(np.diag(Cov))
    Zst = C.ravel() / sd
    sims = rng.multivariate_normal(np.zeros(9), Cov, nsim) / sd
    crit = np.quantile(np.abs(sims).max(1), 0.95)
    return C, Zst, crit


rej_pair, rej_any = [], []
for _ in range(1500):
    ZA, ZB = draw(np.zeros((b, b)))
    C, Zst, crit = analyse(ZA, ZB, 2000)
    rej_pair.append(np.abs(Zst) > 1.96); rej_any.append(np.abs(Zst).max() > crit)
print(f"null: per-pair type-I {np.mean(rej_pair):.3f} (target 0.05); max-T FWER {np.mean(rej_any):.3f} (target 0.05)")

CAB = np.zeros((b, b)); CAB[1, 7] = 0.5 * 8                 # gene A0 x gene B2
true = V.T @ (np.linalg.solve(RA, RA) @ CAB @ np.linalg.solve(RB, RB)) @ V
ests = np.array([analyse(*draw(CAB), 10)[0] for _ in range(6000)])
print("alternative: true C_GH\n", true.round(2), "\nmean C_hat_GH\n", ests.mean(0).round(2),
      "\nMC s.e.\n", (ests.std(0) / np.sqrt(len(ests))).round(2))
