"""Gene-burden conditional test when a gene's causal SNPs are correlated (supplement §S8.2a).

Windows A, B: 9 tags each (AR(1) LD); one gene per window = 3 causal SNPs (w = 30) with within-gene
effect correlation 0.6; cross-gene rho on all 9 SNP pairs (0 under H0, 0.5 under H1).
Burden statistic C_hat = sum_a s_a (g'z_Aa)(h'z_Ba)/||s||^2, g = R_A^-1 v_G, h = R_B^-1 v_H.
Conditional variance (8.2) needs g'G_A g, with G_A = R_A Sigma_A R_A. Choices:
  diag     R_A Diag(w_hat) R_A, w_hat = clipped LRCQ OLS (assumes no local rho; the package default)
  oracle   true G_A (includes within-gene correlation)
  projdiag R_A Diag(diag(R_A^-1 G_hat R_A^-1))_+ R_A: diagonal of the local Sigma_hat (drops within-gene C)
  moment   G_hat = sum_a s_a (z_Aa z_Aa' - c_aa R_A)/||s||^2, PSD-clipped (direct second moment, §S3.9;
           no R^-1, no local-rho assumption). Identical to R Sigma_hat R with Sigma_hat the §S3.9 local estimator.
  boot     moment-studentised statistic, p-value from a conditional parametric bootstrap: Z_A* drawn from
           Gamma (x) R_A + T (x) G_hat_A with Z_B held fixed (B = 200)
"""
import numpy as np
from common import ar1_ld

rng = np.random.default_rng(61)
b = 9
RA, RB = ar1_ld(b, 0.5), ar1_ld(b, 0.5)
gene = np.arange(3, 6)
v = np.zeros(b); v[gene] = 1
g, h = np.linalg.solve(RA, v), np.linalg.solve(RB, v)
W_GENE = 30.0
w = np.ones(b) * 0.3; w[gene] = W_GENE
U0 = np.eye(b); U0[np.ix_(gene, gene)] = 0.6; np.fill_diagonal(U0, 1)
SigA = np.sqrt(np.outer(w, w)) * U0
D_A = RA ** 2; Dinv = np.linalg.inv(D_A)


def setup(q):
    cl = np.repeat(np.arange(4), q // 4)
    rg = np.where(cl[:, None] == cl[None, :], 0.6, 0.1); np.fill_diagonal(rg, 1)
    s = rng.uniform(0.05, 0.2, q)
    T = np.sqrt(np.outer(s, s)) * rg
    Gam = 0.5 * rg + 0.1 * (1 - np.eye(q)); np.fill_diagonal(Gam, 1.0)
    return s, T, Gam


def run(q, rho, reps=2000, do_boot=False):
    s, T, Gam = setup(q); ss = s @ s
    LGa, LT = np.linalg.cholesky(Gam), np.linalg.cholesky(T)
    CAB = np.zeros((b, b)); CAB[np.ix_(gene, gene)] = rho * W_GENE
    Sig = np.block([[SigA, CAB], [CAB.T, SigA]])
    R = np.block([[RA, np.zeros((b, b))], [np.zeros((b, b)), RB]])
    LR, LG = np.linalg.cholesky(R), np.linalg.cholesky(R @ Sig @ R + 1e-9 * np.eye(2 * b))
    GA_true = RA @ SigA @ RA
    rej = {k: 0 for k in ["diag", "projdiag", "oracle", "moment", "boot"]}
    for _ in range(reps):
        Z = LR @ rng.standard_normal((2 * b, q)) @ LGa.T + LG @ rng.standard_normal((2 * b, q)) @ LT.T
        ZA, ZB = Z[:b], Z[b:]
        x, y = g @ ZA, h @ ZB
        C = (s * x) @ y / ss
        ys = s * y
        kR, kT = ys @ Gam @ ys, ys @ T @ ys
        wh = np.clip(Dinv @ ((ZA ** 2 - 1) @ s) / ss, 0, None)
        Gm = (ZA * s) @ ZA.T / ss - RA
        ev, Q = np.linalg.eigh(Gm); Gm = (Q * np.clip(ev, 0, None)) @ Q.T
        Sl = np.linalg.solve(RA, np.linalg.solve(RA, Gm).T)
        GAs = {"diag": RA @ np.diag(wh) @ RA, "projdiag": RA @ np.diag(np.clip(np.diag(Sl), 0, None)) @ RA,
               "oracle": GA_true, "moment": Gm}
        for k, GA in GAs.items():
            var = (g @ RA @ g * kR + g @ GA @ g * kT) / ss ** 2
            rej[k] += abs(C) / np.sqrt(var) > 1.96
        if do_boot:
            zobs = abs(C) / np.sqrt((g @ RA @ g * kR + g @ Gm @ g * kT) / ss ** 2)
            LGm = np.linalg.cholesky(Gm + 1e-9 * np.eye(b)); LRA = np.linalg.cholesky(RA)
            zb = []
            for _ in range(200):
                ZA_ = LRA @ rng.standard_normal((b, q)) @ LGa.T + LGm @ rng.standard_normal((b, q)) @ LT.T
                Gb = (ZA_ * s) @ ZA_.T / ss - RA
                ev, Q = np.linalg.eigh(Gb); Gb = (Q * np.clip(ev, 0, None)) @ Q.T
                Cb = (s * (g @ ZA_)) @ y / ss
                zb.append(abs(Cb) / np.sqrt((g @ RA @ g * kR + g @ Gb @ g * kT) / ss ** 2))
            rej["boot"] += (1 + np.sum(np.array(zb) >= zobs)) / 201 <= 0.05
    return {k: r / reps for k, r in rej.items() if do_boot or k != 'boot'}


for q in [100, 300]:
    for rho in [0.0, 0.5]:
        r = run(q, rho)
        r.update({"boot": run(q, rho, reps=300, do_boot=True)["boot"]})
        print(f"q={q} rho={rho}: " + ", ".join(f"{k} {v:.3f}" for k, v in r.items()))
