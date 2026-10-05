"""Plug-in variance for phenotype-category contrasts with sparse, heavy-tailed w (supplement §S7.4a).

Null: same w for both trait groups. Contrast Delta = 1_c' D^-1 Y2 u_con / |c| (and per-SNP versions).
Variance (7.1) needs G = R W R. Plug-ins compared:
  true      W known
  naive     W_hat = Diag(w_hat_OLS) from the same data
  indep     w_hat from an independent replicate data set (oracle for "independent G")
  debiased  naive, with tau_GG corrected for E[w_s w_t] = w_s w_t + Cov(w_s, w_t)
  clipped   naive with w_hat floored at 0
  boot      parametric bootstrap of the studentised (debiased) contrast under H0, simulating from the
            pooled, PSD-clipped fit (run with --boot; ~10 min)
Results (seed 5): see supplement §S7.4a.
"""
import sys
import numpy as np
from common import ar1_ld, block_diag

rng = np.random.default_rng(5)
nb, bs = 8, 30
R = block_diag(*[ar1_ld(bs, rng.uniform(0.5, 0.85)) for _ in range(nb)]); m = nb * bs
D = R ** 2; Dinv = np.linalg.inv(D)
w = np.full(m, 0.05); big = rng.choice(m, 6, replace=False); w[big] = rng.uniform(40, 80, 6)
w *= m / w.sum()
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
g1 = cl < 2
u_all = s / (s @ s)
u_con = np.where(g1, s / (s[g1] @ s[g1]), -s / (s[~g1] @ s[~g1]))
GG, GT, TT = u_con @ (Gam * Gam) @ u_con, u_con @ (Gam * T) @ u_con, u_con @ (T * T) @ u_con
gg, gt, tt = u_all @ (Gam * Gam) @ u_all, u_all @ (Gam * T) @ u_all, u_all @ (T * T) @ u_all

annot = np.zeros(m, bool); annot[rng.choice(m, 48, replace=False)] = True; annot[big[:2]] = True
targets = {"category": Dinv @ annot / annot.sum(), "top SNP": Dinv[:, big[0]], "null SNP": Dinv[:, np.argmin(w)]}


def var_contrast(a, Gx, debias_V=None):
    """2[tau_RR GG + 2 tau_RG GT + tau_GG TT] for A = Diag(a); optional debiasing of tau_GG."""
    RA = R * a; tRR = np.sum((R * a) * (R * a).T) if False else np.trace(RA @ RA)
    GA = Gx * a
    tRG, tGG = np.trace(RA @ GA), np.trace(GA @ GA)
    if debias_V is not None:   # E tau_GG(G_hat) = tau_GG(G) + sum_st V_st ((R A R)_st)^2
        RAR = R @ (a[:, None] * R)
        tGG = max(tGG - np.sum(debias_V * RAR ** 2), 0.0)
    return 2 * (tRR * GG + 2 * tRG * GT + tGG * TT)


def cov_what(Gx):
    """Cov(w_hat_OLS) = D^-1 Cov(y_bar) D^-1, Cov(y_bar) = 2[gg R∘R + 2 gt R∘G + tt G∘G]."""
    return Dinv @ (2 * (gg * D + 2 * gt * (R * Gx) + tt * (Gx * Gx))) @ Dinv


def draw():
    Z = LR @ rng.standard_normal((m, q)) @ LGa.T + LG @ rng.standard_normal((m, q)) @ LT.T
    return Z ** 2 - 1


reps = 1500
res = {k: {p: [] for p in ["true", "naive", "indep", "debiased", "clipped"]} for k in targets}
for _ in range(reps):
    Y2 = draw(); Y2b = draw()
    wh = Dinv @ (Y2 @ u_all); whb = Dinv @ (Y2b @ u_all)
    plug = {"true": G, "naive": R @ np.diag(wh) @ R, "indep": R @ np.diag(whb) @ R,
            "clipped": R @ np.diag(np.clip(wh, 0, None)) @ R}
    V = cov_what(plug["naive"])
    for k, a in targets.items():
        delta = a @ (Y2 @ u_con)
        for p, Gx in plug.items():
            res[k][p].append(delta / np.sqrt(max(var_contrast(a, Gx), 1e-12)))
        res[k]["debiased"].append(delta / np.sqrt(max(var_contrast(a, plug["naive"], V), 1e-12)))
for k in targets:
    print(k)
    for p, z in res[k].items():
        z = np.array(z)
        print(f"   {p:9s} SD(z) {z.std():.2f}   type-I 5% {np.mean(np.abs(z) > 1.96):.3f}")


if "--boot" in sys.argv:
    def stat(Y2, a):
        wh = Dinv @ (Y2 @ u_all); Gh = R @ np.diag(wh) @ R
        return (a @ (Y2 @ u_con)) / np.sqrt(max(var_contrast(a, Gh, cov_what(Gh)), 1e-12))

    reps, B = 200, 100
    for k, a in targets.items():
        pv = []
        for _ in range(reps):
            Y2 = draw(); z = stat(Y2, a)
            ev, Q = np.linalg.eigh(R @ np.diag(np.clip(Dinv @ (Y2 @ u_all), 0, None)) @ R)
            L = Q * np.sqrt(np.clip(ev, 0, None))
            zb = [stat((LR @ rng.standard_normal((m, q)) @ LGa.T + L @ rng.standard_normal((m, q)) @ LT.T) ** 2 - 1, a)
                  for _ in range(B)]
            pv.append((1 + np.sum(np.abs(zb) >= abs(z))) / (B + 1))
        pv = np.array(pv)
        print(f"boot {k:9s} type-I 5% {np.mean(pv <= .05):.3f}  10% {np.mean(pv <= .10):.3f}")
