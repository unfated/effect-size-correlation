"""Tight LD: ridge vs tag-set LRCQ, and tag-projected LRCP (supplement §S3.11, §S4.7).

LD with clusters of near-perfectly linked SNPs (r^2 ~ 0.9-1). Checks
1. ridge   w_l = (D + l I)^{-1} y_bar          E = (D + l I)^{-1} D w
2. tag-set w_T = (D_.T' D_.T)^{-1} D_.T' y_bar  E = w_T + Pi w_notT,  Pi = (D_.T'D_.T)^{-1} D_.T' D_.notT
3. LRCP across two windows, tags P: C_hat = R_A[P,P]^-1 Ybar[P,P] R_B[P,P]^-1 and tag-level variances
   V_hat_T = R_TT^-1 (Z_T S Z_T' - sum c s R_TT) R_TT^-1 / ||s||^2  (diag); rho_T = C / sqrt(V V)
   targets the correlation of R-projected tag effects; using D-based w_hat_T in the denominator does not.
"""
import numpy as np

rng = np.random.default_rng(41)


def clustered_ld(n_clusters, size, r_in=0.97, phi=0.5, n_ref=4000):
    """Clusters of tightly linked SNPs; clusters AR(1)-correlated."""
    k = n_clusters
    Lc = np.linalg.cholesky(phi ** np.abs(np.subtract.outer(np.arange(k), np.arange(k))))
    F = rng.standard_normal((n_ref, k)) @ Lc.T
    X = np.repeat(F, size, axis=1) * r_in + np.sqrt(1 - r_in ** 2) * rng.standard_normal((n_ref, k * size))
    return np.corrcoef(X, rowvar=False)


def greedy_tags(R, thr=0.5):
    tags, covered = [], np.zeros(len(R), bool)
    for i in np.argsort(-(R ** 2).sum(0)):
        if not covered[i]:
            tags.append(i); covered |= R[i] ** 2 >= thr
    return np.array(sorted(tags))


R = clustered_ld(8, 4); m = len(R); D = R ** 2
print(f"cond(D) = {np.linalg.cond(D):.1e}; median diag(D^-1) = {np.median(np.diag(np.linalg.inv(D))):.0f}")
w = np.zeros(m); w[[1, 9, 22]] = [12, 6, 9]                  # causal SNPs inside clusters 0, 2, 5
w = w + 0.2
T = greedy_tags(R); nT = np.setdiff1d(np.arange(m), T)
q = 200; s = rng.uniform(0.03, 0.12, q)
L = np.linalg.cholesky(R)


def sim():
    beta = np.sqrt(w)[:, None] * rng.standard_normal((m, q)) * np.sqrt(s)
    Z = R @ beta + L @ rng.standard_normal((m, q))
    return Z, (Z ** 2 - 1) @ s / (s @ s)


lam = 0.05
Rl = np.linalg.solve(D + lam * np.eye(m), D)
DT = D[:, T]; LT_ = np.linalg.solve(DT.T @ DT, DT.T)
Pi = LT_ @ D[:, nT]
ests_r, ests_t = [], []
for _ in range(1500):
    Z, yb = sim()
    ests_r.append(np.linalg.solve(D + lam * np.eye(m), yb)); ests_t.append(LT_ @ yb)
er, et = np.mean(ests_r, 0), np.mean(ests_t, 0)
print("1. ridge: max|mean - (D+lI)^-1 D w| =", np.abs(er - Rl @ w).max().round(3),
      "| cluster sums true vs ridge-mean:", [(w[c*4:(c+1)*4].sum().round(1), er[c*4:(c+1)*4].sum().round(1)) for c in (0, 2, 5)])
tgt = w[T] + Pi @ w[nT]
print("2. tag-set: max|mean - (w_T + Pi w_notT)| =", np.abs(et - tgt).max().round(3),
      "| tags:", T.tolist(), "| tag estimand", tgt.round(2).tolist(),
      "| column sums of Pi (1 = full allocation):", Pi.sum(0).round(2).tolist()[:6], "...")
print("   SD per SNP, ridge median %.2f vs tag-set median %.2f" % (np.median(np.std(ests_r, 0)), np.median(np.std(ests_t, 0))))

# 3. LRCP with tight LD in both windows
RA = clustered_ld(6, 4, r_in=0.9); RB = clustered_ld(6, 4, r_in=0.9); b = len(RA)
TA, TB = greedy_tags(RA), greedy_tags(RB)
wA = np.full(b, 0.2); wB = np.full(b, 0.2); i, j = 5, 13; wA[i] = 15; wB[j] = 15
Sig = np.block([[np.diag(wA), np.zeros((b, b))], [np.zeros((b, b)), np.diag(wB)]]); Sig[i, b + j] = Sig[b + j, i] = 0.6 * 15
Ls = np.linalg.cholesky(Sig); LA, LB = np.linalg.cholesky(RA), np.linalg.cholesky(RB)
PA_ = np.linalg.solve(RA[np.ix_(TA, TA)], RA[TA, :]); PB_ = np.linalg.solve(RB[np.ix_(TB, TB)], RB[TB, :])
Vt = np.block([[PA_, np.zeros((len(TA), b))], [np.zeros((len(TB), b)), PB_]]) @ Sig @ \
     np.block([[PA_, np.zeros((len(TA), b))], [np.zeros((len(TB), b)), PB_]]).T
rho_t = Vt[:len(TA), len(TA):] / np.sqrt(np.outer(np.diag(Vt)[:len(TA)], np.diag(Vt)[len(TA):]))
ta, tb = np.argmax(np.abs(PA_[:, i])), np.argmax(np.abs(PB_[:, j]))
Cs, Vs, Ws = [], [], []
DA, DB = RA ** 2, RB ** 2
for _ in range(1500):
    beta = (Ls @ rng.standard_normal((2 * b, q))) * np.sqrt(s)
    ZA = RA @ beta[:b] + LA @ rng.standard_normal((b, q)); ZB = RB @ beta[b:] + LB @ rng.standard_normal((b, q))
    Yb = (ZA * s) @ ZB.T / (s @ s)
    C = np.linalg.solve(RA[np.ix_(TA, TA)], Yb[np.ix_(TA, TB)]) @ np.linalg.inv(RB[np.ix_(TB, TB)])
    vA = np.diag(np.linalg.solve(RA[np.ix_(TA, TA)], np.linalg.solve(RA[np.ix_(TA, TA)], (ZA[TA] * s) @ ZA[TA].T - s.sum() * RA[np.ix_(TA, TA)]).T)) / (s @ s)
    vB = np.diag(np.linalg.solve(RB[np.ix_(TB, TB)], np.linalg.solve(RB[np.ix_(TB, TB)], (ZB[TB] * s) @ ZB[TB].T - s.sum() * RB[np.ix_(TB, TB)]).T)) / (s @ s)
    wTA = np.linalg.solve(DA[:, TA].T @ DA[:, TA], DA[:, TA].T @ ((ZA ** 2 - 1) @ s / (s @ s)))
    wTB = np.linalg.solve(DB[:, TB].T @ DB[:, TB], DB[:, TB].T @ ((ZB ** 2 - 1) @ s / (s @ s)))
    Cs.append(C[ta, tb]); Vs.append((vA[ta], vB[tb])); Ws.append((wTA[ta], wTB[tb]))
Cm, Vm, Wm = np.mean(Cs), np.mean(Vs, 0), np.mean(Ws, 0)
print(f"3. tag pair: true projected rho {rho_t[ta, tb]:.3f}; C/sqrt(V V) with R-projected tag variances "
      f"{Cm / np.sqrt(Vm.prod()):.3f}; with D-based tag w {Cm / np.sqrt(Wm.prod()):.3f}; SNP-level rho 0.600")
