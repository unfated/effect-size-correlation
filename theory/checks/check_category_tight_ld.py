"""Category enrichment E_c when the category splits tight LD clumps (supplement §S7.2a).

16 clusters of 4 SNPs at r ~ 0.99. Category c covers a fraction f_k in {0,1/4,...,1} of
cluster k, so it splits most clumps. True w follows w = tau0 + tau1 * 1_c (3-fold in c).
Estimators of the in-category mean enrichment E_c = 1_c'w/|c|:
  exact   1_c' D^+ y_bar / |c|                  unbiased; variance ~ 1_c'D^-1 1_c explodes
  ridge   1_c'(D + l I)^-1 y_bar / |c|          estimand smooths w to clump means
  clump   sum_k f_k W_k / sum_k f_k g_k         W_k = tag-set clump totals (3.6); estimand E_c^clump
  annreg  tau0 + tau1, tau = (A'D^2A)^-1 A'D y  A = [1, 1_c]; exact if w = A tau
Run twice: (a) w in the annotation span, (b) w concentrated on one SNP per clump (model wrong).
"""
import numpy as np

rng = np.random.default_rng(7)
K, g = 16, 4


def clustered_ld(k, size, r_in=0.99, phi=0.5, n_ref=4000):
    Lc = np.linalg.cholesky(phi ** np.abs(np.subtract.outer(np.arange(k), np.arange(k))))
    F = rng.standard_normal((n_ref, k)) @ Lc.T
    X = np.repeat(F, size, axis=1) * r_in + np.sqrt(1 - r_in ** 2) * rng.standard_normal((n_ref, k * size))
    return np.corrcoef(X, rowvar=False)


R = clustered_ld(K, g); m = len(R)
ev, V = np.linalg.eigh(R ** 2); D = (V * np.clip(ev, 0, None)) @ V.T
clump = np.arange(m) // g
n_in = np.tile([0, 1, 2, 3, 4, 2, 1, 3], 2)                         # SNPs of each clump inside c
inc = np.concatenate([np.arange(g) < n_in[k] for k in range(K)])
f = n_in / g
T = np.array([k * g for k in range(K)])                               # one tag per clump (r^2 < 0.5 between clumps)
LT = np.linalg.solve(D[:, T].T @ D[:, T], D[:, T].T)
A = np.column_stack([np.ones(m), inc.astype(float)])
P_ann = np.linalg.solve(A.T @ D @ D @ A, A.T @ D)
lam = 0.1
M = {"exact": np.linalg.pinv(D, rcond=1e-12), "ridge": np.linalg.inv(D + lam * np.eye(m))}
q = 200; s = rng.uniform(0.03, 0.12, q); ss = s @ s
L = np.linalg.cholesky(R + 1e-10 * np.eye(m))
print(f"cond(D) = {np.linalg.cond(D):.1e}; 1_c'D^-1 1_c/|c| = {inc @ M['exact'] @ inc / inc.sum():.1f}")

scen = {"(a) w = tau0 + tau1*1_c": np.where(inc, 1.8, 0.6)}
wb = np.full(m, 0.2); wb[np.flatnonzero(inc)[::2]] += 2.0; scen["(b) w on single SNPs (model wrong)"] = wb
for name, w in scen.items():
    w = w * m / w.sum()
    E = w[inc].mean(); W = np.bincount(clump, w); E_clump = f @ W / (f @ np.full(K, g))
    wr = M["ridge"] @ D @ w; tstar = P_ann @ D @ w
    print(f"\n{name}: E_c = {E:.3f}, E_c^clump = {E_clump:.3f}, ridge estimand = {wr[inc].mean():.3f}, "
          f"annreg estimand = {tstar.sum():.3f}")
    res = {k: [] for k in ["exact", "ridge", "clump", "annreg"]}
    sw = np.sqrt(w)[:, None]
    for _ in range(3000):
        Z = L @ rng.standard_normal((m, q)) + R @ (sw * rng.standard_normal((m, q))) * np.sqrt(s)
        y = ((Z ** 2 - 1) @ s) / ss
        for k in M:
            res[k].append((M[k] @ y)[inc].mean())
        res["clump"].append(f @ (LT @ y) / (f @ np.full(K, g)))
        res["annreg"].append((P_ann @ y).sum())
    for k, v in res.items():
        print(f"  {k:7s} mean {np.mean(v):6.3f}  SD {np.std(v):6.3f}")
