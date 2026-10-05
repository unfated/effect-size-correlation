"""Prototype (numpy) of the distal-pair LRCP estimator, used to check the
closed form before the package's stage 4 lands.

For SNP sets A and B in different LD windows (r_kl = 0 for k in A, l in B),
theory S1.5 gives  E[z_ka z_la] = t_aa (R_A Σ̃_AB R_B)_kl.  With
Ybar = Σ_a z_.a z_.a' / Σ_a t_aa  (|A| × |B|),  E[Ybar] = R_A Σ̃_AB R_B.
Restricting Σ̃ to prescreened SNPs P_A, P_B, the GLS solution under the
Kronecker null covariance (R_B ⊗ R_A) is
    Σ̂ = R_A[P,P]^{-1} Ybar[P_A,P_B] R_B[P,P]^{-1}
and ρ̂_ij = Σ̂_ij / sqrt(w_i w_j).  SEs: jackknife over trait clusters.
"""
import numpy as np, pandas as pd, sys, itertools

def ar1(p, r):
    i = np.arange(p); return r ** np.abs(i[:, None] - i[None, :])

def simulate(rng, q, p=100, r_ld=0.6, rho=0.5, n=3e5, M=1e6, n_clusters=None, rg=0.0,
             overlap_r=0.0, prop=0.1):
    RA, RB = ar1(p, r_ld), ar1(p, r_ld)
    m = 2 * p
    w = np.zeros(m); nz = rng.choice(m, int(prop * m), replace=False)
    w[nz] = rng.lognormal(0, 1, nz.size); w *= m * prop / w.sum() / prop  # mean over nonzero ~ 1/prop
    w = w / w.mean()
    i = rng.choice(np.where(w[:p] > 0)[0]); j = p + rng.choice(np.where(w[p:] > 0)[0])
    Rb = np.eye(m); Rb[i, j] = Rb[j, i] = rho
    Sig = np.sqrt(w)[:, None] * Rb * np.sqrt(w)[None, :]     # Σ̃ (per-unit), h2/M scaling below
    h2 = rng.beta(2, 8, q)
    lab = np.arange(q) % (n_clusters or q)
    G = np.where(lab[:, None] == lab[None, :], rg, 0.0); np.fill_diagonal(G, 1.0)
    V = np.sqrt(h2)[:, None] * G * np.sqrt(h2)[None, :] / M  # trait genetic covariance per unit
    C = np.where(lab[:, None] == lab[None, :], overlap_r, 0.0); np.fill_diagonal(C, 1.0)
    # B ~ MN(0, Σ̃, V)
    LS = np.linalg.cholesky(Sig + 1e-12 * np.eye(m)); LV = np.linalg.cholesky(V)
    B = LS @ rng.standard_normal((m, q)) @ LV.T
    Rfull = np.zeros((m, m)); Rfull[:p, :p] = RA; Rfull[p:, p:] = RB
    LR = np.linalg.cholesky(Rfull + 1e-10 * np.eye(m)); LC = np.linalg.cholesky(C)
    E = LR @ rng.standard_normal((m, q)) @ LC.T
    Z = np.sqrt(n) * Rfull @ B + E
    t = n * h2 / M
    return dict(Z=Z, RA=RA, RB=RB, w=w, i=i, j=j, t=t, lab=lab, B=B, p=p)

def lrcp_distal(Z, RA, RB, t, PA, PB, lab, wA, wB):
    p = RA.shape[0]
    ZA, ZB = Z[:p], Z[p:]
    def est(idx):
        Y = (ZA[np.ix_(PA, idx)] @ ZB[np.ix_(PB, idx)].T) / t[idx].sum()
        S = np.linalg.solve(RA[np.ix_(PA, PA)], Y) @ np.linalg.inv(RB[np.ix_(PB, PB)])
        return S / np.sqrt(np.outer(wA, wB))
    full = est(np.arange(Z.shape[1]))
    cl = np.unique(lab)
    jk = np.array([est(np.where(lab != c)[0]) for c in cl])
    g = len(cl)
    se = np.sqrt((g - 1) / g * ((jk - jk.mean(0)) ** 2).sum(0))
    return full, se

def one(rng, **kw):
    s = simulate(rng, **kw)
    p, i, j = s["p"], s["i"], s["j"]
    PA = np.where(s["w"][:p] > 0)[0]; PB = np.where(s["w"][p:] > 0)[0]  # oracle prescreen (true w>0)
    est, se = lrcp_distal(s["Z"], s["RA"], s["RB"], s["t"], PA, PB, s["lab"] if kw.get("n_clusters") else np.arange(len(s["t"])) // max(1, len(s["t"]) // 30),
                          s["w"][PA], s["w"][p + PB])
    ii, jj = np.where(PA == i)[0][0], np.where(PB == j - p)[0][0]
    zc = lambda a, b: np.corrcoef(a, b)[0, 1]
    q = s["Z"].shape[1]
    # null pairs: all other prescreened cross pairs
    mask = np.ones(est.shape, bool); mask[ii, jj] = False
    naive = np.array([[zc(s["Z"][a], s["Z"][p + b]) for b in PB] for a in PA])
    return dict(rho_hat=est[ii, jj], se=se[ii, jj], naive=naive[ii, jj],
                oracle=zc(s["B"][i], s["B"][j]),
                null_z_lrcp=(est[mask] / se[mask]), null_naive_z=np.arctanh(naive[mask]) * np.sqrt(q - 3))

if __name__ == "__main__":
    rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 1)
    rows = []
    for q, (nc, rg, ov), rho in itertools.product((30, 100, 300), ((None, 0, 0), (10, 0.6, 0.3)), (0.0, 0.5, 0.8)):
        res = [one(rng, q=q, rho=rho, n_clusters=nc, rg=rg, overlap_r=ov, M=2e5) for _ in range(int(sys.argv[3]) if len(sys.argv) > 3 else 100)]
        rh = np.array([r["rho_hat"] for r in res]); se = np.array([r["se"] for r in res])
        nl = np.concatenate([r["null_z_lrcp"] for r in res]); nn = np.concatenate([r["null_naive_z"] for r in res])
        rows.append(dict(q=q, traits="indep" if nc is None else "10clus_rg0.6_ov0.3", rho=rho,
                         mean_rho_hat=rh.mean(), median_rho_hat=np.median(rh), sd_rho_hat=rh.std(), mean_se=se.mean(),
                         coverage=np.mean(np.abs(rh - rho) < 1.96 * se),
                         mean_naive=np.mean([r["naive"] for r in res]), mean_oracle=np.mean([r["oracle"] for r in res]),
                         typeI_lrcp=np.mean(np.abs(nl) > 1.96), typeI_naive=np.mean(np.abs(nn) > 1.96)))
        print(rows[-1], flush=True)
    pd.DataFrame(rows).to_csv(sys.argv[1], sep="\t", index=False, float_format="%.4g")
