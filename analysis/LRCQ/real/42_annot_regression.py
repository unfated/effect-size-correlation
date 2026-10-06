#!/usr/bin/env python3
"""LRCQ annotation regression on the trait-pooled response (methods, S7.3).

ybar = (D A) tau + e, fitted by least squares over all core SNPs, optionally with
LDSC-style weights 1/l2. Fold enrichment of annotation c is
  E_c = mean_{k in c} (A tau)_k / mean_k (A tau)_k
(for continuous annotations, the A-weighted mean). SEs: 200-block jackknife over
contiguous SNP blocks (valid for annotation-level quantities).
Two fits: joint (all baseline-LF annotations) and marginal (base + one annotation).

Usage: 42_annot_regression.py <annot_ld.npz> <out_tsv> [weights: none|l2]
"""
import sys

import numpy as np
import pandas as pd
from scipy import stats

npz, out = sys.argv[1:3]
wmode = sys.argv[3] if len(sys.argv) > 3 else "none"
z = np.load(npz, allow_pickle=True)
y = z["ybar"]; X = z["DA"].astype(np.float64); A = z["A"].astype(np.float64); l2 = z["l2"]
an = list(z["annot_names"])
ok = np.isfinite(y) & (y < np.quantile(y[np.isfinite(y)], 0.999))   # drop extreme outliers (cf. chi2 > 80 in LDSC)
y, X, A, l2 = y[ok], X[:, :], A, l2
X, A, l2 = X[ok], A[ok], l2[ok]
wt = 1 / np.maximum(l2, 1) if wmode == "l2" else np.ones_like(y)
keep = np.where(A.std(0) > 0)[0]
keep = np.r_[0, keep[keep != 0]]
nb = 200
blk = np.minimum((np.arange(len(y)) * nb) // len(y), nb - 1)


def fit(cols):
    Xc = X[:, cols] * np.sqrt(wt)[:, None]; yc = y * np.sqrt(wt)
    XtX_b = np.stack([Xc[blk == b].T @ Xc[blk == b] for b in range(nb)])
    Xty_b = np.stack([Xc[blk == b].T @ yc[blk == b] for b in range(nb)])
    XtX, Xty = XtX_b.sum(0), Xty_b.sum(0)
    Ac = A[:, cols]
    Acnt = np.stack([Ac[blk == b].sum(0) for b in range(nb)])
    AtA_b = np.stack([Ac[blk == b].T @ Ac[blk == b] for b in range(nb)])   # sum_k A_kc A_kd

    def enr(tau, cnt, AtA):
        tot = cnt[0]                                      # number of SNPs
        h_all = (cnt @ tau) / tot                         # mean (A tau)_k over SNPs
        h_c = (AtA @ tau) / np.maximum(cnt, 1e-12)        # A-weighted mean within c
        return h_c / h_all
    tau = np.linalg.lstsq(XtX, Xty, rcond=None)[0]
    E = enr(tau, Acnt.sum(0), AtA_b.sum(0))
    jt, jE = [], []
    for b in range(nb):
        tb = np.linalg.lstsq(XtX - XtX_b[b], Xty - Xty_b[b], rcond=None)[0]
        jt.append(tb); jE.append(enr(tb, Acnt.sum(0) - Acnt[b], AtA_b.sum(0) - AtA_b[b]))
    jt, jE = np.array(jt), np.array(jE)
    f = (nb - 1) / nb
    return tau, np.sqrt(f * ((jt - jt.mean(0)) ** 2).sum(0)), E, np.sqrt(f * ((jE - jE.mean(0)) ** 2).sum(0))


rows = []
tau, tse, E, Ese = fit(keep)
for i, c in enumerate(keep):
    rows.append(dict(annotation=an[c], model="joint", tau=tau[i], tau_se=tse[i], E=E[i], E_se=Ese[i]))
for c in keep[1:]:
    tau, tse, E, Ese = fit(np.array([0, c]))
    rows.append(dict(annotation=an[c], model="marginal", tau=tau[1], tau_se=tse[1], E=E[1], E_se=Ese[1]))
r = pd.DataFrame(rows)
r["prop"] = [A[:, an.index(a)].mean() for a in r["annotation"]]
r["z_tau"] = r["tau"] / r["tau_se"]; r["p_tau"] = 2 * stats.norm.sf(np.abs(r["z_tau"]))
r["z_E"] = (r["E"] - 1) / r["E_se"]; r["p_E"] = 2 * stats.norm.sf(np.abs(r["z_E"]))
r["weights"] = wmode
r.to_csv(out, sep="\t", index=False)
print(r[r["model"] == "marginal"].sort_values("p_E").head(20).to_string())
