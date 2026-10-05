#!/usr/bin/env python3
"""Motivating analysis: are some SNPs consistently more associated across traits
than their LD score predicts?

Per SNP and trait, y_ka = (z_ka^2 - c_aa) / s_a estimates (D w)_k for that trait
(E = l_k when w = 1). Pooled over a set of traits, ybar_k = sum s_a^2 y_ka / sum s_a^2.
  1. ybar vs LD score (binned): the LDSC expectation is the identity line.
  2. Split-half replication: traits are split at random into two halves; per-SNP
     residuals of ybar on l_k are correlated between halves. Under w = 1 for all
     SNPs, the halves share no signal beyond LD score and the correlation is ~0
     (sampling covariance comes only from shared samples, removed via intercepts in
     expectation). Repeated over random splits; also at 1-Mb block resolution.

Usage: 43_motivation.py <zprefix> <ldsc_prefix> <out_prefix> [n_splits]
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

zp, lp, outp = sys.argv[1:4]
nsp = int(sys.argv[4]) if len(sys.argv) > 4 else 20
rng = np.random.default_rng(1)
uni = pd.read_csv("/home/user/data/ref/snp_universe.tsv", sep="\t")
m = len(uni)
M = float(open("/home/user/data/ref/UKBB.EUR.l2.M_5_50").read().split()[0])
tr = pd.read_csv(zp + ".traits.tsv", sep="\t")
q = len(tr)
h2 = pd.read_csv(lp + ".h2.tsv", sep="\t").set_index("trait_id").loc[tr["trait_id"]]
s = tr["N"].to_numpy() * h2["h2_panukb_ldsc"].to_numpy() / M
c = h2["intercept_panukb"].to_numpy()
Z = np.memmap(zp + ".Z.f32", dtype=np.float32, mode="r", shape=(m, q))
l2 = uni["L2_UKB_EUR"].to_numpy()
# trait-pair intercepts and genetic covariances (as in 30_stage3_lrcq.R) for the null expectation
rd = lambda f: pd.read_csv(f, sep="\t", index_col=0).loc[tr["trait_id"], tr["trait_id"]].to_numpy()
Cm = rd(lp + ".intercept.tsv"); np.fill_diagonal(Cm, c)
Gm = rd(lp + ".gcov.tsv")
d_old = np.sqrt(np.maximum(np.diag(Gm), 1e-6)); d_new = np.sqrt(np.maximum(h2["h2_panukb_ldsc"].to_numpy(), 1e-6))
Gm = Gm / np.outer(d_old, d_old) * np.outer(d_new, d_new)
Tm = np.sqrt(np.outer(tr["N"], tr["N"])) * Gm / M          # Cov(z_a, z_b) = c_ab + T_ab l_k under w = 1


def null_cov(u, v):
    """sum_ab u_a v_b 2 (c_ab + T_ab l)^2 = k0 + k1 l + k2 l^2 (Gaussian z, w = 1)."""
    return 2 * (u @ Cm ** 2 @ v), 4 * (u @ (Cm * Tm) @ v), 2 * (u @ Tm ** 2 @ v)
mhc = (uni["CHR"] == 6) & (uni["BP"] > 25e6) & (uni["BP"] < 34e6)

# per-SNP sums over traits, accumulated in chunks: S_a-weighted numerator and denominator
def pooled(cols):
    num = np.zeros(m); den = np.zeros(m)
    for i0 in range(0, m, 100000):
        zc = np.asarray(Z[i0:i0 + 100000][:, cols], dtype=np.float64)
        ok = ~np.isnan(zc)
        num[i0:i0 + 100000] = np.where(ok, (zc ** 2 - c[cols]) * s[cols], 0).sum(1)
        den[i0:i0 + 100000] = (ok * s[cols] ** 2).sum(1)
    return num / den

ok = (~mhc) & np.isfinite(l2)
yall = pooled(np.arange(q))
ok &= np.isfinite(yall)
cap = np.quantile(yall[ok], 0.999)
ok &= yall < cap
blk = (uni["CHR"].astype(np.int64) * 1000 + (uni["BP"] // 1e6).astype(np.int64)).to_numpy()


def resid(y, mk):
    X = np.c_[np.ones(mk.sum()), l2[mk]]
    b = np.linalg.lstsq(X, y[mk], rcond=None)[0]
    return y[mk] - X @ b


res = []
for sidx in range(nsp):
    perm = rng.permutation(q); h1, h2_ = np.sort(perm[: q // 2]), np.sort(perm[q // 2:])
    y1, y2 = pooled(h1), pooled(h2_)
    mk = ok & np.isfinite(y1) & np.isfinite(y2)
    r1, r2 = resid(y1, mk), resid(y2, mk)
    snp_r = np.corrcoef(r1, r2)[0, 1]
    u = np.zeros(q); u[h1] = s[h1] / (s[h1] ** 2).sum()
    v = np.zeros(q); v[h2_] = s[h2_] / (s[h2_] ** 2).sum()
    L = l2[mk]
    ev = lambda k: (k[0] + k[1] * L + k[2] * L ** 2).mean()
    null_r = ev(null_cov(u, v)) / np.sqrt(ev(null_cov(u, u)) * ev(null_cov(v, v)))
    b = pd.DataFrame(dict(b=blk[mk], r1=r1, r2=r2)).groupby("b").mean()
    res.append(dict(split=sidx, snp_r=snp_r, snp_r_null=null_r, block_r=np.corrcoef(b["r1"], b["r2"])[0, 1]))
    print(res[-1], flush=True)
res = pd.DataFrame(res)
res.to_csv(outp + ".splithalf.tsv", sep="\t", index=False)

# figure
fig, ax = plt.subplots(1, 3, figsize=(13, 4))
bins = np.quantile(l2[ok], np.linspace(0, 1, 51))
bi = np.clip(np.digitize(l2[ok], bins) - 1, 0, 49)
dfb = pd.DataFrame(dict(b=bi, l=l2[ok], y=yall[ok])).groupby("b").agg(l=("l", "mean"), y=("y", "mean"), ys=("y", "sem"))
ax[0].errorbar(dfb["l"], dfb["y"], yerr=1.96 * dfb["ys"], fmt="o", ms=3)
lim = [0, dfb["l"].max() * 1.05]
ax[0].plot(lim, lim, "k--", lw=1)
ax[0].set_xlabel("LD score"); ax[0].set_ylabel("pooled (z² − c)/s across traits")
ax[0].set_title(f"{q} traits")
ax[1].hist(yall[ok] / l2[ok], bins=200, range=(-2, 8), color="grey")
ax[1].set_yscale("log"); ax[1].set_xlabel("pooled signal / LD score"); ax[1].set_ylabel("SNPs")
perm = rng.permutation(q); y1, y2 = pooled(np.sort(perm[: q // 2])), pooled(np.sort(perm[q // 2:]))
mk = ok & np.isfinite(y1) & np.isfinite(y2)
b = pd.DataFrame(dict(b=blk[mk], r1=resid(y1, mk), r2=resid(y2, mk))).groupby("b").mean()
ax[2].scatter(b["r1"], b["r2"], s=3, alpha=0.4)
ax[2].set_xlabel("1-Mb block residual, trait half A"); ax[2].set_ylabel("trait half B")
ax[2].set_title(f"split-half r: SNPs {res['snp_r'].mean():.2f} (null {res['snp_r_null'].mean():.2f}), 1-Mb {res['block_r'].mean():.2f}", fontsize=9)
fig.tight_layout(); fig.savefig(outp + ".png", dpi=150)
pd.DataFrame(dict(l=dfb["l"], y=dfb["y"], ysem=dfb["ys"])).to_csv(outp + ".binned.tsv", sep="\t", index=False)
