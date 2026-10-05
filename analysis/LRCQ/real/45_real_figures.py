#!/usr/bin/env python3
"""Pan-UKB figures and tables: omni-locus Manhattan plot, top omni-loci table with
nearest genes, LOEUF gradient, LRCQ vs mean chi2.

Usage: 45_real_figures.py <results_real_dir> <fig_dir> <gnomad_lof_metrics.bgz> [label]
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

rd, fd, gf = sys.argv[1:4]
lab = sys.argv[4] if len(sys.argv) > 4 else "lrcq_all_ols"
t = pd.read_csv(os.path.join(rd, lab + ".tags.tsv.gz"), sep="\t")
g = pd.read_csv(gf, sep="\t", compression="gzip", low_memory=False, usecols=["gene", "chromosome", "start_position", "end_position", "oe_lof_upper"])
g = g[g["chromosome"].astype(str).str.fullmatch(r"\d+")]
g["chr"] = g["chromosome"].astype(int)

# nearest protein-coding gene (distance 0 if inside)
def nearest(c, bp):
    gc = g[g["chr"] == c]
    d = np.maximum(0, np.maximum(gc["start_position"] - bp, bp - gc["end_position"]))
    i = int(np.argmin(d.to_numpy()))
    return gc["gene"].iloc[i], int(d.iloc[i])

# top omni-loci: FDR 5%, collapse to one row per 1-Mb locus (best tag)
sig = t[t["q_gt1"] < 0.05].sort_values("p_gt1").copy()
sig["locus"] = sig["CHR"].astype(str) + ":" + (sig["BP"] // 1_000_000).astype(str)
loci = sig.drop_duplicates("locus").head(40).copy()
nn = [nearest(c, b) for c, b in zip(loci["CHR"], loci["BP"])]
loci["nearest_gene"] = [x[0] for x in nn]; loci["gene_dist_kb"] = [x[1] / 1000 for x in nn]
cols = ["CHR", "BP", "RSID", "nearest_gene", "gene_dist_kb", "clump_n", "W", "W_se", "enr", "z_gt1", "q_gt1",
        "mean_chi2", "n_sig"]
loci[cols].round(3).to_csv(os.path.join(rd, lab + ".top_loci.tsv"), sep="\t", index=False)
print(f"{len(sig)} tags at FDR 5% in {sig['locus'].nunique()} 1-Mb loci")
print(loci[cols].head(25).round(2).to_string())

# Fig 5: Manhattan of one-sided -log10 p for H0: mean enrichment <= 1
t = t.sort_values(["CHR", "BP"])
off = t.groupby("CHR")["BP"].max().cumsum().shift(fill_value=0)
x = t["BP"] + t["CHR"].map(off)
y = -np.log10(np.clip(t["p_gt1"], 1e-300, 1))
fig, ax = plt.subplots(2, 1, figsize=(13, 7), gridspec_kw=dict(height_ratios=[2, 1.3]))
col = np.where(t["CHR"] % 2 == 0, "#4c72b0", "#8fa8d0")
ax[0].scatter(x, np.minimum(y, 40), s=2, c=col, rasterized=True)
thr = t.loc[t["q_gt1"] < 0.05, "p_gt1"].max()
if pd.notna(thr):
    ax[0].axhline(-np.log10(thr), c="r", lw=0.8, ls="--")
for _, r in loci.head(15).iterrows():
    xi = r["BP"] + off[r["CHR"]]
    ax[0].annotate(r["nearest_gene"], (xi, min(-np.log10(max(r["p_gt1"], 1e-300)), 40)), fontsize=6, rotation=45)
mid = t.groupby("CHR").apply(lambda d: (d["BP"] + off[d.name]).median())
ax[0].set_xticks(mid); ax[0].set_xticklabels(mid.index, fontsize=7)
ax[0].set_ylabel("−log10 P (clump enrichment > 1)"); ax[0].set_title("LRCQ omni-loci across 444 Pan-UKB traits (capped at 40)")
# panel b: LRCQ enrichment vs mean chi2 (binned)
q = pd.qcut(t["mean_chi2"].rank(method="first"), 50, labels=False)
b = t.groupby(q).agg(mc=("mean_chi2", "mean"), e=("enr", "mean"), es=("enr", "sem"), l2=("l2_window", "mean"))
ax[1].errorbar(b["mc"], b["e"], yerr=1.96 * b["es"], fmt="o", ms=3)
ax[1].set_xscale("symlog"); ax[1].axhline(1, c="k", lw=0.6)
ax[1].set_xlabel("mean χ² − c across traits (tag, binned)"); ax[1].set_ylabel("LRCQ clump-mean enrichment")
fig.tight_layout(); fig.savefig(os.path.join(fd, "fig5_omniloci.png"), dpi=150); plt.close(fig)

# LOEUF gradient figure (if gene results exist)
dp = os.path.join(rd, "genes_all.loeuf_deciles.tsv")
if os.path.exists(dp):
    d = pd.read_csv(dp, sep="\t")
    fig, ax = plt.subplots(figsize=(5, 3.5))
    ax.errorbar(d["decile"], d["mean_E"], yerr=1.96 * d["sem_genes"], fmt="o-", label="mean")
    ax.plot(d["decile"], d["median_E"], "s--", label="median")
    ax.axhline(1, c="k", lw=0.6)
    ax.set_xlabel("gnomAD LOEUF decile (1 = most constrained)"); ax.set_ylabel("gene phenome-wide enrichment")
    ax.legend(fontsize=8); fig.tight_layout(); fig.savefig(os.path.join(fd, "fig7_constraint.png"), dpi=150)
