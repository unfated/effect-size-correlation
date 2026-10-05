#!/usr/bin/env python3
"""Simulation figures and the main simulation table for the LRCQ manuscript.

Usage: sim_figures.py <results_sim_dir> <fig_dir>
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

rd, fd = sys.argv[1:3]
os.makedirs(fd, exist_ok=True)
s = pd.read_csv(os.path.join(rd, "sim_summary.tsv"), sep="\t").set_index("scenario")
pr = pd.read_csv(os.path.join(rd, "sim_per_replicate.tsv"), sep="\t")
an = pd.read_csv(os.path.join(rd, "sim_annotation.tsv"), sep="\t")
order = ["base", "q30", "q100", "q1000", "sparse", "dense", "inf", "ukb", "strat", "strat_miss",
         "rb_local", "rb_cross", "het", "ref500", "tag08", "tag02"]
lab = {"base": "baseline", "q30": "q = 30", "q100": "q = 100", "q1000": "q = 1000", "sparse": "0.5% non-null",
       "dense": "10% non-null", "inf": "all non-null", "ukb": "UKB-like traits", "strat": "intercept 1.1",
       "strat_miss": "intercept missed", "rb_local": "local ρ = 0.5", "rb_cross": "cross-block ρ", "het": "heterogeneous w",
       "ref500": "LD from 500", "tag08": "tags r² < 0.8", "tag02": "tags r² < 0.2"}
s = s.loc[order]
x = np.arange(len(order))
se = lambda col: pr.groupby("scenario")[col].sem().reindex(order)

fig, ax = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
for i, (col, name) in enumerate([("pearson_tag_ols", "LRCQ OLS"), ("pearson_tag", "LRCQ WLS"), ("pearson_tag_eq", "LRCQ equal-weight")]):
    ax[0].errorbar(x + (i - 1) * 0.2, s[col], yerr=1.96 * se(col), fmt="o", ms=4, label=name)
ax[0].set_ylabel("Pearson r with tag estimand"); ax[0].legend(fontsize=8); ax[0].set_ylim(0.6, 1.0)
comps = [("auc10_tag_ols", "LRCQ OLS"), ("auc10_ridge_clump", "ridge LRCQ (clump sum)"), ("auc10_mean_chi2", "mean χ²"),
         ("auc10_ldnorm_chi2", "LD-normalised χ²"), ("auc10_omnibus", "omnibus Wald"), ("auc10_n_sig", "# traits P<5e-8")]
for i, (col, name) in enumerate(comps):
    ax[1].errorbar(x + (i - 2.5) * 0.12, s[col], yerr=1.96 * se(col), fmt="o", ms=3, label=name)
ax[1].set_ylabel("AUC, top 10% enriched tags"); ax[1].legend(fontsize=7, ncol=3); ax[1].set_ylim(0.5, 1.0)
ax[1].set_xticks(x); ax[1].set_xticklabels([lab[o] for o in order], rotation=45, ha="right")
fig.tight_layout(); fig.savefig(os.path.join(fd, "fig2_sim_accuracy.png"), dpi=150); plt.close(fig)

fig, ax = plt.subplots(1, 3, figsize=(13, 4))
for i, (e, name) in enumerate([("tag_ols", "OLS"), ("tag", "WLS"), ("tag_eq", "equal-weight")]):
    ax[0].plot(x + (i - 1) * 0.2, s["zres_sd_" + e], "o", label=name)
    ax[1].plot(x + (i - 1) * 0.2, s["cover95_" + e], "o", label=name)
    ax[2].plot(x + (i - 1) * 0.2, s["slope_" + e] if e != "tag_eq" else s["slope_tag_eq_vs_eqtarget"], "o", label=name)
ax[0].axhline(1, c="k", lw=0.8); ax[0].set_title("SD of (ŵ − w)/SE")
ax[1].axhline(0.95, c="k", lw=0.8); ax[1].set_title("95% CI coverage")
ax[2].axhline(1, c="k", lw=0.8); ax[2].set_title("slope of ŵ on its estimand")
for a in ax:
    a.set_xticks(x); a.set_xticklabels([lab[o] for o in order], rotation=60, ha="right", fontsize=7)
ax[0].legend(fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(fd, "fig3_sim_calibration.png"), dpi=150); plt.close(fig)

sc = ["base", "q30", "ukb", "dense", "inf"]
fig, ax = plt.subplots(1, len(sc), figsize=(15, 3.5), sharey=True)
for j, c in enumerate(sc):
    d = an[an["scenario"] == c].sort_values("annot")
    xx = d["annot"].to_numpy()
    ax[j].plot(xx, d["true"], "k_", ms=18, mew=2, label="true E_c")
    ax[j].plot(xx, d["clump_target"], "_", c="grey", ms=18, mew=2, label="clump-level target")
    for k, (col, name, off) in enumerate([("lrcq_clump", "LRCQ clump E_c", -0.15), ("lrcq_ridge", "ridge LRCQ", 0),
                                          ("sldsc_pooled", "annotation regression", 0.15)]):
        ax[j].errorbar(xx + off, d[col], yerr=d["sd_" + col], fmt="o", ms=4, label=name)
    ax[j].set_title(lab[c]); ax[j].set_xticks(xx); ax[j].set_xlabel("annotation level")
ax[0].set_ylabel("fold enrichment"); ax[0].legend(fontsize=7)
fig.tight_layout(); fig.savefig(os.path.join(fd, "fig4_sim_annotation.png"), dpi=150); plt.close(fig)

cols = ["pearson_tag_ols", "slope_tag_ols", "bias_tag_ols", "rmse_tag_ols", "pearson_tag", "slope_tag", "rmse_tag",
        "pearson_tag_eq", "zres_sd_tag_ols", "zres_sd_tag_eq", "typeI_null_tag_ols", "cover95_tag_ols",
        "auc10_tag_ols", "auc10_mean_chi2", "auc10_ldnorm_chi2", "auc10_omnibus", "auc10_n_sig",
        "spearman_nonnull_tag_ols", "spearman_nonnull_mean_chi2", "methodD_n_disc", "methodD_fdp", "methodD_power", "n_rep"]
tab = s[[c for c in cols if c in s.columns]].copy()
tab.index = [lab[o] for o in order]
tab.round(3).to_csv(os.path.join(rd, "table_sim_main.tsv"), sep="\t")
print(tab.round(3).to_string())
