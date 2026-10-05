#!/usr/bin/env python3
"""Post-process stage-3 tag-level LRCQ output.

  * renormalise so that the clump totals sum to the number of core SNPs (mean w = 1)
  * per-tag clump mean enrichment w_hat / g, omni-locus tests H0: W <= g (mean
    enrichment <= 1) and H0: W <= 0, BH q-values
  * method-of-moments spread of true clump enrichment: var(w_hat/g) - mean(se^2/g^2)
  * model-free clump enrichment E_c = sum f_kc W_k / sum f_kc g_k for annotations
    (needs the 41_annot_ld.py npz for the annotation matrix), block-jackknife SEs

Usage: 40_postprocess.py <s3_prefix> <out_dir> [annot_npz]
  s3_prefix: 30_stage3_lrcq.R output without .tsv (needs .tsv and .snp2tag.tsv)
"""
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

sp, od = sys.argv[1:3]
annot = sys.argv[3] if len(sys.argv) > 3 else None
os.makedirs(od, exist_ok=True)
lab = os.path.basename(sp)

t = pd.read_csv(sp + ".tsv", sep="\t")
s2t = pd.read_csv(sp + ".snp2tag.tsv", sep="\t")
n_core = len(s2t)
scale = n_core / t["w_raw"].sum()
t["W"] = t["w_raw"] * scale
t["W_se"] = t["se"] * scale
# trait-cluster jackknife (30_stage3_lrcq.R, ols/equal): renormalise every replicate
jkf = sp + ".jk.f32"
K = None
if os.path.exists(jkf):
    jk = np.fromfile(jkf, dtype=np.float32).astype(np.float64)
    K = jk.size // len(t)
    jk = jk.reshape(len(t), K)
    jk = jk * (n_core / jk.sum(0))
    t["W_se_model"] = t["W_se"]
    t["W_se"] = np.sqrt((K - 1) / K * ((jk - jk.mean(1, keepdims=True)) ** 2).sum(1))
t["enr"] = t["W"] / t["clump_n"]
t["enr_se"] = t["W_se"] / t["clump_n"]
dfree = (K - 1) if K else np.inf


def bh(p):
    p = np.asarray(p); n = len(p); o = np.argsort(p)
    q = np.empty(n); q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1)


t["z_gt1"] = (t["W"] - t["clump_n"]) / t["W_se"]
t["p_gt1"] = stats.t.sf(t["z_gt1"], dfree)
t["q_gt1"] = bh(t["p_gt1"])
t["z_gt0"] = t["W"] / t["W_se"]
t["p_gt0"] = stats.t.sf(t["z_gt0"], dfree)
t["q_gt0"] = bh(t["p_gt0"])

summ = dict(label=lab, se_type="trait-cluster jackknife, K=%s" % K if K else "model (Theorem 3.10)", n_core_snps=int(n_core), n_tags=int(len(t)), scale=float(scale),
            median_clump_n=float(t["clump_n"].median()),
            median_se_enr_singleton=float(t.loc[t["clump_n"] == 1, "enr_se"].median()),
            median_se_W=float(t["W_se"].median()),
            mean_z_gt1=float(t["z_gt1"].mean()), sd_z_gt1=float(t["z_gt1"].std()),
            frac_W_neg=float((t["W"] < 0).mean()),
            n_omni_q05=int((t["q_gt1"] < 0.05).sum()), n_omni_q10=int((t["q_gt1"] < 0.10).sum()),
            n_nonzero_q05=int((t["q_gt0"] < 0.05).sum()))
# spread of true clump-mean enrichment (weights g so that large clumps count by size)
g = t["clump_n"].to_numpy(); e = t["enr"].to_numpy(); v = (t["enr_se"] ** 2).to_numpy()
mu = np.average(e, weights=g)
summ["var_true_enr_mom"] = float(np.average((e - mu) ** 2, weights=g) - np.average(v, weights=g))
summ["var_obs_enr"] = float(np.average((e - mu) ** 2, weights=g))
# heterogeneity: sum of squared z under H0 w = 1 (approximate, ignores within-window covariance)
summ["chi2_flat"] = float((t["z_gt1"] ** 2).sum()); summ["df_flat"] = int(len(t))
for cmp in ["mean_chi2", "ldnorm_chi2", "n_sig", "omnibus"]:
    summ["spearman_enr_" + cmp] = float(stats.spearmanr(t["enr"], t[cmp]).correlation)
    summ["spearman_W_" + cmp] = float(stats.spearmanr(t["W"], t[cmp]).correlation)

# ---- model-free clump enrichment of annotations
if annot:
    z = np.load(annot, allow_pickle=True)
    an = list(z["annot_names"]); A = pd.DataFrame(z["A"], index=z["ID"], columns=an)
    A = A[~A.index.duplicated()]
    s2 = s2t[s2t["ID"].isin(A.index)]
    F = A.loc[s2["ID"]].groupby(s2["tag_ID"].to_numpy()).sum()   # SNP counts per clump and annotation
    F = F.reindex(t["ID"]).fillna(0.0)
    G = F.to_numpy()
    W = t["W"].to_numpy()
    nb = 200
    blk = np.minimum((np.arange(len(t)) * nb) // len(t), nb - 1)
    num_b = np.vstack([(G[blk == b] * W[blk == b, None]).sum(0) for b in range(nb)])
    # E_c = sum_k F_kc W_k / g_k  /  sum_k F_kc  (F_kc / g_k = f_kc, sum f_kc g_k = sum F_kc)
    Wg = W / t["clump_n"].to_numpy()
    num_b = np.vstack([(G[blk == b] * Wg[blk == b, None]).sum(0) for b in range(nb)])
    den_b = np.vstack([G[blk == b].sum(0) for b in range(nb)])
    num, den = num_b.sum(0), den_b.sum(0)
    E = num / den
    jk = np.vstack([(num - num_b[b]) / (den - den_b[b]) for b in range(nb)])
    se = np.sqrt((nb - 1) / nb * ((jk - jk.mean(0)) ** 2).sum(0))
    prop = den / den[an.index("base")]
    res = pd.DataFrame(dict(annotation=an, prop_snps=prop, E_clump=E, E_clump_se=se,
                            z_vs1=(E - 1) / se, p_vs1=2 * stats.norm.sf(np.abs((E - 1) / se))))
    res.to_csv(os.path.join(od, lab + ".annot_clump.tsv"), sep="\t", index=False)

out_cols = ["window", "CHR", "BP", "ID", "RSID", "clump_n", "W", "W_se"] + (["W_se_model"] if K else []) + [ "enr", "enr_se", "z_gt1", "p_gt1",
            "q_gt1", "z_gt0", "p_gt0", "q_gt0", "l2_window", "n_missing", "mean_chi2", "ldnorm_chi2", "n_sig", "omnibus"]
t[out_cols].to_csv(os.path.join(od, lab + ".tags.tsv.gz"), sep="\t", index=False, compression="gzip")
json.dump(summ, open(os.path.join(od, lab + ".summary.json"), "w"), indent=1)
print(json.dumps(summ, indent=1))
