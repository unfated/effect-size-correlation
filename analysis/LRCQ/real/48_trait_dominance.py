#!/usr/bin/env python3
"""Single-trait dominance at omni-loci: for each tag in <lab>.top_loci.tsv, the share of
the trait-pooled OLS response sum_a s_a (z_ka^2 - c_aa) contributed by its largest trait,
and the effective number of contributing traits 1/sum(p_a^2) (p_a = positive shares).
Usage: 48_trait_dominance.py <zprefix> <ldsc_prefix> <results_real_dir> [label]
"""
import os
import sys

import numpy as np
import pandas as pd

zp, lp, rd = sys.argv[1:4]
lab = sys.argv[4] if len(sys.argv) > 4 else "lrcq_all_ols"
uni = pd.read_csv("/home/user/data/ref/snp_universe.tsv", sep="\t")
tr = pd.read_csv(zp + ".traits.tsv", sep="\t")
# M = number of SNPs whose LD enters R (the 1,094,844 HapMap3 universe SNPs), so that
# n h2 l / M matches E[chi2] - c with l from the same reference (lrcpq::check_scale:
# 1.13 here, 7.0 with M_5_50). Renormalising w-hat to mean 1 makes estimates and
# model SEs exactly invariant to M; M_5_50 was used for the published run.
M = 1094844.0
h2 = pd.read_csv(lp + ".h2.tsv", sep="\t").set_index("trait_id").loc[tr["trait_id"]]
s = tr["N"].to_numpy() * h2["h2_panukb_ldsc"].to_numpy() / M
c = h2["intercept_panukb"].to_numpy()
Z = np.memmap(zp + ".Z.f32", dtype=np.float32, mode="r", shape=(len(uni), len(tr)))
row = pd.Series(np.arange(len(uni)), index=uni["ID"])
top = pd.read_csv(os.path.join(rd, lab + ".top_loci.tsv"), sep="\t")
tags = pd.read_csv(os.path.join(rd, lab + ".tags.tsv.gz"), sep="\t", usecols=["CHR", "BP", "ID"])
top = top.merge(tags, on=["CHR", "BP"], how="left")
out = []
for _, r in top.iterrows():
    z = np.nan_to_num(np.asarray(Z[row[r["ID"]]], dtype=np.float64))
    contrib = s * (z ** 2 - c)
    pos = np.clip(contrib, 0, None)
    p = pos / pos.sum()
    j = int(np.argmax(contrib))
    out.append(dict(ID=r["ID"], nearest_gene=r["nearest_gene"], top_trait=tr["description"].iloc[j] if "description" in tr else tr["trait_id"].iloc[j],
                    top_trait_z=round(float(z[j]), 1), top_share=round(float(contrib[j] / contrib.sum()), 3),
                    eff_traits=round(float(1 / (p ** 2).sum()), 1)))
o = pd.DataFrame(out)
o.to_csv(os.path.join(rd, lab + ".top_loci_dominance.tsv"), sep="\t", index=False)
print(o.to_string())
print("median top-trait share", o["top_share"].median(), "; loci with share > 0.5:", int((o["top_share"] > 0.5).sum()), "of", len(o))
