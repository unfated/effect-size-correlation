#!/usr/bin/env python3
"""Inputs for the LRCQ annotation regression (methods: tau = (A'D^2A)^-1 A'D ybar).

For every core SNP of every LD window this writes
  ybar_k = sum_a s_a (z_ka^2 - c_aa) / sum_a s_a^2   (trait-pooled response, observed traits only)
  l2_k   = (D 1)_k                                  (window LD score, bias-corrected r^2)
  DA_k   = (D A)_k                                  (annotation-stratified window LD scores)
so that E[ybar] = D w and, if w = A tau, E[ybar] = (DA) tau.

Annotations: baseline-LF v2.2 UKB (Gazal et al. 2019; Weissbrod et al. 2020). The
_lowfreq/_common pairs are disjoint by MAF, so their sum is the plain annotation.

Usage: 41_annot_ld.py <ld_dir> <annot_dir> <zprefix> <ldsc_prefix> <N_ref> <out_prefix> [trait_ids_file]
"""
import glob
import os
import re
import sys

import numpy as np
import pandas as pd

ld_dir, annot_dir, zp, lp, N_ref, outp = sys.argv[1:7]
N_ref = float(N_ref)
subset = open(sys.argv[7]).read().split() if len(sys.argv) > 7 else None

uni = pd.read_csv("/home/user/data/ref/snp_universe.tsv", sep="\t")
m_all = len(uni)
# M = number of SNPs whose LD enters R (the 1,094,844 HapMap3 universe SNPs), so that
# n h2 l / M matches E[chi2] - c with l from the same reference (lrcpq::check_scale:
# 1.13 here, 7.0 with M_5_50). Renormalising w-hat to mean 1 makes estimates and
# model SEs exactly invariant to M; M_5_50 was used for the published run.
M = 1094844.0

# ---- annotations aligned to the universe
ukeys = set(uni["CHR"].astype(str) + ":" + uni["BP"].astype(str))
parts = []
for f in sorted(glob.glob(os.path.join(annot_dir, "*.annot.gz"))):
    for ch in pd.read_csv(f, sep="\t", chunksize=200000):
        k = ch["CHR"].astype(str) + ":" + ch["BP"].astype(str)
        ch = ch[k.isin(ukeys).to_numpy()]
        if len(ch):
            num = ch.columns.difference(["CHR", "BP", "SNP", "CM"])
            ch[num] = ch[num].astype(np.float32)
            parts.append(ch)
an = pd.concat(parts, ignore_index=True)
cols = [c for c in an.columns if c not in ("CHR", "BP", "SNP", "CM")]
base = {}
for c in cols:
    b = re.sub(r"_(lowfreq|common)$", "", c)
    base.setdefault(b, []).append(c)
A_df = pd.DataFrame({b: an[cs].sum(axis=1) for b, cs in base.items()})
A_df.insert(0, "base", 1.0)
A_df["key"] = an["CHR"].astype(str) + ":" + an["BP"].astype(str)
A_df = A_df.drop_duplicates("key").set_index("key")
del an
ukey = uni["CHR"].astype(str) + ":" + uni["BP"].astype(str)
A = A_df.reindex(ukey).to_numpy(np.float32)
has_annot = ~np.isnan(A[:, 0])
A = np.nan_to_num(A)
anames = list(A_df.columns)
print(f"{len(anames)} annotations; {has_annot.sum()} of {m_all} universe SNPs annotated", flush=True)

# ---- trait scales and intercepts
tr = pd.read_csv(zp + ".traits.tsv", sep="\t")
q_all = len(tr)
h2 = pd.read_csv(lp + ".h2.tsv", sep="\t").set_index("trait_id").loc[tr["trait_id"]]
cols_t = np.arange(q_all) if subset is None else np.where(tr["trait_id"].isin(subset))[0]
s = (tr["N"].to_numpy() * h2["h2_panukb_ldsc"].to_numpy() / M)[cols_t]
c = h2["intercept_panukb"].to_numpy()[cols_t]
Z = np.memmap(zp + ".Z.f32", dtype=np.float32, mode="r", shape=(m_all, q_all))

# ---- windows
uidx = pd.Series(np.arange(m_all), index=uni["ID"])
wins = []
for f in glob.glob(os.path.join(ld_dir, "*.snps.tsv")):
    nm = os.path.basename(f)[: -len(".snps.tsv")]
    ch, st, en = re.match(r"chr(\d+)_(\d+)_(\d+)", nm).groups()
    wins.append((int(ch), int(st) - 1, int(en) - 1, nm))
wins.sort()
w = pd.DataFrame(wins, columns=["chr", "start", "end", "name"])
w["lo"] = w["start"] + 0.5e6
w["hi"] = w["start"] + 2.5e6
w.loc[~w["chr"].duplicated(), "lo"] = 0
w.loc[~w["chr"].duplicated(keep="last"), "hi"] = np.inf
w = w[~((w["chr"] == 6) & (w["end"] > 25e6) & (w["start"] < 34e6))]

out_id, out_y, out_l2, out_DA, out_win, out_nobs = [], [], [], [], [], []
for _, r in w.iterrows():
    sn = pd.read_csv(os.path.join(ld_dir, r["name"] + ".snps.tsv"), sep="\t")
    p = len(sn)
    R = np.fromfile(os.path.join(ld_dir, r["name"] + ".R.f32"), dtype=np.float32).reshape(p, p).astype(np.float64)
    D = R * R
    D = D - (1 - D) / (N_ref - 2)
    np.fill_diagonal(D, 1.0)
    ui = uidx.loc[sn["ID"]].to_numpy()
    core = np.where((sn["BP"] >= r["lo"]) & (sn["BP"] < r["hi"]))[0]
    if len(core) == 0:
        continue
    Aw = A[ui]
    DA = D[core] @ Aw
    l2 = D[core].sum(axis=1)
    zc = np.asarray(Z[ui[core]][:, cols_t], dtype=np.float64)
    ok = ~np.isnan(zc)
    resp = np.where(ok, (zc ** 2 - c) * s, 0.0).sum(axis=1)
    den = (ok * s ** 2).sum(axis=1)
    out_id.append(sn["ID"].to_numpy()[core]); out_y.append(resp / den); out_l2.append(l2)
    out_DA.append(DA.astype(np.float32)); out_win.append(np.repeat(r["name"], len(core)))
    out_nobs.append(ok.sum(axis=1))
    print(r["name"], len(core), flush=True)

np.savez_compressed(outp + ".npz", ID=np.concatenate(out_id), ybar=np.concatenate(out_y),
                    l2=np.concatenate(out_l2), DA=np.concatenate(out_DA), window=np.concatenate(out_win),
                    n_obs=np.concatenate(out_nobs), annot_names=np.array(anames),
                    A=A[uidx.loc[np.concatenate(out_id)].to_numpy()].astype(np.float32))
print("done", flush=True)
