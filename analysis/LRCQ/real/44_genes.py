#!/usr/bin/env python3
"""Gene-level phenome-wide enrichment and its relation to LoF constraint.

Per core SNP k, e_k = W_t / g_t for its tag t (clump-mean enrichment). Gene enrichment
E_g = mean of e_k over HapMap3 SNPs within [start - 10 kb, end + 10 kb] (GRCh37,
gnomAD v2.1.1 gene table, protein-coding genes). E_g = sum_t (n_tg / n_g) W_t / g_t,
so its SE is built from the tag model SEs (Theorem 3.10), ignoring covariance between
tag estimates (tags are pruned at r2 < 0.5). The trait-cluster jackknife SE is written
as E_se_jk, a diagnostic only (invalid under sample overlap; theory S3.6).
Also: mean E_g by gnomAD LOEUF decile, and the genes of the user's original
observation table (observed data hints/table 100+100 genes.xlsx) plus TP53/MDM2,
against their counts of significant GWAS.

Usage: 44_genes.py <s3_prefix> <gnomad_lof_metrics.bgz> <out_prefix> [hints.xlsx]
"""
import sys

import numpy as np
import pandas as pd
from scipy import stats

sp, gf, outp = sys.argv[1:4]
hints = sys.argv[4] if len(sys.argv) > 4 else None
FLANK = 10_000

t = pd.read_csv(sp + ".tsv", sep="\t")
s2t = pd.read_csv(sp + ".snp2tag.tsv", sep="\t")
n_core = len(s2t)
jk = np.fromfile(sp + ".jk.f32", dtype=np.float32).astype(np.float64).reshape(len(t), -1)
K = jk.shape[1]
W = np.c_[t["w_raw"].to_numpy(), jk]
W = W * (n_core / W.sum(0))                      # column 0 = full data, 1..K = replicates
E_tag = W / t["clump_n"].to_numpy()[:, None]
V_tag = (t["se"].to_numpy() * n_core / t["w_raw"].sum() / t["clump_n"].to_numpy()) ** 2
tag_row = pd.Series(np.arange(len(t)), index=t["ID"])
snp_tag = tag_row.reindex(s2t["tag_ID"]).to_numpy()
ok = ~np.isnan(snp_tag)
ids = s2t["ID"].to_numpy()[ok]
e_snp = E_tag[snp_tag[ok].astype(int)]          # n_snp x (K+1)
t_snp = snp_tag[ok].astype(int)
chrom = np.array([int(x.split(":")[0]) for x in ids]); bp = np.array([int(x.split(":")[1]) for x in ids])
o = np.lexsort((bp, chrom)); chrom, bp, e_snp, t_snp = chrom[o], bp[o], e_snp[o], t_snp[o]
cs = np.vstack([np.zeros((1, e_snp.shape[1])), np.cumsum(e_snp, 0)])

g = pd.read_csv(gf, sep="\t", compression="gzip", low_memory=False)
g = g[g["chromosome"].astype(str).str.fullmatch(r"\d+")].copy()
g["chr"] = g["chromosome"].astype(int)
rows = []
for c, gc in g.groupby("chr"):
    idx = np.where(chrom == c)[0]
    if len(idx) == 0:
        continue
    b = bp[idx]
    lo = idx[0] + np.searchsorted(b, gc["start_position"].to_numpy() - FLANK)
    hi = idx[0] + np.searchsorted(b, gc["end_position"].to_numpy() + FLANK, side="right")
    n = hi - lo
    sums = cs[hi] - cs[lo]
    with np.errstate(invalid="ignore", divide="ignore"):
        means = sums / n[:, None]
    for i, (_, r) in enumerate(gc.iterrows()):
        if n[i] == 0:
            continue
        rep = means[i, 1:]
        se_jk = np.sqrt((K - 1) / K * ((rep - rep.mean()) ** 2).sum())
        tt, cnt = np.unique(t_snp[lo[i]:hi[i]], return_counts=True)
        se = np.sqrt(((cnt / n[i]) ** 2 * V_tag[tt]).sum())
        rows.append(dict(gene=r["gene"], chr=c, start=r["start_position"], end=r["end_position"], n_snps=int(n[i]),
                         n_tags=len(tt), E=means[i, 0], E_se=se, E_se_jk=se_jk, loeuf=r["oe_lof_upper"],
                         pLI=r.get("pLI", np.nan)))
gr = pd.DataFrame(rows)
gr["z"] = (gr["E"] - 1) / gr["E_se"]
gr["p_gt1"] = stats.norm.sf(gr["z"])
gr = gr.sort_values("p_gt1")
gr.to_csv(outp + ".genes.tsv.gz", sep="\t", index=False, compression="gzip")

# LOEUF deciles: mean of E_g within decile (each gene weighted equally); jackknife SE
gl = gr.dropna(subset=["loeuf"]).copy()
gl["decile"] = pd.qcut(gl["loeuf"], 10, labels=False) + 1
dec_rows = []
for d, gd in gl.groupby("decile"):
    dec_rows.append(dict(decile=d, loeuf_min=gd["loeuf"].min(), loeuf_max=gd["loeuf"].max(), n_genes=len(gd),
                         mean_E=gd["E"].mean(), sem_genes=gd["E"].std() / np.sqrt(len(gd)),
                         median_E=gd["E"].median()))
dec = pd.DataFrame(dec_rows)
rho = stats.spearmanr(gl["loeuf"], gl["E"])
dec.to_csv(outp + ".loeuf_deciles.tsv", sep="\t", index=False)
print(dec.round(3).to_string())
print(f"Spearman(LOEUF, E_g) = {rho.correlation:.3f}, P = {rho.pvalue:.2g}, genes = {len(gl)}")

if hints:
    x = pd.read_excel(hints, sheet_name=0)
    a = x[["Chr1 gene", "significant GWASs"]].set_axis(["gene", "n_gwas"], axis=1)
    b = x[["Chr2 gene", "significant GWASs.1"]].set_axis(["gene", "n_gwas"], axis=1)
    h = pd.concat([a, b]).dropna().drop_duplicates("gene")
    h = pd.concat([h, pd.DataFrame(dict(gene=["TP53", "MDM2"], n_gwas=[np.nan, np.nan]))])
    hm = h.merge(gr, on="gene", how="left")
    hm["genome_pct"] = [100 * (gr["E"] < e).mean() if pd.notna(e) else np.nan for e in hm["E"]]
    hm.to_csv(outp + ".hint_genes.tsv", sep="\t", index=False)
    v = hm.dropna(subset=["E", "n_gwas"])
    rr = stats.spearmanr(v["n_gwas"], v["E"])
    print(f"hint genes: {len(v)} matched; median genome percentile of E_g {hm['genome_pct'].median():.1f}; "
          f"Spearman(n_gwas, E_g) = {rr.correlation:.3f} (P = {rr.pvalue:.2g})")
    print(hm[["gene", "n_gwas", "n_snps", "E", "E_se", "genome_pct"]].round(2).to_string())
