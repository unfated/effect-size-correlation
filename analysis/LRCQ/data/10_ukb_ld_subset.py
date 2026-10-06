"""Subset one UKB in-sample LD window (Weissbrod et al. 2020 PolyFun release,
s3://broad-alkesgroup-ukbb-ld/UKBB_LD/<chr>_<start>_<end>.{gz,npz}) to the
LRCQ SNP universe (Pan-UKB EUR HapMap3 SNPs), matching on chr:pos and alleles.

Usage: 10_ukb_ld_subset.py <prefix> <snp_universe.tsv> <out_prefix>
Writes <out_prefix>.snps.tsv (universe rows kept, with sign flip vs UKB A1/A2)
and <out_prefix>.R.f32 (dense float32 m*m, row-major), plus a summary line.
"""
import sys
import numpy as np
import pandas as pd

prefix, uni_path, out = sys.argv[1:4]
snps = pd.read_csv(prefix + ".gz", sep="\t")
d = np.load(prefix + ".npz")
uni = pd.read_csv(uni_path, sep="\t")
chrom = int(snps.chromosome.iloc[0])
lo, hi = snps.position.min(), snps.position.max()
u = uni[(uni.CHR == chrom) & (uni.BP >= lo) & (uni.BP <= hi)].copy()
ref = u.ID.str.split(":").str[2]
alt = u.ID.str.split(":").str[3]
u["ref"], u["alt"] = ref.values, alt.values
snps["i"] = np.arange(len(snps))
mg = u.merge(snps, left_on="BP", right_on="position", how="inner")
same = (mg.ref == mg.allele1) & (mg.alt == mg.allele2)
flip = (mg.ref == mg.allele2) & (mg.alt == mg.allele1)
mg = mg[same | flip].copy()
mg["sign"] = np.where((mg.ref == mg.allele1) & (mg.alt == mg.allele2), 1, -1)
mg = mg.drop_duplicates("ID").sort_values("BP")
idx = mg.i.values
# Triplets are sorted by row, so only the row ranges of kept SNPs are read
# (much faster than building a sparse matrix of the whole window).
row_all = d["row"]
assert np.all(np.diff(row_all[::1000]) >= 0), "expected row-sorted triplets"
lo_i = np.searchsorted(row_all, idx, "left")
hi_i = np.searchsorted(row_all, idx, "right")
sel = np.concatenate([np.arange(a, b) for a, b in zip(lo_i, hi_i)])
newi = np.full(int(d["shape"][0]), -1, dtype=np.int64)
newi[idx] = np.arange(len(idx))
row, col = newi[row_all[sel]], newi[d["col"][sel]]
keep = col >= 0
sub = np.zeros((len(idx), len(idx)))
sub[row[keep], col[keep]] = d["data"][sel][keep]
# file stores one triangle; symmetrise
sub = sub + sub.T - np.diag(np.diag(sub))
np.fill_diagonal(sub, 1.0)
sub = sub * np.outer(mg.sign.values, mg.sign.values)
mg[["CHR", "BP", "ID", "RSID", "L2_UKB_EUR", "sign"]].to_csv(out + ".snps.tsv", sep="\t", index=False)
sub.astype(np.float32).tofile(out + ".R.f32")
print(f"{prefix}: {len(u)} universe SNPs in range, {len(mg)} matched; "
      f"min eig {np.linalg.eigvalsh(sub).min():.3g}")
