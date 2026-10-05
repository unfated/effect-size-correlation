"""Build 1000 Genomes EUR LD windows in the same format as the UKB windows.

Windows: 3 Mb, one every 2 Mb (starts 1, 2000001, 4000001, ... as in the UKB
release), so the stage-3 driver can use either reference unchanged.
Reference: 503 EUR samples, HapMap3 SNPs, PLINK (built by the package thread,
/mnt/project-files/papers/software/ref/eur_hm3_chr*.{bed,bim,fam}). SNPs are
matched to the Pan-UKB universe by chr:pos and alleles; R is signed relative to
the universe ALT allele (ID = chr:pos:ref:alt). Strand-ambiguous SNPs are kept
when alleles match exactly as listed (no strand flips are attempted).

Usage: 12_kg_ld_windows.py <ref_dir> <snp_universe.tsv> <out_dir> [chr ...]
"""
import sys, os
import numpy as np
import pandas as pd

ref_dir, uni_path, out = sys.argv[1:4]
chroms = [int(c) for c in sys.argv[4:]] or list(range(1, 23))
os.makedirs(out, exist_ok=True)
uni = pd.read_csv(uni_path, sep="\t")
parts = uni.ID.str.split(":", expand=True)
uni["ref"], uni["alt"] = parts[2], parts[3]

def read_bed(prefix, n, cols):
    nb = (n + 3) // 4
    raw = np.fromfile(prefix + ".bed", dtype=np.uint8)[3:].reshape(-1, nb)[cols]
    g = np.unpackbits(raw[:, :, None], axis=2, bitorder="little").reshape(len(cols), nb, 4, 2)
    lo, hi = g[..., 0], g[..., 1]
    code = (lo + 2 * hi).reshape(len(cols), nb * 4)[:, :n]       # 0=hom A1, 1=missing, 2=het, 3=hom A2
    lut = np.array([2.0, np.nan, 1.0, 0.0])
    return lut[code].T                                            # n x snps, counts of A1

for c in chroms:
    pre = os.path.join(ref_dir, f"eur_hm3_chr{c}")
    if not os.path.exists(pre + ".fam"):
        print(f"chr{c}: reference not ready, skipped"); continue
    n = sum(1 for _ in open(pre + ".fam"))
    bim = pd.read_csv(pre + ".bim", sep=r"\s+", header=None, names=["chr", "rs", "cm", "pos", "a1", "a2"])
    bim["row"] = np.arange(len(bim))
    u = uni[uni.CHR == c]
    mg = u.merge(bim, left_on="BP", right_on="pos")
    same = (mg.alt == mg.a1) & (mg.ref == mg.a2)     # A1 count = ALT count
    flip = (mg.alt == mg.a2) & (mg.ref == mg.a1)
    mg = mg[same | flip].copy()
    mg["sign"] = np.where(same[same | flip], 1.0, -1.0)
    mg = mg.drop_duplicates("ID").sort_values("BP").reset_index(drop=True)
    G = read_bed(pre, n, mg.row.values)
    mu = np.nanmean(G, 0); G = np.where(np.isnan(G), mu, G)
    sd = G.std(0); ok = sd > 0
    mg, G = mg[ok].reset_index(drop=True), G[:, ok]
    Gs = (G - G.mean(0)) / G.std(0) * mg.sign.values
    maxbp = mg.BP.max()
    start = 1
    nwin = 0
    while start <= maxbp:
        end = start + 3_000_000
        sel = np.where((mg.BP.values >= start) & (mg.BP.values < end))[0]
        if len(sel) > 1:
            R = (Gs[:, sel].T @ Gs[:, sel]) / n
            np.fill_diagonal(R, 1.0)
            name = f"chr{c}_{start}_{end}"
            mg.loc[sel, ["CHR", "BP", "ID", "RSID", "L2_UKB_EUR", "sign"]].to_csv(
                os.path.join(out, name + ".snps.tsv"), sep="\t", index=False)
            R.astype(np.float32).tofile(os.path.join(out, name + ".R.f32"))
            nwin += 1
        start += 2_000_000
    print(f"chr{c}: {len(mg)} of {len(u)} universe SNPs matched in 1000G EUR; {nwin} windows")
