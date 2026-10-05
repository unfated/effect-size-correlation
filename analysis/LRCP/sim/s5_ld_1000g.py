"""S5 input: 1000G EUR (503) LD for the SNPs of the targeted UKB LD blocks, signed to the
Pan-UKB alt allele, restricted to SNPs present in both references.
usage: s5_ld_1000g.py <ukb_ld_dir> <ref_dir> <out_dir>"""
import sys, glob, os, numpy as np, pandas as pd
LD, REF, OUT = sys.argv[1:4]
for f in sorted(glob.glob(f"{LD}/*.snps.tsv")):
    nm = os.path.basename(f)[:-9]; c = int(nm.split("_")[0][3:])
    u = pd.read_csv(f, sep="\t"); m = len(u)
    Ru = np.fromfile(f.replace(".snps.tsv", ".R.f32"), dtype=np.float32).reshape(m, m)
    bim = pd.read_csv(f"{REF}/eur_hm3_chr{c}.bim", sep=r"\s+", header=None, names=["chr", "id", "cm", "pos", "a1", "a2"])
    nfam = sum(1 for _ in open(f"{REF}/eur_hm3_chr{c}.fam"))
    u["ref"] = u.ID.str.split(":").str[2]; u["alt"] = u.ID.str.split(":").str[3]
    mg = u.reset_index().merge(bim.reset_index(), left_on="BP", right_on="pos", suffixes=("", "_b"))
    same = (mg.alt == mg.a1) & (mg.ref == mg.a2); flip = (mg.alt == mg.a2) & (mg.ref == mg.a1)
    mg = mg[same | flip].drop_duplicates("ID"); sgn = np.where((mg.alt == mg.a1) & (mg.ref == mg.a2), 1.0, -1.0)[(same | flip)[mg.index]] if False else np.where(mg.alt == mg.a1, 1.0, -1.0)
    nb = (nfam + 3) // 4
    with open(f"{REF}/eur_hm3_chr{c}.bed", "rb") as fh:
        G = np.empty((nfam, len(mg)))
        for j, bi in enumerate(mg.index_b.values):
            fh.seek(3 + int(bi) * nb); raw = np.frombuffer(fh.read(nb), dtype=np.uint8)
            codes = np.stack([(raw >> s) & 3 for s in (0, 2, 4, 6)], 1).ravel()[:nfam]
            g = np.select([codes == 0, codes == 2, codes == 3], [2.0, 1.0, 0.0], np.nan)  # count of A1
            g[np.isnan(g)] = np.nanmean(g); G[:, j] = g * sgn[j]
    sd = G.std(0); ok = sd > 0
    G = G[:, ok]; mg = mg[ok]
    R1 = np.corrcoef(G.T)
    idx = mg["index"].values
    mg[["CHR", "BP", "ID", "RSID"]].to_csv(f"{OUT}/{nm}.snps.tsv", sep="\t", index=False)
    R1.astype(np.float32).tofile(f"{OUT}/{nm}.R1kg.f32")
    Ru[np.ix_(idx, idx)].astype(np.float32).tofile(f"{OUT}/{nm}.Rukb.f32")
    d = np.abs(R1 - Ru[np.ix_(idx, idx)])[np.triu_indices(len(idx), 1)]
    print(nm, "UKB", m, "matched", len(idx), "mean|dR|", d.mean().round(4), "max", d.max().round(3))
