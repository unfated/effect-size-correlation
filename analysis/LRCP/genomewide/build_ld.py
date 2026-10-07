"""Genome-wide LRCP (work plan step 13): LD among the LRCQ candidate tags of each
LRCQ window from the 1000G EUR HapMap3 reference, signed to the Pan-UKB alt allele.
Tags are LRCQ clump representatives (r2 < 0.5 within window), so within-window LD
is partial; the conditional test is calibrated with a 1000G reference (S5).
usage: build_ld.py <candidate_tags.tsv> <ref_dir> <out_dir>
Writes out_dir/w<window>.snps.tsv + .R.f32 (lrcpq::ld_from_windows format) and
out_dir/ld_summary.tsv (window, n_cand, n_matched, max_r2)."""
import sys, os, numpy as np, pandas as pd
CAND, REF, OUT = sys.argv[1:4]
os.makedirs(OUT, exist_ok=True)
c = pd.read_csv(CAND, sep="\t")
c["ref"] = c.ID.str.split(":").str[2]; c["alt"] = c.ID.str.split(":").str[3]
rows = []
for chrom, cc in c.groupby("CHR"):
    bim = pd.read_csv(f"{REF}/eur_hm3_chr{chrom}.bim", sep=r"\s+", header=None,
                      names=["chr", "id", "cm", "pos", "a1", "a2"])
    nfam = sum(1 for _ in open(f"{REF}/eur_hm3_chr{chrom}.fam")); nb = (nfam + 3) // 4
    mg = cc.merge(bim.reset_index(), left_on="BP", right_on="pos")
    ok = ((mg.alt == mg.a1) & (mg.ref == mg.a2)) | ((mg.alt == mg.a2) & (mg.ref == mg.a1))
    mg = mg[ok].drop_duplicates("ID").copy()
    mg["sgn"] = np.where(mg.alt == mg.a1, 1.0, -1.0)
    with open(f"{REF}/eur_hm3_chr{chrom}.bed", "rb") as fh:
        for win, g in mg.groupby("window"):
            G = np.empty((nfam, len(g)))
            for j, (bi, s) in enumerate(zip(g["index"].values, g.sgn.values)):
                fh.seek(3 + int(bi) * nb); raw = np.frombuffer(fh.read(nb), dtype=np.uint8)
                codes = np.stack([(raw >> k) & 3 for k in (0, 2, 4, 6)], 1).ravel()[:nfam]
                x = np.select([codes == 0, codes == 2, codes == 3], [2.0, 1.0, 0.0], np.nan)
                x[np.isnan(x)] = np.nanmean(x); G[:, j] = x * s
            keep = G.std(0) > 0; G = G[:, keep]; g = g[keep]
            if len(g) == 0: continue
            R = np.atleast_2d(np.corrcoef(G.T)) if len(g) > 1 else np.ones((1, 1))
            nm = f"w{win}"
            out = g[["CHR", "BP", "ID", "RSID"]].copy(); out["sign"] = 1
            out.to_csv(f"{OUT}/{nm}.snps.tsv", sep="\t", index=False)
            R.astype(np.float32).tofile(f"{OUT}/{nm}.R.f32")
            r2 = (R ** 2)[np.triu_indices(len(g), 1)]
            rows.append((win, chrom, int((cc.window == win).sum()), len(g), float(r2.max()) if len(r2) else 0.0))
    print("chr", chrom, "matched", len(mg), "of", len(cc), flush=True)
pd.DataFrame(rows, columns=["window", "CHR", "n_cand", "n_matched", "max_r2"]).to_csv(
    f"{OUT}/ld_summary.tsv", sep="\t", index=False)
