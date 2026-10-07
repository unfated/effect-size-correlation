"""Assemble inputs for the targeted LRCP analysis (work plan step 12).

Joins the per-trait block Z scores (fetch_panukb_blocks.py) with the two
named motivating SNPs (TP53 rs78378222, MDM2 rs3730556; not in HapMap3, taken
from the motivating region pull), restricts the UKB in-sample LD windows to
each Berisa-Pickrell block, and computes preliminary stage 1-2 plug-ins.

usage: build_inputs.py <trait-list.tsv> <blocks_dir> <regions_dir> <ld_dir> <out_dir>
Writes to out_dir:
  Z.tsv.gz        SNP x trait Z (rows ordered by chr, pos; ids chr:pos:ref:alt)
  ld/<chrC_S_E>.snps.tsv + .R.f32   one LD "window" per block (lrcpq::ld_from_windows format)
  traits.tsv      key, description, N_eff, h2_obs, ldsc_intercept, indep
  C_hat.tsv.gz    trait x trait Z correlation at near-null SNPs (random windows)
"""
import os, sys, glob, numpy as np, pandas as pd

tl_path, BLK, REG, LD, OUT = sys.argv[1:6]
os.makedirs(f"{OUT}/ld", exist_ok=True)
BLOCKS = {"TP53": (17, 7317398, 8306425), "MDM2": (12, 67909729, 69826542),
          "PCSK9": (1, 54226262, 56413117), "LDLR": (19, 9238393, 11284028),
          "HMGCR": (5, 73759326, 75798866), "APOB": (2, 21050490, 23341383)}
EXTRA = {"17:7571752:T:G": "TP53_rs78378222", "12:69216521:T:G": "MDM2_rs3730556"}

tl = pd.read_csv(tl_path, sep="\t")
tl["key"] = tl.aws_path.str.split("/").str[-1].str.replace(".tsv.bgz", "", regex=False)

Zs, Zn = {}, {}
for k in tl.key:
    d = pd.read_csv(f"{BLK}/{k}.tsv", sep="\t").drop_duplicates("id").set_index("id")
    z = d.beta / d.se
    r = pd.read_csv(f"{REG}/{k}.tsv", sep="\t", dtype={"chr": str})
    r["id"] = r.chr + ":" + r.pos.astype(str) + ":" + r.ref + ":" + r.alt
    x = r[r.id.isin(EXTRA)].drop_duplicates("id").set_index("id")
    z = pd.concat([z, x.beta / x.se])
    Zs[k] = z
    # near-null SNPs for the trait-level null correlation (as in panukb_motivating.py)
    rr = r[r.region.str.startswith("rand") & r.af.between(0.05, 0.95)].drop_duplicates("id").set_index("id")
    Zn[k] = rr.beta / rr.se
Z = pd.DataFrame(Zs)[list(tl.key)]
ok = Z.notna().mean(1) > 0.98
print("block SNPs", Z.shape[0], "kept (>=98% non-missing)", int(ok.sum()))
Z = Z[ok].fillna(0.0)

N = pd.DataFrame(Zn)[list(tl.key)]
N = N[N.notna().mean(1) > 0.98].fillna(0.0)
N = N[(N ** 2).mean(1) < 1.5]
C_hat = np.corrcoef(N.values.T)
q = C_hat.shape[0]
print("near-null SNPs", N.shape[0], "q_eff", round(q ** 2 / np.sum(C_hat ** 2), 1))
pd.DataFrame(C_hat, index=tl.key, columns=tl.key).to_csv(f"{OUT}/C_hat.tsv.gz", sep="\t", float_format="%.5g")

# LD per block, restricted to SNPs with Z
chrom = Z.index.str.split(":").str[0].astype(int); pos = Z.index.str.split(":").str[1].astype(int)
keep_ids = []
for name, (c, s, e) in BLOCKS.items():
    f = [g for g in glob.glob(f"{LD}/chr{c}_*.snps.tsv")]
    assert len(f) == 1, f
    snps = pd.read_csv(f[0], sep="\t")
    m = len(snps)
    R = np.fromfile(f[0].replace(".snps.tsv", ".R.f32"), dtype=np.float32).reshape(m, m)
    sel = np.where((snps.BP >= s) & (snps.BP <= e) & snps.ID.isin(Z.index))[0]
    snps.iloc[sel].to_csv(f"{OUT}/ld/chr{c}_{s}_{e}.snps.tsv", sep="\t", index=False)
    R[np.ix_(sel, sel)].astype(np.float32).tofile(f"{OUT}/ld/chr{c}_{s}_{e}.R.f32")
    keep_ids += list(snps.ID.iloc[sel])
    print(name, f"chr{c}:{s}-{e}", "SNPs with Z and LD:", len(sel))
Z = Z.loc[keep_ids]
o = np.lexsort((Z.index.str.split(":").str[1].astype(int), Z.index.str.split(":").str[0].astype(int)))
Z = Z.iloc[o]
Z.to_csv(f"{OUT}/Z.tsv.gz", sep="\t", float_format="%.5g")
tl[["key", "description", "category_top", "N_eff", "h2_obs", "ldsc_intercept", "indep"]].to_csv(
    f"{OUT}/traits.tsv", sep="\t", index=False)
print("Z", Z.shape, "targets present:", [i for i in EXTRA if i in Z.index])
